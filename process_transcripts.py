import pandas as pd
import geopandas as gpd
import numpy as np
import scipy.signal as sps
from scipy.ndimage import label, maximum_filter
from shapely.geometry import Point, Polygon
from shapely import union_all
import matplotlib.pyplot as plt
import argparse
from skimage.draw import line
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import networkx as nx

def find_density_peaks(df_for_peaks, voxel_size, density_radius, min_transcripts_per_peak, background_sd_multiplier):
    """
    Finds high-density peaks in transcript locations using a 2D histogram and convolution.
    Uses 'x_location' and 'y_location' columns.
    Returns peak locations, the full density map, and a coordinate mapping dictionary.
    """
    print("Finding density peaks...")
    xmin, ymin = df_for_peaks[["x_location", "y_location"]].min().values
    xmax, ymax = df_for_peaks[["x_location", "y_location"]].max().values
    nx = int(np.ceil((xmax - xmin) / voxel_size))
    ny = int(np.ceil((ymax - ymin) / voxel_size))
    H, _, _ = np.histogram2d(
        df_for_peaks["y_location"], df_for_peaks["x_location"],
        bins=[ny, nx],
        range=[[ymin, ymax], [xmin, xmax]]
    )
    r_pix = int(round(density_radius / voxel_size))
    kernel = np.zeros((2 * r_pix + 1, 2 * r_pix + 1), dtype=np.uint8)
    yy, xx = np.ogrid[-r_pix:r_pix + 1, -r_pix:r_pix + 1]
    kernel[xx**2 + yy**2 <= r_pix**2] = 1

    # --- Automatic threshold calculation ---
    if min_transcripts_per_peak == 0:
        print("Calculating background density using convex hull...")
        # Create a convex hull around all transcripts to define the tissue area
        points = df_for_peaks[['x_location', 'y_location']].values
        hull = ConvexHull(points)
        hull_path = Path(points[hull.vertices])

        # Create a grid of voxel center coordinates
        vx, vy = np.meshgrid(
            xmin + voxel_size * (np.arange(nx) + 0.5),
            ymin + voxel_size * (np.arange(ny) + 0.5)
        )
        voxel_coords = np.vstack([vx.ravel(), vy.ravel()]).T

        # Create a mask for voxels inside the convex hull
        mask_in_hull = hull_path.contains_points(voxel_coords).reshape(ny, nx)
        
        voxels_in_hull = H[mask_in_hull]

        if len(voxels_in_hull) > 0:
            # Use mean of voxels within the hull for background estimation
            mean_voxel_count = np.mean(voxels_in_hull)
            std_voxel_count = np.std(voxels_in_hull)
            thresh_voxel_count = mean_voxel_count + (background_sd_multiplier * std_voxel_count)
            num_voxels_in_kernel = np.sum(kernel)
            expected_transcripts_in_radius = thresh_voxel_count * num_voxels_in_kernel
            min_transcripts_per_peak = expected_transcripts_in_radius
            print(f"Mean background transcripts per voxel: {mean_voxel_count:.2f}")
            print(f"SD background transcripts per voxel: {std_voxel_count:.2f}")
            print(f"Threshold per voxel (mean + {background_sd_multiplier}*SD): {thresh_voxel_count:.2f}")
            print(f"Voxels in density radius: {num_voxels_in_kernel}")
            print(f"Automatically calculated MIN_TRANSCRIPTS_PER_PEAK: {min_transcripts_per_peak:.2f}")
        else:
            print("Warning: No voxels found within the convex hull. Using a fallback.")
            min_transcripts_per_peak = 100  # Fallback value

    print("Convolving histogram with kernel to create density map...")
    density = sps.convolve(H, kernel, mode="same")

    mask = density >= min_transcripts_per_peak
    local_max = density == maximum_filter(density, footprint=kernel, mode="nearest")
    peaks = mask & local_max
    lbl, n_lbl = label(peaks, structure=np.ones((3, 3)))
    coord_map = {'xmin': xmin, 'ymin': ymin, 'voxel_size': voxel_size, 'nx': nx, 'ny': ny}

    if n_lbl == 0:
        print("No peaks found meeting the criteria.")
        return None, density, coord_map

    peak_coords = np.column_stack(np.where(lbl > 0))
    peak_xy = np.c_[xmin + peak_coords[:, 1] * voxel_size + voxel_size / 2,
                    ymin + peak_coords[:, 0] * voxel_size + voxel_size / 2]
    print(f"Found {n_lbl} initial density peaks.")
    return peak_xy, density, coord_map

def merge_overlapping_nuclei(initial_nuclei, peak_xy):
    """
    Stage 1: Geometrically merge any nuclei that physically overlap.
    Returns a list of merged polygons and a mapping from the new polygons back to their original constituent peaks.
    """
    print("Stage 1: Merging overlapping nuclei...")
    if not initial_nuclei:
        return [], {}

    gdf = gpd.GeoDataFrame(geometry=initial_nuclei)
    overlaps_df = gpd.sjoin(gdf, gdf, how="inner", predicate="intersects")

    G = nx.Graph()
    G.add_nodes_from(range(len(initial_nuclei)))
    edge_list = overlaps_df[overlaps_df.index != overlaps_df.index_right][['index_right']].reset_index().values
    G.add_edges_from(edge_list)

    components = list(nx.connected_components(G))
    
    merged_polygons = []
    # This dictionary will map the index of a merged polygon to the indices of the original peaks it contains
    merged_to_original_peaks = {}

    for i, comp in enumerate(components):
        # Get the actual polygon objects for the component
        polygons_to_merge = [initial_nuclei[j] for j in comp]
        # Use union_all for efficient merging
        merged_poly = union_all(polygons_to_merge)
        merged_polygons.append(merged_poly)
        # Store the mapping from this new polygon's index to the original peak indices
        merged_to_original_peaks[i] = list(comp)

    print(f"Merged {len(initial_nuclei)} initial nuclei into {len(merged_polygons)} non-overlapping groups.")
    return merged_polygons, merged_to_original_peaks

def filter_nuclei_by_density(merged_nuclei, merged_to_original_peaks, peak_xy, density_map, coord_map):
    """
    Stage 2: Filters the merged groups based on the density of the path between their peaks.
    """
    print("Stage 2: Filtering merged groups by inter-peak density...")
    if not merged_nuclei or len(merged_nuclei) < 2:
        if not merged_nuclei:
            return None
        # If there's only one group, no filtering is needed.
        nuclei_gdf = gpd.GeoDataFrame(geometry=merged_nuclei)
        nuclei_gdf['cell_id'] = [f'anucleated_{i+1}' for i in range(len(nuclei_gdf))]
        return nuclei_gdf

    from scipy.spatial import cKDTree

    # Helper to convert real-world coords to voxel indices
    def world_to_voxel(xy):
        vx = int((xy[0] - coord_map['xmin']) / coord_map['voxel_size'])
        vy = int((xy[1] - coord_map['ymin']) / coord_map['voxel_size'])
        return vy, vx

    # Get original peak densities
    peak_voxels_all = [world_to_voxel(p) for p in peak_xy]
    peak_densities_all = [density_map[vy, vx] for vy, vx in peak_voxels_all]

    # We need to find the representative peak for each merged group (the densest one)
    representative_peaks = []
    for i in range(len(merged_nuclei)):
        original_peak_indices = merged_to_original_peaks[i]
        if not original_peak_indices: continue
        
        densest_original_peak_idx = max(original_peak_indices, key=lambda idx: peak_densities_all[idx])
        representative_peaks.append({
            'merged_idx': i,
            'peak_coord': peak_xy[densest_original_peak_idx],
            'peak_density': peak_densities_all[densest_original_peak_idx]
        })

    if not representative_peaks:
        return None

    rep_peak_coords = np.array([p['peak_coord'] for p in representative_peaks])
    rep_peak_densities = [p['peak_density'] for p in representative_peaks]
    
    G = nx.Graph()
    G.add_nodes_from(range(len(representative_peaks)))

    search_radius = 2 * (int(round(coord_map.get('density_radius', 12.0) / coord_map['voxel_size'])) * coord_map['voxel_size'])
    tree = cKDTree(rep_peak_coords)
    pairs = tree.query_pairs(r=search_radius)

    print(f"Checking {len(pairs)} potential connections between merged groups...")
    for i, j in pairs:
        y0, x0 = world_to_voxel(rep_peak_coords[i])
        y1, x1 = world_to_voxel(rep_peak_coords[j])
        rr, cc = line(y0, x0, y1, x1)
        
        valid_indices = (rr >= 0) & (rr < density_map.shape[0]) & (cc >= 0) & (cc < density_map.shape[1])
        rr, cc = rr[valid_indices], cc[valid_indices]

        if len(rr) < 3: continue

        path_voxels = density_map[rr[1:-1], cc[1:-1]]
        min_peak_density = min(rep_peak_densities[i], rep_peak_densities[j])
        connection_threshold = 0.5 * min_peak_density

        if np.all(path_voxels > connection_threshold):
            G.add_edge(i, j)

    components = list(nx.connected_components(G))
    final_polygons = []
    print(f"Found {len(components)} final components from {len(representative_peaks)} merged groups.")

    for comp in components:
        if len(comp) > 1:
            densest_rep_peak_idx = max(comp, key=lambda idx: rep_peak_densities[idx])
            original_merged_idx = representative_peaks[densest_rep_peak_idx]['merged_idx']
            final_polygons.append(merged_nuclei[original_merged_idx])
        else:
            original_merged_idx = representative_peaks[list(comp)[0]]['merged_idx']
            final_polygons.append(merged_nuclei[original_merged_idx])

    nuclei_gdf = gpd.GeoDataFrame(geometry=final_polygons)
    nuclei_gdf['cell_id'] = [f'anucleated_{i+1}' for i in range(len(nuclei_gdf))]
    print(f"Created {len(nuclei_gdf)} final filtered nuclei.")
    return nuclei_gdf

def merge_and_filter_nuclei(initial_nuclei, peak_xy, density_map, coord_map):
    """
    Orchestrates the two-stage merging and filtering process.
    """
    if not initial_nuclei or peak_xy is None:
        return None
    
    # Stage 1: Merge overlapping polygons
    merged_nuclei, merged_to_original_peaks = merge_overlapping_nuclei(initial_nuclei, peak_xy)

    # Stage 2: Filter the merged groups based on density
    final_nuclei_gdf = filter_nuclei_by_density(merged_nuclei, merged_to_original_peaks, peak_xy, density_map, coord_map)
    
    return final_nuclei_gdf

def create_nuclei_from_peaks(peak_xy, transcripts_df, density_radius, nucleus_radius):
    """
    For each peak, finds the true center of density and creates a circular nucleus.
    """
    print("Creating initial nuclei from true density centers...")
    if peak_xy is None:
        return []
    from scipy.spatial import cKDTree
    transcripts_gdf = gpd.GeoDataFrame(
        transcripts_df,
        geometry=gpd.points_from_xy(transcripts_df.x_location, transcripts_df.y_location)
    )
    initial_nuclei = []
    for peak in peak_xy:
        # Step 1: Find the local cloud of transcripts around the coarse peak
        density_circle = Point(peak).buffer(density_radius)
        local_cloud = transcripts_gdf[transcripts_gdf.within(density_circle)]
        if len(local_cloud) < 3:
            continue
        # Step 2: Find the point with the most neighbors in a small radius (true density center)
        tree = cKDTree(local_cloud[['x_location', 'y_location']].values)
        # Count neighbors within a 2um radius for each point
        counts = tree.query_ball_point(local_cloud[['x_location', 'y_location']].values, r=2.0, return_length=True)
        if counts.max() > 0:
            # Get the coordinate of the point with the most neighbors
            true_center_coord = local_cloud.iloc[counts.argmax()][['x_location', 'y_location']].values
            # Step 3: Create the 3um nucleus around this true center
            nucleus = Point(true_center_coord).buffer(nucleus_radius)
            initial_nuclei.append(nucleus)
    print(f"Created {len(initial_nuclei)} initial nuclei before merging.")
    return initial_nuclei

def update_all_transcripts(transcripts_df, nuclei_gdf):
    """
    Updates the provided transcript DataFrame with new cell_id and distance info.
    This function now performs a two-step update as requested.
    """
    print("Updating transcript information...")
    if nuclei_gdf is None or nuclei_gdf.empty:
        print("No nuclei found, returning original transcript data.")
        return transcripts_df
    # Create a GeoDataFrame from our working set, using the correct column names
    transcripts_gdf = gpd.GeoDataFrame(
        transcripts_df,
        geometry=gpd.points_from_xy(transcripts_df.x_location, transcripts_df.y_location)
    )
    # Step A: Find transcripts INSIDE nuclei
    print("Step A: Identifying transcripts inside nuclei...")
    transcripts_inside = gpd.sjoin(transcripts_gdf, nuclei_gdf, how="inner", predicate='within')
    if not transcripts_inside.empty:
        print(f"Found {len(transcripts_inside)} transcripts inside nuclei.")
        # Use .loc with the index of transcripts_inside to update the main DataFrame
        transcripts_df.loc[transcripts_inside.index, 'cell_id'] = transcripts_inside['cell_id_right']
        transcripts_df.loc[transcripts_inside.index, 'overlaps_nucleus'] = 1
        transcripts_df.loc[transcripts_inside.index, 'nucleus_distance'] = 0.0
    else:
        print("No transcripts found inside any nucleus.")
    # Step B: Find distance for transcripts OUTSIDE nuclei
    print("Step B: Calculating distance for transcripts outside nuclei...")
    outside_mask = ~transcripts_df.index.isin(transcripts_inside.index)
    if outside_mask.any():
        outside_gdf = transcripts_gdf.loc[outside_mask]
        print(f"Calculating distances for {len(outside_gdf)} outside transcripts using sjoin_nearest...")

        # Use sjoin_nearest for efficient distance calculation.
        # This performs a spatial join and calculates the distance to the nearest nucleus in one step.
        # We use a left join to ensure all 'outside' transcripts are kept.
        outside_with_distances = gpd.sjoin_nearest(outside_gdf, nuclei_gdf, how='left', distance_col='nucleus_distance')

        # The join might create duplicate indices if a transcript is equidistant
        # from multiple nuclei. We'll keep the first match.
        outside_with_distances = outside_with_distances[~outside_with_distances.index.duplicated(keep='first')]

        # Update the main DataFrame using the index from the result
        transcripts_df.loc[outside_with_distances.index, 'nucleus_distance'] = outside_with_distances['nucleus_distance']
    else:
        print("No transcripts found outside any nucleus.")
    print("Finished updating transcripts.")
    return transcripts_df

def generate_plot(df_to_plot, nuclei_gdf, output_plot_file):
    """
    Generates and saves a plot of background transcripts and the final nuclei.
    """
    print(f"Generating plot and saving to {output_plot_file}...")
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Plot background transcripts
    ax.scatter(df_to_plot['x_location'], df_to_plot['y_location'], s=0.005,edgecolors='none', c='gray', label='Background Transcripts', alpha=1)
    
    # Plot nuclei
    if nuclei_gdf is not None and not nuclei_gdf.empty:
        nuclei_gdf.plot(ax=ax, edgecolor='red', facecolor='none', linewidth=0.5)
    ax.set_title('Background Transcripts and Final Anucleated Nuclei')
    ax.set_xlabel('X coordinate (µm)')
    ax.set_ylabel('Y coordinate (µm)')
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(output_plot_file, dpi=1200, bbox_inches='tight')
    print("Plot saved.")

def main():
    """
    Main function to execute the transcript processing pipeline.
    """
    parser = argparse.ArgumentParser(description="Find and process anucleated nuclei from transcript data.")
    parser.add_argument('--transcript-metadata-file', type=str, default='transcript-metadata.csv.gz', help='Path to the transcript metadata file.')
    parser.add_argument('--transcripts-file', type=str, default='transcripts.parquet', help='Path to the transcripts parquet file.')
    parser.add_argument('--output-parquet-file', type=str, default='transcripts_anuc.parquet', help='Path to save the output parquet file.')
    parser.add_argument('--output-plot-file', type=str, default='anucleated_nuclei_plot.png', help='Path to save the output plot.')
    parser.add_argument('--voxel-size', type=float, default=3.0, help='Voxel size in µm for density calculation.')
    parser.add_argument('--density-radius', type=float, default=12.0, help='Radius in µm for density calculation.')
    parser.add_argument('--min-transcripts-per-peak', type=int, default=0, help='Minimum number of transcripts to call a peak. Set to 0 for automatic calculation.')
    parser.add_argument('--nucleus-radius', type=float, default=3.0, help='Radius in µm for creating nuclei.')
    parser.add_argument('--background-sd-multiplier', type=float, default=1.0, help='Multiplier for the standard deviation in automatic threshold calculation.')

    args = parser.parse_args()

    print("--- Starting Anucleated Nuclei Identification Pipeline ---")

    # 1. Load and Filter Data for Peak Finding
    print(f"Loading transcript metadata from {args.transcript_metadata_file}...")
    try:
        transcript_meta_df = pd.read_csv(args.transcript_metadata_file) #note the locations in this are scrambled due to allowable transcript diffusion in the proseg model!
    except FileNotFoundError:
        print(f"Error: {args.transcript_metadata_file} not found.")
        return

    background_transcripts = transcript_meta_df[transcript_meta_df['background'] == 1].copy()
    if background_transcripts.empty:
        print("No background transcripts found. Exiting.")
        return

    # 4. Load Full Transcript Dataset and Subset It
    print(f"Loading full transcript data from {args.transcripts_file}...")
    try:
        all_transcripts_df = pd.read_parquet(args.transcripts_file)
    except FileNotFoundError:
        print(f"Error: {args.transcripts_file} not found.")
        return
    
    # Get the set of IDs to keep
    background_ids = set(background_transcripts['transcript_id'])
    
    # Create the working copy by subsetting the original parquet file
    working_df = all_transcripts_df[all_transcripts_df['transcript_id'].isin(background_ids)].copy()
    print(f"Subsetted to {len(working_df)} background transcripts, preserving original schema.")

    # Reset existing cell IDs to ensure a clean slate
    print("Resetting existing cell IDs to 'UNASSIGNED'...")
    working_df['cell_id'] = 'UNASSIGNED'

    # Find Density Peaks using the consistent working dataframe
    coarse_peaks, density_map, coord_map = find_density_peaks(working_df, args.voxel_size, args.density_radius, args.min_transcripts_per_peak, args.background_sd_multiplier)

    # Create initial nuclei centered on the true point of highest density
    initial_nuclei = create_nuclei_from_peaks(coarse_peaks, working_df, args.density_radius, args.nucleus_radius)

    # Merge nuclei based on the new density-aware logic
    if coarse_peaks is not None:
        # Pass the necessary data to the new merging function
        coord_map['density_radius'] = args.density_radius # Pass this for the search radius calculation
        nuclei_gdf = merge_and_filter_nuclei(initial_nuclei, coarse_peaks, density_map, coord_map)
    else:
        nuclei_gdf = None # No peaks, so no nuclei

    # Update The Subsetted DataFrame
    updated_df = update_all_transcripts(working_df, nuclei_gdf)

    # Save Final Output
    print(f"Saving updated transcript data to {args.output_parquet_file}...")
    # Ensure data types are consistent before saving
    updated_df['nucleus_distance'] = updated_df['nucleus_distance'].astype(np.float32)
    
    # # Reorder columns to match the original file exactly, just in case
    # updated_df = updated_df[all_transcripts_df.columns]
    updated_df.to_parquet(args.output_parquet_file, index=False)

    # 7. Generate Plot
    # The plot should show the locations from the final data, which use x_location/y_location
    generate_plot(updated_df, nuclei_gdf, args.output_plot_file)

    print("--- Pipeline Finished Successfully ---")

if __name__ == '__main__':
    main()
import argparse
import pathlib
import gzip
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    """
    Main function to process and merge proseg outputs.
    """
    parser = argparse.ArgumentParser(description="Process and merge proseg outputs.")
    parser.add_argument("--original_dir", type=pathlib.Path, required=True, help="Path to the original proseg output directory.")
    parser.add_argument("--anucleated_dir", type=pathlib.Path, required=True, help="Path to the anucleated proseg output directory.")
    parser.add_argument("--output_dir", type=pathlib.Path, required=True, help="Path to the output directory.")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Load Data ---
    print("Loading data...")
    # Load cell polygons
    with gzip.open(args.original_dir / "union-cell-polygons.geojson.gz", "rt") as f:
        orig_polys = gpd.read_file(f)
    with gzip.open(args.anucleated_dir / "union-cell-polygons.geojson.gz", "rt") as f:
        anuc_polys = gpd.read_file(f)

    # Set CRS to None to avoid geographic projection issues
    orig_polys.crs = None
    anuc_polys.crs = None

    # Load transcript metadata
    orig_tx = pd.read_csv(args.original_dir / "transcript-metadata.csv.gz", compression='gzip')
    anuc_tx = pd.read_csv(args.anucleated_dir / "transcript-metadata.csv.gz", compression='gzip')

    # Load cell metadata
    orig_meta = pd.read_csv(args.original_dir / "cell-metadata.csv.gz", compression='gzip')
    anuc_meta = pd.read_csv(args.anucleated_dir / "cell-metadata.csv.gz", compression='gzip')

    # --- Sanity Check ---
    print("Performing sanity check on cell metadata and polygons...")
    from shapely.geometry import Point
    
    # Check that the x,y from metadata is within the polygon at the same index
    orig_mismatches = sum(
        not poly.contains(Point(row['centroid_x'], row['centroid_y']))
        for i, (poly, (_, row)) in enumerate(zip(orig_polys.geometry, orig_meta.iterrows()))
    )
    if orig_mismatches > 0:
        print(f"WARNING: Found {orig_mismatches} mismatched coordinates in the original dataset.")

    anuc_mismatches = sum(
        not poly.contains(Point(row['centroid_x'], row['centroid_y']))
        for i, (poly, (_, row)) in enumerate(zip(anuc_polys.geometry, anuc_meta.iterrows()))
    )
    if anuc_mismatches > 0:
        print(f"WARNING: Found {anuc_mismatches} mismatched coordinates in the anucleated dataset.")

    # --- 2. Filter and Plot ---
    print("Generating scatter plot...")
    # Filter original transcripts
    orig_tx_filtered = orig_tx[orig_tx['background'] == 0]

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(anuc_tx[anuc_tx['background'] == 1]['x'], anuc_tx[anuc_tx['background'] == 1]['y'], s=0.05, linewidths=0, color='grey', label='Background')
    ax.scatter(orig_tx_filtered['x'], orig_tx_filtered['y'], s=0.05, linewidths=0, color='cornflowerblue', label='Nucleated')
    ax.scatter(anuc_tx[anuc_tx['background'] == 0]['x'], anuc_tx[anuc_tx['background'] == 0]['y'], s=0.05, linewidths=0, color='forestgreen', label='Anucleated')
    
    ax.set_aspect("equal")
    orig_polys.boundary.plot(ax=ax, linewidth=0.1, color="blue", label="Nucleated Cells")
    anuc_polys.boundary.plot(ax=ax, linewidth=0.1, color="lime", label="Anucleated Cells")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.legend(markerscale=10, frameon=False, loc="upper right")
    plt.tight_layout()
    plt.savefig(args.output_dir / "summary_scatterplot.png", dpi=4800)
    print(f"Saved scatter plot to {args.output_dir / 'summary_scatterplot.png'}")

    # --- 3. Calculate Morphology Features ---
    print("Calculating morphology features...")

    def calculate_morphology_features(gdf):
        """Calculates morphology features for a GeoDataFrame."""
        features = pd.DataFrame(index=gdf.index)
        features['area'] = gdf.geometry.area
        features['perimeter'] = gdf.geometry.length
        # Calculate solidity, handling potential division by zero for invalid geometries
        convex_hull_area = gdf.geometry.convex_hull.area
        features['solidity'] = features['area'] / convex_hull_area.where(convex_hull_area > 0, np.nan)
        # Calculate circularity, handling potential division by zero
        features['circularity'] = (4 * np.pi * features['area']) / (features['perimeter']**2).where(features['perimeter'] > 0, np.nan)
        return features

    orig_features = calculate_morphology_features(orig_polys)
    anuc_features = calculate_morphology_features(anuc_polys)

    # --- 4. Load and Combine All Data ---
    print("Loading and combining all data sources...")
    orig_counts = pd.read_csv(args.original_dir / "expected-counts.csv.gz", compression='gzip')
    anuc_counts = pd.read_csv(args.anucleated_dir / "expected-counts.csv.gz", compression='gzip')
    # --- Gene list consistency check ---
    print("Checking for gene list consistency...")
    orig_genes = set(orig_counts.columns)
    anuc_genes = set(anuc_counts.columns)

    if orig_genes != anuc_genes:
        missing_in_anuc = orig_genes - anuc_genes
        extra_in_anuc = anuc_genes - orig_genes
        error_message = "Gene lists do not match between original and anucleated datasets.\n"
        if missing_in_anuc:
            error_message += f"Genes in original but not in anucleated: {missing_in_anuc}\n"
        if extra_in_anuc:
            error_message += f"Genes in anucleated but not in original: {extra_in_anuc}\n"
        raise ValueError(error_message)

    # Ensure the column order of anuc_counts matches orig_counts
    anuc_counts = anuc_counts[orig_counts.columns]

    # Concatenate all data sources by index (side-by-side)
    # Order: metadata, counts, morphology features
    orig_merged = pd.concat([orig_meta, orig_counts, orig_features], axis=1)
    anuc_merged = pd.concat([anuc_meta, anuc_counts, anuc_features], axis=1)

    # Add 'Nucleated' column
    orig_merged['Nucleated'] = 'Nucleated'
    anuc_merged['Nucleated'] = 'Anucleated'

    # --- 5. Combine and Export ---
    print("Combining and exporting final data...")
    final_df = pd.concat([orig_merged, anuc_merged], ignore_index=True)
    output_path = args.output_dir / "expected-counts-complete.csv"
    final_df.to_csv(output_path, index=False)
    print(f"Saved final merged data to {output_path}")

if __name__ == "__main__":
    main()
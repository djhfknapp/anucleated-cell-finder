# Anucleated Nuclei Identification Pipeline

This script is designed to solve a specific challenge in 3D spatial transcriptomics (e.g., 10x Xenium): identifying cells whose nuclei were not captured in the imaging plane but whose cell bodies (and thus, transcript clouds) are still present. These "anucleated" cells are often missed by standard segmentation workflows that rely on a DAPI-stained nucleus as a seed.

This pipeline identifies these cells by finding high-density "transcript shadows" in regions previously classified as background. It takes the outputs from a standard Xenium Ranger and [proseg](https://github.com/dcjones/proseg/tree/main) run, identifies new cell candidates, and generates a new transcript file. This new file can then be used as input for a secondary `proseg` run to fully segment these previously missed cells.

The pipeline performs the following key steps:
1.  **Density Map Generation:** Creates a 2D density map from transcript locations.
2.  **Automatic Peak Thresholding:** Intelligently calculates a minimum density threshold to distinguish true peaks from background noise.
3.  **Peak Identification:** Finds all local density maxima that exceed the calculated threshold.
4.  **Nucleus Creation:** Creates initial circular nucleus candidates around the true center of density for each peak.
5.  **Two-Stage Merging & Filtering:** A robust, two-stage process first merges any physically overlapping nuclei and then filters the resulting groups based on the density of the space between them, ensuring only the most prominent peak in a connected region is kept.
6.  **Transcript Assignment:** Assigns transcripts to the newly identified nuclei and calculates their distance to the nearest nucleus.
7.  **Output Generation:** Saves the updated transcript data to a new Parquet file and generates a plot visualizing the final nuclei.

## Analysis Pipeline Workflow

The following diagram illustrates the data flow and decision-making process within the script:

```mermaid
graph TD
    A[Start] --> B{Load Transcript Data};
    B --> C{Generate Convex Hull of All Transcripts};
    B --> D{Create 2D Density Map};
    C & D --> E{Calculate Mean Density of Voxels within Hull};
    E --> F{Set MIN_TRANSCRIPTS_PER_PEAK};
    D & F --> G{Find Density Peaks};
    G --> H{Create Initial Nuclei from Peaks};
    
    subgraph Merging [Two-Stage Merging]
      direction LR
      H --> M1{Stage 1: Merge Overlapping Nuclei with union_all};
      M1 --> M2{Stage 2: Filter Merged Groups by Inter-Peak Density};
    end

    Merging --> K{Assign Transcripts to Final Nuclei};
    K --> L{Save Output Parquet & Plot};
    L --> M[End];

    style E fill:#f9f,stroke:#333,stroke-width:2px
    style Merging fill:#bbf,stroke:#333,stroke-width:2px
```

## Workflow Context

The primary goal of this tool is to prepare data for a second `proseg` run.

*   **2D Analysis**: The current implementation works in 2D (`x` and `y` coordinates). Therefore, the subsequent `proseg` run should be performed using the `--ignore-z-coord` option to ensure compatibility.
*   **Input**: The script requires the `transcripts.parquet` and `transcript-metadata.csv.gz` files from a previous segmentation run.
*   **Output**: The output is a modified `transcripts_anuc.parquet` file where the newly identified cell candidates are assigned a `cell_id`. This file is ready to be used as input for a new `proseg` run.

## Installation

To set up the necessary environment, install the required Python packages using the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
```

## Usage

The script is executed from the command line. All parameters are configurable via command-line arguments.

```bash
python process_transcripts.py [OPTIONS]
```

### Parameters

#### File Paths
*   `--transcript-metadata-file`: Path to the transcript metadata CSV file (e.g., `transcript-metadata.csv.gz`). Default: `transcript-metadata.csv.gz`.
*   `--transcripts-file`: Path to the main transcripts Parquet file. Default: `transcripts.parquet`.
*   `--output-parquet-file`: Path to save the updated transcripts data. Default: `transcripts_anuc.parquet`.
*   `--output-plot-file`: Path to save the output visualization plot. Default: `anucleated_nuclei_plot.png`.

#### Core Algorithm Parameters
*   `--voxel-size`: The size of each square voxel in the 2D histogram, in micrometers (µm). This parameter controls the resolution of the density map. Default: `3.0`.
*   `--density-radius`: The radius around each voxel used to calculate the local transcript density, in micrometers (µm). Default: `12.0`.
*   `--nucleus-radius`: The radius of the final circular nuclei created around the center of density, in micrometers (µm). Default: `3.0`.
*   `--min-transcripts-per-peak`: The minimum number of transcripts required within the `density-radius` to consider a voxel a potential peak. **Set to 0 for automatic calculation (recommended)**. Default: `0`.
*   `--background-sd-multiplier`: A multiplier for the standard deviation of the background density. Used only during automatic threshold calculation. Default: `1.0`.

### Parameter Interactions and Tuning

The most critical and interactive parameters are `--voxel-size`, `--density-radius`, and `--background-sd-multiplier`. Understanding how they work together is key to tuning the algorithm.

1.  **`--voxel-size`**: This sets the fundamental resolution of your analysis.
    *   A **smaller** value (e.g., 1.0) creates a higher-resolution density map. This can be good for separating closely spaced peaks but may increase noise and computation time.
    *   A **larger** value (e.g., 5.0) creates a lower-resolution, smoother map, which can be faster but may merge distinct nearby peaks.

2.  **`--density-radius`**: This defines the "neighborhood" for counting transcripts. It is converted from µm into pixels based on the `voxel-size`. For example, a `density-radius` of 12.0 with a `voxel-size` of 3.0 means the kernel will have a radius of 4 pixels.

3.  **`--background-sd-multiplier`**: This directly controls the sensitivity of the peak detection during automatic thresholding. The threshold is calculated as:
    `Threshold = (Mean + (Multiplier * SD)) * NumVoxelsInKernel`
    *   `Mean` and `SD` are the mean and standard deviation of transcript counts in voxels within the tissue's convex hull.
    *   A **higher** multiplier (e.g., 2.0, 3.0) will result in a higher threshold, leading to the detection of only the most prominent density peaks. Use this if you are getting too many false positives.
    *   A **lower** multiplier (e.g., 0.5, 1.0) will lower the threshold, making the algorithm more sensitive and allowing it to detect weaker peaks. Use this if the algorithm is missing real nuclei.

**How they interact:** The automatic threshold (`min_transcripts_per_peak`) is a product of the per-voxel background estimate and the number of voxels in the density kernel.
*   If you **decrease** `voxel-size`, the number of voxels within the `density-radius` increases quadratically. The algorithm compensates for this by calculating the expected number of transcripts within this larger kernel, making the threshold robust to changes in resolution.
*   The `background-sd-multiplier` is the ultimate lever for sensitivity. After setting reasonable `voxel-size` and `density-radius` values for your tissue, this multiplier should be your primary tool for tuning the number of nuclei detected.

## Future Development

A future version of this pipeline is planned to incorporate the Z-dimension for a full 3D analysis, which will eliminate the need for the `--ignore-z-coord` flag in the secondary `proseg` run.
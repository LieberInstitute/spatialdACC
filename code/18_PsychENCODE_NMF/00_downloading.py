import anndata
import pandas as pd

def export_raw_X_to_csv(dataset: str, file_path: str) -> None:
    """
    Export the raw X counts matrix from an h5ad file to a CSV file.
    
    Args:
    - dataset: Name of the dataset to be used in the output CSV filename.
    - file_path: Path to the .h5ad file.
    """
    adata = anndata.read_h5ad(file_path)
    raw_X = adata.raw.X
    raw_X_df = pd.DataFrame(raw_X.toarray())
    csv_file_path = f'/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/18_PsychENCODE_NMF/raw_{dataset}.csv'
    raw_X_df.to_csv(csv_file_path, index=False)
    print(f"Exported raw X counts matrix to {csv_file_path}")

# List of dataset files
data_list = [
    "CMC_annotated.h5ad", "DevBrain-snRNAseq_annotated.h5ad",
    "IsoHuB_annotated.h5ad", "MultiomeBrain-DLPFC_annotated.h5ad",
    "PTSDBrainomics_annotated.h5ad", "SZBDMulti-Seq_annotated.h5ad",
    "UCLA-ASD_annotated_mismatches_removed.h5ad"
]

# File path base
file_path_base = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/"

# Loop to process each dataset
for dataset_file in data_list:
    dataset_name = dataset_file.split('_')[0]  # Extract dataset name for CSV filename
    file_path = file_path_base + dataset_file
    export_raw_X_to_csv(dataset_name, file_path)


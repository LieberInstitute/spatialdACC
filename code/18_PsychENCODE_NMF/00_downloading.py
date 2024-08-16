import anndata
import pandas as pd

def export_raw_X_to_csv(file_path: str, output_dir: str, dataset_name: str):
    # Read the .h5ad file
    adata = anndata.read_h5ad(file_path)
    
    # Access the raw X matrix
    raw_X = adata.raw.X
    
    # Convert to DataFrame
    raw_X_df = pd.DataFrame(raw_X.toarray())
    
    # Set the output file path
    output_file_path = f'{output_dir}/raw_{dataset_name}.csv'
    
    # Export to CSV
    raw_X_df.to_csv(output_file_path, index=False)
    
    print(f"Exported {dataset_name} to {output_file_path}")

# List of file paths and dataset names
data_list = [
    "CMC_annotated.h5ad", 
    "DevBrain-snRNAseq_annotated.h5ad",
    "IsoHuB_annotated.h5ad",
    "MultiomeBrain-DLPFC_annotated.h5ad", 
    "PTSDBrainomics_annotated.h5ad",
    "SZBDMulti-Seq_annotated.h5ad", 
    "UCLA-ASD_annotated_mismatches_removed.h5ad"
]

# Directory where CSV files will be saved
output_dir = '/processed-data/18_PsychENCODE'

# Run the function for each dataset
for dataset in data_list:
    file_path = f"/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/{dataset}"
    export_raw_X_to_csv(file_path, output_dir, dataset.split('.')[0])

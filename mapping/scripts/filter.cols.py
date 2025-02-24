import pandas as pd

def filter_tsv(input_file, output_file, column_name):
    # Load the TSV file
    df = pd.read_csv(input_file, sep='\t', dtype=str)
    
    # Filter out rows where the column starts with 'XR'
    filtered_df = df[~df[column_name].str.startswith('XR', na=False)]
    
    # Save the filtered data to a new TSV file
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Filtered file saved as: {output_file}")

if __name__ == "__main__":
    input_file = "mapping/Triticum_aestivum.IWGSC.60.refseq.tsv"
    output_file = "mapping/Triticum_aestivum.IWGSC.60.refseq.XM.tsv"
    column_name = "xref"
    
    filter_tsv(input_file, output_file, column_name)
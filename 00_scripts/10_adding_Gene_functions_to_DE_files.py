import pandas as pd
import sys
import os

def merge_dfs(input_file, filtered_file):
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Read the filtered CSV file
    filtered_df = pd.read_csv(filtered_file)
    
    # Merge the DataFrames on the 'Gene' column
    combined_df = pd.merge(df, filtered_df, on="Gene", how="left")
    
    # Define the output file name
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_interpro_functions{ext}"
    
    # Save the combined DataFrame to the new CSV file
    combined_df.to_csv(output_file, index=False)
    
    print(f"Combined data saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 09_adding_Gene_functions_to_DE_files.py input_file.csv filtered_file.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    filtered_file = sys.argv[2]
    merge_dfs(input_file, filtered_file)

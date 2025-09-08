import os
import multiprocessing as mp
from Bio import SeqIO
from functools import partial

# Path to your MAG files folder
mags_folder_path = '/usr/home/qgg/zexi/INCOME/meta_pipeline/kraken_database/version4/mags_data/mags/'
# Path to your taxid.map file
taxid_map_path = '/usr/home/qgg/zexi/INCOME/meta_pipeline/kraken_database/version4/mags_data/taxid.map'
# Output folder where renamed files will be saved
output_folder = '/usr/home/qgg/zexi/INCOME/meta_pipeline/kraken_database/version4/mags_data/mags_corrected_header/'

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Load the taxid.map file into a dictionary
def load_taxid_map(map_path):
    genome_to_taxid = {}
    
    with open(map_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                genome_id, taxid = parts
                genome_to_taxid[genome_id] = taxid
    
    return genome_to_taxid

# Function to rename headers based on the filename and taxid map
def rename_fasta_headers(fasta_file, mags_folder_path, output_folder, genome_to_taxid):
    input_file_path = os.path.join(mags_folder_path, fasta_file)
    output_file_path = os.path.join(output_folder, fasta_file)
    
    # Extract the genome name from the filename (remove extension)
    genome_name = fasta_file.replace('.fa', '')
    
    # Look up the taxid for this genome
    taxid = genome_to_taxid.get(genome_name, 'UnknownTaxID')
    
    try:
        with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
            # Initialize counter for contig numbers
            contig_counter = 1
            
            for record in SeqIO.parse(input_handle, 'fasta'):
                # Create the new standardized header with sequential contig numbers
                new_id = f"{genome_name}_{contig_counter}|kraken:taxid|{taxid}"
                contig_counter += 1
                
                # Update the record
                record.id = new_id
                record.description = ''
                
                # Write the renamed record
                SeqIO.write(record, output_handle, 'fasta')
        
        return f"Successfully processed {fasta_file} - renamed {contig_counter-1} contigs"
    except Exception as e:
        return f"Error processing {fasta_file}: {str(e)}"

def main():
    # Load taxid map
    genome_to_taxid = load_taxid_map(taxid_map_path)
    
    # Get list of FASTA files
    fasta_files = [f for f in os.listdir(mags_folder_path) if f.endswith('.fa')]
    
    # Set up the process pool
    num_cpus = mp.cpu_count()  # Use all available CPUs
    print(f"Using {num_cpus} CPUs for processing")
    
    # Create a partial function with the fixed arguments
    process_file = partial(
        rename_fasta_headers,
        mags_folder_path=mags_folder_path,
        output_folder=output_folder,
        genome_to_taxid=genome_to_taxid
    )
    
    # Process files in parallel
    with mp.Pool(processes=num_cpus) as pool:
        results = pool.map(process_file, fasta_files)
    
    # Print results
    for result in results:
        print(result)
    
    print("Renaming completed. The renamed files are saved in:", output_folder)

if __name__ == "__main__":
    main()
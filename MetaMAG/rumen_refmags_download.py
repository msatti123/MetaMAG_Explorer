import os
import shutil
import requests
import threading
import subprocess
import traceback
import time
import argparse
from bs4 import BeautifulSoup

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
        print(f"[INFO] Created directory: {directory}")

def download_dataset(output_dir, script_dir=None):
    """
    Download reference MAGs from three sources to the specified output directory.
    """
    if script_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    
    dataset_dir = output_dir
    ensure_directory_exists(dataset_dir)
    
    # Find dataset list file
    dataset_list_file = os.path.join(script_dir, "dataset_file_list.txt")
    if not os.path.exists(dataset_list_file):
        print(f"[ERROR] Dataset list file not found: {dataset_list_file}")
        return False

    # Read dataset file names as a list
    with open(dataset_list_file, 'r') as file:
        expected_files = [line.strip() for line in file if line.strip()]

    existing_files = set(os.listdir(dataset_dir)) if os.path.exists(dataset_dir) else set()
    missing_files = [file for file in expected_files if file not in existing_files]

    print(f"[INFO] Total expected files: {len(expected_files)}")
    print(f"[INFO] Existing files: {len(existing_files)}")
    print(f"[INFO] Missing files: {len(missing_files)}")

    # Categorize existing files by their prefix
    rgmgc_existing_files = [f for f in existing_files if f.startswith("GC")]
    mg_existing_files = [f for f in existing_files if f.startswith("MGYG")]
    rug_existing_files = [f for f in existing_files if f.startswith("RUG")]

    # Categorize missing files by their prefix
    rgmgc_missing_files = [f for f in missing_files if f.startswith("GC")]
    mg_missing_files = [f for f in missing_files if f.startswith("MGYG")]
    rug_missing_files = [f for f in missing_files if f.startswith("RUG")]

    # Check if all required files are already present
    if len(rgmgc_existing_files) == 10373 and len(rug_existing_files) == 4941 and len(mg_existing_files) == 2729:
        print("[INFO] All dataset files have already been downloaded. Skipping download phase.")
    else:
        threads = []
        
        if rgmgc_missing_files:
            print("[DEBUG] Starting RGMGC download thread...")
            t1 = threading.Thread(target=download_rgmgc, args=(rgmgc_missing_files, dataset_dir, script_dir))
            threads.append(t1)
        
        if mg_missing_files:
            print("[DEBUG] Starting MG download thread...")
            t2 = threading.Thread(target=download_mgnify, args=(dataset_dir,))
            threads.append(t2)

        if rug_missing_files:
            print("[DEBUG] Starting RUG download thread...")
            t3 = threading.Thread(target=download_rug, args=(dataset_dir,))
            threads.append(t3)
        else:
            print(f"[INFO] RUG dataset files already exist in {dataset_dir}. Skipping download.")

        # Start download threads if required
        for t in threads:
            print(f"[DEBUG] Starting thread {t.name}")
            t.start()

        # Wait for all downloads to finish
        for t in threads:
            t.join()
            print(f"[DEBUG] Thread {t.name} has finished.")

        print("[INFO] Dataset download process completed.")
    
    # Step 1: Delete unwanted files
    delete_unwanted_files(dataset_dir)

    # Step 2: Extract .gz files
    extract_gz_files(dataset_dir)

    # Step 3: Rename all files to .fa
    rename_all_to_fa(dataset_dir)

    print("[INFO] Dataset cleanup, extraction, and renaming completed.")
    return True

def download_rgmgc(files, download_folder, script_dir):
    """Downloads missing RGMGC files from ENA browser."""
    try:
        print(f"[DEBUG] RGMGC download started with {len(files)} files.")

        ena_tsv = os.path.join(script_dir, "ena_mags_list.tsv")
        if not os.path.exists(ena_tsv):
            print(f"[ERROR] RGMGC TSV file not found: {ena_tsv}")
            return

        # Ensure the download folder exists
        os.makedirs(download_folder, exist_ok=True)

        downloaded_count = 0

        with open(ena_tsv, 'r') as file:
            for line in file:
                accession_number = line.strip().split("\t")[0].replace('"', '').replace("'", "")
                
                # Construct the download URL
                url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession_number}?download=true&gzip=true"
                output_file_path = os.path.join(download_folder, f"{accession_number}.fasta.gz")

                # Skip already downloaded files
                if os.path.exists(output_file_path):
                    print(f"[INFO] File already exists: {output_file_path}. Skipping.")
                    continue

                print(f"[INFO] Downloading: {url} ? {output_file_path}")

                try:
                    response = requests.get(url, stream=True)
                    response.raise_for_status()

                    with open(output_file_path, 'wb') as fasta_file:
                        for chunk in response.iter_content(chunk_size=8192):
                            fasta_file.write(chunk)

                    print(f"[INFO] Download complete: {output_file_path}")
                    downloaded_count += 1

                except requests.exceptions.RequestException as e:
                    print(f"[ERROR] Failed to download {url}: {e}")
                    if os.path.exists(output_file_path):
                        os.remove(output_file_path)  # Remove incomplete file

        print(f"[INFO] RGMGC downloads completed. Total downloaded: {downloaded_count}")

    except Exception as e:
        print(f"[ERROR] Exception in download_rgmgc: {e}")
        traceback.print_exc()

def download_rug(dataset_dir):
    """Downloads and verifies the RUG dataset before extracting."""
    rug_url = "https://datashare.ed.ac.uk/bitstream/handle/10283/3224/rug2.tar.gz"
    rug_tar_path = os.path.join(dataset_dir, "rug2.tar.gz")

    os.makedirs(dataset_dir, exist_ok=True)

    # Check if files already exist
    existing_files = [f for f in os.listdir(dataset_dir) if f.startswith("RUG")]
    if existing_files:
        print(f"[INFO] RUG dataset files already exist in {dataset_dir}. Skipping download.")
        return

    # Step 1: Download the file if it doesn't exist
    if not os.path.exists(rug_tar_path):
        print(f"[INFO] Downloading RUG dataset to {rug_tar_path}...")
        process = subprocess.run(f'wget -O "{rug_tar_path}" "{rug_url}"', shell=True, check=False, stderr=subprocess.PIPE)

        if process.returncode != 0:
            print(f"[ERROR] Failed to download RUG dataset. Error: {process.stderr.decode()}")
            return
    else:
        print(f"[INFO] RUG dataset archive already exists: {rug_tar_path}")

    # Step 2: Verify the downloaded file before extraction
    print("[INFO] Verifying the integrity of the downloaded RUG dataset...")
    check_cmd = f'tar -tzf "{rug_tar_path}"'
    process = subprocess.run(check_cmd, shell=True, check=False, stderr=subprocess.PIPE)

    if process.returncode != 0:
        print(f"[ERROR] The RUG dataset archive appears to be corrupted. Re-downloading...")
        os.remove(rug_tar_path)  # Delete corrupted file
        download_rug(dataset_dir)  # Retry download
        return

    # Step 3: Extract only the required files
    print(f"[INFO] Extracting RUG dataset into {dataset_dir}...")
    extract_cmd = f'tar -xz -C "{dataset_dir}" --strip-components=2 --wildcards "*/mags/*" -f "{rug_tar_path}"'
    process = subprocess.run(extract_cmd, shell=True, check=False, stderr=subprocess.PIPE)

    if process.returncode == 0:
        print(f"[INFO] Successfully extracted RUG dataset into {dataset_dir}.")
    else:
        print(f"[ERROR] Failed to extract RUG dataset. Error: {process.stderr.decode()}")

MAX_RETRIES = 5   # Number of retries before skipping a file
RETRY_DELAYS = [5, 10, 20, 40, 60]  # Increasing delays (seconds)

def download_mgnify(download_folder):
    """Downloads missing MGnify files from EBI FTP server."""
    base_url = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/cow-rumen/v1.0.1/species_catalogue/"
    print(f"[DEBUG] Checking MGnify files from {base_url}")

    list_files_mgnify(base_url, ".fna", download_folder)

def list_files_mgnify(url, extension, download_folder):
    """Recursively finds and downloads missing MGnify .fna files."""
    try:
        print(f"[DEBUG] Searching for files at: {url}")

        # GET request with increased timeout to prevent failure
        req = requests.get(url, timeout=60)
        req.raise_for_status()  # Ensure no HTTP errors

        soup = BeautifulSoup(req.text, 'html.parser')

        for link in soup.find_all('a'):
            href = link.get('href')

            # Ensure properly formatted URLs
            full_url = url.rstrip("/") + "/" + href.lstrip("/")

            if href.endswith('/'):  # Recursively process directories
                print(f"[DEBUG] Entering directory: {full_url}")
                time.sleep(1)  # Slow down to avoid being blocked
                list_files_mgnify(full_url, extension, download_folder)

            elif href.endswith(extension) and 'pan-genome.fna' not in href:
                file_name = href.split('/')[-1]
                file_path = os.path.join(download_folder, file_name)

                if not os.path.exists(file_path):
                    print(f"[INFO] Downloading MGnify file: {file_name}")
                    download_file(full_url, file_path)

    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Error accessing {url}: {e}")

def download_file(url, output_file):
    """Downloads a file with retry mechanism."""
    for attempt in range(MAX_RETRIES):
        try:
            print(f"[INFO] Downloading: {url} ? {output_file} (Attempt {attempt+1}/{MAX_RETRIES})")

            with requests.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(output_file, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            print(f"[INFO] Download complete: {output_file}")
            return  # Success! Exit function

        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Failed attempt {attempt+1}: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAYS[attempt])  # Wait before retrying
            else:
                print(f"[ERROR] Max retries reached for {url}. Skipping.")

def delete_unwanted_files(dataset_dir):
    """ Deletes all files that do not start with RUG, GC, or MGYG. """
    allowed_prefixes = ("RUG", "GC", "MGYG")

    for file in os.listdir(dataset_dir):
        file_path = os.path.join(dataset_dir, file)
        if not file.startswith(allowed_prefixes):
            print(f"[INFO] Deleting unwanted file: {file}")
            os.remove(file_path)

def extract_gz_files(dataset_dir):
    """ Extracts all .gz files and removes the original compressed files. """
    for file in os.listdir(dataset_dir):
        if file.endswith(".gz"):
            gz_file = os.path.join(dataset_dir, file)
            extracted_file = gz_file.rstrip(".gz")  # Remove .gz extension

            print(f"[INFO] Extracting {gz_file} ? {extracted_file}")

            try:
                subprocess.run(["gunzip", "-c", gz_file], stdout=open(extracted_file, "wb"), check=True)
                os.remove(gz_file)  # Delete original .gz file
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] Failed to extract {gz_file}: {e}")

def rename_all_to_fa(dataset_dir):
    """ Renames all files to .fa. """
    for file in os.listdir(dataset_dir):
        old_path = os.path.join(dataset_dir, file)
        new_path = os.path.join(dataset_dir, f"{file}.fa")

        if not file.endswith(".fa"):
            print(f"[INFO] Renaming {old_path} ? {new_path}")
            os.rename(old_path, new_path)

def run(output_dir, script_dir=None):
    """Main function to run the dataset download process."""
    print("[INFO] Starting rumen reference MAGs download process")
    success = download_dataset(output_dir, script_dir)
    if success:
        print("[INFO] Rumen reference MAGs download completed successfully")
        return 0
    else:
        print("[ERROR] Rumen reference MAGs download failed")
        return 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download rumen reference MAGs datasets.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save downloaded datasets")
    parser.add_argument("--script_dir", type=str, help="Directory containing supporting files (default: script location)")
    
    args = parser.parse_args()
    exit_code = run(args.output_dir, args.script_dir)
    exit(exit_code)
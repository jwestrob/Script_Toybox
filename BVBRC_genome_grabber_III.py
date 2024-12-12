#!/Users/jacob/.pyenv/shims/python

import os
import sys
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import subprocess
import io
import logging
from typing import Tuple, List, Optional, Dict, Union
from enum import Enum, auto
import requests
from urllib.parse import quote
import json
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('bvbrc_download.log')
    ]
)
logger = logging.getLogger(__name__)

class SearchType(Enum):
    EXACT = '--eq'
    SUBSTRING = '--in'

class SequenceType(Enum):
    NUCLEOTIDE = auto()
    PROTEIN = auto()
    BOTH = auto()

def fetch_genomes_by_lineage(taxonomy_id: str, batch_size: int = 25000, start: int = 0) -> Union[List[Dict], str]:
    """
    Fetch a batch of genome IDs using the BV-BRC API.
    """
    base_url = "https://www.bv-brc.org/api/genome/"
    field = quote("taxon_lineage_names")  # Changed from taxon_lineage_ids
    value = quote(str(taxonomy_id))
    query = f"in({field},({value}))&limit({batch_size},{start})&select(genome_id,genome_name)"
    
    url = f"{base_url}?{query}"
    
    try:
        response = requests.get(url, timeout=300)
        if response.status_code == 200:
            return response.json()
        else:
            return f"Error {response.status_code}: {response.text}"
    except requests.exceptions.RequestException as e:
        return f"Request failed: {str(e)}"

def fetch_all_genomes_api(taxonomy_id: str, batch_size: int = 25000, max_records: Optional[int] = None, threads: int = 1) -> pd.DataFrame:
    """
    Fetch all genome IDs using pagination through the BV-BRC API.
    """
    all_genome_records = []
    start = 0
    concurrent_requests = 100  # Use at most 10 concurrent requests
    empty_batches = 0  # Counter for empty results
    
    logger.info("Starting genome retrieval...")
    
    while True:
        # Create list of start positions for this small batch of parallel requests
        start_positions = []
        for _ in range(concurrent_requests):
            if max_records and start >= max_records:
                break
            start_positions.append(start)
            start += batch_size
            
        if not start_positions:
            break
            
        # Process this small batch of positions
        with ThreadPoolExecutor(max_workers=concurrent_requests) as executor:
            futures = []
            for start_pos in start_positions:
                futures.append(executor.submit(fetch_genomes_by_lineage, taxonomy_id, batch_size, start_pos))
            
            batch_results = []
            with tqdm(total=len(futures), desc=f"Fetching genomes {start-len(start_positions)*batch_size:,}-{start:,}") as pbar:
                for future in concurrent.futures.as_completed(futures):
                    try:
                        batch = future.result()
                        if not isinstance(batch, str) and batch:
                            batch_results.extend(batch)
                            empty_batches = 0  # Reset counter on successful batch
                        else:
                            empty_batches += 1
                    except Exception as e:
                        logger.error(f"Error in batch: {e}")
                        empty_batches += 1
                    pbar.update(1)
        
        # If we got no results in several consecutive batches, we're probably done
        if empty_batches >= concurrent_requests:
            logger.info("Received multiple empty batches, assuming completion")
            break
            
        if batch_results:
            all_genome_records.extend(batch_results)
            logger.info(f"Retrieved {len(batch_results):,} records in this batch. Total so far: {len(all_genome_records):,}")
        
        # Small sleep between batches
        time.sleep(1)
    
    # Convert to DataFrame
    if all_genome_records:
        df = pd.DataFrame(all_genome_records)
        # Rename columns to match expected format
        df.columns = [f'genome.{col}' for col in df.columns]
        logger.info(f"Total records retrieved: {len(df):,}")
        return df
    return pd.DataFrame()

def build_download_cmd(row: pd.Series, outdir: str, seq_type: SequenceType) -> List[Tuple[List[str], str]]:
    """
    Build command(s) to download genome sequences using ID.
    Skip if output file already exists.
    
    Args:
        row: Pandas Series containing genome information
        outdir: Output directory path
        seq_type: Type of sequence to download (nucleotide, protein, or both)
    
    Returns:
        List of tuples (command list, output path)
    """
    name = row['genome.genome_name'].replace('/','__')
    genome_id = row['genome.genome_id']
    commands = []

    if seq_type in (SequenceType.NUCLEOTIDE, SequenceType.BOTH):
        nuc_path = os.path.join(outdir, f"{name}.fna")
        if not os.path.exists(nuc_path):
            nuc_cmd = ["p3-genome-fasta", genome_id]
            commands.append((nuc_cmd, nuc_path))
        
    if seq_type in (SequenceType.PROTEIN, SequenceType.BOTH):
        prot_path = os.path.join(outdir, f"{name}.faa")
        if not os.path.exists(prot_path):
            prot_cmd = ["p3-genome-fasta", "--protein", genome_id]
            commands.append((prot_cmd, prot_path))
        
    return commands

def make_p3_cmd(filters: Optional[List[str]], search_type: SearchType, attributes: List[str]) -> List[str]:
    """
    Create command for p3-all-genomes with filters and specified search type.
    
    Args:
        filters: Optional list of filter strings in format "field,value"
        search_type: SearchType.EXACT for exact matches or SearchType.SUBSTRING for substring matches
        attributes: List of attributes to retrieve
    
    Returns:
        Command as list of strings
    """
    logger.info(f"Creating p3-all-genomes command with {search_type.name} search")
    if filters:
        logger.info(f"Filters: {filters}")
    logger.info(f"Attributes: {attributes}")
    
    cmd = ['p3-all-genomes']
    for attr in attributes:
        cmd.extend(['--attr', attr])
    
    if filters:
        for filter in filters:
            try:
                field, value = filter.split(',')
                cmd.extend([search_type.value, ','.join((field.strip(), value.strip()))])
            except ValueError:
                logger.error(f"Invalid filter format: {filter}. Expected format: field,value")
                raise ValueError(f"Invalid filter format: {filter}. Expected format: field,value")
    
    return cmd

def run_cmd(cmd: List[str], output_path: Optional[str] = None, max_retries: int = 3) -> Tuple[Optional[bytes], int]:
    """
    Run a command with retries and proper error handling.
    
    Args:
        cmd: Command to run as list of strings
        output_path: Optional path to save output
        max_retries: Maximum number of retry attempts
    
    Returns:
        Tuple of (stdout bytes, return code)
    """
    cmd_str = ' '.join(cmd)
    logger.debug(f"Running command: {cmd_str}")
    
    for attempt in range(max_retries):
        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate(timeout=300)

            if process.returncode != 0:
                logger.warning(f"Attempt {attempt + 1}/{max_retries} failed. Return code: {process.returncode}")
                logger.warning(f"stderr: {stderr.decode('utf-8')}")
                if attempt < max_retries - 1:
                    logger.info("Retrying...")
                continue

            if not stdout:
                logger.warning(f"Attempt {attempt + 1}/{max_retries} returned empty output")
                if attempt < max_retries - 1:
                    logger.info("Retrying...")
                continue

            if output_path and stdout:
                try:
                    with open(output_path, 'wb') as f:
                        f.write(stdout)
                    logger.debug(f"Successfully wrote output to {output_path}")
                except IOError as e:
                    logger.error(f"Failed to write output to {output_path}: {e}")
                    return None, 1

            return stdout, process.returncode

        except subprocess.TimeoutExpired:
            logger.warning(f"Attempt {attempt + 1}/{max_retries} timed out after 300 seconds")
            process.kill()
            if attempt < max_retries - 1:
                logger.info("Retrying...")
            continue
        
        except Exception as e:
            logger.error(f"Unexpected error running command {cmd_str}: {e}")
            return None, 1

    logger.error(f"All {max_retries} attempts failed for command: {cmd_str}")
    return None, 1

def get_genome_info(filters: Optional[List[str]], search_type: SearchType, attributes: List[str], use_api: bool = False, threads: int = 1) -> Optional[pd.DataFrame]:
    """
    Get genome information using either CLI or API approach.
    """
    # If --big flag is set, go straight to API approach
    if use_api and filters and len(filters) == 1 and 'taxon_lineage_names' in filters[0]:
        try:
            taxonomy_id = filters[0].split(',')[1].strip()
            names = fetch_all_genomes_api(taxonomy_id, threads=threads)
            if not names.empty:
                logger.info(f"Successfully retrieved {len(names)} genomes using API")
                return names
        except Exception as e:
            logger.error(f"API approach failed: {str(e)}")
            
    # Try CLI approach if API wasn't used or failed
    try:
        p3_cmd = make_p3_cmd(filters, search_type, attributes)
        logger.info("Attempting to fetch genome information using CLI...")
        stdout, return_code = run_cmd(p3_cmd)
        
        if return_code == 0 and stdout is not None:
            dtypes = {
                'genome.genome_name': str, 
                'genome.genome_id': str,
                'genome.taxon_lineage_ids': str,
                'genome.taxon_lineage_names': str
            }
            data = io.StringIO(stdout.decode())
            names = pd.read_csv(data, sep='\t', dtype=dtypes)
            if not names.empty:
                logger.info(f"Successfully retrieved {len(names)} genomes using CLI")
                return names
    except Exception as e:
        logger.warning(f"CLI approach failed: {str(e)}")

    # If CLI fails or returns no results, try API approach if we haven't already
    if not use_api and filters and len(filters) == 1 and 'taxon_lineage_names' in filters[0]:
        try:
            taxonomy_id = filters[0].split(',')[1].strip()
            names = fetch_all_genomes_api(taxonomy_id, threads=threads)
            if not names.empty:
                logger.info(f"Successfully retrieved {len(names)} genomes using API")
                return names
        except Exception as e:
            logger.error(f"API approach failed: {str(e)}")
    
    return None

def main(args: argparse.Namespace) -> None:
    """
    Main function to coordinate genome downloads.
    
    Args:
        args: Parsed command line arguments
    """
    filters = args.filter if hasattr(args, 'filter') else None
    threads = args.threads
    outdir = args.outdir
    search_type = SearchType.SUBSTRING if args.substring_search else SearchType.EXACT
    seq_type = args.sequence_type
    
    if search_type == SearchType.EXACT and not filters:
        logger.error("Filters must be provided when using exact search")
        sys.exit(1)
    
    attributes = ['genome_name', 'genome_id']
    if args.include_taxonomy:
        attributes.extend(['taxon_lineage_ids', 'taxon_lineage_names'])

    logger.info(f"Starting genome search with {search_type.name} matching")
    logger.info(f"Output directory: {outdir}")
    logger.info(f"Sequence type: {seq_type.name}")

    try:
        os.makedirs(outdir, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directory {outdir}: {e}")
        sys.exit(1)

    # Get genome information using either CLI or API approach
    names = get_genome_info(filters, search_type, attributes, use_api=args.big, threads=threads)
    
    if names is None or names.empty:
        logger.error("Failed to retrieve genome information using both CLI and API approaches")
        sys.exit(1)

    logger.info(f"Found {len(names)} genomes matching the filters")
    
    names['genome.genome_name'] = names['genome.genome_name'].str.replace(r'[^\w\s-]', '').str.replace(' ', '_')

    # Generate all download commands and output paths
    all_downloads = []
    skipped_files = 0
    for _, row in names.iterrows():
        commands = build_download_cmd(row, outdir, seq_type)
        all_downloads.extend(commands)
        # If we expected files but got no commands, they must have been skipped
        expected_files = 2 if seq_type == SequenceType.BOTH else 1
        skipped_files += expected_files - len(commands)

    logger.info(f"Found {skipped_files} existing files that will be skipped")
    logger.info(f"Proceeding with {len(all_downloads)} downloads")

    # Multithreaded download with progress bar
    failed_downloads = []
    successful_downloads = 0
    
    with ThreadPoolExecutor(max_workers=100) as executor:
        futures = {executor.submit(run_cmd, cmd, output_path): (cmd, output_path) 
                  for cmd, output_path in all_downloads}
        
        with tqdm(total=len(futures), desc="Downloading sequences") as pbar:
            for future in concurrent.futures.as_completed(futures):
                cmd, output_path = futures[future]
                try:
                    stdout, return_code = future.result()
                    if return_code != 0 or stdout is None:
                        failed_downloads.append((cmd, output_path))
                    else:
                        successful_downloads += 1
                except Exception as e:
                    logger.error(f"Error downloading {output_path}: {e}")
                    failed_downloads.append((cmd, output_path))
                pbar.update(1)

    logger.info(f"Download complete. Successfully downloaded {successful_downloads} sequences")
    if failed_downloads:
        logger.warning(f"Failed to download {len(failed_downloads)} sequences")
        with open('failed_downloads.txt', 'w') as f:
            for cmd, output_path in failed_downloads:
                f.write(f"Command: {' '.join(cmd)}\nOutput path: {output_path}\n\n")
        logger.warning("Failed downloads have been logged to failed_downloads.txt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download genomes from BV-BRC database.')
    parser.add_argument('--filter', nargs='+', required=False, 
                      help='Filter criteria in the form field_name,value. Example: --filter taxon_lineage_names,Alphacoronavirus')
    parser.add_argument('--threads', type=int, default=1,
                      help='The number of threads to be used in multithreading.')
    parser.add_argument('--outdir', required=True,
                      help='The directory where results will be stored.')
    parser.add_argument('--substring_search', action='store_true',
                      help='Use substring search (--in) instead of exact match (--eq)')
    parser.add_argument('--include_taxonomy', action='store_true',
                      help='Include taxonomic lineage information in output')
    parser.add_argument('--sequence_type', type=str, choices=['nucleotide', 'protein', 'both'],
                      default='nucleotide', help='Type of sequence to download')
    parser.add_argument('--debug', action='store_true',
                      help='Enable debug logging'),
    parser.add_argument('--big', action='store_true',
                  help='Use API approach directly for large taxonomy groups')


    args = parser.parse_args()
    
    # Convert sequence_type string to enum
    args.sequence_type = SequenceType[args.sequence_type.upper()]
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
        
    main(args)
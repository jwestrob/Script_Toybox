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
from typing import Tuple, List, Optional, Dict
from enum import Enum, auto

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

def build_download_cmd(row: pd.Series, outdir: str, seq_type: SequenceType) -> List[Tuple[List[str], str]]:
    """
    Build command(s) to download genome sequences using ID.
    
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
        nuc_cmd = ["p3-genome-fasta", genome_id]
        commands.append((nuc_cmd, nuc_path))
        
    if seq_type in (SequenceType.PROTEIN, SequenceType.BOTH):
        prot_path = os.path.join(outdir, f"{name}.faa")
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

    try:
        p3_cmd = make_p3_cmd(filters, search_type, attributes)
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)

    logger.info("Fetching genome information...")
    stdout, return_code = run_cmd(p3_cmd)
    
    if return_code != 0 or stdout is None:
        logger.error("Failed to fetch genome information")
        sys.exit(1)

    try:
        dtypes = {
            'genome.genome_name': str, 
            'genome.genome_id': str,
            'genome.taxon_lineage_ids': str,
            'genome.taxon_lineage_names': str
        }
        data = io.StringIO(stdout.decode())
        names = pd.read_csv(data, sep='\t', dtype=dtypes)
        
        if names.empty:
            logger.error("No genomes found matching the specified filters")
            sys.exit(1)
            
        logger.info(f"Found {len(names)} genomes matching the filters")
            
    except Exception as e:
        logger.error(f"Failed to parse genome information: {e}")
        sys.exit(1)

    names['genome.genome_name'] = names['genome.genome_name'].str.replace(r'[^\w\s-]', '').str.replace(' ', '_')

    # Generate all download commands and output paths
    all_downloads = []
    for _, row in names.iterrows():
        all_downloads.extend(build_download_cmd(row, outdir, seq_type))

    # Multithreaded download with progress bar
    failed_downloads = []
    successful_downloads = 0
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
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
                      help='Enable debug logging')

    args = parser.parse_args()
    
    # Convert sequence_type string to enum
    args.sequence_type = SequenceType[args.sequence_type.upper()]
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
        
    main(args)

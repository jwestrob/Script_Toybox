
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Create a logger
logging.basicConfig(filename='hmm_search.log', level=logging.INFO)

# Parse command line arguments
parser = argparse.ArgumentParser(description="Scan protein FASTA files with a HMM and extract regions with hits.")
parser.add_argument('--input', type=str, required=True, help='Input protein FASTA file or directory.')
parser.add_argument('--output_filename', type=str, required=True, help='Output file name.')
parser.add_argument('--hmm_file', type=str, required=True, help='Path to HMM file or directory.')
parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
parser.add_argument('--len_label', action='store_true', help='Append the length of the original ORF to the FASTA headers.')
args = parser.parse_args()

# File paths
hmm_file_path = args.hmm_file
input_path = args.input

# Check if the HMM file path is a directory or a file
if os.path.isdir(hmm_file_path):
    try:
        hmm_files = [os.path.join(hmm_file_path, f) for f in os.listdir(hmm_file_path) if f.endswith('.hmm')]
    except FileNotFoundError:
        logging.error(f'Could not find the directory: {hmm_file_path}')
        raise
elif os.path.isfile(hmm_file_path) and hmm_file_path.endswith('.hmm'):
    hmm_files = [hmm_file_path]
else:
    logging.error(f'Invalid input: {hmm_file_path}. Please provide a valid HMM file or directory containing HMM files.')
    raise ValueError(f'Invalid input: {hmm_file_path}. Please provide a valid HMM file or directory containing HMM files.')

# Check if the input is a directory or a file
if os.path.isdir(input_path):
    try:
        input_files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith('.faa')]
    except FileNotFoundError:
        logging.error(f'Could not find the directory: {input_path}')
        raise
elif os.path.isfile(input_path) and input_path.endswith('.faa'):
    input_files = [input_path]
else:
    logging.error(f'Invalid input: {input_path}. Please provide a valid FASTA file or directory containing FASTA files.')
    raise ValueError(f'Invalid input: {input_path}. Please provide a valid FASTA file or directory containing FASTA files.')

# Create a list of tasks
tasks = [(hmm_file, fasta_file) for hmm_file in hmm_files for fasta_file in input_files]

def process_task(task):
    hmm_file, fasta_file = task

    logging.info(f'Processing HMM file: {hmm_file}')
    logging.info(f'Processing FASTA file: {fasta_file}')

    # Construct full paths for input files
    fasta_path = fasta_file
    hmmsearch_output = 'hmmsearch_output_' + os.path.basename(fasta_file) + '_' + os.path.basename(hmm_file) + '.txt'

    # Run HMMsearch with STDOUT directed to /dev/null
    with open(os.devnull, 'w') as devnull:
        subprocess.run(["hmmsearch", "--domtblout", hmmsearch_output, hmm_file, fasta_path], stdout=devnull)

    # Extract hit regions from HMMsearch output
    hit_regions = {}
    with open(hmmsearch_output, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # Ignore comment lines
                fields = line.split()
                seq_id = fields[0]
                e_value = float(fields[12])
                if e_value < 1e-5:  # Check e-value
                    start = int(fields[19]) - 1  # Python is 0-indexed
                    end = int(fields[20])  # End is inclusive
                    if seq_id not in hit_regions or (start, end) not in hit_regions[seq_id]:
                        if seq_id not in hit_regions:
                            hit_regions[seq_id] = []
                        hit_regions[seq_id].append((start, end))

    # Filter original FASTA file for hit regions
    filtered_records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in hit_regions:
            for start, end in hit_regions[record.id]:
                hit_seq = record.seq[start:end]
                if args.len_label:
                    new_record = SeqRecord(hit_seq, id=record.id + "_" + str(start) + "_" + str(end) + "__len" + str(len(record.seq)), description="")
                else:
                    new_record = SeqRecord(hit_seq, id=record.id + "_" + str(start) + "_" + str(end), description="")
                filtered_records.append(new_record)

    # Remove the HMMsearch output file
    os.remove(hmmsearch_output)

    return filtered_records

# Use a ProcessPoolExecutor to process the tasks
with ProcessPoolExecutor(max_workers=args.threads) as executor:
    futures = [executor.submit(process_task, task) for task in tasks]
    
    all_filtered_records = []
    for future in tqdm(as_completed(futures), total=len(futures), desc='Processing tasks'):
        all_filtered_records.extend(future.result())

# Write all filtered sequences to a single output FASTA file
try:
    SeqIO.write(all_filtered_records, args.output_filename, "fasta")
except Exception as e:
    logging.error(f'Error writing output file: {e}')
    raise

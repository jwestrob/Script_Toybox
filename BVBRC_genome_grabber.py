import os, sys, pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import argparse
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import subprocess
import io

def build_download_cmd(row, outdir):
    # builds command to download genome fasta using ID
    # saves to proper name instead of genome ID
    name = row['genome.genome_name']
    genome_id = row['genome.genome_id']
    output_path = os.path.join(outdir, f"{name}.fna")
    cmd = ["p3-genome-fasta", genome_id]
    return cmd, output_path



def make_p3_cmd(filters):
    cmd = ['p3-all-genomes', '--attr', 'genome_name']
    for filter in filters:
        field, value = filter.split(',')
        cmd.extend(['--eq', ','.join((field.strip(), value.strip()))])
    return cmd
    
def run_cmd(cmd, output_path=None, max_retries=3):
    for _ in range(max_retries):
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # If an output path was specified in the command, write stdout to it
        if output_path and stdout:
            with open(output_path, 'wb') as f:
                f.write(stdout)

        # retry if process fails or if download is empty
        if process.returncode != 0 or not stdout:
            print(f"stderr: {stderr.decode('utf-8')}")
            print(f"Download failed. Retrying...")
        else:
            return stdout, process.returncode

    # if we've reached here, all retries have failed
    print(f"All retries failed for command: {' '.join(cmd)}")
    return None, 1


def main():
    filters = args.filter
    threads = args.threads
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    # Generate command to obtain BVBRC IDs with filters
    p3_cmd = make_p3_cmd(filters)

    # Execute command to download info
    stdout, _ = run_cmd(p3_cmd)

    dtypes = {'genome.genome_name': str, 'genome.genome_id': str}
    data = io.StringIO(stdout.decode())
    names = pd.read_csv(data, sep='\t', dtype=dtypes)

    # Replace spaces with underscores in genome names
    names['genome.genome_name'] = names['genome.genome_name'].apply(lambda x: '_'.join(x.split(' ')))

    # Generate download commands and output paths
    download_cmds_and_paths = names.apply(build_download_cmd, outdir=outdir, axis=1)

    # Multithreaded download
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_cmd, cmd, output_path) for cmd, output_path in download_cmds_and_paths}
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download genomes from BV-BRC database.')
    parser.add_argument('--filter', nargs='+', required=True, help='Filter criteria in the form field_name,value. Example: --filter host_common_name,Human genome_status,Complete')
    parser.add_argument('--threads', type=int, default=1, help='The number of threads to be used in multithreading.')
    parser.add_argument('--outdir', required=True, help='The directory where results will be stored.')

    args = parser.parse_args()
    main()
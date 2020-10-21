import os, sys, pandas as pd
from Bio import Entrez
from Bio import SeqIO
import urllib.request
import argparse
from functools import reduce
from tqdm import tqdm


#Make a folder for the genomes and for the metadata

#Make a log and a record of biosample IDs

def check_biosample_summary(biosample, lines, outdir):
    try:
        biosample_line = list(filter(lambda x: 'BioSample' in x, lines))[0]
    except:
        print("NCBI Entrez returned no information for BioSample ID " + biosample)

    biosample_ID = biosample_line.split('BioSample: ')[1].strip(';')
    if biosample_ID != biosample:
        print("NCBI Entrez returned bad information for BioSample ID " + biosample)
    return

def downloadingAssemblySummary(assembly_summary_filename) :
    url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
    if not os.path.exists(assembly_summary_filename):
        print('downloading assembly_summary_genbank.txt from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS')
        liste = urllib.request.urlretrieve(url, assembly_summary_filename)
        #print(liste)

def entrez_grab_ids(tax_searchterm):
    Entrez.email = "jwestrob@berkeley.edu"
    handle = Entrez.esearch(db="assembly", term=tax_searchterm, retmax='5000')
    record = Entrez.read(handle)
    ids = record['IdList']
    return ids

def get_accession_and_biosample_id(id):
    """Get Accession number and Biosample ID for assembly db ID"""
    Entrez.email = "jwestrob@berkeley.edu"
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'], esummary_record['DocumentSummarySet']['DocumentSummary'][0]['BioSampleId']


def download_genome_by_accession(accession, outdir, ftp_path):
    genome_id_waccession = ftp_path.split('/')[-1]
    output_genome_filename = os.path.join(outdir, genome_id_waccession + '.fna.gz')
    genome_filename = genome_id_waccession + '_genomic.fna.gz'
    if os.path.exists(output_genome_filename):
        return
    try:
        list = urllib.request.urlretrieve(ftp_path + '/' + genome_filename, output_genome_filename)
    except (HTTPError, URLError) as error:
        print("Data of " + accession + " not retrieved because " + error + " at " + ftp_path)
    return list


def grab_biosample_info(accession, biosample, outdir):
    Entrez.email = "jwestrob@berkeley.edu"
    if os.path.exists(os.path.join(outdir, accession + '_BS_' + biosample + '_info.txt')):
        return
    handle = Entrez.efetch(db='BioSample', id=biosample, retmode='text')
    lines = [x.rstrip() for x in handle.readlines()]
    if not check_biosample_summary(biosample, lines):
        bad_ids = open('bad_biosample_requests.tsv','a')
        bad_ids.write(accession + '\t' + biosample + '\n')
        return
    with open(os.path.join(outdir, accession + '_BS_' + biosample + '_info.txt'), 'w') as outfile:
        for element in lines:
            outfile.writelines(element + '\n')
    return lines

def main():
    ap = argparse.ArgumentParser(description='Fetches sequence metadata from NCBI, genomes, and associated metadata.')

    ap.add_argument('-taxfile',
                  type=str,
                  required=False,
                  help="If you'd like to download all the hits for a specific taxonomy; will download the IDs with Biopython Entrez")

    ap.add_argument('-assembly_summary',
                  type=str,
                  required=False,
                  help="Summary file from NCBI with genbank genomes info; point me to it if you already downloaded it. Otherwise it will be downloaded and named assembly_summary_genbank_[DATE].txt")

    ap.add_argument('-accession_list',
                  type=str,
                  required=False,
                  help="File with GCA genbank accessions on each line please.")

    ap.add_argument('-threads',
                  type=int,
                  required=False,
                  default=1,
                  help="Number of threads to use.")

    ap.add_argument('-outdir',
                  type=str,
                  required=True,
                  help="Name of directory to output files. Required.")

    ap.add_argument('-no_genomes',
                  default=False,
                  action='store_true',
                  help="Don't download genomes, only biosample metadata.")

    ap.add_argument('-no_metadata',
                  default=False,
                  action='store_true',
                  help="Don't download BioSample metadata.")

    args = ap.parse_args()

    threads = args.threads

    outdir = args.outdir

    no_genomes = args.no_genomes
    no_metadata = args.no_metadata

    if os.path.exists(outdir):
        print("Warning: output directory already exists. This run may overwrite files in that directory.")
        print("Please exit now if you do not wish this to occur.")
    else:
        os.system('mkdir ' + outdir)

    assemblyfile_names =  ['assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name',
       'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
       'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
       'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
       'excluded_from_refseq', 'relation_to_type_material']



    if args.assembly_summary == None:
        if not os.path.exists('assembly_summary_genbank.txt'):
            print("Downloading assembly summary:")
            downloadingAssemblySummary('assembly_summary_genbank.txt')
        print("Loading assembly summary...")
        assembly_summary = pd.read_csv('assembly_summary_genbank.txt', sep='\t', skiprows=2, names=assemblyfile_names, low_memory=False)
    else:
        print("Loading assembly summary...")
        assembly_summary = pd.read_csv(args.assembly_summary, sep='\t', skiprows=2, names=assemblyfile_names, low_memory=False)

    accession2biosample = dict(zip(assembly_summary.assembly_accession.tolist(), assembly_summary.biosample.tolist()))
    accession2ftp = dict(zip(assembly_summary.assembly_accession.tolist(), assembly_summary.ftp_path.tolist()))

    if args.accession_list == None:
        if args.taxfile == None:
            print("You need to provide something to search for or download! Specify -taxfile or -accession_list, pleeease.")
            sys.exit(420)
    else:
        with open(args.accession_list, 'r') as infile:
            accession_ids = [x.rstrip() for x in infile.readlines()]

        biosample_ids = [accession2biosample[x] for x in accession_ids]

    if args.taxfile != None:
        with open(args.taxfile, 'r') as infile:
            search_terms = [x.rstrip() for x in infile.readlines()]

        all_ids = set()
        print("Grabbing IDs for search terms...")
        all_ids_lists = []
        for x in tqdm(search_terms):
            all_ids_lists.append(entrez_grab_ids(x))
        #all_ids_lists = [entrez_grab_ids(x) for x in search_terms]


        flatten = lambda l: [item for sublist in l for item in sublist]

        all_ids_list = flatten(all_ids_lists)
        print(str(len(all_ids_list)), "IDs grabbed!")

        accessions_tax = []
        biosample_ids_tax = []
        print("Retrieving biosample IDs...")
        for id in tqdm(all_ids_list):
            accession, biosample_id = get_accession_and_biosample_id(id)
            accessions_tax.append(accession)
            biosample_ids_tax.append(biosample_id)


    accessions = pd.Series(accession_ids + accessions_tax).unique().tolist()

    #Function to find the row in the assembly summary corresponding to an accession
    #if NCBI fucked up the accession in some subtle way so you can't do direct
    #lookup because they're a bunch of... well you know.
    def find_equivalent(messed_up_accession):
        #Thanks Eugene
        true_format_accessions = accession2biosample.keys()
        if messed_up_accession in true_format_accessions:
            return messed_up_accession
        id_without_GCA_and_decimal = messed_up_accession.split('_')[1].split('.')[0]
        matching_id = list(filter(lambda x: id_without_GCA_and_decimal in x, true_format_accessions))[0]
        return matching_id

    accessions_fixed = list(map(find_equivalent, accessions))

    biosamples = [accession2biosample[find_equivalent(x)] for x in accessions]

    with open(os.path.join(outdir, 'biosample_list.txt'), 'w') as outfile:
        for index in range(len(biosamples)):
            outfile.writelines(accessions_fixed[index] + '\t' + biosamples[index] + '\n')

    downloaded_genomes = []
    if not no_genomes:
        print("Downloading genomes:")
        for accession in tqdm(accessions_fixed):
            if accession not in accession2ftp.keys():
                #this means NCBI is trying to  pull some shit... find the right accession ID
                try:
                    genomes_outdir = os.path.join(outdir, 'genomes')
                    if not os.path.exists(genomes_outdir):
                        os.system('mkdir ' + genomes_outdir)
                    download_genome_by_accession(accession, outdir, accession2ftp[accession])
                except:
                    print("Bad accession! Please investigate: ", accession)
            else:
                download_genome_by_accession(accession, outdir, accession2ftp[accession])
    #downloaded_genomes = list(map(lambda accession: download_genome_by_accession(accession, outdir, accession2ftp[accession]), accessions))
    """
    if not no_metadata:
        biosample_list = []
        print("Downloading biosample information:")
        for accession in tqdm(accessions_fixed):
            if accession not in accession2ftp.keys():
                try:

                    biosample_list.append(accession2ftp[accession])
                    #grab_biosample_info(accession, accession2ftp[accession], outdir)

                except:
                    print("Bad accession! Please investigate: ", accession)
            else:
                biosample_list.append(accession2biosample[accession])
                #grab_biosample_info(accession, accession2biosample[accession], outdir)
    #biosample_info = list(map(lambda accession: grab_biosample_info(accession, accession2biosample[accession], outdir),  accessions))
        with open(os.path.join(outdir, 'biosample_list.txt'), 'w') as outfile:
            for element in biosample_list:
                outfile.writelines(element + '\n')
    """

    print("Complete!")









if __name__ == '__main__':
    main()

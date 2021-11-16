import os, sys, pandas
from Bio import SearchIO, SeqIO
import argparse, subprocess, glob

def filter_contigs_notshit(assembled_contigs, outdir):
    ''' Calls pullseq to filter assembled contigs, doesn't shit the bed while doing so '''

    filtered_fasta = os.path.join(outdir, assembled_contigs.split('/')[-1].split('.')[0] + '_min5000.fa')
    print("Filtered fasta: ", filtered_fasta)
    if os.path.exists(filtered_fasta):
        print("You already filtered your fasta! Not doing this again.")
        return filtered_fasta
    pullseq_command = 'pullseq -m 5000 -i ' + assembled_contigs + ' > ' + filtered_fasta
    print("Filtering contigs by length with pullseq...")
    print(pullseq_command)
    send = subprocess.Popen(pullseq_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print(stderr)
        print("Pullseq step failed!")
        sys.exit()
    print("Contigs filtered!")
    return os.path.abspath(filtered_fasta)

def run_prodigal(filtered_fasta, outdir):
    ''' Calls meta-prodigal to predict proteins on min5000 contigs '''

    prodigal_out = filtered_fasta.split('.fa')[0] + '_genes.faa'
    if os.path.exists(prodigal_out):
        print("You already did prodigal! Not repeating this...")
        return prodigal_out

    prodigal_command = "prodigal -p meta -i " + filtered_fasta + " -a " + prodigal_out
    print(prodigal_command)
    send = subprocess.Popen(prodigal_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print(stderr)
        print("Prodigal gene prediction step failed!")
        sys.exit()
    return prodigal_out

def pfam_hmmsearch(protfile, threads, outdir):
    ''' Just searches with the PFAM hmm on provided proteins. '''
    pfam_db = '/groups/banfield/users/jwestrob/Viral_Identification/earth_virome_pipeline/reference_files/Pfam-A.hmm'
    pfam_out = os.path.join(outdir, 'pfam.hmmout')
    if os.path.exists(pfam_out):
        print("You already did your pfam search! Not repeating.")
        return pfam_out
    #num_cpu = multiprocessing.cpu_count()
    hmmsearch_command = "hmmsearch --cpu " + str(threads) + " --cut_tc --domtblout " + pfam_out + ' ' + pfam_db + ' ' + protfile
    #final_command = "prodigal -i {0} -a {1} && hmmsearch --cpu {2} --cut_tc --domtblout {3} {4} {1}".format(filtered_contigs, prodigal_ORF_out, threads, pfam_out, pfam_db)

    print('Annotating proteins with pfam...')
    print(hmmsearch_command)
    send = subprocess.Popen(hmmsearch_command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print("Pfam annotation failed!")
        print(stderr)
        sys.exit()

    return pfam_out

def format_prodigal_out_headers(protfile, outdir):
    '''Takes a prodigal predicted protein fasta file and renames the headers such that the format is ">CONTIG|GENE"
    where CONTIG is the contig ID and GENE is the gene ID for a gene in the CONTIG'''
    reformatted_file = os.path.join(outdir, protfile.replace('.fa', '.reformatted.fa'))

    print('Refomatting protein fasta headers...')

    new_recs = []
    for rec in SeqIO.parse(protfile, 'fasta'):
        rec_id = rec.id
        gene = rec_id
        scaffold = '_'.join(gene.split('_')[:-1])
        new_rec = rec
        new_rec.id = scaffold + '|' + gene
        new_recs.append(new_rec)
    SeqIO.write(new_recs, reformatted_file, 'fasta')

    print('Done! Formatted prodigal output at: {0}\n'.format(reformatted_file))
    return reformatted_file

def parse_pfam_out(pfam_out, outdir):
    '''Takes in pfam_out file in domtblout format and returns a custom trimmed version as outlined in the Earth Virome pipeline'''
    df = pandas.read_csv(pfam_out, skiprows=[0,1,2], sep = '\s+', header = None).dropna()
    important_col = [0, 4, 21, 15, 16, 17, 18, 6, 7, 19, 20]
    out = df.iloc[:, important_col].copy()
    out.columns = ['GeneID', 'pfamID','pident','query_start','query_end','subject_start','subject_end','e-value','bit_score','align_start','align_end']
    out['Alignment_Length'] = out['align_end'] - out['align_start'] + 1
    out['pident'] = (out['pident']*100).apply(int)
    out = out.drop(columns=['align_start','align_end'])

    pfam_final = os.path.join(outdir, pfam_out.split('.hmmout')[0] + '.txt')
    out.to_csv(pfam_final, index = None, sep = '\t')
    print('Done! Parsed Pfam out at: {0}\n'.format(pfam_final))
    return pfam_final

def viral_hmmsearch(prodigal_formatted, threads, outdir):
    hmm_out = os.path.join(outdir, prodigal_formatted.split('.fa')[0] + '_hits_to_vHMMs.hmmout')
    if os.path.exists(hmm_out):
        print("You already performed viral HMMsearch, not repeating...")
        return hmm_out
    viral_hmm = '/groups/banfield/users/jwestrob/Viral_Identification/earth_virome_pipeline/reference_files/viral_reference_model.hmm'

    final_command = 'hmmsearch --cpu {0} -E 1.0e-05 --tblout {1} {2} {3}'.format(str(threads), hmm_out, viral_hmm, prodigal_formatted)

    print('Identifying viral proteins using hmmsearch...')
    send = subprocess.Popen(final_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print("Viral hmmsearch failed!")
        print(stderr)
        sys.exit()

    print('Done! Hits to viral protein can be found at: {0}\n'.format(hmm_out))
    return hmm_out

def get_hits_to_VPFs(hmmout_file):
    '''Takes a HMMER3 hmmsearch tab output file as an input and
    returns a dictionary mapping each scaffold with the number of unique genes that match a protein family

    Input:
        - hmmout_file (str): path to HMMER3 hmmsearch out file in tab format

    Returns:
        - hits_to_VPFs (dict): dictionary where key are scaffold IDs and values are number of unique genes that matched a protein family
    '''
    hits_to_VPFs = {}
    with open(hmmout_file, 'r') as input:
        for qresult in SearchIO.parse(input, 'hmmer3-tab'):
            hits = qresult.hits
            num_hits = len(hits)
            if num_hits > 0:
                for i in range(0,num_hits):
                    query_seq_id = hits[i].id
                    scaffold, gene = query_seq_id.split('|')
                    hits_to_VPFs[scaffold] = hits_to_VPFs.get(scaffold, set([])).union([gene])

    for key, value in iter(hits_to_VPFs.items()):
        hits_to_VPFs[key] = len(value)
    return hits_to_VPFs

def get_num_of_genes(fasta_file):
    '''Takes a FASTA file as an input where each FASTA header is defined as: Scaffold_ID|Gene_ID
    and returns a dictionary mapping each scaffold with the number of unique genes in that scaffold

    Input:
        - fasta_file (str): path to fasta file in which headers are formatted Scaffold_ID|Gene_ID

    Returns:
        - gene_count (dict): dictionary where key are scaffold IDs and values are number of unique genes that are part of the scaffold
        - gene_scaffold_map (dict): dictionary that maps each gene ID back to its respective scaffold ID
    '''

    gene_count = {}
    gene_scaffold_map = {}
    with open(fasta_file) as fasta:
        for line in fasta:
            if line[0] == '>':
                scaffold, gene = line[1:].strip().split('|')
                gene_count[scaffold] = gene_count.get(scaffold, set([])).union([gene])
                gene_scaffold_map[gene] = scaffold
    for key, value in iter(gene_count.items()):
        gene_count[key] = len(value)
    return gene_count, gene_scaffold_map

def blastn(viral_contigs, threads, outdir):

    reference_database = '/groups/banfield/users/jwestrob/Viral_Identification/earth_virome_pipeline/reference_files/mVCs_PaezEspino_Nature.fna'

    blastn_out = os.path.join(outdir, 'blastn_mVCs.blout')
    #num_cpu = threads
    blastn_command = "blastn -query {0} -db {1} -outfmt '6 std qlen slen' -out {2} -evalue 1.0e-50 -perc_identity 80 -num_threads {3}".format(viral_contigs, reference_database, blastn_out, threads)

    print('Performing blastn on viral contigs')
    print(blastn_command)
    send = subprocess.Popen(blastn_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print("Blastn failed!")
        print(stderr)
        sys.exit()
    print('blastn complete, output can be found at: {0}'.format(blastn_out))
    return blastn_out

def write_master_table(master_table_file, VPF_hits, num_gene, pfam_gene_count):
    with open(master_table_file, 'w') as master:
        master.write('\t'.join(['Scaffold_ID', 'hits_to_VPFs', '#_of_genes', '%covered_VPFs', 'genes_with_pfams', '%genesPfams']) + '\n')
        for scaffold in VPF_hits.keys():
            out = []
            out.append(scaffold)
            out.append(VPF_hits.get(scaffold, 0))
            out.append(num_gene.get(scaffold, 0))
            out.append(VPF_hits.get(scaffold, 0) * 100.0 / float(num_gene.get(scaffold, 0)))
            out.append(pfam_gene_count.get(scaffold, 0))
            out.append(pfam_gene_count.get(scaffold, 0) * 100.0 / float(num_gene.get(scaffold, 0)))
            master.write('\t'.join([str(col) for col in out]) + '\n')
    return master_table_file


def parse_cluster_blast(blastn_out, outdir):
    # Remove self hits
    self_hits_out = out_folder + blastn_out.split('/')[-1].split('.')[0] + 'noSelfHits.blout'
    rm_self_hits = "cat {0} | awk '$1 != $2' > {1}".format(blastn_out, self_hits_out)

    print("Removing self hits...")
    send = subprocess.Popen(rm_self_hits, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    if stderr is not None:
        print("Something's funky with the blast output, and we've run into a snag. Please investigate.")
        print(stderr)
        sys.exit()


    # Parse unique BLAST results:
    parsed_out = os.path.join(outdir, 'blast_parsed.blout')
    parse_blast = "java Parse_BLAST {0} > {1}".format(self_hits_out, parsed_out)

    # Single linkage clustering
    slc_out = out_folder + parsed_out.split('/')[-1].split('.')[0] + '.slc'
    slc = "perl SLC.pl {0} {1}".format(parsed_out, slc_out)

    merged_command = ' && '.join([rm_self_hits, parse_blast, slc])

    print('Removing self hits, parsing blastn output and performing single linkage clustering')
    print(merged_command)
    send = subprocess.Popen(merged_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print(stderr)

    print('Operations complete! Output for single linkage clustering can be found at: {0}'.format(slc_out))
    return slc

def get_pfam_genes(pfam_file, gene_scaffold_map):
    '''Takes a pfam file as an input and returns a dictionary mapping each scaffold with the number o$
    
    Input:
        - pfam_file (str): path to pfam out file
    
    Returns:
        - scaffold_pfam_count (dict): dictionary where key are scaffold IDs and values are number of $
    '''
    
    pfam = pandas.read_table(pfam_file)
    scaffold_pfam_count = {}
    for gene in pfam['GeneID'].values:
        if gene_scaffold_map.get(gene, None) != None:
            scaffold_pfam_count[gene_scaffold_map[gene]] = scaffold_pfam_count.get(gene_scaffold_map[gene], 0) + 1
    return scaffold_pfam_count


def filter_and_extract_sequences(out_folder, master_table, gene_fasta_file, assembly_fasta_file):
    '''Takes a master table as an input and filters sequences, keeping only those that are inferred as viral'''
    filter1 = os.path.join(out_folder,'Filter1.out')
    filter2 = os.path.join(out_folder,'Filter2.out')
    filter3 = os.path.join(out_folder,'Filter3.out')
    out_filtered_contigs = os.path.join(out_folder,gene_fasta_file.split('/')[-1].split('genes')[0] + 'filtered_viral_contigs.txt')
    out_filtered_contigs_fa = os.path.join(out_folder,gene_fasta_file.split('/')[-1].split('genes')[0] + 'filtered_viral_contigs.fa')
    
    filter_command_1 = "cat {0} | awk '$2 >= 5' | awk '$6 <= 40' | awk '$4 >= 10' | cut -f 1 > {1}".format(master_table, filter1)
    print(filter_command_1)
    send = subprocess.Popen(filter_command_1, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print("Stdout, stderr: ", stdout, stderr)

    filter_command_2 = "cat {0} | awk '$2 >= 5' | awk '$2 >= $5' | cut -f 1 > {1}".format(master_table, filter2)
    print(filter_command_2)
    send = subprocess.Popen(filter_command_2, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print("Stdout, stderr: ", stdout, stderr)

    filter_command_3 = "cat {0} | awk '$2 >= 5' | awk '$4 >= 60' | cut -f 1 > {1}".format(master_table, filter3)
    print(filter_command_3)
    send = subprocess.Popen(filter_command_3, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print("Stdout, stderr: ", stdout, stderr)

    merge_command = "cat {0} {1} {2} | sort | uniq > {3}".format(filter1, filter2, filter3, out_filtered_contigs)
    print(merge_command)
    send = subprocess.Popen(merge_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print("Stdout, stderr: ", stdout, stderr)

    extract_command = "seqtk subseq {0} {1} > {2}".format(assembly_fasta_file, out_filtered_contigs, out_filtered_contigs_fa)
    
    #final_command = " && ".join([filter_command_1, filter_command_2, filter_command_3, merge_command, extract_command])
    print(extract_command)
    send = subprocess.Popen(extract_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (stdout, stderr) = send.communicate()
    print(stdout)
    print(stderr)
    
    return out_filtered_contigs_fa


def main():
    parser = argparse.ArgumentParser(description='Wanna find some viruses? Let us find some viruses.')
    parser.add_argument('-contigs',
                        help='path to fasta file containing the contigs from assembly')
    parser.add_argument('-outdir', required=True,
                        help = 'path to output folder [REQUIRED]')
    parser.add_argument('-threads', default = 6, help="Number of threads to use")
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    print("Outdir: ", outdir)
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
    contigs = os.path.abspath(args.contigs)
    if not os.path.exists(contigs):
        print("Please provide a valid path to your contigs file. You have failed in this.")
        sys.exit()

    threads = args.threads
    #Functions from annotate_assembled_contigs.py
    filtered_contigs = filter_contigs_notshit(contigs, outdir)
    print("Filtered contigs: ", filtered_contigs)
    protfile = run_prodigal(filtered_contigs, outdir)
    print("Protfile: ", protfile)
    pfam_out = pfam_hmmsearch(protfile, threads, outdir)
    reformatted_proteins = format_prodigal_out_headers(protfile, outdir)
    pfam_final = parse_pfam_out(pfam_out, outdir)
    viral_hmm_out = viral_hmmsearch(reformatted_proteins, threads, outdir)


    #Functions from filter_viral_contigs_master_table.py

    master_table_file = os.path.join(outdir, 'Viral_master_table.txt')

    VPF_hits = get_hits_to_VPFs(viral_hmm_out)
    num_gene, gene2scaffold = get_num_of_genes(reformatted_proteins)
    pfam_gene_count = get_pfam_genes(pfam_final, gene2scaffold)

    write_master_table(master_table_file, VPF_hits, num_gene, pfam_gene_count)
    print('Final master table can be found at: {0}'.format(master_table_file))

    viral_contigs = filter_and_extract_sequences(outdir, master_table_file, reformatted_proteins, filtered_contigs)
    print('Final list of viral contigs can be found at: {0}'.format(viral_contigs))


    #Functions from viral_contig_clustering.py
    cluster_dir = os.path.join(outdir, 'contig_clustering')
    if not os.path.exists(cluster_dir):
        os.system('mkdir ' + cluster_dir)


    #Clustering not yet implemented because it calls an external java thing and it's a piece of shit!!!
    #blast_out = blastn(viral_contigs, viral_contigs_db, threads, cluster_dir)
    #parse_cluster_blast(blast_out, cluster_dir)

    print("Complete!!")

if __name__ == "__main__":
    main()

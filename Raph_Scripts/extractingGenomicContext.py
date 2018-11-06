#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
import argparse

def kNearestNeighbors(orfName,geneList,k) :
    genomicContextList = list()
    sortedList = sorted(geneList,key=lambda x:x[1])
    for i in range( len(sortedList) ) :
        if sortedList[i][3] != orfName :
            continue
        else:
            if i-k < 0 :
                start = 0
            else :
                start = i-k

            if i+k+1 > len(sortedList) :
                end = len(sortedList)
            else :
                end = i+k+1

            for j in range(start,end) :
                centroid = sortedList[i][3]
                if j == i :
                    genomicContextList.append(sortedList[j])
#                    print('centroid\t'+str(sortedList[i]) )
                else :
#                    print('neighbor\t'+str(sortedList[j]) )
                    neighbor = sortedList[j][2]
                    genomicContextList.append(sortedList[j])

        return genomicContextList



def genome2orfOrder(feature_filename) :
    orf2genome = dict()
    orf2scaffold = dict()
    genome2scaffold2orfList = dict()
    file = open(feature_filename,'r')
    for line in file :
        if line.startswith('patric_id'):
            continue
        line = line.rstrip()
        orfName,genome,scaffold,start,end,strand = line.split("\t")
        orfName = orfName.rstrip()
        genome = genome.rstrip()
        orf2genome[orfName] = genome
        orf2scaffold[orfName] = scaffold
        if genome not in genome2scaffold2orfList :
            genome2scaffold2orfList[genome] = dict()

        if scaffold not in genome2scaffold2orfList[genome] :
            genome2scaffold2orfList[genome][scaffold] = list()

        genome2scaffold2orfList[genome][scaffold].append([ int(start) , int(end) , strand , orfName ])
    file.close()

    return genome2scaffold2orfList,orf2genome,orf2scaffold




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the k up and down neighboor orfs from a list of orfs')
    parser.add_argument('orfList_filename', help='the path of the ORFLIST_FILE, this file contains the list of ORFs you want to extract the neighborood (one ORF per line)')
    parser.add_argument('feature_filename',help='the path of the FEATURE_FILE, this file is a gff-like file, it contains the locations of the ORFs for every scaffolds')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('-k',type=int,default=5,help='the number of neighboors up and down you want to extract (default: 5)')

    args = parser.parse_args()

    if os.path.exists(args.orfList_filename) :
        orfList_filename = os.path.abspath(args.orfList_filename)
    else:
        sys.exit(args.orfList_filename+' does not exist, exit')

    if os.path.exists(args.feature_filename) :
        feature_filename = os.path.abspath(args.feature_filename)
    else:
        sys.exit(args.feature_filename+' does not exist, exit')



    output_filename = os.path.abspath(args.output_filename)

    k = args.k
    print(k)
    orfSet = set()
    file = open(orfList_filename,'r')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        orfSet.add(liste[0])
    file.close()

    genome2scaffold2orfList,orf2genome,orf2scaffold = genome2orfOrder(feature_filename)



    orf2genomicContext = dict()
    for orf in orfSet :
#        print(orf)
        try:
            genome = orf2genome[orf]
        except:
            continue
        scaffold = orf2scaffold[orf]
        orf2genomicContext[orf] = kNearestNeighbors(orf,genome2scaffold2orfList[genome][scaffold],k)


    output = open(output_filename,'w')
    output.write('centroid'+'\t'+'genome'+'\t'+'scaffold'+'\t'+'start'+'\t'+'end'+'\t'+'strand'+'\t'+'orfname'+'\n')
    for orf in orfSet :
#        print(orf)
        try:
            genome = orf2genome[orf]
        except:
            continue
        scaffold = orf2scaffold[orf]
        for elt in orf2genomicContext[orf] :
            output.write(orf+'\t'+genome+'\t'+scaffold+'\t'+str(elt[0])+'\t'+str(elt[1])+'\t'+elt[2]+'\t'+elt[3]+'\n')

    output.close()

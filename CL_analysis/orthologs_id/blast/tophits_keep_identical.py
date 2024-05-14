#!/usr/bin/python2.6
# -*- coding: utf-8 -*-

''' Get RECIPROCAL ORTHOLOGS
This script takes the outputs blast files of a reciprocal blast. 
First, identifies tophit for each blast. Minimum 30 pidentity, then picks the tophit with greatest bitscore. 
If bitscores are identical then the blast hit with greatest pidentity is picked as the tophit. 
If multiple sequenced have the same bitscore and pidentity, the ortholog is discarded.
Second, finds 1-1 reciprocal orthologs. The blastout output format supported is:
"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"
'''
#==============================================================================
import argparse
import sys
import os
from itertools import combinations
import cPickle as pickle
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="A folder containing blastout files")
parser.add_argument("outfile", type=str,
                    help="The prefix of the outfile [script adds .pkl and .txt]")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".bla")]

def find_pairs_ensembl(infiles):
    '''Finds files in a list that have the form a <--> b and b <--> a
    Where a and b are different species. Each pair will be used to determine
    the reciprocal_besthit later. Returns list of pairs
    '''
    pairs = []
    for a, b in combinations(infiles, 2):
        abase = os.path.basename(a).split(".")
        bbase = os.path.basename(b).split(".")
        a_db = abase[0].split("_")[0]
        a_query = abase[1].split("_")[0]
        b_db = bbase[0].split("_")[0]
        b_query = bbase[1].split("_")[0]
        if a_db == b_query and a_query == b_db:
            pairs.append((a,b))
    return pairs

def read_blastout(source):
    countdeleted = 0
    try:
        with open(source, "r") as infile:
            sbase = os.path.basename(source).split(".")
            s_db = sbase[0].split("_")[0]
            s_query = sbase[1].split("_")[0]
            count =0
            blasthits = {}
            blastlist = []
            currentscaffold = None
            for line in infile:
                count +=1
                qseqid = line.rstrip().split()[0]
                if currentscaffold is None:
                    currentscaffold = qseqid
                    blastlist.append(line)
                else:
                    if currentscaffold == qseqid:
                        blastlist.append(line)
                    else:
                        blasthits, countdel = read_blastout_bitscore(blastlist, blasthits, sbase)
                        countdeleted = countdeleted + countdel
                        blastlist = []
                        currentscaffold = qseqid
                        blastlist.append(line)
            #final loop for final scaffold
            blasthits, countdel = read_blastout_bitscore(blastlist, blasthits, sbase)
            countdeleted = countdeleted + countdel
            print source
            print "Number of lines =", count
            print "Number of blast hits filtered because no one tophit =", countdeleted
        return blasthits
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % source
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def read_blastout_bitscore(blastlist, blasthits, sbase):
    ''' The blastout output format supported is:
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"
    For each hit only the besthit is stored, rest is discarded.
    Minimum 30% identity, then pick tophit with greatest bitscore. If bitscores are identical
    pick the tophit with greatest %identity. If equal %identity then ortholog discarded.
    '''
    if len(blastlist) == 0:
        line = line.rstrip().split()
        s_db = sbase[0].split("_")[0]
        s_query = sbase[1].split("_")[0]
        qseqid = s_query+"_"+line[0]
        sseqid = s_db+"_"+line[1]
        pidentity = float(line[2])
        bitscore = float(line[11])
        blasthits[qseqid] = (sseqid, bitscore, pidentity)
    else:
        countdel = 0
        check = []
        for line in blastlist:
            line = line.rstrip().split()
            s_db = sbase[0].split("_")[0]
            s_query = sbase[1].split("_")[0]
            qseqid = s_query+"_"+line[0]
            sseqid = s_db+"_"+line[1]
            pidentity = float(line[2])
            bitscore = float(line[11])
            if pidentity > 30: # Check for min 30% identity
                if qseqid in blasthits:
                    if blasthits[qseqid][1] < bitscore:
                        blasthits[qseqid] = (sseqid, bitscore, pidentity)
                    elif blasthits[qseqid][1] == bitscore:
                        if blasthits[qseqid][2] < pidentity:
                            blasthits[qseqid] = (sseqid, bitscore, pidentity)
                        elif blasthits[qseqid][2] == pidentity:
                            info = [bitscore, pidentity]
                            check.append(info)
                else:
                    blasthits[qseqid] = (sseqid, bitscore, pidentity)
        #check that there arent two hits which have the highest bitscore and pidentity but are identical
    return blasthits, countdel

def reciprocal_besthit(a,b):
    '''Finds the tophit in each file and determines if they are the same. 
    Identifies and returns list of 1-1 reciprocal orthologs.
    '''
    recip_besthit = []
    a_blastout = read_blastout(a)
    b_blastout = read_blastout(b)
    for g in a_blastout:
        query = g
        query_besthit = a_blastout[g][0]
        if query_besthit in b_blastout:
            target_besthit = b_blastout[query_besthit][0] 
            if query == target_besthit:
                recip_besthit.append((query, query_besthit))
    return recip_besthit

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = list_folder(args.infolder)
    pairs = find_pairs_ensembl(infiles)
    print len(pairs)
    print pairs
    recip_besthit_all = []
    for a, b in pairs:
        #pick tophit
        recip_besthit = reciprocal_besthit(a, b)
        #combine all tophit dictionaries into one
        recip_besthit_all.append(recip_besthit)
    print "Number of pairs = ", len(recip_besthit_all)

    print "Starting to pickle..."
    with open(args.outfile+".pkl", "w") as outfile:
        pickle.dump(recip_besthit_all, outfile, -1)
    
    print "Writing to file ..."
    count = 0
    with open(args.outfile+".txt", "w") as outtext:
        for gene in recip_besthit_all:
            for g in gene:
                count += 1
                outtext.write(str(g))
                outtext.write("\n")
    print "Number of reciprocal tophits =", count
        
if __name__ == '__main__':
    main()

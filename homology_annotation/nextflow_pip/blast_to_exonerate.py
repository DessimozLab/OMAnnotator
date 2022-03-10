#!/usr/bin/env python

from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

__author__='Yannis Nevers'

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BLAST output and generate a FASTA file with only the subsequence of the genomes matched this way, on which Exonerate will be run.')
    parser.add_argument('--input', '-i', required=True,
                    help='The input BLAST file to be read')
    parser.add_argument('--genome', '-g', required=True,
                    help='The input genome file to extract the sequence from.')
    parser.add_argument('--output', '-o', required=True,
                    help = 'The output FASTA file')
    args = parser.parse_args()
    
    return args


def read_blast_output(filename):
    qresult = SearchIO.read(filename, 'blast-tab',  comments=True)
    return qresult

def extract_loci(blast_output, extension=25000):
    all_loci = dict()
    for hit in blast_output:
        t_chr = hit.id
        if t_chr not in all_loci:
            all_loci[t_chr] = []
        for hsp in hit:
            start = hsp.hit_start+1
            end = hsp.hit_end
            if start > end:
                tp_start = start
                start = end
                end = tp_start
            if len(all_loci[t_chr])==0:
                all_loci[t_chr] = [[start-extension, end+extension]]
            else:
                overlap = False
                for cur_loc in all_loci[t_chr]:
                    if (start>=cur_loc[0] and start<=cur_loc[1]) or (end>=cur_loc[0] and end<=cur_loc[1]):
                            new_start = max(min(start-extension, cur_loc[0]),1)
                            new_end = max(end+extension, cur_loc[1])
                            cur_loc[0] = new_start
                            cur_loc[1] = new_end
                            overlap = True
                if not overlap :
                        all_loci[t_chr].append([max(start-extension,1),end+extension])
    for k,v in all_loci.items():
        fused_loci = list()
        for loc in sorted(v):
            if fused_loci and fused_loci[-1][1] >= loc[0]:
                fused_loci[-1][1] = loc[1]
            else:
                fused_loci.append(loc)
        all_loci[k] = fused_loci
    return all_loci


def extract_g_coordinate(genome_file, all_loci):
    all_subseq = []
    with open(genome_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in all_loci:
                wh_seq = str(record.seq)
                max_pos = len(record.seq)
                for pos in all_loci[record.id]:
                    sub_seq = wh_seq[pos[0]-1:min(pos[1],max_pos)]
                    all_subseq.append(SeqRecord(Seq(sub_seq),
                                                id=record.id+"_"+str(pos[0]-1)+"_"+str(min(pos[1],max_pos)-1)))
    return all_subseq

def write_genome_file(output, fcoord):
    with open(output, "w") as whandle:
        SeqIO.write(fcoord, whandle, "fasta")
        

if __name__=='__main__':
    args = parse_arguments()
    filename = args.input
    genome_file = args.genome
    output = args.output

    blast_output = read_blast_output(filename)
    loci = extract_loci(blast_output)
    final_seq = extract_g_coordinate(genome_file, loci)
    write_genome_file(output, final_seq)


#!/usr/bin/env python3
#-*-coding:utf8-*-

import sys, time, os
import pandas as pd 
import numpy as np
import argparse
import re
from  multiprocessing import Pool


def disp(txt, nt=0):
    if not args.quiet: sys.stderr.write('[Bam2Bed] @%s \t%s \n' %(time.asctime(), txt))


def markBed(bed,cglist):
    '''
    mark bed file base on the cglist
    '''
    regiondict = {}
    with os.popen(f"bedtools intersect -a {bed} -b {cglist} -wa -wb") as f1:
        for line in f1:
            line = line.strip().split("\t")
            region = line[0] + ":" + str(line[1]) +"-" + str(line[2])
            cgpos  = line[3] + ":" + str(line[4])
            if region in regiondict:
                regiondict[region].append(cgpos)
            else:
                regiondict[region] = [cgpos]
        return regiondict   
    

def MarkRead2Bed(bam,region,regioncglist):
    '''abstract information from reads'''
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T'),'f':('C','T','G','A'),'r':('G','A','C','T')}
    outlines=[]
    for line in bam:
        pattern  = ""
        col = line.strip().split('\t')
        if len(col) <10:continue
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', abs(int(col[8]))
        if "ZS:Z:" in line:
            strand_index = line.find('ZS:Z:')
        elif "YD:Z:" in line:
            strand_index = line.find('YD:Z:')
        else:
            disp('missing strand information "ZS:Z:xx" or "YD:Z:x"')
            exit()
        strand  = line[strand_index+5:strand_index+7]
        gap_pos, gap_size = 0, 0
        while 'I' in cigar or 'D' in cigar:
            for sep in 'MID':
                try: gap_size = int(cigar.split(sep, 1)[0])
                except ValueError: continue
                break
            if sep == 'M': gap_pos += gap_size
            elif sep == 'I': seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            elif sep == 'D':
                seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
                gap_pos += gap_size
            cigar = cigar[cigar.index(sep)+1:]
        strand = strand[0]
        pos2 = pos + len(seq)
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        for cglist in regioncglist:
            index = int(cglist.split(":")[1])
            if index > pos and index < pos2:
                seq_index = index - pos
                if (seq[seq_index] == match or seq[seq_index] == rc_match):
                    pattern += "C"
                elif (seq[seq_index] == convert or seq[seq_index] == rc_convert):
                    pattern += "T"
                else:
                    pattern += "N"
            else:
                pattern += "N"
        outline = "\t".join([region,pattern,col[0],str(insert),','.join([x.split(":")[1] for x in regioncglist])])
        outlines.append(outline)
    return outlines

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='quantifying methylation level from aligned files')
    parser.add_argument('-i','--infiles',dest="infiles",help="input aligned files",required=True)
    parser.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
    parser.add_argument("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
    parser.add_argument('-r',"--region",dest="bedfile",metavar="FILE",help="the region of intersect",required=True)
    parser.add_argument('-cg',"--cglist",dest="cglist",metavar="FILE",help="the cglist file",default="/ehpcdata/Database/hg38/hg38_cpg_postion.bed")
    parser.add_argument("-b",'--bed-path',dest="bed_path",metavar="PATH",help="path to bedtools ,[default:none]",default="") 
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
    parser.add_argument("-t", "--thread", dest="thread", help="number of thread use for this program", default=12)
    
    args = parser.parse_args()

    if len(args.sam_path) > 0:
        if args.sam_path[-1] != '/': args.sam_path += '/'

    if len(args.outfile) == 0: disp("Missing output file name, write to STDOUT.")
    bam    = args.infiles
    bed    = args.bedfile
    cglist = args.cglist
    out    = args.outfile
  ## pool 
    disp("Marking Bed")
    regiondic = markBed(bed,cglist) 
    disp("Start Convert reads to Pattern")
    pool = Pool(int(args.thread))
    with open(out,'w') as f2,open(bed,'r') as f1:
        f2.write("region\tPattern\treadid\tInsert\tcglist\n")
        for  line in f1:
            line = line.strip().split("\t")
            region = line[0] +":" + str(line[1]) + "-" + str(line[2])
            ## read bamfile
            try:
                outbam = os.popen(f"samtools view -F 3844 {bam} {region}").read()
                outbam = outbam.split("\n")
            except IOError as e:
                print(f"No such bam file: {bam}" + e)

            if (region not in regiondic) or (len(regiondic[region]) <10):continue
            if int(args.thread) == 1:
                disp("Signle CPU mode")
                lines = MarkRead2Bed(outbam,region,regiondic[region])
                f2.write("\n".join(lines)+"\n")
            else:
                pool.apply_async(MarkRead2Bed,args=(outbam,region,regiondic[region],),callback=lambda x:f2.write("\n".join(x) +"\n"))
        pool.close()
        pool.join()

    disp("done") 

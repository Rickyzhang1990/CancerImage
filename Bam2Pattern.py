#!/usr/bin/env python3
#-*-coding:utf8-*-

import sys, time, os
import pandas as pd 
import numpy as np
import argparse
import re
import pysam
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
        
def read_bam(bam,bed):
    '''
     use generator to export bam data
    '''
    reads = []
    with open(bed,'r') as f1,pysam.AlignmentFile(bam,'rb') as bam:
        disp("start Reading Bed file")
        for line in f1:
            line = line.strip().split("\t")
            region = line[0] +":" + str(line[1]) + "-" + str(line[2]) 
            try:
                outbam = bam.fetch(line[0],int(line[1]),int(line[2]))
            except IOError as e:
                print(f"No such bam file: {bam}" + e)
            finally:
                yield region,outbam
        yield region,outbam

def MarkRead2Bed(bam,region,regiondict):
    '''abstract information from reads'''
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T'),'f':('C','T','G','A'),'r':('G','A','C','T')}
#    disp("start Marking reads from bed")
    outlines=[]
    for line in bam:
        pattern  = ""
#        col = line.strip().split('\t')
    #    print(len(col))
#        if len(col) < 10: continue
    #    flag = col[1]
   #     if len(col) <10:    print(line)
#        try:
        if line.is_qcfail:continue
        cr, pos, cigar, seq, strand, insert = line.reference_name, int(line.reference_start), line.cigarstring, line.query, '', abs(line.template_length) 
        if len(line.tags) == 3:
            strand_index = line.tags[2][1]
#        assert strand_index >= 0, disp('missing strand information "ZS:Z:xx"')
#        strand = line[strand_index+5:strand_index+7]
        elif len(line.tags) == 7:
            strand_index = line.tags[5][1]
#            strand = line[strand_index+5:strand_index+7]
        else:
            disp('missing strand information "ZS:Z:xx" or "YD:Z:x"')
            exit()
       
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
#        if pos + len(seq) >= len(ref[cr]): continue
        strand = strand[0]
        pos2 = int(line.reference_end)
        match, convert, rc_match, rc_convert = BS_conversion[strand]
#        for i in re.finditer(match,refseq[pos:pos2]):
#        print(regiondict[region])
        for cglist in regiondict[region]:
            index = int(cglist.split(":")[1])
            if index > pos and index < pos2:
                seq_index = index - pos
#                print(cglist,pos,seq_index,len(seq))
#            index = int(i.span()[0])
                if (seq[seq_index] == match or seq[seq_index] == rc_match):
                    pattern += "C"
                elif (seq[seq_index] == convert or seq[seq_index] == rc_convert):
                    pattern += "T"
            else:
#            elif index <= pos or index >= pos2:
                pattern += "N"
#            else:pass
        outline = "\t".join([region,pattern,col[0],str(abs(insert))])
#        print(outline)
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

    bed    = args.bedfile
    cglist = args.cglist
  ## pool 
    regiondic = markBed(bed,cglist) 
    pool = Pool(int(args.thread))
    with open(args.outfile,'w') as f2:
        for region,reads in read_bam(args.infiles,bed):
            print(regiondic['chr1:1041093-1041213'])
            if region not in regiondic:continue
            reads = reads.split("\n")
#            print(len(reads))
#            lines = MarkRead2Bed(reads,region,regiondic)
#            f2.write("\n".join(lines)+"\n")
            pool.apply_async(MarkRead2Bed,(reads,region,regiondic,),callback=lambda x:f2.write("\n".join(x) +"\n"))

    disp("done") 

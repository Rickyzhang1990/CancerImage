#!/usr/nin/env python3
#-coding:utf8-*-


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
 
def MarkRead2Bed(bam,region,regioncglist):
    '''abstract information from reads'''
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T'),'f':('C','T','G','A'),'r':('G','A','C','T')}
    disp("start Marking reads from bed")
    outlines=[]
    for line in bam:
        pattern  = ""
        if line.is_qcfail:continue
        cr, pos, cigar, seq, strand, insert = line.reference_name, int(line.reference_start), line.cigarstring, line.query, '', abs(line.template_length) 
        if len(line.tags) == 3:
            strand_index = line.tags[2][1]
        elif len(line.tags) == 7:
            strand_index = line.tags[5][1]
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
        strand = strand_index[0]
        pos2 = int(line.reference_end)
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
        outline = "\t".join([region,pattern,line.query_name,str(insert)])
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
    out    = args.outfile
  ## pool 
    regiondic = markBed(bed,cglist)
    print(len(regiondic)) 
    pool = Pool(int(args.thread))
    with open(out,'w') as f2:
        for region,reads in read_bam(args.infiles,bed):
#            print(region,reads)
#            print(len(regiondic))
            if (region not in regiondic) or len(regiondic[region]) <10 :
#                print(region)
                continue
#            else:
#                if len(regiondic[region]) <10:continue
#                else:pass
#            reads = [read for read in reads]
            if int(args.thread) == 1:
#                disp("Signle CPU mode")
                lines = MarkRead2Bed(reads,region,regiondic[region])
                f2.write("\n".join(lines)+"\n")
#            else:
#                pool.apply_async(MarkRead2Bed,args=(reads,region,regiondic[region],),callback=lambda x:f2.write("\n".join(x) +"\n"))
#        pool.close()
#        pool.join()

    disp("done") 

#!/usr/bin/env python3
#-*-coding:utf8-*-

import sys, time, os, array
from  multiprocessing import Pool
import pandas as pd 
import numpy as np
import _pickle  as cpickle
import argparse
import re



def disp(txt, nt=0):
    if not args.quiet: sys.stderr.write('[methratio] @%s \t%s' %(time.asctime(), txt))

def read_bam(bam,chunksize=1000):
    '''
     use generator to export bam data
    '''
    reads = []
    with os.popen(f"samtools view -F 3844  {bam}") as f1:
        for line in f1:
            line = line.strip()
            if not len(reads) >= chunksize:
                reads.append(line)
            else:
                yield reads
        yield reads               

def get_alignment(line):
    '''extract information from bamline'''

    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = col[1]
        if 'u' in flag: return []
        if args.unique and 's' in flag: return []
        if args.pair and 'P' not in flag: return []
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
        if cr not in args.chroms: return []
        strand_index = line.find('ZS:Z:')
        assert strand_index >= 0, 'missing strand information "ZS:Z:xx"'
        strand = line[strand_index+5:strand_index+7]
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
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if args.unique and flag != 'UM': return []
        if args.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in args.chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if args.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if args.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-args.trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[args.trim_fillin:], pos+args.trim_fillin
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (seq, strand[0], cr, pos)

def MarkRead2Bed(reads,ref,refmark):
    '''abstract information from reads'''
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
    outlines = []
    for line in reads:
        pattern,index = "",[]
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info
        pos2 = pos + len(seq)
        refseq, refmarkcr = ref[cr], refmark[cr]
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        for i in re.finditer(match,refseq[pos:pos2]):
            index = int(i.span()[0])
            if (seq[index] == match or seq[index] == rc_match) and refmarkcr[pos+index] in seq_context:
                pattern += "C"
                index.append(str(pos+index))
            elif (seq[index] == convert or seq[index] == rc_convert) and refmarkcr[pos+index] in seq_context:
                pattern += "T"
                index.append(str(pos+index))
            else:pass
        outline = "\t".join([cr,str(pos),str(pos2),pattern,",".join(index)])
        outlines.append(outline)
    return outlines

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='quantifying methylation level from aligned files')
    parser.add_argument('infiles',help="input aligned files")
    parser.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
    parser.add_argument("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
    parser.add_argument("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
    parser.add_argument("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
    parser.add_argument("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
    parser.add_argument("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
    parser.add_argument("-t", "--trim-fillin", dest="trim_fillin", type=int, metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
    parser.add_argument("-T", "--thread", dest="thread", help="number of thread use for this program", default=12)
    parser.add_argument("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
    
    args = parser.parse_args()

    if len(args.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
    if len(args.infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.")
    if len(args.chroms) > 0: args.chroms = set(args.chroms.split(','))
    if args.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
    seq_context_str, CG, CHG, CHH = ['CG','CHG','CHH'], 1, 2, 3
    if len(args.context) > 0:
        args.context = set(args.context.upper().split(','))
        try: seq_context = set([seq_context_str.index(cc)+1 for cc in args.context])
        except ValueError: parser.error('Invalid -x value, must be one or more of "CG", "CHG", or "CHH"')
    else: seq_context = set([1, 2, 3])

    if len(args.sam_path) > 0:
        if args.sam_path[-1] != '/': args.sam_path += '/'

    if len(args.outfile) == 0: disp("Missing output file name, write to STDOUT.")
    # Read in chromosomes
    ref, cr, seq = {}, '', ''
    disp('loading reference genome file: %s ...' % args.reffile)
    for line in open(args.reffile):
        if line[0] == '>':
            if len(cr) > 0:
                if len(args.chroms) == 0 or cr in args.chroms: ref[cr] = seq.upper()
            cr, seq = line[1:-1].split()[0], ''
        else: seq += line.strip()

    if len(args.chroms) == 0 or cr in args.chroms: ref[cr] = seq.upper()
    del seq
    args.chroms = set(ref.keys())
    disp('marking reference genome ...')
    if not os.path.isfile(args.reffile+".methmark"):
        refmark = {}
        for cr in ref:
            refmark[cr] = array.array('b', [0]) * len(ref[cr])   ##kaka
            refcr, refmarkcr = ref[cr], refmark[cr]
            index = refcr.find('C', 0, len(refcr)-2)
            while index >= 0:
                if refcr[index+1] == 'G': refmarkcr[index] = CG
                elif refcr[index+2] == 'G': refmarkcr[index] = CHG
                else: refmarkcr[index] = CHH
                index = refcr.find('C', index+1, len(refcr)-2)
            index = refcr.find('G', 2, len(refcr))
            while index >= 0:
                if refcr[index-1] == 'C': refmarkcr[index] = CG
                elif refcr[index-2] == 'C': refmarkcr[index] = CHG
                else: refmarkcr[index] = CHH
                index = refcr.find('G', index+1, len(refcr))
        if os.access(os.path.dirname(os.path.abspath(args.reffile)),os.W_OK):
            with open(args.reffile+".methmark",'wb') as dicref:
                cpickle.dump(refmark,dicref)
        else:
            print("Please ensure that you have the write access to reference file foler")
            pass
    else:
        with open(args.reffile+".methmark",'rb') as dicref:
            refmark = cpickle.load(dicref)
  ## pool 
    pool = Pool(int(args.thread))
    with open(args.outfile,'w') as f2:
        for reads in read_bam(args.infile):
            pool.apply_async(MarkRead2Bed,args=(reads,ref,refmark,),callback=lambda x:f2.write("\n".join(x)))

    print("done") 
        


    

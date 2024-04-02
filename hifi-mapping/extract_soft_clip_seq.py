#!/bin/env python
import sys
import pysam

inbam = pysam.AlignmentFile(sys.argv[1],'r')
#outbam = pysam.AlignmentFile(sys.argv[2],'wb',template=inbam)
    

def is_soft_cipped(contig):
    is_sc = False
    if contig.query_alignment_start > 1000:
        is_sc = True
    if (contig.query_length - contig.query_alignment_end) > 1000:
        is_sc = True
    return is_sc
    
for contig in inbam:
    if contig.is_secondary:
        continue
    if contig.is_supplementary:
        continue
    if not is_soft_cipped(contig):
        out = [[contig.query_name,contig.query_sequence]]
        #outbam.write(contig)        
    else:
        left_soft_clip_end = contig.query_alignment_start
        right_soft_clip_start = contig.query_alignment_end
        out = [
            ["%s_left" % contig.query_name,contig.query_sequence[0:left_soft_clip_end]],
            ["%s_right" % contig.query_name,contig.query_sequence[right_soft_clip_start:]],
            [contig.query_name,contig.query_alignment_sequence]
            ]
    for seq_name,seq in out:
        if len(seq) != 0:
            print(">%s\n%s" % (seq_name,seq))

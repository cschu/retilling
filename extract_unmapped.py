#!/usr/bin/env python
import sys
import csv

from ktoolu_io import readFastq, openFile

# call with samtools view -f2 or on bamfile with only -f2-mapping reads
#mappedReads = se()
#for row in csv.reader(sys.stdin, delimiter='\t'):
#    mappedReads.add(row[0])
mappedReads = set(line.strip().split('\t')[0] for line in sys.stdin)

fq1 = readFastq(sys.argv[1])
fq2 = readFastq(sys.argv[2])

with openFile(sys.argv[1].replace('.fq.gz', '.unmapped_iwgsc10.fq.gz'), fmt='gz', mode='wt') as fo1, openFile(sys.argv[2].replace('.fq.gz', '.unmapped_iwgsc10.fq.gz'), fmt='gz', mode='wt') as fo2:
    while 1:
        try:
            _id1,_seq1,_qual1 = next(fq1)
        except:
            break
        _id2,_seq2,_qual2 = next(fq2)
        
        if _id1[1:] not in mappedReads and _id2[1:] not in mappedReads:
            fo1.write('{}\n{}\n+\n{}\n'.format(_id1, _seq1, _qual1))
            fo2.write('{}\n{}\n+\n{}\n'.format(_id2, _seq2, _qual2))
       




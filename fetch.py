#!/usr/bin/python env
import sys
import re

fa = sys.argv[1]
target = sys.argv[2]

seq_dic = {}
with open (fa, 'r') as fi:
   for each in fi:
        if each.startswith('>'):
            chr_id = each.split()[0].replace('>', '')
            seq_dic.setdefault(chr_id,[])
        else:
            seq_dic.setdefault(chr_id,[]).append(each.strip())

#for chr_id in range(1,22) + ['X','Y']:
for chr_id, seq in seq_dic.items():
    seq = seq_dic[chr_id]
    for a in re.finditer(target, ''.join(seq)):
        start, end = a.span()
        start = start + 1
        end = end +1
        print ("{}\t{}\t{}".format(chr_id, start, end))

from snapgene_file_to_seq_record import snapgene_file_to_seq_record

from Bio import SeqIO

import os
import posixpath

TEST_DIR = './test_samples'

files = os.listdir(TEST_DIR)
for fname in files:
    #filter files only end with 'dna'
    _, ext = posixpath.splitext(fname)
    if ext != '.dna':
        continue
    fpath = posixpath.join(TEST_DIR, fname)
    j = snapgene_file_to_seq_record(fpath)
    with open(posixpath.join(TEST_DIR, fname+'.gb'), 'w') as f:
        SeqIO.write([j,], f, 'genbank')
    print(j)
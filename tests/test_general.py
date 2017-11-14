from SnapGeneFileReader import snapgene_file_to_seq_record
from Bio import SeqIO
import posixpath
import os

TEST_DIR = './test_samples'

def test_parse(tmpdir):
    files = os.listdir(TEST_DIR)
    print('start')
    i = 0
    for fname in files:
        print(fname)
        #filter files only end with 'dna'
        _, ext = posixpath.splitext(fname)
        if ext != '.dna':
            continue
        i += 1
        print(fname)
        fpath = posixpath.join(TEST_DIR, fname)
        j = snapgene_file_to_seq_record(fpath)
        print(j)
        assert j
        assert j.name
        assert j.seq
        with open(posixpath.join(tmpdir, fname+'.gb'), 'w') as f:
            SeqIO.write([j,], f, 'genbank')
    print('processed {} files'.format(i))
    assert i > 0
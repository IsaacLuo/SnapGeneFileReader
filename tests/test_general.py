from snapgene_reader import snapgene_file_to_seqrecord, snapgene_file_to_gbk
from Bio import SeqIO
import os

TEST_DIR = os.path.join('tests', 'test_samples')


def test_snapgene_file_to_seqrecord(tmpdir):
    all_files = [f for f in os.listdir(TEST_DIR) if f.endswith('.dna')]
    assert len(all_files)
    for fname in all_files:
        fpath = os.path.join(TEST_DIR, fname)
        record = snapgene_file_to_seqrecord(fpath)
        assert len(record.seq) > 10
        with open(os.path.join(str(tmpdir), fname + '.gb'), 'w') as f:
            SeqIO.write([record, ], f, 'genbank')

def test_snapgene_file_to_gbk(tmpdir):
    all_files = [f for f in os.listdir(TEST_DIR) if f.endswith('.dna')]
    assert len(all_files)
    for fname in all_files:
        print (fname)
        fpath = os.path.join(TEST_DIR, fname)
        target = os.path.join(str(tmpdir), 'testfile.gbk')
        with open(fpath, 'r') as fsource:
            with open(target, 'w') as ftarget:
                snapgene_file_to_gbk(fsource, ftarget)
        verification_reparsing = SeqIO.read(target, 'genbank')

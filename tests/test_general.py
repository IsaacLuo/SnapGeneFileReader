from snapgene_reader import snapgene_file_to_seqrecord
from Bio import SeqIO
import os

TEST_DIR = os.path.join('tests', 'test_samples')


def test_parse(tmpdir):
    all_files = [f for f in os.listdir(TEST_DIR) if f.endswith('.dna')]
    assert len(all_files)
    for fname in all_files:
        fpath = os.path.join(TEST_DIR, fname)
        record = snapgene_file_to_seqrecord(fpath)
        assert len(record.seq) > 10
        with open(os.path.join(str(tmpdir), fname + '.gb'), 'w') as f:
            SeqIO.write([record, ], f, 'genbank')

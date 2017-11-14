# SnapGeneFileReader
a Python project to read and write Snapgene *.dna

## install
```bash
pip install SnapGeneFileReader
```
## test:
```bash
python runTest.py
```
## usage:

### to read a SnapGene file into a dict
```python
from SnapGeneFileReader import SnapGeneFileReader

file_path = './snap_gene_file.dna'
dna_dict = SnapGeneFileReader.read_file(file_path)
```

### to convert dna file to genbank file:
```python
from Bio import SeqIO
from SnapGeneFileReader import snapgene_file_to_seq_record

file_path = './snap_gene_file.dna'

seqObject = snapgene_file_to_seq_record(file_path)
with open('genbank_file.gb', 'w') as f:
    SeqIO.write([seqObject,], f, 'genbank')
```

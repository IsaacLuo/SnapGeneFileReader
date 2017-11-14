# SnapGeneFileReader
a Python project to read and write Snapgene *.dna

## requirement:
```python
pip install xmltodict
pip install biopython
```
## test:
```bash
python runTest.py
```
## usage:
### to convert dna file to genbank file:
```python
from Bio import SeqIO
from snapgene_file_to_seq_record import snapgene_file_to_seq_record
seqObject = snapgene_file_to_seq_record(fpath)
with open('new_file_name.gb', 'w') as f:
    SeqIO.write([seqObject,], f, 'genbank')
```

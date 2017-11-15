SnapGene Reader
===============

SnapGene Reader is a Python library to parse Snapgene ``*.dna`` files into dictionnaries or Biopython SeqRecords:

.. code:: python

  from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord

  file_path = './snap_gene_file.dna'
  dictionnary = snapgene_file_to_dict(filepath)
  seqrecord = snapgene_file_to_seqrecord(filepath)

Installation
------------

Install with PIP:

.. code:: bash

    pip install snapgene_reader

Test with Pytest:

.. code:: bash

    python -m pytest
    (or simply "pytest")

Licence = MIT
-------------

SnapGene Reader is an open-source software originally written by `Isaac Luo <https://github.com/IsaacLuo>`_ at the Edinburgh Genome Foundry and released on Github <https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader>`_ under the MIT licence (¢ Edinburg Genome Foundry).

```

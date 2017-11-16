SnapGene Reader
===============

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/SnapGeneReader.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/SnapGeneReader
   :alt: Travis CI build status

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

SnapGene Reader is an open-source software originally written by `Isaac Luo <https://github.com/IsaacLuo>`_ at the Cai Lab. This fork is released on `Github <https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader>`_ under the MIT licence (¢ Isaac Luo and Cai Lab) and maintained by the Edinburgh Genome Foundry. Everyone is welcome to contribute.

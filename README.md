# GeneAbacus Python

This repository provides Python code for [GeneAbacus](https://github.com/vejnar/geneabacus). Using the `profileio` module, you can import profiles exported from GeneAbacus using the *binary* format. Profiles are loaded into [Numpy](https://numpy.org) arrays, usable for analysis.

From high-throughput sequencing **mapped reads** ([SAM/BAM](https://samtools.github.io/hts-specs/)), **GeneAbacus**:
* Creates **profiles** representing coverage depth per nucleotide,
* **Counts** reads mapped within user selected features such as chromosomes or genes.

## Download

See [tags](/../../tags) page.

## Install

```bash
pip3 install geneabacus
```

## Reading profiles from Python

```python
import geneabacus.profileio
profiles = geneabacus.profileio.pfopen('profiles.bin.lz4', 'danrer_cdna_protein_coding_ensembl104.fon1.json')
```

To get a transcript profile:
```python
profiles['ENSDART00000000486']
# will return
array([0., 0., 21., ..., 0., 3., 0.], dtype=float32)
```

## License

*GeneAbacus* is distributed under the Mozilla Public License Version 2.0 (see /LICENSE).

Copyright (C) 2015-2021 Charles E. Vejnar

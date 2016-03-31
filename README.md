[![Build Status](https://travis-ci.org/cdiener/corda.svg?branch=master)](https://travis-ci.org/cdiener/corda)
[![codecov.io](https://codecov.io/github/cdiener/corda/coverage.svg?branch=master)](https://codecov.io/github/cdiener/corda?branch=master)
[![Code Health](https://landscape.io/github/cdiener/corda/master/landscape.svg?style=flat)](https://landscape.io/github/cdiener/corda/master)

# CORDA for Python

This is a Python implementation based on the paper of Schultz et. al.

[Reconstruction of Tissue-Specific Metabolic Networks Using CORDA](http://journals.plos.org/ploscompbiol/article/authors?id=10.1371%2Fjournal.pcbi.1004808)

## What does it do?

CORDA, short for Cost Optimization Reaction Dependency Assessment is a method
for the reconstruction of metabolic networks from a given reference model
(a database of all known reactions) and a confidence mapping for reactions.
It allows you to reconstruct metabolic models for tissues, patients or specific
experimental conditions from a set of transcription or proteome measurements.

## How do I install it

CORDA for Python works only for Python 3.4+ and requires
[cobrapy](http://github.com/opencobra/cobrapy) to work. We recommend
installation via the [anaconda or miniconda](https://www.continuum.io/downloads)
distribution. After installing anaconda or miniconda you can install cobrapy
from the bioconda repository

```bash
conda install -c bioconda cobra
```
After that you can install CORDA using the pip from conda

```bash
pip install https://github.com/cdiener/corda/archive/master.zip
```

After CORDA for Python comes out of its infancy I will prepare a conda package
as well. For now the master branch is usually working and tested whereas all
new stuff is kept in the devel branch.

## What do I need to run it?

CORDA requires a base model including all reactions that could possibly included
such as Recon 1/2 or HMR. You will also need gene expression or proteome data
for our tissue/patient/experimental setting. This data has to be translated into
5 distinct classes: unknown (0), not expressed/present (-1), low confidence (1),
medium confidence (2) and high confidence (3). CORDA will then ensure to include
as many high confidence reactions as possible while minimizing the inclusion
of absent (-1) reactions while maintaining a set of metabolic requirements.

## What's the advantage over other reconstruction algorithms

I would say there are two major advantages:

1. It does not require any commercial solvers, in fact it works fastest with
   the free glpk solver that already comes together with cobrapy. For instance
   for the small central metabolism model included in data/ the reconstruction
   uses the following times:
   - cglpk: 0.04 s
   - cplex: 0.53 s
   - gurobi: 0.21 s
2. It's fast. CORDA for Python uses a strategy similar to FastFVA, where a previous
   solution basis is recycled repeatedly (speed-up of ~4-10 times). A normal
   reconstruction for Recon 1 with mCADRE can take several hours. With the original
   Matlab implementation of CORDA this takes about 40 minutes and with CORDA
   for Python it takes less than 5 minutes. A Recon 2 reconstruction can be
   achieved in less than 30 minutes.

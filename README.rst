|Build Status| |codecov.io| |Code Health|

CORDA for Python
================

This is a Python implementation based on the paper of Schultz et. al.

`Reconstruction of Tissue-Specific Metabolic Networks Using
CORDA <http://journals.plos.org/ploscompbiol/article/authors?id=10.1371%2Fjournal.pcbi.1004808>`__

What does it do?
----------------

CORDA, short for Cost Optimization Reaction Dependency Assessment is a
method for the reconstruction of metabolic networks from a given
reference model (a database of all known reactions) and a confidence
mapping for reactions. It allows you to reconstruct metabolic models for
tissues, patients or specific experimental conditions from a set of
transcription or proteome measurements.

How do I install it
-------------------

CORDA for Python works only for Python 3.4+ and requires
`cobrapy <http://github.com/opencobra/cobrapy>`__ to work. After having
a working Python installation with pip (Anaconda or Miniconda works fine
here as well) you can install corda with pip

.. code:: bash

    pip install corda

This will download and install cobrapy as well. I recommend using a
version of pip that supports manylinux builds for faster installation
(pip>=8.1).

For now the master branch is usually working and tested whereas all new
features are kept in its own branch. To install from the master branch
directly use

.. code:: bash

    pip install https://github.com/cdiener/corda/archive/master.zip

What do I need to run it?
-------------------------

CORDA requires a base model including all reactions that could possibly
included such as Recon 1/2 or HMR. You will also need gene expression or
proteome data for our tissue/patient/experimental setting. This data has
to be translated into 5 distinct classes: unknown (0), not
expressed/present (-1), low confidence (1), medium confidence (2) and
high confidence (3). CORDA will then ensure to include as many high
confidence reactions as possible while minimizing the inclusion of
absent (-1) reactions while maintaining a set of metabolic requirements.

How do I use it?
----------------

A small tutorial is found at https://cdiener.github.io/corda.

What's the advantage over other reconstruction algorithms
---------------------------------------------------------

I would say there are two major advantages:

1. It does not require any commercial solvers, in fact it works fastest
   with the free glpk solver that already comes together with cobrapy.
   For instance for the small central metabolism model (101 irreversible
   reactions) included together with CORDA the reconstruction uses the
   following times:

   cglpk: 0.02 s
   cplex: 0.30 s
   gurobi: 0.12 s
   mosek: 0.23 s

2. It's fast. CORDA for Python uses a strategy similar to FastFVA, where
   a previous solution basis is recycled repeatedly (speed-up of ~4-10
   times). A normal reconstruction for Recon 1 with mCADRE can take
   several hours. With the original Matlab implementation of CORDA this
   takes about 40 minutes and with CORDA for Python it takes less than 5
   minutes. A Recon 2 reconstruction can be achieved in less than 30
   minutes.

.. |Build Status| image:: https://travis-ci.org/cdiener/corda.svg?branch=master
   :target: https://travis-ci.org/cdiener/corda
.. |codecov.io| image:: https://codecov.io/github/cdiener/corda/coverage.svg?branch=master
   :target: https://codecov.io/github/cdiener/corda?branch=master
.. |Code Health| image:: https://landscape.io/github/cdiener/corda/master/landscape.svg?style=flat
   :target: https://landscape.io/github/cdiener/corda/master

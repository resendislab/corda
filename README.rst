|actions| |codecov.io|

CORDA for Python
================

This is a Python implementation based on the papers of Schultz et. al. with
some added optimizations. It is based on the following two publiactions:

- `Reconstruction of Tissue-Specific Metabolic Networks Using
  CORDA <http://journals.plos.org/ploscompbiol/article/authors?id=10.1371%2Fjournal.pcbi.1004808>`_
- `IDENTIFYING CANCER SPECIFIC METABOLIC SIGNATURES USING CONSTRAINT-BASED MODELS
  <http://dx.doi.org/10.1142/9789813207813_0045>`_

This Python package is developed in the
`Human Systems Biology Group <https://resendislab.github.io>`_ of
Prof. Osbaldo Resendis Antonio at the `National Institute of Genomic
Medicine Mexico <https://inmegen.gob.mx>`_ and includes recent updates to
the method (*CORDA 2*).


How to cite?
------------

This particular implementation of CORDA has not been published so far. In the
meantime you should if you cite the respective publications for the method
mentioned above and provide a link to this GitHub repository.

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

    pip install https://github.com/resendislab/corda/archive/master.zip

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

You can follow the [introduction](docs/index.ipynb).

What's the advantage over other reconstruction algorithms?
----------------------------------------------------------

No commercial solver needed
***************************

It does not require any commercial solvers, in fact it works fastest
with the free glpk solver that already comes together with cobrapy.
For instance for the small central metabolism model (101 irreversible
reactions) included together with CORDA the glpk version is a bout 3 times
faster than the fastest tested commercial solver (cplex).

Fast reconstructions
********************

CORDA for Python uses a strategy similar to FastFVA, where
a previous solution basis is recycled repeatedly.

Some reference times for reconstructing the minimal growing models for
iJO1366 (*E. coli*) and Recon3:

.. code::

    Python 3.10.8 (main, Oct 24 2022, 10:07:16) [GCC 12.2.0]
    Type 'copyright', 'credits' or 'license' for more information
    IPython 8.4.0 -- An enhanced Interactive Python. Type '?' for help.

    In [1]: from cobra.io import load_model

    In [2]: from corda import benchmark

    In [3]: ecoli = load_model("iJO1366")
    Restricted license - for non-production use only - expires 2023-10-25

    In [4]: opt = benchmark(ecoli)
    Running setup for model `iJO1366`.
    Running CORDA setup... ✔ [0.479 s]
    Running CORDA build... ✔ [7.44 s]
    Running validation on reduced model... ✔ [0.448 s]

    In [5]: print(opt)
    build status: reconstruction complete
    Inc. reactions: 447/2583
    - unclear: 0/0
    - exclude: 446/2582
    - low and medium: 0/0
    - high: 1/1


    In [6]: recon3 = load_model("Recon3D")

    In [7]: opt = benchmark(recon3)
    Running setup for model `Recon3D`.
    Running CORDA setup... ✔ [2 s]
    Running CORDA build... ✔ [13.7 s]
    Running validation on reduced model... ✔ [1.68 s]

    In [8]: print(opt)
    build status: reconstruction complete
    Inc. reactions: 114/10600
    - unclear: 0/0
    - exclude: 113/10599
    - low and medium: 0/0
    - high: 1/1

.. |actions| image:: https://github.com/resendislab/corda/actions/workflows/pythonpackage.yml/badge.svg
   :target: https://github.com/resendislab/corda/actions/workflows/pythonpackage.yml
.. |codecov.io| image:: https://codecov.io/github/resendislab/corda/coverage.svg?branch=main
   :target: https://codecov.io/github/resendislab/corda?branch=main

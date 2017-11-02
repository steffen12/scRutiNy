.. scrutiny documentation master file, created by
   sphinx-quickstart on Mon Oct 30 10:28:00 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scRutiNy's documentation!
====================================

This document describes how to use the Python package scRutiNy_ to analyze single-cell RNA-sequencing (scRNA-seq) data.

.. toctree::
   :maxdepth: 1

   README
   example

The package contains two modules

- :class:`RNAscrutiny.MatriSeq` : generate scRNA-seq data from a gene regulatory network correlation matrix.

- :class:`RNAscrutiny.RegNetInference` : Infer the gene regulatory network correlation matrix from scRNA-seq data.

.. _scRutiNy: https://pypi.python.org/pypi/RNAscrutiny

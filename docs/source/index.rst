tcrdist_rs documentation
=========================

tcrdist_rs is a Rust package with Rust bindings for Python for fast, parallelized CPU TCRdist, hamming, and Levenshtein distance routines, with the latter two based on the `triple_accel <https://github.com/Daniel-Liu-c0deb0t/triple_accel>`_ Rust package. 
Previous versions of TCRdist have been written in `python/numba <https://tcrdist3.readthedocs.io/en/latest/index.html>`_ and `C++ <https://github.com/phbradley/conga/tree/master/tcrdist_cpp>`_.
All code is freely available at `<https://github.com/zacmon/tcrdist_rs>`_.

At the moment, all functions associated with ``tcrdist_allele`` and ``tcrdist_gene`` operate assuming V genes are human V genes.
Future work will be done to allow for custom V gene databases.

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   install

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

.. toctree::
   :maxdepth: 0 
   :caption: API Docmentation

   api

.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   developer
   todo

References
==========
- Dash et al. (2017) "Quantifiable predictive features define epitope-specific T cell receptor repertoires." Nature 547(89-93). `<https://doi.org/10.1038/nature22383>`_
- Mayer-Blackwell et al. (2021) "TCR meta-clonotypes for biomarker discovery with tcrdist3 enabled identification of public, HLA-restricted clusters of SARS-CoV-2 TCRs." eLife 10:68605. `<https://doi.org/10.7554/eLife.68605>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

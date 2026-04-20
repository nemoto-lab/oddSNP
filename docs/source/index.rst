oddSNP documentation
====================

Donor multiplexing is a powerful strategy to increase scale, lower the costs, 
and reduce batch effects in single-cell RNA sequencing (scRNAseq), but clear 
guidelines for experimental design are lacking, forcing researchers to risk 
costly demultiplexing failures. To address this, we introduce SNP-Information 
Content (SNP-IC) and cell-paird SNP-Information Content (cpSNP-IC), quantitative 
metrics that can be computed from simple, unpooled pilot data and that 
accurately predict the success of demultiplexing. **oddSNP** is an open-source 
framework for computing these metrics, enabling in-silico titration of 
sequencing depth and donor complexity to optimize experimental design before 
committing to large-scale studies. 


Details on these metrics and the implementation of the tool are available in the 
manuscript entitled: **OddSNP: a predictive framework for optimizing multiplexed 
single-cell RNA-seq**. 

**oddSNP** is developed at the Nemoto-lab, The University of Osaka.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   notebooks/tutorial
   api/index


How to Cite
-----------

If you use *oddSNP* in your research, please cite:

| **OddSNP: a predictive framework for optimizing multiplexed single-cell RNA sequencing**
| Allendes Osorio, R.S., Nishimura, T., Shigihara, Y., Kimura, M., Takebe, T. and Nemoto, T. (2025)  
| https://www.biorxiv.org/content/10.64898/2025.12.08.692882v1

BibTeX
++++++

.. code-block:: bibtex

    @article{osorio2025oddsnp,
      title   = {OddSNP: a predictive framework for optimizing multiplexed single-cell RNA sequencing},
      author  = {Allendes Osorio, R.S. and Nishimura, T. and Shigihara, Y. and Kimura, M. and Takebe, T. and Nemoto, T.},
      journal = {bioRxiv},
      year    = {2025},
      doi     = {10.64898/2025.12.08.692882},
      url     = {https://www.biorxiv.org/content/10.64898/2025.12.08.692882v1},
      note    = {Preprint}
    }

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
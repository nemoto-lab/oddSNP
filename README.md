<p align="center">
  <img src="docs/source/_static/banner.png" width=300px>
</p>

Donor multiplexing is a powerful strategy to increase scale, lower the costs, and reduce batch effects in single-cell RNA sequencing (scRNAseq), but clear guidelines for experimental design are lacking, forcing researchers to risk costly demultiplexing failures. To address this, we introduce SNP-Information Content (SNP-IC) and cell-paird SNP-Information Content (cpSNP-IC), quantitative metrics that can be computed from simple, unpooled pilot data and that accurately predict the success of demultiplexing. *oddSNP* is an open-source framework for computing these metrics, enabling in-silico titration of sequencing depth and donor complexity to optimize experimental design before committing to large-scale studies. 

Details on these metrics and the implementation of the tool are available in the manuscript entitled: *OddSNP: a predictive framework for optimizing multiplexed single-cell RNA-seq* ([https://doi.org/10.64898/2025.12.08.692882](https://doi.org/10.64898/2025.12.08.692882)).

*oddSNP* is developed at the Nemoto-lab, The University of Osaka.

The full documentation of *oddSNP* is available at:
https://nemoto-lab.github.io/oddSNP/
and at 
https://oddsnp.readthedocs.io/

## Installation

The recommended way to install *oddSNP* is by using a virtual environments manager such as [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) (or [venv](https://docs.python.org/3/library/venv.html)). 

### Using Bioconda:

We create a new conda environment and directly install *oddSNP* from its bioconda source.

```bash
:~$ conda create --name oddsnp
:~$ conda activate oddsnp
(oddsnp):~$ conda install -c bioconda oddsnp
```

### Using PyPI:

Still, we recommend to install *oddSNP* inside a virtual environment. In this case, we need to make sure to also install `pip` to the created environment to avoid interfering with system libraries.

```bash
:~$ conda create --name oddsnp
:~$ conda activate oddsnp
(oddsnp):~$ conda install pip
(oddsnp):~$ pip install oddsnp
```

### From source:

Details on how to install *oddSNP* from source are given in the [Tutorial notebook](notebooks/tutorial.ipynb).

### After installation

An installation of [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/index.html) is required to perform pileup calculations within oddSNP. To install it, use the following command inside your activated conda environment:

```bash
(oddsnp):~$ conda install -c bioconda cellsnp-lite 
```

**NOTE** Other installation methods for `cellsnp-lite` are described in their [original website](https://cellsnp-lite.readthedocs.io/en/latest/index.html).



To check the installation finished properly, we can try and run *oddSNP* from the command line without any sub-commands. The output should be as follows:

```bash
(oddsnp):~$ oddSNP 
Usage: oddSNP [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  cpsnpic
  downsample
  genotype
  snpic
  utils
```

## Using oddSNP

For details on how to use *oddSNP* please refer to the accompanying tutorial notebook: [tutorial.ipynb](notebooks/tutorial.ipynb)

## 📄 How to Cite

If you use this repository in your research, please cite our bioRxiv preprint:

> **OddSNP: a predictive framework for optimizing multiplexed single-cell RNA sequencing**  
> Allendes Osorio, R.S., Nishimura, T., Shigihara, Y., Kimura, M., Takebe, T. and Nemoto, T. (2025)  
> https://www.biorxiv.org/content/10.64898/2025.12.08.692882v1

### BibTeX

```bibtex
@article{osorio2025oddsnp,
  title   = {OddSNP: a predictive framework for optimizing multiplexed single-cell RNA sequencing},
  author  = {Allendes Osorio, R.S. and Nishimura, T. and Shigihara, Y. and Kimura, M. and Takebe, T. and Nemoto, T.},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.64898/2025.12.08.692882},
  url     = {https://www.biorxiv.org/content/10.64898/2025.12.08.692882v1},
  note    = {Preprint}
}
```

Usage
=====

Command-line Interface
----------------------

oddSNP exposes several CLI commands grouped by submodule.

cpSNP-IC
^^^^^^^^

Run the full cpSNP-IC pipeline::

    python -m oddSNP cpsnpic run-all <INPATH> <OUPATH>

Run individual steps::

    python -m oddSNP cpsnpic calculate-cpsnpic <INPATH> <OUPATH>
    python -m oddSNP cpsnpic generate-histogram <CPFILE> <OUPATH>
    python -m oddSNP cpsnpic save-cpsnpic-plot <HISTOFILE> <OUTPUT>

SNP-IC
^^^^^^

::

    python -m oddSNP snpic [OPTIONS] COMMAND [ARGS]...

Python API
----------

You can also call oddSNP functions directly from Python::

    from oddSNP import cpsnpic

    cpsnpic.call_run_all(
        inpath='path/to/pileup/',
        oupath='path/to/output/',
        nproc=20,
        batch_size=100000,
        force=False
    )
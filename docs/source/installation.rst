Installation
============

.. note::
   The following instructions are for installing *oddSNP* directly from source code.
   On publication, *oddSNP* will be available for installation as a normal package
   from the Python Package Index (PyPI).

The recommended way to install *oddSNP* is by using a virtual environment manager
such as `Conda <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_
or `venv <https://docs.python.org/3/library/venv.html>`_.

Setting up a virtual environment
---------------------------------

Using Conda
^^^^^^^^^^^

.. code-block:: bash

   ~$ conda create --name oddsnp
   ~$ conda activate oddsnp
   (oddsnp) ~$ conda install pip

.. note::
   pip is not installed by default in new Conda environments, so the last command above is required.

Using venv
^^^^^^^^^^

.. code-block:: bash

   ~$ python -m venv oddsnp
   ~$ source oddsnp/bin/activate
   (oddsnp) ~$

.. note::
   venv installs both Python and pip by default, so no extra commands are required.

Installing oddSNP
-----------------

With your virtual environment activated, navigate to the oddSNP source directory
and install with pip:

.. code-block:: bash

   (oddsnp) ~$ cd path/to/oddsnp/
   (oddsnp) ~/path/to/oddsnp$ pip install .

This will also install all of oddSNP's dependencies into the currently activated environment.

Installing cellsnp-lite
-----------------------

An installation of `cellsnp-lite <https://cellsnp-lite.readthedocs.io/en/latest/index.html>`_
is required to perform pileup calculations within oddSNP. Install it inside your
activated Conda environment with:

.. code-block:: bash

   (oddsnp) ~$ conda install -c bioconda cellsnp-lite

Verifying the installation
--------------------------

To confirm the installation completed successfully, run *oddSNP* from the command
line without any sub-commands. The expected output is:

.. code-block:: bash

   (oddsnp) ~$ oddSNP
   Usage: oddSNP [OPTIONS] COMMAND [ARGS]...

   Options:
     --help  Show this message and exit.

   Commands:
     cpsnpic
     downsample
     genotype
     snpic
     utils
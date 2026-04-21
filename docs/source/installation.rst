Installation
============

The recommended way to install *oddSNP* is by using a virtual environment manager
such as `Conda <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_
or `venv <https://docs.python.org/3/library/venv.html>`_.

Using Bioconda:
---------------

We create a new conda environment and directly install *oddSNP* from its 
bioconda source.

.. code-block:: bash

  :~$ conda create --name oddsnp python=3.12
  :~$ conda activate oddsnp
  (oddsnp):~$ conda install -c bioconda oddsnp

Using PyPI:
-----------

Still, we recommend to install *oddSNP* in a virtual environment. In  this case,
we need to make sure to also install `pip` to the created environment to avoid
interfering with system libraries.

.. code-block:: bash

  :~$ conda create --name oddsnp python=3.12
  :~$ conda activate oddsnp
  (oddsnp):~$ conda install pip
  (oddsnp):~$ pip install oddsnp


From source:
------------

It is also possible to install *oddSNP* directly from source to a previously 
created environment.

For this, we first clone the contents of the GitHub repository to local folder 
(i.e `path/to/oddsnp`), then, in order to install *oddSNP* as a command line 
tool in the current environment, we simply use the following:

.. code-block:: bash
  
  :~$ conda create --name oddsnp python=3.12
  :~$ conda activate oddsnp
  (oddsnp):~$ conda install pip
  (oddsnp):~$ cd path/to/oddsnp/
  (oddsnp):~/path/to/oddsnp$ pip install .


In this case, make sure that the `pip` version used for installation is the one 
associated to the environment and not a system version (i.e. `$ which pip` 
should point to an environment directory and not your system's pip).

Notice that by running this command, you will also install to the currently 
activated environment all of the tool's dependencies.

After installation
------------------

An installation of `cellsnp-lite <https://cellsnp-lite.readthedocs.io/en/latest/index.html>`_
is required to perform pileup calculations within oddSNP. Install it inside your
activated Conda environment with:

.. code-block:: bash

   (oddsnp):~$ conda install -c bioconda cellsnp-lite

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
# Copyright 2025-2026 Nemoto-lab, The University of Osaka
# email: rodolfo.allendes.prime@osaka-u.ac.jp
# email: nemoto.takahiro.prime@osaka-u.ac.jp

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ---

import click
import traceback
import subprocess
import sys
import os
import pysam

from pathlib import Path
import pandas as pd

def assert_cellsnplite():
  """ Asserts that cellsnp-slite has been installed and executable.
  """
  try:
    # Try to run the command to check it exists.
    # stdout=subprocess.PIPE and stderr=subprocess.PIPE prevent output to console.
    # check=True raises CalledProcessError if the command returns a non-zero exit code.
    result = subprocess.run(['cellsnp-lite', '--version'], 
      text=True,
      stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE
    )
    print(f"cellsnp-lite found: {result.stdout}")
  
  # catch exceptions when the cellsnp-lite is not found/not working 
  except subprocess.CalledProcessError as e:
    print("cellsnp-lite failed to execute.")
    print(e.stderr.decode().strip())
    sys.exit(1)
  except Exception as e:
    print("An unexpected error occurred while checking command cellsnp-lite.")
    print(f"{e}")
    sys.exit(1)
    
def assert_vireo():
  """ Asserts that vireo has been installed and executable.
  """
  try:
    # Try to run the command to check it exists.
    # stdout=subprocess.PIPE and stderr=subprocess.PIPE prevent output to console.
    result = subprocess.run(['vireo'],
      text=True,
      stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE
    )
    print(f"Vireo found: {result.stdout.splitlines()[0]}")
  except FileNotFoundError:
    raise AssertionError("vireo not found. Please ensure it installed in your system's PATH.")
  except Exception as e:
    raise AssertionError('An unexpected error occurred while checking command vireo: {}'.format(e))

 
def generate_bam_index(bam, nproc):
    """Generate an index for the input file using pysam.samtools

    Args:
        bam: The BAM file to index
        nproc: Number of parallel processes to use for indexing
    """
    pysam.samtools.index('--threads', str(nproc), bam)

def generate_vcf_index(vcf, nproc):
    """Generate an index for the input VCF file using pysam.bcftools

    Args: 
        vcf: The VCF file to index
        nproc: Number of parallel processes to use for indexing
    """
    pysam.bcftools.index('--threads', str(nproc), vcf)

# Third-party tools that can be used directly from within oddSNP command-line 
# interface, including cell-snp and vireo for pileup and demultiplexation 
# respectively
@click.group()
def utils():
  pass


@utils.command()
@click.argument('bam', type=click.Path(exists=True, dir_okay=False))#: Input BAM file
@click.argument('reference', type=click.Path(exists=True, dir_okay=False))#: List of reference variants
@click.argument('barcodes', type=click.Path(exists=True, dir_okay=False))#: List of cell barcodes
@click.argument('oupath', type=str)
@click.option('--celltag', type=str, default='XC', help='Tag used to reference cell barcodes in the BAM.')
@click.option('--umitag', type=str, default='XM', help='Tag used to reference UMI codes in the BAM.')
@click.option('--mincount', type=int, default=20)
@click.option('--nproc', type=int, default=1, help='Number of processes to use')
@click.option('--force', is_flag=True, help='Override previous results')
def cslite_pileup(bam, reference, barcodes, oupath, celltag, umitag, mincount, nproc, force):
    """ Run cellsnp-lite to generate pileup files
    
    Arguments:
    
        BAM: input bam file
    
        REFERENCE: List of reference variants  
    
        BARCODES: List of cell barcodes
    
        OUPATH: Path to directory to store results.
    """
    # make sure that cellsnp-lite is installed and working
    try:
      assert_cellsnplite()
    except Exception:
      print(traceback.format_exc())
      sys.exit(1)

    # remove possible trailing slashes
    oupath = oupath.rstrip('/') if len(oupath)>1 else oupath

    # run the pileup routine if no results are found or we need to override them
    tgt = Path(f"{oupath}/cellSNP.samples.tsv")
    if tgt.is_file() and not force:
      print(f'Pileup results at {oupath} already exists, skipping...')
      return 
    
    # check that the index for the bam file exits
    idx = Path(f"{bam}.bai")
    if not idx.is_file():
      print(f"Index file for {bam} not found. Generating index.")
      generate_bam_index(bam, nproc)
      
    # create output directory if it does not exist
    os.makedirs(oupath, exist_ok=True)
    
    try:
      proc = subprocess.run([
        'cellsnp-lite',
        '-s', bam,
        '-b', barcodes,
        '-O', oupath,
        '-R', reference,
        '-p', str(nproc),
        '--minMAF', '0.1',
        '--minCOUNT', str(mincount),
        '--cellTAG', celltag,
        '--UMItag', umitag,
        '--gzip'
      ])
    except Exception:
      print(traceback.format_exc())
      sys.exit(1)

    return


def call_cslite_pileup(bam, reference, barcodes, oupath, celltag, umitag, mincount, nproc, force):
    """ Python wrapper for :func:`cslite_pileup`.

    Args:
        bam: The BAM file to process
        reference: List of reference variants
        barcodes: List of cell barcodes
        oupath: Path to directory to store results
        celltag: Tag used to reference cell barcodes in the BAM
        umitag: Tag used to reference UMI codes in the BAM
        mincount: Minimum count threshold for pileup
        nproc: Number of processes to use
        force: If `True` overwrites existing results
    Returns:
        None
    """
    return cslite_pileup.callback(
      bam=bam,
      reference=reference,
      barcodes=barcodes,
      oupath=oupath,
      celltag=celltag,
      umitag=umitag,
      mincount=mincount,
      nproc=nproc,
      force=force
    )

@utils.command()
@click.argument('bam', type=click.Path(exists=True, dir_okay=False))
@click.option('--nproc', type=int, default=1, help='Number of parallel processes to use')
def count_lines(bam, nproc):
    """ Count the number of lines in a BAM file using samtools
    
    Arguments:
        
        BAM: Name of the BAM file
    """
    lines = pysam.samtools.view('-c', '--threads', str(nproc), bam)
    nlines = int(lines.strip())
    return nlines

def call_count_lines(bam, nproc):
    """ Python wrapper for :func:`count_lines`.

    Args:
        bam: The BAM file to process
        nproc: Number of parallel processes to use
    Returns:
        The number of lines in the BAM file as an integer
    """
    return count_lines.callback(
      bam=bam,
      nproc=nproc
    )
  
@utils.command()
@click.argument('inpath', type=str)
@click.argument('oupath', type=str)
@click.option('--genotype', type=click.File('rb'), default=None, help='Path to the donors genotype VCF file.')
@click.option('--genotag', type=str, default='PL', help='Genotype tag to use from the VCF file (GT, GP, or PL). Default is PL.')
@click.option('--ndonor', type=int, default=None, help='Number of donors in the sample. If not provided, vireo will try to estimate it automatically.')
@click.option('--nproc', type=int, default=1, help='Number of parallel processes to use')
@click.option('--force', is_flag=True, help='Override target files.')
def vireo(inpath, oupath, genotype, genotag, ndonor, nproc, force):
    """ Run vireo to demultiplex single-cell data based on genotype information.

    Arguments:

        INPATH: Path to the input file containing the single-cell data to be demultiplexed.
        
        OUPATH: Path to the output directory to store results.
    """
    try:
      utils.assert_vireo()
    except Exception:
      print(traceback.format_exc())
      sys.exit(1)

    # with genotype for all samples
    # vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT_DIR
    if genotype and not ndonor:
      print('Running vireo demultiplexation with genotype...')
      subprocess.run([
        'vireo',
        '-c', inpath,
        '-d', genotype.name,
        '-o', oupath,
        '--genotag', genotag,
        '--nproc', str(nproc)
      ])
    # without any genotype:
    # vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR
    elif ndonor and not genotype:
      print('Running vireo demultiplexation without genotype...')
      subprocess.run([
        'vireo',
        '-c', inpath,
        '-o', oupath,
        '-N', str(ndonor) if ndonor is not None else 'auto',
        '--nproc', str(nproc)
      ])
    # unknown combination of parameters
    else:
      print('Error: Invalid combination of parameters for vireo demultiplexation.')
      print('Please provide either a genotype file (--genotype) or the number of donors (--ndonor), but not both.')
      sys.exit(1)
    
    return

def call_vireo(inpath, oupath, genotype, genotag, ndonor, nproc, force):
    """ Python wrapper for :func:`vireo`

    Args:
        inpath: Path to the input file containing the single-cell data to be demultiplexed.
        oupath: Path to the output directory to store results.
        genotype: Path to the donors genotype VCF file.
        genotag: Genotype tag to use from the VCF file (GT, GP, or PL). Default is PL.
        ndonor: Number of donors in the sample. If not provided, vireo will try to estimate it automatically.
        nproc: Number of parallel processes to use
        force: If `True` overwrites existing results
    """
    return vireo.callback(
      inpath=inpath,
      oupath=oupath,
      genotype=genotype,
      genotag=genotag,
      ndonor=ndonor,
      nproc=nproc,
      force=force
    )

if __name__ == '__main__':
  utils()

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
import pysam, pysam.bcftools

from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np

from . import utils

@click.group()
def downsample():
  pass

@downsample.command()
@click.argument('bam', type=click.Path(exists=True, dir_okay=False))
@click.argument('reads', type=float)
@click.argument('oupath', type=str)
@click.option('--celltag', default='XC', help='Tag used inside the bam file for cell barcodes')
@click.option('--seed', type=int, default=1234, help='A seed to make the selection of reads deterministic')
@click.option('--nproc', type=int, default=1, help='Number of processes to use')
@click.option('--force', is_flag=True, help='Override previous results')
def reads_downsampling(bam, reads, oupath, celltag, seed, nproc, force):
    """ Downsample BAM file by extracting a given percentage of reads.
    It saves both the filtered BAM file and a TSV file with the cell barcodes to 
    the output directory. 

    Arguments:
    
        BAM: The input BAM file
    
        READS: The fraction of reads to keep (between 0 and 1)
        
        OUTDIR: The path to the directory where to store the filtered BAM file
    """
    # remove possible trailing slash
    oupath = oupath.rstrip('/') if len(oupath)>1 else oupath

    # only perform calculation if file doesnt exist or if force flag is set
    tgt = Path('{}/{}_{:.2f}.bam'.format(oupath, Path(bam).stem, reads)) # resulting file
    if tgt.is_file() and not force:
      print('Filtered file {} already exists, skipping...'.format(tgt))
      return tgt

    os.makedirs(oupath, exist_ok=True)
    # get total number of lines in the file
    totallines= utils.call_count_lines(bam=bam, nproc=nproc)
    print('Total number of reads in {}: {:,}'.format(bam, totallines))

    # generate a list of as many true/false values as needed
    keep = round(totallines*reads)
    arr = np.concatenate((
      np.ones(keep, dtype=bool), 
      np.zeros(totallines-keep, dtype=bool)
    ))
    # and shuffle the values to randoly determine the lines to keep
    rng = np.random.default_rng(seed)
    rng.shuffle(arr) 

    # process the file, keeping the predetermined lines and dropping the
    # remaining ones
    inp = pysam.AlignmentFile('{}'.format(bam), 'rb')
    out = pysam.AlignmentFile(
      '{}'.format(tgt),
      'wb', template=inp)
        
    # re-map all reads in the original file
    bcodes = set()
    for i,read in tqdm(enumerate(inp.fetch(until_eof=True)), 
      total=totallines,
      desc='Filtering reads'):
      if arr[i]:
        try:
          out.write(read) # make sure to copy the read first
          bcodes.add(read.get_tag(tag=celltag)) # and then try to fetch the cell bcode
        except KeyError:
          continue
    # close the files after finishing
    inp.close()
    out.close()

    #convert bcodes to DataFrame and save to file
    bcodes_df = pd.DataFrame.from_dict(data={'cells': list(bcodes)})
    bname = '{}/{}_{:.2f}_barcodes.tsv'.format(oupath, Path(bam).stem, reads)
    bcodes_df.to_csv(bname,
      index=False,
      header=False
    )

    return tgt

def call_reads_downsampling(bam, reads, oupath, celltag, seed, nproc, force):
    """ Python wrapper for :func:`reads_downsampling`.

    Args:
        bam: The input BAM file
        reads: The percentage of reads to keep (between 0 and 1)
        oupath: The path to the directory where to store the filtered BAM file
        celltag: The tag used inside the bam file for cell barcodes
        seed: A seed to make the selection of reads deterministic
        nproc: Number of processes to use
        force: Override previous results

    Returns:
        The path to the filtered BAM file
    """
    return reads_downsampling.callback(
      bam=bam,
      reads=reads,
      oupath=oupath,
      celltag=celltag,
      seed=seed,
      nproc=nproc,
      force=force
    )

@downsample.command()
@click.argument('bam', type=click.File('rb'))
@click.argument('barcodes', type=click.File('r'))
@click.argument('oupath', type=str)
@click.option('--celltag', default='XC', help='Tag used inside the bam file for cell barcodes')
@click.option('--nproc', type=int, default=1, help='Number of processes to use')
@click.option('--force', is_flag=True, help='Override previous results')
def barcode_downsampling(bam, barcodes, oupath, celltag, nproc, force):
    """ Downsample BAM file by extracting selected barcodes.
    Given a list of cell barcodes, filter a BAM file in order to only include 
    reads associated to those barcodes.
    
    Arguments:
      
        BAM: Path to the BAM file to filter.
    
        BARCODES: A file with the list of cell barcodes to include.
    
        OUPATH: The path to the directory where to store the filtered BAM file.
    """
    # remove possible trailing slash
    oupath = oupath.rstrip('/') if len(oupath)>1 else oupath

    # only perform filtering if the output file doesnt exist or if force flag is set
    tgt = Path('{}/{}_filtered.bam'.format(oupath, Path(bam.name).stem)) # resulting file
    if tgt.is_file() and not force:
      print('Filtered file {} already exists, skipping...'.format(tgt))
      return tgt
    
    # create output directory if it does not exist
    os.makedirs(oupath, exist_ok=True)

    # process the files in parallel for read extraction
    pysam.samtools.view(
        '-bh', 
        '-D', '{}:{}'.format(celltag, barcodes.name),
        '-o', str(tgt),
        '--threads', str(nproc),
        bam.name,
        catch_stdout=False
      )
    
    return tgt

def call_barcode_downsampling(bam, barcodes, oupath, celltag, nproc, force):
    """ Python wrapper for :func:`barcode_downsampling`.

    Args:
        bam: Path to the BAM file to filter.
        barcodes: A file with the list of cell barcodes to include.
        oupath: The path to the directory where to store the filtered BAM file.
        celltag: The tag used inside the bam file for cell barcodes
        nproc: Number of processes to use
        force: If `True` override previous results
    """
    return barcode_downsampling.callback(
      bam=bam,
      barcodes=barcodes,
      oupath=oupath,
      celltag=celltag,
      nproc=nproc,
      force=force
    )

@downsample.command()
@click.argument('vcf', type=click.Path(exists=True, dir_okay=False))
@click.argument('regions', type=click.Path(exists=True, dir_okay=False))
@click.argument('out', type=str)
@click.option('--nproc', type=int, default=1, help='Number of processes to use')
@click.option('--force', is_flag=True, help='Override previous results')
def vcf_downsampling(vcf, regions, out, nproc, force):
    """ Filter a reference VCF file to only include variants in given regions.

    Arguments:
        VCF: Path to the VCF file to filter.
      
        REGIONS: A file with the list of regions to include.
        
        OUT: The path to the directory where to store the filtered VCF file.
    """
    # only perform calculation if file doesnt exist or if force flag is set
    tgt = Path(out) # resulting file
    if tgt.is_file() and not force:
      print('Filtered file {} already exists, skipping...'.format(tgt))
      return tgt

    # make sure the vcf file has an index
    if not os.path.isfile(f'{vcf}.csi'):
      print(f"Index file for {vcf} not found. Generating index.")
      utils.generate_vcf_index(vcf, nproc)
    
    # create output directory if it does not exist
    os.makedirs(Path(out).parent, exist_ok=True)

    pysam.bcftools.view(
      "-R", regions,
      "-Oz", "-o", out,
      "--threads", str(nproc),
      vcf,
      catch_stdout=False
    )

    return tgt
  
def call_vcf_downsampling(vcf, regions, out, nproc, force):
    """ Python wrapper for :func:`vcf_downsampling`.

    Args:
        vcf: Path to the VCF file to filter.
        regions: A file with the list of regions to include.
        out: The path to the directory where to store the filtered VCF file.
        nproc: Number of processes to use
        force: If `True` override previous results
    """
    return vcf_downsampling.callback(
      vcf=vcf,
      regions=regions,
      out=out,
      nproc=nproc,
      force=force
    )

if __name__ == '__main__':
  downsample()
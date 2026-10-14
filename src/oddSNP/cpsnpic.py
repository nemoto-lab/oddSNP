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
import math
import tqdm

from pathlib import Path
import pandas as pd
import numpy as np
import scipy.sparse as sp

from . import genotype
from . import plotting

@click.group()
def cpsnpic():
  pass

@cpsnpic.command
@click.argument('inpath', type=click.Path(exists=True))
@click.argument('oupath', type=str)
@click.option('--nproc', type=int, default=20, help='Deprecated, ignored.')
@click.option('--batch_size', type=int, default=100000, help='Deprecated, ignored.')
@click.option('--force', is_flag=True, help='Override target file.')
def calculate_cpsnpic(inpath, oupath, nproc, batch_size, force):
  """ Calculate the cpSNP-IC for pairs of cells. Returns the name of the file 
      where the DataFrame with the cpSNP-IC counts per cell pair is stored.
  
  Arguments:

    INPATH: Path to the folder containing pileup related files. The folder should
            contain at least the following files: ``cellSNP.samples.tsv``, 
            ``cellSNP.tag.AD/DP`` and ``OTH`` matrices
    
    OUPATH: Path to the output folder where results will be saved.
  """
  # only regenerate if the file does not exist or if the force flag is set
  tgt = Path('{}/cpsnpic_counts.pkl.gz'.format(oupath))
  if tgt.is_file() and not force:
    print('cpSNP-IC counts found at: {}'.format(tgt))
    return tgt

  # calculate or fetch the aggregated genotype file
  ag = genotype.call_aggregate_pileup(
    inpath=inpath, 
    oupath=oupath,
    force=force
    )
  # cells x SNPs matrix of read counts (DP+OTH); a SNP not covered by a cell is 0
  ag = ag.reset_index()  # columns: SNP, Cell, AD, DP, OTH, bcode
  counts = (ag['DP'] + ag['OTH']).astype(np.int64).values
  ag, counts = ag.loc[counts > 0], counts[counts > 0]
  cells = pd.Index(sorted(ag.bcode.unique()))
  snps = pd.Index(ag.SNP.unique())
  row = cells.get_indexer(ag.bcode)
  col = snps.get_indexer(ag.SNP)
  ncells = len(cells)

  # For integer counts, min(x_a, x_b) is the number of thresholds t = 1, 2, ...
  # reached by both cells. With B_t the cells x SNPs matrix holding a 1 where a
  # cell has >= t reads, min_sum for every pair is sum_t B_t . B_t^T, and SNPs
  # covered by only one of the cells add nothing. B_t only changes at the
  # distinct count values, so each one is weighted by the gap to the previous.
  # The first product holds nearly all pairs; later ones are small, so they are
  # accumulated apart and added once, instead of rebuilding the large matrix.
  print('Calculating cpSNP-IC for {} cells'.format(ncells))
  minsum = sp.csr_matrix((ncells, ncells), dtype=np.int64)
  extra = sp.csr_matrix((ncells, ncells), dtype=np.int64)
  prev = 0
  for i, t in enumerate(tqdm.tqdm(np.unique(counts), desc='Read thresholds')):
    keep = counts >= t
    Bt = sp.csr_matrix(
      (np.ones(keep.sum(), dtype=np.int64), (row[keep], col[keep])),
      shape=(ncells, len(snps))
    )
    # each pair once (cell1 < cell2), no self-pairs
    prod = sp.triu(Bt @ Bt.T, k=1, format='csr') * int(t - prev)
    prev = t
    if i == 0:
      minsum = prod
    else:
      extra = extra + prod
  minsum = (minsum + extra).tocoo()

  # pairs that share at least one SNP (non-zero entries)
  df = pd.DataFrame({
    'cell1': cells[minsum.row],
    'cell2': cells[minsum.col],
    'min_sum': minsum.data,
  })
  # save the resulting file
  print('Saving results to {}'.format(tgt))
  Path(oupath).mkdir(parents=True, exist_ok=True)
  df.to_pickle(tgt, compression='gzip')
  return tgt

def call_calculate_cpsnpic(inpath,oupath,nproc,batch_size,force):
    """Python API wrapper for :func:`calculate_cpsnpic`.

    Args:
        inpath: Path to the folder containing pileup related files. The folder should
            contain at least the following files: ``cellSNP.samples.tsv``,
            ``cellSNP.tag.AD/DP`` and ``OTH`` matrices
        oupath: Path to the output folder where results will be saved.
        nproc: Number of parallel processes to use
        batch_size: Number of cell pairs to process in each batch.
        force: If ``True``, override the target file.
    """
    return calculate_cpsnpic.callback(
      inpath=inpath,
      oupath=oupath,
      nproc=nproc,
      batch_size=batch_size,
      force=force
    )

@cpsnpic.command()
@click.argument('cpfile', type=str)
@click.argument('oupath', type=str)
@click.option('--force', is_flag=True, help='Override target file.')
def generate_histogram(cpfile, oupath, force):
    """ Generate a histogram of cpSNP-IC values from a cpSNP-IC counts file. 
    Returns the name of the file where the histogram DataFrame is stored.
    
    Arguments:

        CPFILE: Path to the cpSNP-IC counts file (pkl.gz) generated by calculate_cpsnpic.
    
        OUPATH: Path to the folder where to store the histogram file.
    """
    # only regenerate if the file does not exist or if the force flag is set
    tgt = Path('{}/cpsnpic_histogram.pkl'.format(oupath))
    if tgt.is_file() and not force:
      print('cpSNP-IC histogram file found at: {}'.format(tgt))
      return tgt
    
    # load the cpSNP-IC counts file
    print('Loading cpSNP-IC counts from: {}'.format(cpfile))
    agdf = pd.read_pickle(cpfile)

    # determine the maximum value to set the histogram bins
    m = math.ceil(agdf.min_sum.max())
    nb = 100 if m>100 else m # fix the maximum number of bins
    step = m//nb
    bin_edges = np.arange(start=0,stop=m+step, step=step) # make sure the bins have int-edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # midpoint of each bin
    
    bin_counts, _ = np.histogram(
      agdf.min_sum.dropna().values,
      bins=bin_edges,
    )
    
    hist = pd.DataFrame({
      'centers': bin_centers,
      'counts': bin_counts
      })
    
    # save the resulting file
    print('Saving cpSNP-IC histogram to: {}'.format(tgt))
    Path(oupath).mkdir(parents=True, exist_ok=True)
    hist.to_pickle(tgt)

    # return the name of the file where the histogram is stored
    return tgt

def call_generate_histogram(cpfile, oupath, force):
    """ Python API wrapper for :func:`generate_histogram`.

    Args:
        cpfile: Path to the cpSNP-IC counts file (pkl.gz)
        oupath: Path to the folder where to store the histogram file
        force: If ``True``, override the target file.
    """
    return generate_histogram.callback(
      cpfile=cpfile,
      oupath=oupath,
      force=force
    )

@cpsnpic.command()
@click.argument('histofile', type=str)
@click.argument('output', type=str)
@click.option('--force', is_flag=True, help='Override target file.')
def save_cpsnpic_plot(histofile, output, force=False):
    """ Save a cpSNP-IC histogram figure

    Arguments:
    
        HISTOFILE The file where histogram data is stored
    
        OUTPUT The output file where to save the figure (with extension .html)
    """
    tgt = Path(output)
    if tgt.is_file() and not force:
      print('Figure file {} already exists, skipping...'.format(tgt))
      return
    
    fig = plotting.generate_cpsnpic_plot(histofile)
    print('Saving cpSNP-IC figure to: {}'.format(tgt))
    plotting.save_plot(
      figure=fig,
      output=output,
      force=force
    )

def call_save_cpsnpic_plot(histofile, output, force):
    """ Python API wrapper for :func:`save_cpsnpic_plot`.

    Args:
        histofile: The file where histogram data is stored
        output: The output file where to save the figure (with extension .html)
        force: If ``True``, override the target file.
    """
    return save_cpsnpic_plot.callback(
      histofile=histofile,
      output=output,
      force=force
    )

@cpsnpic.command()
@click.argument('inpath', type=str)
@click.argument('oupath', type=str)
@click.option('--nproc', type=int, default=20, help='Number of parallel processes to use')
@click.option('--batch_size', type=int, default=100000, help='Number of cell pairs to process in each batch.')
@click.option('--force', is_flag=True, help='Override target files.')
def run_all(inpath, oupath, nproc, batch_size, force):
    """ Run all cpSNP-IC calculation steps.
    
    Arguments:
    
        INPATH: Path to the folder containing pileup related files. The folder should
            contain at least the following files: cellSNP.samples.tsv, 
            cellSNP.tag.AD/DP and OTH matrices
    
        OUPATH: Path to the output folder where results will be saved.
    """
    # step 1: calculate cpSNP-IC counts
    cpfile = calculate_cpsnpic.callback(
      inpath=inpath,
      oupath=oupath,
      nproc=nproc,
      batch_size=batch_size,
      force=force
    )
    # step 2: generate histogram
    histfile = generate_histogram.callback(
      cpfile=cpfile,
      oupath=oupath,
      force=force
    )
    # step 3: plot histogram
    fig = plotting.generate_cpsnpic_plot(histfile)
    plotting.save_plot(
      figure=fig,
      output='{}/cpsnpic_figure.html'.format(oupath),
      force=force
    )

def call_run_all(inpath, oupath, nproc, batch_size, force):
    """ Python API wrapper for :func:`run_all`.

    Args:
        inpath: Path to the folder containing pileup related files.
        oupath: Path to the output folder where results will be saved.
        nproc: Number of parallel processes to use.
        batch_size: Number of cell pairs to process in each batch.
        force: If ``True``, override target files.
    """
    return run_all.callback(
      inpath=inpath,
      oupath=oupath,
      nproc=nproc,
      batch_size=batch_size,
      force=force
    )

if __name__ == '__main__':
  cpsnpic()
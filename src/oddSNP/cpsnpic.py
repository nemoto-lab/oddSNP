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
import pickle, gzip
import itertools, more_itertools
import joblib
import tqdm 

from pathlib import Path
import pandas as pd
import numpy as np

from . import genotype
from . import plotting

@click.group()
def cpsnpic():
  pass

def matchingSNPs(snps1, snps2):
  """ Calculate the sum of min OTH+DP for matching SNPs between two cells.
  
  Args:
    snps1: DataFrame with SNPs for cell 1
    snps2: DataFrame with SNPs for cell 2
  
  Returns:
    A list with [cell1_bcode, cell2_bcode, sum_min_OTH+DP]
  """
  tmp = pd.concat([snps1,snps2], axis=1, keys=[1,2]).dropna()
  if tmp.shape[0]==0: #empty intersection means no matching SNPs
    return None
  
  # flat the column names and calculate the min OTH+DP for each SNP in the intersection
  tmp.columns = [ '{}_{}'.format(b,a) for (a,b) in tmp.columns.to_flat_index() ]
  tmp['sum_1'] = tmp['DP_1']+tmp['OTH_1']
  tmp['sum_2'] = tmp['DP_2']+tmp['OTH_2']
  tmp['min_sum'] = tmp[['sum_1','sum_2']].min(axis=1)

  return [
    snps1.bcode.iat[0],
    snps2.bcode.iat[0],
    tmp.min_sum.sum() # sum of minimum across matching SNPs
  ]

@cpsnpic.command
@click.argument('inpath', type=click.Path(exists=True))
@click.argument('oupath', type=str)
@click.option('--nproc', type=int, default=20, help='Number of parallel processes to use')
@click.option('--batch_size', type=int, default=100000, help='Number of cell pairs to process in each batch.')
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
  # remove cellid
  ag.reset_index('Cell', drop=True, inplace=True)
  
  # calculate or fetch the groups of SNPs by cell barcode 
  # as this is not a DF, we need to save/load it using pickle + gzip
  spc = Path('{}/snps_per_cell.pkl.gz'.format(oupath))
  if not spc.is_file() or force:
    grps = ag.groupby('bcode')
    if not Path(oupath).exists(): Path(oupath).mkdir(parents=True)
    print('Saving SNPs per cell to {}'.format(spc))
    with gzip.open(spc, 'wb') as f:
      pickle.dump(grps, f)
  else:
    print('Loading SNPs per cell from {}'.format(spc))
    with gzip.open(spc, 'rb') as f:
      grps = pickle.load(f)
  
  # generate all cell pairs
  cells = list(grps.groups.keys())
  print('Generating cell pairs for {} cells'.format(len(cells)))
  pairs = itertools.combinations(cells, 2)
  npairs = (len(cells) * (len(cells) - 1)) // 2

  # Perform the acutal matching of snps and counts calculation
  # Given that the cells share SNPs, we can obtain a series of statistics 
  # associated to them, in particular, what is the minimum Allele frequency
  # or the number of times (DP+DP) that the SNP was found for each cell
  full_results = []
  total_batches= npairs // batch_size
  print('total pairs {}, total batches {}'.format(npairs, total_batches))

  for i,batch in enumerate(tqdm.tqdm(
    more_itertools.ichunked(pairs, batch_size), 
    total=total_batches+1, 
    desc="Batch progress")):

      results = joblib.Parallel(
        n_jobs=nproc
      )(
        tqdm.tqdm([
          joblib.delayed(matchingSNPs)(grps.get_group(a), grps.get_group(b)
          ) for a,b in batch
        ],
          total=batch_size if i!=total_batches else npairs % batch_size,
          desc="Parallel jobs"
        )
      )
      results = [r for r in results if r is not None] # filter out None results
      full_results.extend(results)
  
  # Create a DataFrame from the list of results
  df = pd.DataFrame(results, columns=['cell1','cell2','min_sum'])
  # save the resulting file
  print('Saving results to {}'.format(tgt))
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
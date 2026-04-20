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
import sys, traceback
import pandas as pd

import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from pathlib import Path

def generate_snpic_plot(histofile):
    """ Generate a SNP-IC histogram figure and return the generated figure object.
    
    Args:
        histofile: File that contains histogram data, a list of bin 
            centers and counts.
    
    Returns:
        plotly.graph_objects.Figure: The generated figure object. 
    """
    histogram = pd.read_pickle(histofile)
    
    fig = go.Figure(
      layout = dict(
        template='plotly_white',
        barmode='stack',
        height=400,
        width=800,
        margin=dict(l=20,r=25,t=40,b=20),
        uniformtext_minsize=12,
        uniformtext_mode='hide',
        xaxis_title_text='SNP-IC values',
        yaxis_title_text='Number of cells',
      )
    ).update_xaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True,
      rangemode='tozero',
      range=[0,4.5],
      type='log',
    ).update_yaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True,
    )

    # add the bars
    fig.add_trace(go.Bar(
      x=histogram.centers,
      y=histogram.counts,
      marker_color = plotly.colors.qualitative.Set2[2], # manuscript unassigned
      showlegend=False,
    ))

    # threshold line
    fig.add_shape(go.layout.Shape( 
      type='line', 
      line=dict(
        color=plotly.colors.qualitative.Set2[7],
        dash='dot',),
      xref='x',x0=50, x1=50, 
      yref='paper',y0=0, y1=1
    ))
    return fig
  
def generate_cpsnpic_plot(histofile):
    """ Generate a cpSNP-IC histogram figure
    
    Args:
        histofile: File that contains histogram data, a list of bin centers and 
            counts.
    
    Returns:
        plotly.graph_objects.Figure: The generated figure object.
    """
    histogram = pd.read_pickle(histofile)
    
    fig = go.Figure(
      layout = dict(
        template='plotly_white',
        barmode='stack',
        height=400,
        width=800,
        margin=dict(l=20,r=25,t=40,b=20),
        uniformtext_minsize=12,
        uniformtext_mode='hide',
        xaxis_title_text='cpSNP-IC values',
        yaxis_title_text='Number of cells',
      )
    ).update_xaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True,
      rangemode='tozero',
      range=[0,4.5],
      type='log',
    ).update_yaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True,
    )

    # add the bars
    fig.add_trace(go.Bar(
      x=histogram.centers,
      y=histogram.counts,
      marker_color = plotly.colors.qualitative.Antique[1], # manuscript other
      showlegend=False,
    ))

    # threshold line
    fig.add_shape(go.layout.Shape( 
      type='line', 
      line=dict(
        color=plotly.colors.qualitative.Set2[7],
        dash='dot',),
      xref='x',x0=3000, x1=3000, 
      yref='paper',y0=0, y1=1
    ))
    return fig
    
def save_plot(figure, output, force=False):
    """ Save a SNP-IC or cpSNP-IC histogram figure
    
    Args:
        figure: The figure object to save (with extension .html)
        output: The output file where to save the figure
        force: If `True` ovewrites existing files
    """
    out = '{}.html'.format(output) if Path(output).suffix.lower() != '.html' else output
    tgt = Path(out)
    if tgt.is_file() and not force:
      print('Figure file {} already exists, skipping...'.format(tgt))
      return  
    
    try:
      figure.write_html(out)  
      print('Figure saved to: {}'.format(tgt))
      return 
    except Exception:
      print('Image could not be saved')
      print(traceback.format_exc())
      sys.exit(1)  

  
def generate_donor_distribution_plot(barsfile):
    """ Generate a stacked bar plot showing the distribution of correct and 
    incorrect classifications across donors.

    Args:
        barsfile: File that contains the data for the plot, a DataFrame with 
            columns 'donor' and 'correct'.
    Returns:
        plotly.graph_objects.Figure: The generated figure object.
    """    
    bars = pd.read_pickle(barsfile)
    x = bars.donor.unique()
    correct = []
    incorrect = []
    n = len(bars)

    for i in x:
      correct.append(len(bars.loc[(bars.donor==i) & (bars.correct==True)]))
      incorrect.append(len(bars.loc[(bars.donor==i) & (bars.correct==False)]))

    fig = go.Figure(
      layout = dict(
        template='plotly_white',
        barmode='stack',
        height=300,
        margin=dict(l=20,r=25,t=40,b=20),
        xaxis=dict(
          ticks='outside', 
          tickvals=x,
          ticktext=[f"D-{i}" for i in x] ),
      )
    ).update_xaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True
    ).update_yaxes(
      showline=True, linecolor='black', showgrid=False, mirror=True,
      # range=[0,1]
    )


    fig.add_trace(go.Bar(
      x=x,
      y=correct,
      marker_color = plotly.colors.qualitative.Set2[0],
      name='Correct'
    ))

    fig.add_trace(go.Bar(
      x=x,
      y=incorrect,
      marker_color = plotly.colors.qualitative.Set2[2],
      name='Incorrect'
    ))

    return fig
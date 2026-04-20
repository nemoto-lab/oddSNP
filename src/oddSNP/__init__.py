"""oddSNP is an open-source framework for the calculation of SNP-Information 
Content (SNP-IC), a quantitative metric computable from simple unpooled pilot 
data that accurately predicts the success of genotype-based demultiplexing; 
and its equivalent for genotype-free approaches, cell-paired SNP-Information 
Content (cpSNP-IC)."""

__version__ = '0.0.2'
__author__ = 'Rodolfo Allendes'
__email__ = 'rodolfo.allendes.prime@osaka-u.ac.jp'

# make oddSNP available if someone imports
from .__main__ import oddSNP as main

__all__ = ['main']
# Copyright (C) 2025 Nemoto-lab, Osaka University
# email: rodolfo.allendes.prime@osaka-u.ac.jp
# email: nemoto.takahiro.prime@osaka-u.ac.jp
#
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 2 of the License, or (at your option) any later 
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
# Place, Suite 330, Boston, MA 02111-1307 USA
import click

from . import utils
from . import downsampling
from . import genotype
from . import snpic
from . import cpsnpic

@click.group()
@click.version_option(version='0.1.0', prog_name='oddSNP')
def oddSNP():
  pass

oddSNP.add_command(utils.utils)
oddSNP.add_command(downsampling.downsample)
oddSNP.add_command(genotype.genotype)
oddSNP.add_command(snpic.snpic)
oddSNP.add_command(cpsnpic.cpsnpic)
################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2018
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2018
#
#########################################################################
#
import os
import sys
from argparse import ArgumentParser

import ADFRcc
from ADFR import ADFRRigidReceptor
from ADFR.GA import GA
from ADFR.utils.maps import MapsFile
from MolKit2 import Read
from MolKit2.molecule import getAtomIndicesPerType
from mglutil.math.rmsd import HungarianMatchingRMSD_prody


parser = ArgumentParser(description="Randomize a ligand and minimize its internal energy to avoid clashes")
parser.add_argument('-l', "--ligand", dest="ligandFile",
                    metavar=('ligand.pdbqt',),
                    required=True,
                    help="ligand to be docked (in PDBQT format)",)
parser.add_argument('-t', "--target", dest="receptorMapsFile",
                    metavar=('targetFile.trg',),
                    required=True,
                    help="file produced by AGFR and describing the receptor (see http://agfr.scripps.edu)",)
parser.add_argument('-o', "--output", dest="outputFilename",
                    metavar=('randomLigand.pdbqt',),
                    required=False,
                    help="filename for randomized ligand",)
parser.add_argument('-s', "--silent", dest="silent",
                    action="store_true", default=False,
                    help="specify to supress printing to stdout",)

advancedGroup = parser.add_argument_group(
    'ADVANCED OPTIONS:',
    'the following options are less commonly used')
# BFGS specific
advancedGroup.add_argument(
    "--bfgs_maxSteps", dest="bfgs_max_steps",
    metavar=('bfgs_maxSteps',), default=1000, type=int,
    help="maximum number of steps per BFGS optimiations. Use -1 for unlimited. Default 100",)

advancedGroup.add_argument(
    "--bfgs_ftol_abs", dest="bfgs_ftol_abs",
    metavar=('bfgs_ftol_abs',), default=-1.0, type=float,
    help="stop BFGS optimization if score changes by less than this value. Use -1 to disable",)

advancedGroup.add_argument(
    "--bfgs_ftol_rel", dest="bfgs_ftol_rel",
    metavar=('bfgs_ftol_rel',), default=-1.0, type=float,
    help="stop BFGS optimization if score changed by less than this value*abs(score). Use -1 to disable",)

advancedGroup.add_argument(
    "--bfgs_xtol_abs", dest="bfgs_xtol_abs",
    metavar=('bfgs_xtol_abs',), default=-1.0, type=float,
    help="stop BFGS optimization if every variable changes by less than this value. Use -1 to disable",)

advancedGroup.add_argument(
    "--bfgs_xtol_rel", dest="bfgs_xtol_rel",
    metavar=('bfgs_xtol_rel',), default=-1.0, type=float,
    help="stop BFGS optimization if every variable changes by less than this value*abs(variable). Use -1 to disable",)

options = parser.parse_args()

if options.silent:
    sys.stdout = open(os.devnull, 'w')

ligand = Read(options.ligandFile)
maps = MapsFile(options.receptorMapsFile)
maps.unzipMaps()

adfr = ADFRRigidReceptor(ligand, maps.folder)

# disable lRR scoring to noly optimize internal energy during minimization
adfr.scorer.scoreLRR(False)

ind = adfr.createPopulation(1)[0]
s0 = ind.score()

ga = GA([ind], ligand, localSearchMode='BFGS')

adfr.scorer.resetNumEvals()
miniOpt = { 'aMaxSteps':options.bfgs_max_steps,
            'xtol_abs':options.bfgs_xtol_abs,
            'xtol_rel':options.bfgs_xtol_rel,
            'ftol_abs':options.bfgs_ftol_abs,
            'ftol_rel':options.bfgs_ftol_rel , 'trajFilename':None}
sc, nb = ga.minimize(ind, **miniOpt)

print 'energy: %.2f evals: %d'%(-sc, adfr.scorer.getNumEvals())

if options.outputFilename is None:
    outputName = os.path.splitext(options.ligandFile)[0]+'_randomzied'
else:
    outputName = os.path.splitext(options.outputFilename)[0]

adfr.writeSolution(solution=ind,
                   torTree=adfr.ligandFT.torTree,
                   filename=outputName+'.pdbqt',
                   solNum=0, rmsdCalculators=[],
                   status='ligand randomized by ADFR')

if options.silent:
    sys.stdout = sys.__stdout__

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
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2015
#
# Updated in 2018 by Antoine P. SANNER
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/minimize.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: minimize.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
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


parser = ArgumentParser(description="Minimizes a ligand using either the BFGS or Solis Wets and stores the result.")
parser.add_argument('-l', "--ligand", dest="ligandFile",
                    metavar=('ligand.pdbqt',),
                    required=True,
                    help="ligand to be docked (in PDBQT format)",)
parser.add_argument('-t', "--target", dest="receptorMapsFile",
                    metavar=('targetFile.trg',),
                    required=True,
                    help="file produced by AGFR and describing the receptor (see http://agfr.scripps.edu)",)
parser.add_argument('-r', "--reference", dest="reference",
                    metavar=('referenceLigForRMSD.pdbqt',),
                    help="filename of reference ligand used for RMSD reporting",)
parser.add_argument('-L', "--localSearch", dest="localSearchMode", type=str,
                    metavar='local search mode can be "SolisWets" or "BFGS"',
                    choices=["SolisWets", "BFGS", ], default="BFGS",
                    help='select local search method, default is "BFGS"')
parser.add_argument('-d', "--debug", dest="debug",
                    action="store_true", default=False,
                    help='If True allows the printing of information during the file reading, default is False')
parser.add_argument('-o', "--logFile", dest="logfile",
                    metavar=("logFilename",), default=None,
                    help="summary docking. A .dlg extension will be added if missing")
parser.add_argument('-s', "--score", dest="score",
                    action="store_true", default=False,
                    help='If True this script doesn\'t minimize the individual and only scores it.')

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

if not options.debug:
    sys.stdout = open(os.devnull, 'w')

ligand = Read(options.ligandFile)
maps = MapsFile(options.receptorMapsFile)
maps.unzipMaps()
adfr = ADFRRigidReceptor(ligand, maps.folder)
ind = adfr.createPopulation(1)[0]
ind.setGenes(ind.genomePy.getIdentityGenesPy())

# ind.genes= ind.genome.identityGenes()
s0 = ind.score()

atomsIn = adfr.ligand.select()
d1 = getAtomIndicesPerType(atomsIn)
# RMSD calculator to measure ligand changes upon minimization
rmsdCalcIn = HungarianMatchingRMSD_prody(atomsIn.getCoords(), d1, d1)

if options.reference:
    refLig = Read(options.reference)
    atomsRef = refLig.select()
    d2 = getAtomIndicesPerType(atomsRef)
    rmsdCalcRef = HungarianMatchingRMSD_prody(atomsRef.getCoords(), d2, d1)
else:
    rmsdCalcRef = None

ga = GA([ind], ligand, localSearchMode=options.localSearchMode)

if options.localSearchMode=='SolisWets':
    miniOpt = { 'nbSteps': 10, 'noImproveStop': 3, 'max_steps': 300,
                'MAX_FAIL': 4, 'searchRate': 0.05, 'trajFilename':None}
elif options.localSearchMode=='BFGS':
    miniOpt = { 'aMaxSteps':options.bfgs_max_steps,
                'xtol_abs':options.bfgs_xtol_abs,
                'xtol_rel':options.bfgs_xtol_rel,
                'ftol_abs':options.bfgs_ftol_abs,
                'ftol_rel':options.bfgs_ftol_rel , 'trajFilename':None}

if options.score:
    sc, nb = ind.score(), 1
    evals = 1
else:
    adfr.scorer.resetNumEvals()
    print miniOpt
    sc, nb = ga.minimize(ind, **miniOpt)
    evals = adfr.scorer.getNumEvals()
    ## scores = []
    ## name = os.path.split(options.ligandFile)[0].split(os.sep)[-1]
    ## import numpy, pickle
    ## _scoresStats = []
    ## _evalsStats = []
    ## if options.localSearchMode=='SolisWets':
    ##     for i in range(1, 200):
    ##         _scores = []
    ##         _evals = []
    ##         miniOpt['aMaxEvals'] = i
    ##         for i in range(100):
    ##             indc = ind.clone()
    ##             adfr.scorer.resetNumEvals()
    ##             sc, nb = ga.minimize(indc, **miniOpt)
    ##             _evals.append(adfr.scorer.getNumEvals())
    ##             _scores.append(sc)
    ##             indc.score()
    ##         #import pdb; pdb.set_trace()
    ##         score = numpy.average(_scores)
    ##         scoredev = numpy.std(_scores)
    ##         scoremin = min(_scores)
    ##         scoremax = max(_scores)

    ##         evals = numpy.average(_evals)
    ##         evalsdev = numpy.std(_evals)
    ##         evalsmin = min(_evals)
    ##         evalsmax = max(_evals)
    ##         _scoresStats.append([score,scoredev,scoremin,scoremax])
    ##         _evalsStats.append([evals,evalsdev,evalsmin,evalsmax])
    ##     f = open('%s_scoresSW.pkl'%name, 'w')
    ##     pickle.dump([_scoresStats,_evalsStats], f)
    ##     f.close()
    ##     print ' file: %s evals: %5d dev %5d min %5d max %5d | startE: %f endE: avg %.3f std %.3f min %.3f max %.3f ediff: %f rmsdI: %.2f rmsdR: %.2f'%(
    ## os.path.split(options.ligandFile)[1], evals, evalsdev, evalsmin, evalsmax, s0, score, scoredev, scoremin, scoremax, score-s0, rmsdi, rmsdr),
    ## else:
    ##     for i in range(1, 200):
    ##         indc = ind.clone()
    ##         miniOpt['aMaxSteps'] = i
    ##         adfr.scorer.resetNumEvals()
    ##         sc, nb = ga.minimize(indc, **miniOpt)
    ##         _evalsStats.append(adfr.scorer.getNumEvals())
    ##         _scoresStats.append(sc)
    ##     #import pdb; pdb.set_trace()
    ##     name = os.path.split(options.ligandFile)[0].split(os.sep)[-1]
    ##     f = open('%s_scoresBFGS.pkl'%name, 'w')
    ##     pickle.dump((_scoresStats, _evalsStats),f)
    ##     f.close()
    
rmsdi = rmsdCalcIn.computeRMSD(ind.phenotype[1])
if rmsdCalcRef:
    rmsdr = rmsdCalcRef.computeRMSD(ind.phenotype[1])
else:
    rmsdr = -1.

if not options.debug:
    sys.stdout = sys.__stdout__

print ' file: %s evals: %5d startE: %f endE: %.3f ediff: %f rmsdI: %.2f rmsdR: %.2f, LL: %9.2f RRL: %9.2f FRL: %9.2f FRFR: %9.2f RRFR: %9.2f wRR: %4.2f FEB: %9.2f'%(
    os.path.split(options.ligandFile)[1], evals, -s0, -sc, sc-s0, rmsdi, rmsdr, ind.energies['LL'], ind.energies['RRL'],
        ind.energies['FRL'], ind.energies['FRFR'], ind.energies['RRFR'], ind.energies['wRR'],
        ind.energies['RRL'] + ind.energies['FRL'] + adfr.ligandFT.torTree.torsdof*ADFRcc._parameters.feCoeffTors),

path, file_name = os.path.split(os.path.abspath(options.ligandFile))
new_file_name = file_name.split('.')[0] + '_minimized.' + file_name.split('.')[1]
path = options.logfile if vars(options)['logfile'] is not None else path

if not options.score:
    adfr.writeSolution(solution=ind,
                       torTree=adfr.ligandFT.torTree,
                       filename=path + '/' + new_file_name,
                       solNum=0,
                       rmsdCalculators=[rmsdCalcRef] if rmsdCalcRef is not None else [],
                       status='Locally minimized using ' + options.localSearchMode)

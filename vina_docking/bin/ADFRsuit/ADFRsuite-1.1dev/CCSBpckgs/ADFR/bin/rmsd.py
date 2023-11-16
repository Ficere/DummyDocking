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
# Copyright: M. Sanner and TSRI 2018
#
#########################################################################
#

import os
import sys
from argparse import ArgumentParser
from MolKit2 import Read
from MolKit2.molecule import getAtomIndicesPerType
from mglutil.math.rmsd import RMSDCalculator, HungarianMatchingRMSD_prody, MorpicRMSD

parser = ArgumentParser(description="Compute rmsd between 2 conformations of the same molecule. Atom matching can be 1to1, Hngarian or automophism.")
parser.add_argument('-r', "--reference", dest="ref",
                    metavar=('refLig.pdb',),
                    required=True,
                    help="reference ligand",)
parser.add_argument('-m', "--moving", dest="mov",
                    metavar=('movLig.pdb',),
                    required=True,
                    help="moving ligand with atoms in same order as refLig",)
parser.add_argument('-s', "--atoms", dest="selstr",
                    metavar=('"name N CA C and bb or name CB"',),
                    required=False, default=None,
                    help="ProDy selection string. Default is None i.e. using all atoms",)
parser.add_argument('-o', "--1to1", dest="oneToOne",
                    action="store_false", #default="False",
                    help="do not report 1to1 RMSD",)
parser.add_argument('-a', "--autoMorphic", dest="autoMorph",
                    action="store_true", default=False,
                    help="Compute automorphisms and report lowest RMSD",)
parser.add_argument('-H', "--hungarianMatching", dest="hungarian",
                    action="store_true", default=False,
                    help="Compute RMSD using Hungarian matching",)
parser.add_argument('-M', "--models", dest="models",
                    type=str, default='first',
                    help='specify which models to use (if more than one). Use "all" or space separated list of integers',)

options = parser.parse_args()

ref = Read(options.ref)
mov = Read(options.mov)

if options.selstr is None:
    atref = ref._ag
    atmov = mov._ag
    selstr = 'all atoms'
else:
    atref = ref._ag.select(options.selstr)
    atmov = mov._ag.select(options.selstr)
    selstr = options.selstr

assert len(atref)==len(atmov) and len(atref)>0

if options.models=='first':
    modelNums = [0]
else:
    nbModels = mov._ag.numCoordsets()
    if options.models.lower()=='all':
        modelNums = range(nbModels)
    else:
        modelNums = [int(x)-1 for x in options.models.split()]

if options.autoMorph:
    autos = ref.getAutomorphisms()
    if options.selstr is not None:
        ll = []
        atRefInd = set(atref.getIndices())
        atMovInd = set(atmov.getIndices())
        for auto in autos:
            l = []
            for i,j in auto:
                if i in atRefInd and j in atMovInd:
                    l.append((i,j))
            if len(l)==len(atRefInd):
                ll.append(l)
            #else:
            #    print len(auto), len(atRefInd)
                # sometimes we get weird automorphisms e.g. 1BC5_FH_ligand.pdbqt
        autos = ll

#import pdb; pdb.set_trace()
for modelNum in modelNums:
    mov._ag.setACSIndex(modelNum)
    atmov.setACSIndex(modelNum)
    print '%3d %s:"%s" ModelRMSD: '%(modelNum+1, options.mov, selstr),

    if options.oneToOne:
        rmsdCalc1 = RMSDCalculator(atref.getCoords())
        rmsd1 = rmsdCalc1.computeRMSD(atmov.getCoords())
        print '1to1: %5.2f |'%rmsd1,

    if options.hungarian:
        d1 = getAtomIndicesPerType(atref)
        d2 = getAtomIndicesPerType(atmov)
        rmsdCalcH = HungarianMatchingRMSD_prody(atref.getCoords(), d1, d2)
        rmsdH = rmsdCalcH.computeRMSD(atmov.getCoords())
        print 'hung: %5.2f |'%rmsdH,

    if options.autoMorph:
        rmsdMCalc = MorpicRMSD(autos, ref._ag.getCoords())
        rmsdM = rmsdMCalc.computeRMSD(mov._ag.getCoords())
        print 'autoMorph: %5.2f'%rmsdM,

    print

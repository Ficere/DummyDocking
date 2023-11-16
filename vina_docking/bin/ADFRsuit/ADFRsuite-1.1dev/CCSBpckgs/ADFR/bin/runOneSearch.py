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
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/runOneGA.py,v 1.20.2.1 2017/08/21 22:23:11 annao Exp $
#
# $Id: runOneGA.py,v 1.20.2.1 2017/08/21 22:23:11 annao Exp $
#

import os, sys, traceback

## sys.path.append('/mgl/ms1/people/sanner/python/mglTools2/MGLTools2-1.1/MGLToolsPckgs')
## sys.path.append('/mgl/ms1/people/sanner/python/mglTools2/MGLTools2-1.1/lib/python2.7/site-packages')

from ADFR.utils.optParser import ArgParser
parser = ArgParser('OneGA')
args = parser.parse_args()

from ADFR.MCsearch2 import MCsearch
# use 4.2 by default
#from ADFRcc import setForceFieldVersion
#parameters = setForceFieldVersion('4.2')
from ADFRcc.adfrcc import Parameters
parameters = Parameters.getParameters()
customParams = args['customParams']
if customParams:
    if os.path.exists(customParams):
        parameters.loadFromDatFile(customParams)
    else:
        raise ValueError('runOneGA: custome parameter file %s not found'%customParams)
    
if len(args)<1:
    parser.print_help()
    sys.exit(1)
    
from ADFR import ADFR
import os
from mglutil.math.rmsd import RMSDCalculator, MorpicRMSD

seedValue = args['seedValue']
searchEngine = args['search']

from time import time
from ADFR import getLigandFromFile
ligand, error = getLigandFromFile(args['ligandFile'])

## check ligand for automorphisms
##
autos = ligand.getAutomorphisms()
assert len(autos)>=1
rmsdMCalc = MorpicRMSD(autos)

if error:
    print error
    sys.exit(1)

mapsFolder = None
reference = args['reference']

neighborSearchCutoff = args['neighborSearchCutoff']
neighborSearchGroup = args['neighborSearchGroup']
if neighborSearchGroup is None:
    neighborSearchGroup = "all"
    
flexRes = args['flexRes']
fixedRoot = args['fixedRoot']
if flexRes is not None:
    from ADFR.utils.maps import flexResStr2flexRes
    flexRes = flexResStr2flexRes(flexRes)
    
covalentLigand = args.get('covalentLigand', None)
if args['receptorMapsFile'] is not None:
    from ADFR.utils.maps import MapsFile, flexResStr2flexRes
    mf = MapsFile(args['receptorMapsFile'])
    mf.unzipMaps()
    mapsFolder = mf.getMapsFolder()
    receptorFilename = os.path.join(mf.getMapsFolder(),
                                    mf.getReceptorFilename())
    flexResStr = mf.getFlexResStr()
    covalentRec = mf.getCovalentBond()
    #mapsFolder, receptorFilename, flexRes, msg, amap = unzipMaps(
    #    args['receptorMapsFile'])
    #flexRes = flexResStr2flexRes(flexRes)
    #f = open(os.path.join(mapsFolder, 'data.pkl'))
    #gridData = pickle.load(f)
    #f.close()
    #covalentRec = d.get('covalentBond', None)
    if covalentRec is not None:
        covalentRec += [d.get('covalentBondTorsionAtoms')[-1]]
    else:
        covalentRec = None    
    if mapsFolder is None:
        print msg
        sys.exit(1)
    tPtsFile = os.path.join(mapsFolder, "translationPoints.npy")
    if os.path.exists(tPtsFile):
        tpointsFilename = os.path.join(mapsFolder, "translationPoints.npy")
    else:
        tpointsFilename = None
    mapFilesRoot = 'rigidReceptor'

else:
    mapsFolder = args['mapFilesFolder']
    receptorFilename = args['receptorFile']
    mapFilesRoot = args['mapFilesRoot']
    tpointsFilename = args.get('tpointsFilename', None)
    cr = args.get('covalentReceptor', None)
    covalentRec = args.get('covalentRec', None)

if covalentLigand is not None:
    if covalentRec is None:
        raise RuntimeError("ERROR: covalentLigand specified but no covalentRec")
elif covalentRec is not None:
    if covalentLigand is None:
        raise RuntimeError("ERROR: covalentRec specified but covalentLigand")

if covalentLigand is not None and covalentRec is not None:
    covalentIndices = covalentLigand + covalentRec
else:
    covalentIndices = None

if receptorFilename is not None:
    from MolKit2 import Read
    receptor = Read(receptorFilename)
else:
    receptor = None

useRama = not args['doNotUseRama']
print 'USE RAMA', useRama
print 'search engine', searchEngine
#from ADFRcc import adfrcc
#import pdb; pdb.set_trace()
#ADFRcc.adfrcc.Logger.setLogLevelInfo()
#print parameters.printDebugDescription()
#ADFRcc.adfrcc.Logger.setLogLevelInfo()
t00= time()
try:
    adfr = ADFR(ligand, mapsFolder, logFile=args['logfile'],
                receptor=receptor,
                flexibleResidues=flexRes,
                mapFilesRoot=mapFilesRoot,
                tpointsFilename=tpointsFilename,
                seedValue=seedValue,
                covalentIndices=covalentIndices,
                FTRecsrc=args.get('FTRecsrc', None),
                fixedRoot=fixedRoot,
                localSearchMode = args['localSearchMode'],
                neighborSearchCutoff=neighborSearchCutoff,
                useRama = useRama)

except Exception, e:
    print traceback.format_exc()
    sys.exit(1)
#           torsionInfoFilename=args['torsInfo)

clusterEcutoff = args['eCutoffForSolutions']

import ADFRcc
from ADFRcc import getFFParameters
parameters = getFFParameters()
print "ForceField cutoffs: vdw %4.2f estat %4.2f hbond %4.2f desolv %4.2f\n"%(
    parameters.maxDistVdw, parameters.maxDistEstat,
    parameters.maxDistHbond, parameters.maxDistSolv)

ligTorsions = adfr.ligandFT.torTree.torsdof
recTorsions = len(adfr.FT.allMotions) - 2 - ligTorsions
if searchEngine=='LGA':
    if args['maxEvals']=='auto':
        args['maxEvals']=2500000
    if args['popSize'] is 'auto':
        nLigGenes = 7+adfr.ligandFT.torTree.torsdof
        popSize = 50 + nLigGenes*10
    else:
        popSize = args['popSize']
    maxEvals = int(args['maxEvals'])
    
elif searchEngine=='moca':
    maxTry = 1500 + 100*ligTorsions
    popSize = 1
    if args['maxEvals']=='auto':
        maxEvals =  200000 + 60000*ligTorsions + 20000*recTorsions
    else:
        maxEvals = int(args['maxEvals'])

print 'receptor', receptorFilename
print 'flexRes', flexRes
print 'mapFilesRoot', mapFilesRoot
print 'tpointsFilename', tpointsFilename
print 'covalentIndices', covalentIndices
print 'RMSDMatching', args['RMSDMatching']
print 'MAXGENS', args['maxGens']
print 'ligTorsions', ligTorsions
print 'recTorsions', recTorsions
print 'MAXEVALS', maxEvals
print 'NOIMPROVESTOP', args['noImproveStop']

if neighborSearchCutoff>0:
    adfr.scorer.setNeighborRMSDcutoff(neighborSearchCutoff,
                         ligand._ag.select(neighborSearchGroup).getIndices())
    pop = adfr.createNeighborPopulation(popSize, neighborSearchCutoff)
else:
    pop = adfr.createPopulation(popSize, neighborSearchCutoff)

#adfr.savePopulation(pop, 'test')
#name = os.path.splitext(os.path.basename(args['ligandFile']))[0]
#adfr.savePopulationGenes(pop, 'astex/ligands/%s_genes.pop'%name[:4])
#sys.exit(0)
#import pdb; pdb.set_trace()

#ind = pop[0]
#if args['swBiasCoeff1:
#    sw.biasCoef1 = args['swBiasCoeff1
#if args['swBiasCoeff2:
#    sw.biasCoef2 = args['swBiasCoeff2

if searchEngine=='LGA':
    adfr.createGA(pop, reference, RMSDMatching=args['RMSDMatching'],
                  neighborSearchCutoff=neighborSearchCutoff,
                  clusterEnergyCut=clusterEcutoff,
                  clusteringRMSDCutoff=args['GAclusteringRMSDCutoff'])

    adfr.ga.setDebug(args['debug'])
    stat = adfr.ga.evolve(maxGens=args['maxGens'],
                          maxEvals=maxEvals,
                          noImproveStop=args['noImproveStop'])
    # cluster top solutions using automorphisms if any
    cut = pop[0]._score - clusterEcutoff
    nbMorphs = len(rmsdMCalc.morphisms)
    toCluster = [i for i,x in enumerate(pop) if x._score > cut]
    clusters = adfr.ga.cluster(toCluster, RMSDcalc=rmsdMCalc)
    print "clustered %d individuals with cutoff %.2f and %d morphisms into %d clusters"%(
        len(toCluster), clusterEcutoff, nbMorphs, len(clusters))

else:
    if args['reference']:
        refLig = Read(args['reference'])
    else:
        refLig = None
    adfr.referenceLigand = refLig
    pose, status = MCsearch(pop[0], maxEvals=maxEvals, maxTry= maxTry, debug=args['debug'])
    
# reset neighborSearchCutoff so that solutions can be re-scored to NOT include
# RMSD distance penalty
rescoreNeeded = False
if neighborSearchCutoff>0:
    adfr.scorer.setNeighborRMSDcutoff(-1.0,ligand._ag.select(neighborSearchGroup).getIndices())
    rescoreNeeded = True
    ## for finalpop in pop:
    ##     finalpop.genomePy.neighborSearchCutoff=-1.0
    ##     #sscore=finalpop._score
    ##     finalpop.score()
    ##     #finalpop.genomePy.neighborSearchCutoff=neighborSearchCutoff
    ##     if finalpop._neighborRMSD<neighborSearchCutoff:
    ##         outputpop.append(finalpop)
    ## if len(outputpop)>0:
    ##     #pop=outputpop
    ##     print "Rescore populations with rmsd larger than RMSD cutoff left with %3d populations"%(len(outputpop))
    ## else:
    ##     print "Nothing left"
    ## pop.sort()

# build Hungarian matching RMSD calculator regarless of matching method used
# in GA
if adfr.referenceLigand:
    refAtoms = adfr.referenceLigand.select('not hydrogen')
    refAtInd = refAtoms.getIndices().tolist()
    #refAtoms = adfr.referenceLigand.select()
    #from MolKit2.molecule import getAtomIndicesPerType
    #from mglutil.math.rmsd import HungarianMatchingRMSD_prody
    #d1 = getAtomIndicesPerType(refAtoms)
    #d2 = getAtomIndicesPerType(adfr.ligand.select())
    #d2 = getAtomIndicesPerType(adfr.ligand.select('not hydrogen'))
    #rmsdHCalc = HungarianMatchingRMSD_prody(refAtoms.getCoords(), d1, d2)

    rmsd1Calc = RMSDCalculator(refAtoms.getCoords())
    rmsdMCalc.setRefCoords(adfr.referenceLigand.select().getCoords())

else:
    rmsd1 = -1
    rmsdM = -1

if args.get('logfile', None):
    import os
    folder = os.path.split(args['logfile'])[0]
else:
    folder = '.'

if args.get('jobNumber',None) is None:
    jobNum = adfr.initialSeed
else:
    jobNum = args['jobNumber']

torTree = adfr.ligandFT.torTree

#adfr.savePopulation(pop, adfr.ligandFT.torTree, 'pop')

if searchEngine=='LGA':
    solTosave = args['nbPosesToSave']
    if solTosave==-1:
        nbSol = len(clusters)
    else:
        nbSol = min(len(clusters), solTosave)

    for i in range(nbSol):
        sol = pop[clusters[i][0]]
        if rescoreNeeded:
            sol.score()

        if adfr.referenceLigand:
            #rmsdH = rmsdHCalc.computeRMSD(sol.phenotype[1][refAtInd])
            rmsd1 = rmsd1Calc.computeRMSD(sol.phenotype[1][refAtInd])
            #rmsdH = rmsdHCalc.computeRMSD(sol.phenotype[1])
            #rmsd1 = rmsd1Calc.computeRMSD(sol.phenotype[1])
            rmsdM = rmsdMCalc.computeRMSD(sol.phenotype[1])

        #import pdb; pdb.set_trace()
        print 'END%d bestScore: %.3f RMSD: %5.1f %5.1f Ene: FLL: %9.2f RRL: %9.2f FRL: %9.2f FRFR: %9.2f RRFR: %9.2f wRR: %4.2f FEB: %9.2f'%(
            i+1, -sol._score, rmsd1, rmsdM, sol.energies['LL'], sol.energies['RRL'],
            sol.energies['FRL'], sol.energies['FRFR'], sol.energies['RRFR'], sol.energies['wRR'],
            sol.energies['RRL'] + sol.energies['FRL'] + adfr.ligandFT.torTree.torsdof*ADFRcc._parameters.feCoeffTors)

        if adfr.referenceLigand:
            adfr.writeSolution(
                sol, torTree, os.path.join(folder, '%s%04d_%02d'%(args['jobName'],jobNum,i+1)),
                jobNum, rmsdCalculators=[rmsdMCalc], status=stat, seed=adfr.initialSeed)
        else:
            adfr.writeSolution(
                sol,  torTree, os.path.join(folder, '%s%04d_%02d'%(args['jobName'],jobNum,i+1)),
                jobNum, status=stat, seed=adfr.initialSeed)

    print 'evolution stopped after %d evals in %.2f seconds'%(
        adfr.scorer.getNumEvals(), time()-t00)
    print '  status %s'%stat
    print '  %d clusters'%len(clusters)
    print '  seed %d'%adfr.initialSeed
    print 'GA %d terminated successfully with %d solutions'%(jobNum, nbSol)
else:
    from mglutil.math.rmsd import MorpicRMSD
    _autos = adfr.ligand.getAutomorphisms()
    assert len(_autos)>0
    if adfr.referenceLigand:
        rmsdCalcMorph = MorpicRMSD(_autos)
        rmsdCalcMorph.setRefCoords(adfr.referenceLigand._ag.getCoords())
        rmsd = rmsdCalcMorph.computeRMSD(pose.phenotype[1])
    else:
        rmsd = -1
    print "maxevals %d maxtry %d"%(maxEvals, maxTry)
    print "END bestScore: %9.5f rmsda %5.1f %5.1f Ene: FLL: %9.2f RRL: %9.2f FRL: %9.2f FRFR: %9.2f RRFR: %9.2f wRR: %4.2f FEB: %9.2f"%(
        -pose._score, rmsd, rmsd, pose.energies['LL'], pose.energies['RRL'],
        pose.energies['FRL'], pose.energies['FRFR'], pose.energies['RRFR'], pose.energies['wRR'],
        pose.energies['RRL'] + pose.energies['FRL'] + adfr.ligandFT.torTree.torsdof*ADFRcc._parameters.feCoeffTors)
    print 'evolution stopped after %d evals in %.2f seconds'%(adfr.scorer.getNumEvals(), time()-t00)
    print '  status %s'%status
    print '  1 cluster'
    print '  seed %d'%adfr.initialSeed
    print 'GA %d terminated successfully with 1 solution'%(jobNum,)
    
    if adfr.referenceLigand:
        adfr.writeSolution(
            pose, torTree, os.path.join(folder, '%s%04d_%02d'%(args['jobName'],jobNum,1)),
            jobNum, rmsdCalculators=[rmsdMCalc], status=status, seed=adfr.initialSeed)
    else:
        adfr.writeSolution(
            pose,  torTree, os.path.join(folder, '%s%04d_%02d'%(args['jobName'],jobNum,1)),
            jobNum, status=status, seed=adfr.initialSeed)

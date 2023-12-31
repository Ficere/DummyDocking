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
#
#
import sys
#sys.path.append('/mgl/ms1/people/sanner/python/mglTools2/MGLTools2-1.1/MGLToolsPckgs')
#sys.path.append('/mgl/ms1/people/sanner/python/mglTools2/MGLTools2-1.1/lib/python2.7/site-packages')

import os, sys, ADFR, numpy, platform, datetime, tempfile, shutil, random, \
       tarfile, pickle
from glob import glob
import ADFRcc
import openbabel

from MolKit2 import Read
from ADFR import getLigandFromFile, getTORSDOF
from ADFR.utils.cluster import clusterPoses, oneCluster
#from mglutil.math.rmsd import HungarianMatchingRMSD_prody
from mglutil.math.rmsd import MorpicRMSD
from MolKit2.molecule import getAtomIndicesPerType

def _saveBestInClusterAsModels(folder, clusters, jobNums, solNums, pdbqtBaseName, pdbqtFormat):
    # folder seems to always end with pdbqtBaseName do the name is repeated
    #f = open('%s_%s_out.pdbqt'%(folder, pdbqtBaseName), 'w')
    f = open('%s_out.pdbqt'%(folder,), 'w')
    for modelNumber, cl in enumerate(clusters):
        f1 = open(os.path.join(folder, pdbqtFormat%(pdbqtBaseName, jobNums[cl[0]], solNums[cl[0]])))
        lines = f1.readlines()
        f1.close()
        f.write("MODEL %d\n"%(modelNumber+1))
        [f.write(l) for l in lines]
        f.write("ENDMDL\n")
    f.close()

def packageResults(basename, pdbqtBaseName, ligand, dataDict, fullReceptor=None,
                   overwrite=False, noTargetFileOutput=False):
    """
    create filename.dro file. 
    """
    bname = basename+'_dro'
    # create a temporary folder
    if overwrite and os.path.exists(bname):
        shutil.rmtree(bname)
    os.mkdir(bname)

    ## copy solutions file
    shutil.copy("%s_out.pdbqt"%(basename,), bname)

    ## copy summary file
    shutil.copy("%s_summary.dlg"%basename, bname)

    ## copy  GA solutions/logs folder
    shutil.copytree(basename, os.path.join(bname,'dockingDetails'))

    ## make a folder for input data
    inputDir = os.path.join(bname, 'input') 
    os.mkdir(inputDir)
    if not noTargetFileOutput:
        ## copy maps to inputDir
        shutil.copy(dataDict['receptorMapsFile'], inputDir)
        #for mapfile in glob(os.path.join(mapsFolder, "*.map")):
        #     shutil.copy(mapfile, inputDir)

    ## copy mapsfile data dict to inputDir
    #shutil.copy(os.path.join(mapsFolder, 'data.pkl'), inputDir)
    
    ## copy ligand to inputDir
    dataDict['basename'] = basename
    dataDict['ligandFilename'] = os.path.basename(ligand)
    dataDict['dockingPosesFilename'] = "%s_out.pdbqt"%basename
    
    shutil.copy(ligand, inputDir)
    #if FRAtomsIndices is not None:
    #    dataDict['FRAtomsIndices']= FRAtomsIndices
    if fullReceptor:
        shutil.copy(fullReceptor, inputDir)
        dataDict['fullReceptorFilename'] = os.path.basename(fullReceptor)
    #if tpoints:
    #    shutil.copy(tpoints, inputDir)
    #    dataDict['tpoints'] = os.path.basename(tpoints)

    with open(os.path.join(bname, 'data.pkl'), 'w') as f:
        pickle.dump(dataDict, f)
    #f = open(os.path.join(bname, 'data.npz'), 'wb')
    #numpy.savez(f, **dataDict)
    #f.close()
    #tar = tarfile.open(basename+'.dro', "w:gz")
    tar = tarfile.open(basename+'.dro', "w")
    tar.add(bname)
    tar.close()
    shutil.rmtree(bname)

class runADFR:

    def myprint(self, str, newline=True):
        sys.stdout.write(str)
        self.summaryFP.write(str)
        if newline:
            sys.stdout.write('\n')
            self.summaryFP.write('\n')
            
    def __init__(self):

        import multiprocessing
        self.ncpu = multiprocessing.cpu_count()
        import platform, subprocess
        system_info = platform.uname()
        _platform = system_info[0]
        if _platform == 'Windows':
            self.shell=False
            pythonsh = sys.executable
        else:
            self.shell=True
            pythonsh = sys.executable+'sh'

        cmd = os.path.join(os.path.abspath(ADFR.__path__[0]), 'bin', 'runOneSearch.py')
        self._argv = ['"%s"'%pythonsh, '"%s"'%cmd]

        self.dockedLigand = None
        self.referenceLig = None

        self.completedJobs = 0    
        self.failedJobs = 0    
        self.numberOfJobs = 0
        # per solution
        self._scores =  []
        self._FEB = []
        self._rmsdsRef =  []
        self._energies =  []
        self._jobNum =  [] # 1-based GA evolution number for this solution
        self._solNum =  [] # 1-based solution number in a GA run for this solution
        # per run
        self._seeds =  []
        self._evals =  []
        self._endStatus =  {
            'searching':0,
            '1':0,
            'no':0,
            'ran':0,
            'population':0,
            'other':0}
        self._cmds = [] # holds the commands for runnign each GA
        self.dockedLigand = None
        self.outputBaseName = None # a folder with that name will be crated to store log files and ligands
                                   # a _summary.dlg file will be create with this name too
                                   # is specified using -o on the command line
        self.fullOutput = False # when True the fodler containg the logfile and poses from each GA is kept else it
                                # will be deleted
        self.summaryFP = None
        self.jobName = None
        self.cb_start = self.cbs
        self.cb_end = self.cbe
        self.clusters = None
        self.rmsdCalc = None
        self.rmsdRefCalc = None
        self._autos = None
        
    def getPoseData(self, logFile):
        # get seed, status, evals, [scores], [rmsds], [{energies}]
        f = open(logFile)
        lines = f.readlines()
        f.close()
        nbSol = int(lines[-1].split()[-2])
        seed = int(lines[-2].split()[1])
        status = lines[-4][9:-1]
        evals = int(lines[-5].split()[3])
        scores = []
        rmsds = []
        energies = []
        #import pdb; pdb.set_trace()
        for i in range(nbSol):
            w = lines[-5-nbSol+i].split()
            scores.append(float(w[2]))
            rmsds.append(float(w[5]))
            energies.append({
                'LL': float(w[8]), 'RRL': float(w[10]), 'FRL': float(w[12]), 'FRFR': float(w[14]),
                'RRFR': float(w[16]), 'wRR': float(w[18]), 'FEB': float(w[20])})
        return seed, status, evals, scores, rmsds, energies
    
    def cbs(self, jobNum, logFile):
        pass

    def cbe(self, jobNum, logFile, percent, status, error, search, sortBy):
        if status=='FAILED':
            self._jobStatus[jobNum-1] = 3
            print 'ERROR', error
            self.failedJobs += 1
            if search=='moca':
                percent = float(self.completedJobs+self.failedJobs)/self.numberOfJobs
                sys.stdout.write('%s\r' % ('*'*int(50*percent)))
                sys.stdout.flush()
        elif status=='OK':
            self._jobStatus[jobNum-1] = 2
            self.completedJobs += 1
            if search=='moca':
                percent = float(self.completedJobs+self.failedJobs)/self.numberOfJobs
                sys.stdout.write('%s\r' % ('*'*int(50*percent)))
                sys.stdout.flush()

            # get score
            # endStatus can be 'searching', '1', 'no', 'ran' or 'population'
            seed, endStatus, evals, scores, rmsds, energies = self.getPoseData(logFile)

            if endStatus.split()[0] in self._endStatus:
                self._endStatus[endStatus.split()[0]] += 1
            else:
                self._endStatus['other'] += 1
            nbSol = len(scores)
            self._seeds[jobNum-1] =  seed
            self._evals[jobNum-1] =  evals
            self._scores.extend(scores)
            self._rmsdsRef.extend(rmsds)
            self._energies.extend(energies)
            for i in range(nbSol):
                self._jobNum.append(jobNum)
                self._solNum.append(i+1)
                self._FEB.append(energies[i]['RRL'] + energies[i]['FRL'] +
                                 self.torsdof*ADFRcc._parameters.feCoeffTors)

                # get pose coordinates
                ligandFilename = os.path.join(
                    self.outputBaseName, '%s%04d_%02d.pdbqt'%(self.jobName, jobNum, i+1))
                lig = Read(ligandFilename)
                ag = self.dockedLigand._ag
                poseCoords = lig._ag.select('segment LIG').getCoords()
                ag.addCoordset(poseCoords, 'pose %d'%(i))

            if percent==1.0:
                if search=='LGA':
                    print '%s' % ('*'*51)
                #self.myprint( "done.\n")

                # print termination info
                self.myprint( "%d successful GA evolutions consumed %d evals\n"%(
                    self.completedJobs, numpy.sum(self._evals)), newline=False)
                t = self.nbRuns
                n = self.completedJobs
                f = self.failedJobs
                a = self._endStatus['searching']
                b = self._endStatus['1']
                c = self._endStatus['no']
                d = self._endStatus['ran']
                e = self._endStatus['population']
                self.myprint(" %4d/%4d %4.1f%s runs failed"%\
                             (f, t, 100*f/float(t),'%'))
                self.myprint(" %4d/%4d %4.1f%s runs exhausted their evaluations"%\
                             (a, n, 100*a/float(n),'%'))
                self.myprint(" %4d/%4d %4.1f%s runs stopped converged 1 or 2 clusters"%\
                             (b, n, 100*b/float(n),'%'))

                self.myprint(" %4d/%4d %4.1f%s runs stopped after no improvement in clusters"%\
                             (c, n, 100*c/float(n),'%'))
                self.myprint(" %4d/%4d %4.1f%s runs stopped because GA ran out of choices"%\
                             (d, n, 100*d/float(n),'%'))
                self.myprint(" %4d/%4d %4.1f%s runs stopped because GA population converged"%\
                             (e, n, 100*e/float(n),'%'))

                self.myprint( "\nRefining results ...\n", newline=False)

                ## cluster solutions
                order = []
                scores = [] # list of scores from the self._scores for jobs that have completed
                #build scores list and list of indices of solutions to be clustered
                #for i, sc in enumerate(self._scores):
                if sortBy=='score':
                    best = min(self._scores)
                    for i, sc in enumerate(self._scores):
                        if sc is not None and sc < best+self._eCutoffForSolutions:
                            order.append(i) # because solution coords start at self.dockedLigand._ag._coords[1]
                            scores.append(sc)
                else:
                    best = min(self._FEB)
                    for i, sc in enumerate(self._FEB):
                        if sc is not None and sc < best+self._eCutoffForSolutions:
                            order.append(i) # because solution coords start at self.dockedLigand._ag._coords[1]
                            scores.append(sc)

                # make sure the 'order' list is sorted by score
                oorder = numpy.argsort(scores)
                order = numpy.array(order)[oorder]
                if len(order)>1:
                    self.clusters = clusterPoses(self.dockedLigand._ag._coords[1:], order,
                                                 self.rmsdCalc, cutOff=self.POSEclusteringRMSD)
                else:
                    self.clusters = [[0]]

                print "done.\n"
                self.printClusters()
                self.saveBestInClusterAsModels()

    def saveBestInClusterAsModels(self):
        self.myprint( 'Writing poses to '+ self.outputBaseName+'_out.pdbqt\n')
        _saveBestInClusterAsModels(self.outputBaseName, self.clusters, self._jobNum, self._solNum,
                                   self.jobName, '%s%04d_%02d.pdbqt')
        ## f = open(self.outputBaseName+'_out.pdbqt', 'w')
        ## for modelNumber, cl in enumerate(self.clusters):
        ##     poseIndex = cl[0]
        ##     f1 = open(os.path.join(self.outputBaseName,
        ##                            '%s%04d.pdbqt'%(self.jobName, poseIndex+1)))
        ##     lines = f1.readlines()
        ##     f1.close()
        ##     f.write("MODEL %d\n"%(modelNumber+1))
        ##     [f.write(l) for l in lines]
        ##     f.write("ENDMDL\n")
        ## f.close()
        
    def printClusters(self):
        self.myprint("mode |  affinity  | clust. | ref. | clust. | rmsd | energy | best |    Score  |")
        self.myprint("     | (kcal/mol) | rmsd   | rmsd |  size  | stdv |  stdv  | run  |           |")
        self.myprint("-----+------------+--------+------+--------+------+--------+------+-----------+")

        eneList = []
        seedList = []
        rsmdsList = []
        for i, cl in enumerate(self.clusters):
            # calculate rmsd0 between best in cluster and overal best
            self.rmsdCalc.setRefCoords(self.dockedLigand._ag._coords[self.clusters[0][0]+1])
            rmsd0 = self.rmsdCalc.computeRMSD(self.dockedLigand._ag._coords[cl[0]+1])
            # calculate rmsd with reference ligand for this cluster's best
            if self.referenceLig:
                rmsdRef = self.rmsdRefCalc.computeRMSD(self.dockedLigand._ag._coords[cl[0]+1])
            else:
                rmsdRef = -1
            rmsds = [0.0]
            #ene0 = self._scores[cl[0]]
            ene0 = self._FEB[cl[0]]
            ene = []
            enel = [ene0]
            seeds = [self._seeds[self._jobNum[cl[0]]-1]]
            if len(cl)>1:
                # compute stats over solutions in cluster
                self.rmsdCalc.setRefCoords(self.dockedLigand._ag._coords[cl[0]+1])
                for j in cl[1:]:
                    rmsds.append(self.rmsdCalc.computeRMSD(self.dockedLigand._ag._coords[j+1]))
                    #enel.append(self._scores[j])
                    #ene.append(ene0-self._scores[j])
                    enel.append(self._FEB[j])
                    ene.append(ene0-self._FEB[j])
                    seeds.append(self._seeds[self._jobNum[j]-1])

                self.myprint( "%4d  %11.1f %7.1f %7.1f  %6d   %5.1f %7.1f    %03d   %11.3f "%(
                    i+1, self._FEB[cl[0]], rmsd0, rmsdRef, len(cl), numpy.std(rmsds), numpy.std(ene), self._jobNum[cl[0]], self._scores[cl[0]]))
            else:
                self.myprint( "%4d  %11.1f %7.1f %7.1f  %6d      NA      NA    %03d   %11.3f"%(
                    i+1, self._FEB[cl[0]], rmsd0, rmsdRef, len(cl), self._jobNum[cl[0]], self._scores[cl[0]]))
            rsmdsList.append(rmsds)
            eneList.append(enel)
            seedList.append(seeds)
            
        self.summaryFP.write('\n\nClustering information\n')
        self.summaryFP.write('\nClustering cutoffs: GA: %5.2f final poses: %5.2f\n'%(
            self.GAclusteringRMSD,self.POSEclusteringRMSD))
        for i, cl in enumerate(self.clusters):
            self.summaryFP.write('\ncluster %d\n'%(i+1))
            self.summaryFP.write('  runs %s\n'%' '.join(['%d'%(x+1) for x in cl]))
            self.summaryFP.write('  energies %s\n'%' '.join(['%.2f'% x for x in eneList[i]]))
            self.summaryFP.write('  rmsdbest %s\n'%' '.join(['%.2f'% x for x in rsmdsList[i]]))
            self.summaryFP.write('  seeds %s\n'%' '.join(['%d'% x for x in seedList[i]]))
        
                
    def __call__(self, **kw):
        #
        # run ADFR GAs using the list of command line arguments from the sysargv list
        # 
        import subprocess, datetime
        dataDict = {}
        # get the ligand name
        ligFilename = kw['ligandFile']
        if not os.path.exists(ligFilename):
            print "ERROR: ligand file %s not found"%(ligFilename,)
            sys.exit(1)

        if not os.path.exists(kw['receptorMapsFile']):
            print "ERROR: target file %s not found"%(kw['receptorMapsFile'],)
            sys.exit(1)
            
        ligandName = os.path.splitext(os.path.split(ligFilename)[1])[0]
        dataDict['jobName'] = self.jobName = kw['jobName']
        self.outputBaseName = '%s_%s'%(ligandName, self.jobName)
        #if self.jobName != "NoName":
        #    self.outputBaseName = '%s_%s'%(ligandName, self.jobName)
        #else:
        #    now = datetime.datetime.now()
        #    self.outputBaseName = '%s_%s'%(ligandName, now.isoformat())
        seed = None
        rncpu= None
        unzippedMapsFolder = None
        nbRuns = kw.pop('nbRuns') # 50
        if nbRuns=='auto':
            if kw['search']=='LGA':
                nbRuns=50
            else:
                nbRuns=8
        else:
            nbRuns=int(nbRuns)

        if kw['search']=='LGA':
            if kw['maxEvals']=='auto':
                maxEvals = 2500000
                kw['maxEvals'] = maxEvals
            else:
                maxEvals = int(kw['maxEvals'])
        elif kw['search']=='moca':
            if kw['maxEvals']=='auto':
                f = open(kw['ligandFile'])
                lines = f.readlines()
                f.close()
                for line in lines:
                    if line[:7]=='TORSDOF':
                        torsdof = int(line.split()[1])
                kw['maxEvals'] = maxEvals =  200000 + 60000*torsdof
            else:
                maxEvals = int(kw['maxEvals'])
        skip = False

        # check for -o option so we know the folder used for results
        basename = kw.pop('logfile', None)
        if basename:
            if basename.endswith(os.sep):
                print "ERROR: output file basename is missing"%(self.outputBaseName)
                print "       please correct your command"
                sys.exit(1)
            
            self.outputBaseName = basename

        if os.path.exists(self.outputBaseName+'_out.pdbqt'):
            if kw['overwriteFiles']:
                os.remove(self.outputBaseName+'_out.pdbqt')
            else:
                print "ERROR: output file %s already exists"%(
                    self.outputBaseName+'_out.pdbqt')
                print "       please remove, or use the -o command line to provide unused output name, or specify -O to overwrite"
                sys.exit(1)
        
        if os.path.exists(self.outputBaseName+'_summary.dlg'):
            if kw['overwriteFiles']:
                os.remove(self.outputBaseName+'_summary.dlg')
            else:
                print "ERROR: output file %s already exists"%(
                    self.outputBaseName+'_summary.pdbqt')
                print "       please remove, or use the -o command line to provide unused output name, or specify -O to overwrite"
                sys.exit(1)

        self.summaryFP = open(self.outputBaseName+'_summary.dlg', 'w')
        dataDict['summaryFile'] = self.outputBaseName+'_summary.dlg'

        # create a folder for writing ligand and GA log files
        if os.path.exists(self.outputBaseName):
            if kw['overwriteFiles']:
                shutil.rmtree(self.outputBaseName)
            else:
                print "ERROR: cannot create folder %s file or folder already exists"%self.outputBaseName
                print "       please use the -o command line to provide unused output folder name or specify -O to overwrite"
                sys.exit(1)
	os.mkdir(self.outputBaseName)

        ## if os.path.exists(self.outputBaseName+'_dro'):
        ##     if kw['overwriteFiles']:
        ##         shutil.rmtree(self.outputBaseName+'_dro')
        ##     else:
        ##         print "ERROR: docking result folder %s in in the way"%(
        ##             self.outputBaseName+'_dro',)
        ##         print "       please use the -o command line to provide unused output folder name or specify -O to overwrite"
        ## os.mkdir(self.outputBaseName+'_dro')

        self.myprint( "#################################################################")
        self.myprint( "# If you used ADFR in your work, please cite:                   #")
        self.myprint( "#                                                               #")
        self.myprint( "# P.A. Ravindranath S. Forli, D.S. Goodsell, A.J. Olson and     #")
        self.myprint( "# M.F. Sanner                                                   #")
        self.myprint( "#                                                               #")
        self.myprint( "# AutoDockFR: Advances in Protein-Ligand Docking with           #")
        self.myprint( "# Explicitly Specified Binding Site Flexibility                 #")
        self.myprint( "# PLoS Comput Biol 11(12): e1004586                             #")
        self.myprint( "#                                                               #")
        self.myprint( "# DOI:10.1371/journal.pcbi.1004586                              #")
        self.myprint( "#                                                               #")
        self.myprint( "# Please see http://ccsb.scripps.edu/adfr for more information. #")
        self.myprint( "#################################################################")
        self.myprint( "")
        self.myprint( 'Docking on %s a %s computer'%(platform.node(), platform.platform(), ))
        self.myprint( 'Date %s'%datetime.datetime.now().ctime())

        dataDict['date'] = datetime.datetime.now().ctime()
        dataDict['platform'] = platform.platform()
        dataDict['node'] = platform.node()

        # load the ligand name
        self.myprint( "reading ligand %s"%os.path.abspath(ligFilename))
        mol, error = getLigandFromFile(ligFilename)
        if error is None:
            self.dockedLigand = mol
        else:
            print error
            return

        self.torsdof = getTORSDOF(ligFilename)
        atoms = self.dockedLigand.select()

        ## build RMSD calculator used for clustering of solutions
        ## check ligand for automorphisms
        ##
        self._autos = self.dockedLigand.getAutomorphisms()
        assert len(self._autos)>0
        self.rmsdCalc = MorpicRMSD(self._autos)

        # remove options not needed by runOneGA
        self.fullOutput = kw.pop('fullOutput', None)
        seed = kw.pop('seedValue')
        dataDict['seed'] = seed

        rncpu = kw.pop('maxCores')
        self.GAclusteringRMSD = kw.pop('GAclusteringRMSDCutoff')
        self.POSEclusteringRMSD = kw.pop('poseClusteringRMSDCutoff')
        if rncpu is 0:
            ncores = self.ncpu
            self.myprint( 'Detected %d cores, using %d cores'%(self.ncpu, ncores))
        else:
            assert rncpu > 0, "ERROR: maxCores a positive number, got %d"%rncpu
            ncores = min(self.ncpu, rncpu)
            self.myprint( 'Detected %d cores, request %d cores, using %d cores'%(self.ncpu, rncpu, ncores))

        dataDict['ncores'] = ncores
        dataDict['nbRuns'] = nbRuns
        self._eCutoffForSolutions = kw['eCutoffForSolutions']

        self.nbRuns = nbRuns
        self.numberOfJobs = nbRuns
        self._jobStatus =  [None]*nbRuns
        self._seeds = [None]*nbRuns
        self._evals = [0]*nbRuns
        self._jobNum = []
        self._scores =  []
        self._FEB =  []
        self._rmsdsRef = []
        self._energies = []
        self._endStatus =  {
            'searching':0,
            '1':0,
            'no':0,
            'ran':0,
            'population':0,
            'other':0}

        # build cmdline args for OneGA
        argv = self._argv
        # if target is zip file unzip and replace cmdline arguments
        receptorMapsFile = kw.pop('receptorMapsFile')
        dataDict['receptorMapsFile'] = receptorMapsFile

        covalentLigand = kw.get('covalentLigand', None)
        dataDict['covalentLigand'] = covalentLigand

        # unzip mapsFile
        from ADFR.utils.maps import MapsFile
        self.myprint( "Unpacking maps %s"%os.path.abspath(receptorMapsFile))
        mf = MapsFile(receptorMapsFile)
        mf.unzipMaps()
        unzippedMapsFolder = mf.getMapsFolder()
        receptorFilename = os.path.join(mf.getMapsFolder(),
                                        mf.getReceptorFilename())
        flexResStr = mf.getFlexResStr()
        covalentRec = mf.getCovalentBond()

        if covalentRec is not None:
            covalentRec.insert(0,
                    int(mf.getCovalentBondTorsionAtom().split()[1][1:-1]))

        if unzippedMapsFolder is None:
            self.myprint( msg)
            sys.exit(1)
        else:
            # extend argv as we do not want to pass zip file to OneGA and have the file
            # unzipped by every single GA run
            argv.extend(['-F', '"%s"'%unzippedMapsFolder,
                         '-M', 'rigidReceptor',
                         '-R', '"%s"'%receptorFilename])
            if flexResStr is not None:
                argv.extend(['-X', '"%s"'%flexResStr])
            tPtsFile = os.path.join(unzippedMapsFolder, "translationPoints.npy")
            if os.path.exists(tPtsFile):
                argv.extend(['-T', '"%s"'%tPtsFile])
                tpoints = tPtsFile
            else:
                tpoints = None
        if covalentLigand is not None:
            if covalentRec is None:
                raise RuntimeError("ERROR: covalentLigand specified but no covalentRec")
        elif covalentRec is not None:
            if covalentLigand is None:
                raise RuntimeError("ERROR: covalentRec specified but covalentLigand")

        if flexResStr is not None:
            mol = Read(receptorFilename)
            from ADFR.utils.MakeGrids import splitFlexRes
            from ADFR.utils.maps import flexResStr2flexRes
            flexRes = flexResStr2flexRes(flexResStr)
            receptorAtoms, sideChainAtoms = splitFlexRes(mol, flexRes, exclude='C O')
            self._FRAtomsIndices = sideChainAtoms.getIndices()
            self._FRAtomsIndices.sort()
        else:
            self._FRAtomsIndices = None
        dataDict['FRAtomsIndices'] = self._FRAtomsIndices
        # get reference ligand
        referenceLigName = kw.pop('reference')
        dataDict['referenceLigName'] = referenceLigName
        # load the reference ligand
        if referenceLigName is not None:
            self.myprint( "reading reference ligand %s"%os.path.abspath(referenceLigName))
            mol, error = getLigandFromFile(referenceLigName)
            if error is None:
                self.referenceLig = mol
            else:
                print error
                return
        #build rmsd calculator for ref rmsd
        if self.referenceLig:
            self.rmsdRefCalc = MorpicRMSD(self._autos, self.referenceLig._ag.getCoords())

        maxEvals = maxEvals
        dataDict['maxEvals'] = maxEvals
        self.jobName = kw.get('jobName')
        dataDict['jobName'] = self.jobName
        self.dryRun = kw.pop('dryRun')
        fixedRoot = kw.get('fixedRoot', False)
        dataDict['fixedRoot'] = fixedRoot
        # pass all other arguments:
        skip = False
        for i, word in enumerate(sys.argv[1:]):
            # overwritten arguments of ADFR only
            if word == '-o' or word== '--logFile' or \
               word == "-t" or word == "--target" or \
               word == "-S" or word == "--seed" or \
               word == "-c" or word == "--maxCores" or \
               word == "-n" or word == "--nbRuns":
                skip = True # skip value befind the argument name
            elif word == "-f" or word == "--fullOutput":
                pass
            elif word == "-T" or word == "--noTargetFileOutput":
                pass
            elif word == "--noDro":
                pass
            else:
                if skip:
                    skip = False
                else:
                    argv.append('"%s"'%word)

        if covalentRec is not None:
            argv.append('-V')
            for v in covalentRec:
                argv.append('"%s"'%v)

        # makes sure we have enough coord sets to store poses
        # first job is in coordinate set 1 NOT 0
        ag = self.dockedLigand._ag
        coords = ag.getCoords()
        #if ag.numCoordsets() < nbRuns:
        #    for i in range(ag.numCoordsets(), nbRuns):
        #        self.dockedLigand._ag.addCoordset(coords, 'pose %d'%(i))
        
        # add arguments that will be set during the loop submitting jobs
        # for seed jubNum and outputName
        argv.extend(['-S', '-1', '-j', '0', '-o', 'NoName'])
        jobs = {} # key will be process until process.poll() is not None (i.e. finished)

        from time import time, sleep
        t0 = time()
        runStatus = [None]*(nbRuns)
        procToRun = {}
        nbStart = 0 # number of started runs
        nbDone = 0 # number of completed runs

        self.myprint( "Performing search %d searches with %s ..."%(nbRuns, kw['search']))
        print "0%   10   20   30   40   50   60   70   80   90   100%"
        print "|----|----|----|----|----|----|----|----|----|----|"

        # submit the first set of jobs
        self._cmds = []
        evalsDone = {}
        for j in range(nbRuns):
            evalsDone[j] = 0
        for jobNum in range(1,min(nbRuns,ncores)+1):
            # overwrite seed
            if seed == -1:
                argv[-5] = str(random.randint(1,999999))
            else:
                argv[-5] = str(seed+jobNum-1)
            # overwrite jobNum
            argv[-3] = '%d'%jobNum
            logFile = os.path.join(self.outputBaseName, '%s%04d.dlg'%(self.jobName, jobNum))
            #logFile = '%s%04d.dlg'%(logFileRoot, jobNum)
            # overwrite outputBaseName
            argv[-1] = '"%s"'%logFile
            #print( ' '.join(argv))
            #raise
            #self.myprint( 'STARTING', jobNum, ' '.join(argv))
            if self.cb_start:
                self.cb_start(jobNum, logFile)
            # remove output file in case it exists
            try:
                os.remove(logFile)
            except OSError:
                pass

            if self.dryRun:
                print '\n*************** command ***************************\n'
                print ' '.join(argv)
                print
                sys.exit(0)
            if jobNum==1:
                cmd1 = ' '.join(argv)
            self._cmds.append(' '.join(argv))
            process = subprocess.Popen(' '.join(argv),
                                       stdout=subprocess.PIPE , 
                                       stderr=subprocess.PIPE, 
                                       bufsize = 1, shell=self.shell, cwd=os.getcwd())
            procToRun[process] = jobNum-1
            nbStart += 1
        if kw['search']== 'moca':
            totEvals = nbRuns*maxEvals
        else:
            totEvals = None
        # check for completion and start new runs ntil we are done
        while nbDone < nbRuns:
            # check running processes
            for proc, jnum in procToRun.items():
                if proc.poll() is not None: # process finished
                    if proc.returncode !=0:
                        #import pdb; pdb.set_trace()
                        runStatus[jnum] = ('Error', proc.stderr.readlines())
                        error = '\n'.join(runStatus[jnum][1])
                        status = 'FAILED'
                        #self.myprint( 'job %d ENDED WITH ERROR'%(jnum+1))
                        #self.myprint( '    '+str(self._cmds[jnum]))
                    else:
                        status = 'OK'
                        error = ''
                        runStatus[jnum] = ('OKAY', '%s%04d.dlg'%(self.jobName, jnum+1))
                        #self.myprint( 'job %d ENDED OK'%(jnum+1))
                    nbDone += 1
                    if self.cb_end:
                        self.cb_end(jnum+1,
                                    os.path.join(self.outputBaseName, '%s%04d.dlg'%(self.jobName, jnum+1)),
                                    nbDone/float(nbRuns), status, error, kw['search'], kw['sortBy'])

                    # remove process
                    del procToRun[proc]
                    if nbStart < nbRuns:
                        # start new one
                        jobNum += 1
                        if seed is not None:
                            argv[-5] = '%d'%(seed+jobNum)
                        else:
                            argv[-5] = '-1'
                        argv[-3] = '%d'%jobNum
                        logFile =  os.path.join(self.outputBaseName, '%s%04d.dlg'%(self.jobName, jobNum))
                        argv[-1] = logFile
                        #self.myprint( 'STARTING', jobNum, ' '.join(argv))
                        if self.cb_start:
                            self.cb_start(jobNum, logFile)
                        # remove output file in case it exists
                        try:
                            os.remove(logFile)
                        except OSError:
                            pass
                        process = subprocess.Popen(' '.join(argv),
                                                   stdout=subprocess.PIPE , 
                                                   stderr=subprocess.PIPE, 
                                                   bufsize = 1, shell=self.shell, cwd=os.getcwd())
                        procToRun[process] = jobNum-1
                        nbStart += 1
                        evalsDone[jnum] = maxEvals
                else:
                    logFile = os.path.join(self.outputBaseName, '%s%04d.dlg'%(self.jobName, jnum))
                    if kw['search']=='moca' and os.path.exists(logFile):
                        f = open(logFile)
                        lines = f.readlines()
                        f.close()
                        #output = proc.stdout.readlines()
                        for i in xrange(len(lines)-1, -1, -1):
                            #print jnum, lines[i]
                            w = lines[i].split()
                            if len(w) and w[0]=='LOG':
                                evalsDone[jnum] = int(w[2])
                                #print 'UP', jnum, evalsDone[jnum], logFile
                                break
            if kw['search']== 'moca':
                evals = 0
                for j in range(nbRuns):
                    evals += evalsDone[j]
                percent = float(evals)/totEvals
            else:
                percent = float(self.completedJobs+self.failedJobs)/self.numberOfJobs
            #print evals, totEvals, percent, int(50*percent)
            sys.stdout.write('%s\r' % ('*'*int(50*percent)))
            sys.stdout.flush()
            #print
            sleep(0.1)
        dt = time()-t0
        h,m,s = str(datetime.timedelta(seconds=dt)).split(':')
        self.myprint( 'Docking performed in %.2f seconds, i.e. %s hours %s minutes %s seconds '%(dt, h, m, s))

        if not os.path.exists(self.outputBaseName+'_out.pdbqt'):
            print 'ERROR: not poses generated'
            sys.exit(1)
            
        self.myprint('\n*************** first GA command ***************************')
        self.myprint(cmd1)
        # create .dro Docking Result Object file
        if kw.get("noDro", False) is False:
            self.myprint('\npackaging docking results in to %s.dro '%self.outputBaseName,
                         newline= False)
            from time import time
            t0 = time()
            noTargetFileOutput = kw.get("noTargetFileOutput", False) 
            packageResults(self.outputBaseName, self.jobName, ligFilename, dataDict,
                           fullReceptor=receptorFilename,
                           overwrite=kw['overwriteFiles'],
                           noTargetFileOutput=noTargetFileOutput)
            self.myprint('in %5.2f (s.)'%(time()-t0))
        else:
            self.myprint('results is %s.dlg and %s_out.pdbqt'%(self.outputBaseName, self.outputBaseName))
        if not self.fullOutput:
            shutil.rmtree(self.outputBaseName)
        self.summaryFP.close()

if __name__=='__main__':

    #from ADFR.utils.runADFR import runADFR
    from ADFR.utils.optParser import ArgParser

    parser = ArgParser('ADFR')
    kw = parser.parse_args()

    runner = runADFR()        
    runner(kw)

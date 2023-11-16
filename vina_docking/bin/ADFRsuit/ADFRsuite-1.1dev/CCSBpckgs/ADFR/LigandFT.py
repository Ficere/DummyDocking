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
# $Header: /mnt/raid/services/cvs/ADFR/LigandFT.py,v 1.14 2017/07/12 02:40:08 sanner Exp $
#
# $Id: LigandFT.py,v 1.14 2017/07/12 02:40:08 sanner Exp $
#
import numpy
from math import sqrt

from ADFRcc.adfr import FTBase, AtomSet, AtomSetStatic, \
     FTAtomsNode, FTDeclareAtomSet, FTDiscreteTranslation, \
     FTRotationAboutPointQuat, FTTorsion
from ADFR import _intelContest_

from random import gauss, random
class MotionBase(FTBase):

##     def __init__(self):
##         self._deltaAmplitude = None # original amplitude of deviate
##         self.deltaAmplitude = None  # current amplitued modified by SW
##         self.isCyclic = [True]*self.getNumVariables()# vector of True False for each variable
##         self.mini = [None]*self.getNumVariables()
##         self.maxi = [None]*self.getNumVariables()
##         self.length = [None]*self.getNumVariables()
##         self.active = True # set to False to disable this motion
##         self.idGenes = [0]*self.getNumVariables()
        
    ## def setBounds(self, mini, maxi):
    ##     # set maxi and mini real values for the genes of this motion object
    ##     # Use None for unboudned values
    ##     assert len(mini)==self.getNumVariables()
    ##     assert len(maxi)==self.getNumVariables()
    ##     self.mini = mini
    ##     self.maxi = maxi
    ##     self.length = [float(ma-mi) for ma, mi in zip(maxi, mini)]

    ## def setCyclic(self, isCyclic):
    ##     assert len(isCyclic)==self.getNumVariables()
    ##     self.isCyclic = isCyclic
    
    ## def setValuesFromGenes(self, genes):
    ##     realValues = []
    ##     for i,gene in enumerate(genes):
    ##         realValues.append(self.mini[i] + gene*self.length[i])
    ##     self.setVariables(realValues)
        
    ## def getGenesForValues(self, values):
    ##     genes = []
    ##     for i, value in enumerate(values):
    ##         maxi = self.maxi[i]
    ##         mini = self.mini[i]
    ##         length = self.length[i]
    ##         isCyclic = self.isCyclic[i]
    ##         if maxi is not None and value > maxi:
    ##             if isCyclic:
    ##                 value = mini + (value - maxi)%length
    ##             else:
    ##                 value = maxi
    ##         elif mini is not None and value < mini:
    ##             if isCyclic:
    ##                 value = maxi - (mini - value)%length
    ##             else:
    ##                 value = mini
    ##         genes.append( (value - self.mini[i]) / self.length[i])
    ##     return genes

    ## def scaleDeltaAmplitudeBase(self, factor):
    ##     self._deltaAmplitude = self._deltaAmplitude*factor

    ## def setDeltaAmplitude(self, delta=0.01):
    ##     self._deltaAmplitude = delta
    ##     self.deltaAmplitude = delta

    ## def resetDeltaAmplitude(self):
    ##     self.deltaAmplitude = self._deltaAmplitude

    ## def genDeviate(self, searchRate):
    ##     dev = []
    ##     ct = 0
    ##     for i in range(self.getNumVariables()):
    ##         if searchRate==1.0 or random() < searchRate:
    ##             dev.append(gauss(0.0, self.deltaAmplitude))
    ##             ct += 1
    ##         else:
    ##             dev.append(0)
    ##     return ct, dev

    ## def initialBias(self):
    ##     return [0.]*self.getNumVariables()

    ## def biasedDevs(self, dev, bias):
    ##     return [d+b for d,b in zip(dev, bias)]

    ## def scaleBias(self, bias, factor):
    ##     return [b*factor for b in bias]

    ## def newBias(self, coef1, bias, coef2, dev):
    ##     return [a*coef1 + b*coef2 for a,b in zip(bias, dev)]
    
    ## def scaleUpAmplitude(self, factor):
    ##     self.deltaAmplitude *= factor

    ## def scaleDownAmplitude(self, factor):
    ##     self.deltaAmplitude *= factor
    ##     if self.deltaAmplitude < self._deltaAmplitude:
    ##         return True
    ##     return False

    ## def applyDelta(self, values, delta, direction='+'):
    ##     """ apply delta deviation to set of motion values
    ##     while taking care of gene cyclicity
    ##     """
    ##     newValues = []
    ##     i = 0
    ##     for val, d, isCyclic in zip(values, delta, self.isCyclic):
    ##         if direction=='+':
    ##             nv = val+d
    ##         else:
    ##             nv = val-d
    ##         if nv > 1.0:
    ##             if isCyclic:
    ##                 nv = (nv - 1.)%1.0
    ##             else:
    ##                 nv = 1.0
    ##         elif nv < 0.0:
    ##             if isCyclic:
    ##                 nv = 1.0 - (0. - nv)%1.0
    ##             else:
    ##                 nv = 0.0
    ##         newValues.append(nv)
    ##     return newValues

    def mutate(self, currentGenes, mutationRate, dev=0.2):
        mutated = 0
        newGenes = []
        isCyclic = self.getCyclicPy()
        for i, gene in enumerate(currentGenes):
            if random() < mutationRate:
                mutated += 1
                nv = gauss(gene, dev)
                if nv > 1.0:
                    if isCyclic[i]:
                        nv = (nv - 1.)%1.0
                    else:
                        nv = 1.0
                elif nv < 0.0:
                    if isCyclic[i]:
                        nv = 1.0 - (0. - nv)%1.0
                    else:
                        nv = 0.0
                newGenes.append(nv)
            else:
                newGenes.append(gene)
        self.setGenes(newGenes)
        return mutated

class PyFTRotationAboutPointQuat(FTRotationAboutPointQuat, MotionBase):

    pass
    ## def __init__(self, name):
    ##     FTRotationAboutPointQuat.__init__(self, name)
    ##     MotionBase.__init__(self)
    ##     self.setDeltaAmplitude(0.01) # emulate ADFR VAR
    ##     self.setBounds( [-1.]*4, [1.]*4)
    ##     self.idGenes = self.getGenesForValues([0, 0, 0, 1])

    ## def initialize(self):
    ##     # configure rotation object with identity matrix
    ##     self.setVariables([0.,0.,0.,1.])

    ## def setValuesFromGenes(self, genes):
    ##     # need to sub class to normalize the quaternion
    ##     realValues = []
    ##     for i,gene in enumerate(genes):
    ##         realValues.append(self.mini[i] + gene*self.length[i])
    ##     a,b,c,d = realValues
    ##     n1 = 1./sqrt(a*a +b*b +c*c +d*d)
    ##     self.setVariables([a*n1, b*n1, c*n1, d*n1])

class PyFTDiscreteTranslation(FTDiscreteTranslation, MotionBase):

    pass
    ## def __init__(self, name):
    ##     FTDiscreteTranslation.__init__(self, name)
    ##     MotionBase.__init__(self)
    ##     self.setDeltaAmplitude(0.01)
    ##     # these bounds need to be set to the box corners later
    ##     self.setBounds( [-1.]*3, [1.]*3)
        
    ## def initialize(self):
    ##     self.setVariables([0.,0.,0.])


class PyFTTorsion(FTTorsion, MotionBase):

    pass
    ## def __init__(self, name):
    ##     FTTorsion.__init__(self, name)
    ##     MotionBase.__init__(self)
    ##     self.setDeltaAmplitude(0.01)
    ##     self.setBounds( [0.], [360.])

    ## def setDeltaAmplitude(self, delta=2.):
    ##     self._deltaAmplitude = delta
    ##     self.deltaAmplitude = delta

    ## def initialize(self):
    ##     self.setVariables([self.originalValue])

    ## def setValuesFromGenes(self, genes):
    ##     # need to sub class because ADFR use angle dev rather than absolute angels
    ##     origAngle = self.getOrigAngle()
    ##     angle = genes[0]*360. + origAngle
    ##     if angle >360.:
    ##         angle = 0.0 + (angle - 360.)%360.
    ##     elif angle < 0.:
    ##         angle =  360. - (0. - angle)%360.
    ##     self.setVariables([angle])
    

from mglutil.math.torsion import torsion
from ADFR import getDonorTypes

class LigandFT:
    ## def getDonorTypes(self, mol):
    ##     """
    ##     for all atoms with HD neighbors we need to assign a new type

    ##     N-HD becomes NX-HD
    ##     NA-HD becomes N2-HD
    ##     OA-HD becomes OX-HD
    ##     SA-HD becomes SX-HD
    ##     """
    ##     def hasHD(atom, atypes):
    ##         for neighbor in atom.iterBonded():
    ##             if atypes[neighbor.getIndex()]=='HD':
    ##                 return True
    ##         return False
        
    ##     adTypes = mol._ag.getData("AD_element")
    ##     for i,a in enumerate(mol._ag):
    ##         if adTypes[i]=='N':
    ##             if hasHD(a, adTypes):
    ##                adTypes[i] = 'NX' 
    ##         elif adTypes[i]=='NA':
    ##             if hasHD(a, adTypes):
    ##                 adTypes[i] = 'N2'                    
    ##         elif adTypes[i]=='OA':
    ##             if hasHD(a, adTypes):
    ##                adTypes[i] = 'OX' 
    ##         elif adTypes[i]=='SA':
    ##             if hasHD(a, adTypes):
    ##                adTypes[i] = 'SX' 
    ##     return adTypes
        
    def isPhi(self, a1, a2, a3, a4):
        # return True is a1,a2,a3,a4 for a PHI angle in a peptide i.e.
        #   elementC(previous) -- nameN -- nameCA -- elementC or
        #   elementC -- nameCA -- nameN -- elementC(nrext residue)
        # sw1 and sw2 are: either None when no switching is needed, or pairs of
        # atom indices where the first index is the current atom index defining
        # the torsion and the second index is the atom that should be used. This
        # is needed when the torsion tree picked HA as 4th atom rather the C connected
        # the the CA for instance.
        sw1 = sw2 = None
        for a in [a1,a2,a3,a4]:
            if a.getFlag('protein') is False:
                return False, sw1, sw2
        for a in [a2,a3]:
            if a.getFlag('bb') is False:
                return False, sw1, sw2
        if a2.getName()=='CA' and a3.getName()=='N':
            # check if a1 is C of same residue
            if a1.getElement()!='C':
                for b in a2.iterBonded():
                    if b==a1: continue
                    if b.getElement()=='C':
                        sw1 = [a1.getIndex(), b.getIndex()]
                        break
                if not sw1: return False, sw1, sw2
            # check a4 to be C of previous residue
            if a4.getElement()!='C':# and a4.getResnum()!=self._fstResnum: # we always want C
                for b in a3.iterBonded():
                    if b==a3: continue
                    if b.getElement()=='C' and b.getResnum()!=a3.getResnum():
                        sw2 = [a4.getIndex(), b.getIndex()]
                        break
                if not sw2: return False, sw1, sw2
            return True, sw1, sw2

        if a2.getName()=='N' and a3.getName()=='CA':
            # check if a1 is C of previous residue
            if a1.getElement()!='C':# and a2.getResnum()!=self._fstResnum:
                for b in a2.iterBonded():
                    if b==a1: continue
                    if b.getElement()=='C' and b.getResnum()!=a2.getResnum():
                        sw1 = [a1.getIndex(), b.getIndex()]
                        break
                if not sw1: return False, sw1, sw2
            # check a4 to be C of previous residue
            if a4.getElement()!='C':
                for b in a3.iterBonded():
                    if b==a3: continue
                    if b.getElement()=='C' and b.getResnum()==a3.getResnum():
                        sw2 = [a4.getIndex(), b.getIndex()]
                        break
                if not sw2: return False, sw1, sw2
            return True, sw1, sw2
        return  False, sw1, sw2

    def isPsi(self, a1, a2, a3, a4):
        # will return True if a2 and a3 have names "CA" and "C" and they are flanked
        # by element "C" and element "N"
        # sw1 privide possible needed substituion of frst atom and sw2 provides
        # a possible substituttion for last atom (see def isPhi() )
        sw1 = sw2 = None
        for a in [a1,a2,a3,a4]:
            if a.getFlag('protein') is False:
                return False, sw1, sw2
        for a in [a2,a3]:
            if a.getFlag('bb') is False:
                return False, sw1, sw2

        if a2.getName()=='CA' and a3.getName()=='C':
            # check if a1 is N of same residue
            if a1.getElement()!='N':
                for b in a2.iterBonded():
                    if b==a1: continue
                    if b.getElement()=='N' and b.getResnum()==a2.getResnum():
                        sw1 = [a1.getIndex(), b.getIndex()]
                        break
                if not sw1: return False, sw1, sw2
            # check a4 to be N of next residue
            if a4.getElement()!='N':# and a3.getResnum()!=self._lstResnum:
                for b in a3.iterBonded():
                    if b==a3: continue
                    if b.getElement()=='N' and b.getResnum()!=a3.getResnum():
                        sw2 = [a4.getIndex(), b.getIndex()]
                        break
                if not sw2:
                    return False, sw1, sw2
            return True, sw1, sw2

        #import pdb; pdb.set_trace()
        if a2.getName()=='C':# and a3.getName()=='CA':
            # check if a1 is N of next residue
            if a1.getElement()!='N':# and a2.getResnum()!=self._lstResnum:
                for b in a2.iterBonded():
                    if b==a1: continue
                    if b.getElement()=='N' and b.getResnum()!=a2.getResnum():
                        sw1 = [a1.getIndex(), b.getIndex()]
                        break
                if not sw1: return False, sw1, sw2
            # check a4 to be N of same residue
            if a4.getElement()!='N':
                for b in a3.iterBonded():
                    if b==a3: continue
                    if b.getElement()=='N' and b.getResnum()==a3.getResnum():
                        sw2 = [a4.getIndex(), b.getIndex()]
                        break
                if not sw2: return False, sw1, sw2
            return True, sw1, sw2
        return  False, sw1, sw2

    def __init__(self, mol, torTree, transPoints=None, covalentLigandAnchorAtoms=False, fixedRoot=False, neighborSearchCutoff=-1.0, useRama=True):
        """Build a flexibility tree for a prody molecule"""

        from MolKit2.molecule import Molecule
        assert isinstance(mol, Molecule)
        self.mol = mol
        
        from .torTreeFromPDBQT import TreeNode
        assert isinstance(torTree, TreeNode)
        self.torTree = torTree
        self.atomsAffectedByTorsion = {}
        self._ramaTors = None
        self.useRama = useRama
        
        #import pdb; pdb.set_trace()
        
        # build the cAutoDock atomSetStatic object for the ligand
        ag = mol._ag
        ag.setSegnames(['LIG']*len(ag))
        # use str(mol.name) in case it is unicode which adfrcc does not like
        atomSetStatic = AtomSetStatic(len(ag),  str(mol.name))
        mol._adtypes = getDonorTypes(mol.select())
        atomSetStatic.setAtomTypes(mol._adtypes)
        atomSetStatic.setCharges(ag.getCharges())
        atomSetStatic.setOrigCoords(ag.getCoords())
        # remove 1-2 ansd 1-3 interactions for scorable pairs
        # but keep 1-4 interactions. Some of these will be removed later
        # if they are within a rigid body
        #atomSetStatic.setNoScoreHops(2)
        atomSetStatic.setNoScoreHops(3) # ignore 1-4 interactions
        bonds = mol.select().getBonds()
        bl = []
        for bds in bonds[1:]:
            for bond in bds:
                bl.extend(bond)
        bl.append(-1)
        atomSetStatic.setCovalentBonds(bl)

        self.atomSetStatic = atomSetStatic

        #create a list of indices into mol._ag for atoms in the root of the torsion tree 
        rootAtomIndices = [torTree.serialToIndex[x] for x in torTree.atoms]
        rootAtoms = mol._ag[rootAtomIndices]

        # remove non bonded pairs within rigid body
        from ADFR.torTreeFromPDBQT import weedBonds
        pairs = weedBonds(mol, torTree)
        for i,j in pairs:
            self.atomSetStatic.setPairScorable(i, j, False)

        # make al1 not pairwise scorable as al1 is included in receptor grids
        # make al2 and all its neighbors not grid scorable to avoid 1-2 and 1-3
        # interactions
        if covalentLigandAnchorAtoms:

            # set all root atoms as not grid scorable
            for ind in rootAtomIndices:
                #print 'IND', ind
                atomSetStatic.setScorableGrid(ind, False)
                
            # set all neighbors of al2 as not grid scorable
            for at3 in covalentLigandAnchorAtoms[1].iterBonded():
                #print 'IND3', at3.getIndex()
                atomSetStatic.setScorableGrid(at3.getIndex(), False)

            rat1Ind = covalentLigandAnchorAtoms[0].getIndices()[0]
            #print 'IND4', rat1Ind
            atomSetStatic.setScorablePairwise(rat1Ind, False) # N not pairwise
                    
        # Create an AtomSet for the ligand
        atomSet = AtomSet(atomSetStatic)
        self.atomSet = atomSet

        # compute center of gravity of root node
        if not covalentLigandAnchorAtoms:
            g = numpy.sum( rootAtoms.getCoords(), 0)/len(rootAtoms)
            # MS
            g = mol._ag[0].getCoords() # first atom
            self.rotationCenter = g
            #print 'Rotate about', g

        if _intelContest_:
            f = open('ligandFT.txt', 'w')
            f.write('NB_LIGAND_ATOMS %d\n'%len(mol._adtypes))
            f.write('NAMES ')
            for x in ag.getNames():
                f.write('%s '%x)
            f.write('\n')
            f.write('TYPES ')
            for x in mol._adtypes:
                f.write('%s '%x)
            f.write('\n')

            f.write('CHARGES ')
            for x in ag.getCharges():
                f.write('%.3f '%x)
            f.write('\n')

            f.write('COORDS ')
            for x in ag.getCoords():
                f.write('%.3f %.3f %.3f '%tuple(x))
            f.write('\n')

            f.write('BONDS ')
            for x in bl:
                f.write('%d '%x)
            f.write('\n')

            f.write('PAIRS_NON_CORABLE ')
            for i,j in pairs:
                f.write('%d %d '%(i, j))
            f.write('\n')

            f.write('ROTATE_ABOUT %.3f %.3f %.3f\n'%tuple(g))

            f.write('NODE 0 ROOT ', )
            for x in rootAtomIndices:
                f.write('%d '%x)
            f.write('\n')
 
        FTYPE = 'd'
        self.motions = []
        self.nodes = []

        # create FT root
        ftRoot = FTDeclareAtomSet()
        ftRoot.setAtomSet(atomSet)
        self.ftRoot = ftRoot
        if not covalentLigandAnchorAtoms and not fixedRoot:
            # add translation motion for moving ligand to translation points
            ftTrans = PyFTDiscreteTranslation('Translation')
            ftRoot.addChild(ftTrans)
            self.motions.append(ftTrans)

            # set translation to move ligand center of gravity to box
            # lower left corner
            # ftTrans.setOriginOffset(g[0], g[1], g[2])

            # set preferred translations
            if transPoints and neighborSearchCutoff<0:
                ftTrans.setPreferredPoints(transPoints)
            if neighborSearchCutoff>0:
                neighborPoints=[]
                for i in numpy.arange(-2.,2.,0.2):
                    for j in numpy.arange(-2.,2.,0.2):
                        for k in numpy.arange(-2.,2.,0.2):
                            neighborPoints.append([g[0]+i*neighborSearchCutoff,g[1]+j*neighborSearchCutoff,g[2]+k*neighborSearchCutoff])
                ftTrans.setPreferredPoints(neighborPoints)
                #node = ftRoot
                #import pdb; pdb.set_trace()
            # add ligand rotation
            ftRotQuat = PyFTRotationAboutPointQuat('Rotation')
            ftTrans.addChild(ftRotQuat)
            ftRotQuat.setRotPoint(g[0], g[1], g[2])
            self.motions.append(ftRotQuat)
            node = ftRotQuat
        else:
            node = ftRoot
            
        # create core atoms rigid body
        ftCore = FTAtomsNode()
        if _intelContest_:
            self._num = 0
            self._f = f
            ftCore._num = self._num
            self._num+=1
        self.nodes.append(ftCore)
        node.addChild(ftCore)
        ftCore.setAtomIndexes(rootAtomIndices)
        #print 'ROOT', rootAtomIndices
            
        def _handleChild(torTree, ttchild, ftParent):
            # create child node
            inds = a, b, c, d = ttchild.torsionIndices
            na, nb, nc, nd = ttchild.torsionAtomNames

            print 'TOR_lig', na, nb, nc, nd, ttchild.distanceFromLeaf,
            a1, a2, a3, a4 = [self.mol._ag[x] for x in inds]
            if self.useRama:
                if self.mol._ag[b].getResname()=='PRO':
                    _isphi = _isphi = False
                    print
                else:
                    _isphi, sw1, sw2 = self.isPhi(a1, a2, a3, a4)
                    if _isphi:
                        if sw1:
                            a = sw1[1]
                            na = self.mol._ag[a].getName()
                            print 'RamaMODPhi', na, nb, nc, nd,
                        if sw2:
                            d = sw2[1]
                            nd = self.mol._ag[d].getName()
                            print 'RamaMODPhi', na, nb, nc, nd,
                        print 'PHI', self.mol._ag[b].getResname(), self.mol._ag[b].getResnum()
                    else:
                        _ispsi, sw1, sw2 = self.isPsi(a1, a2, a3, a4)
                        if _ispsi:
                            if sw1:
                                a = sw1[1]
                                na = self.mol._ag[a].getName()
                                print 'RamaMODPsi', na, nb, nc, nd,
                            if sw2:
                                d = sw2[1]
                                nd = self.mol._ag[d].getName()
                                print 'RamaMODPsi', na, nb, nc, nd,
                            print 'PSI', self.mol._ag[b].getResname(), self.mol._ag[b].getResnum()
                        else:
                            print 
            else:
		print

	    ftTorsion = PyFTTorsion('Torsion [%s]-%s--%s-[%s]'%(na,nb,nc,nd))
            #print inds, a,b,c,d
            #import pdb; pdb.set_trace()
            ftTorsion.setDihedralAtomIndexes(a,b,c,d, atomSet)
            #import pdb; pdb.set_trace()
            #a1, a2, a3, a4 = [self.mol._ag[x] for x in [a,b,c,d]]
            #c1, c2, c3, c4 = self.mol._ag[[a,b,c,d]].getCoords() # does not respect a,b,c,d order
            c1, c2, c3, c4 = [self.mol._ag.getCoords()[x] for x in [a,b,c,d]]
            #ftTorsion.originalValue = torsion(c1, c2, c3, c4)
            #ftTorsion.setDeltaAmplitude(
            #    ftTorsion._deltaAmplitude/float(ttchild.distanceFromLeaf+1))
            ftParent.addChild(ftTorsion)
            self.motions.append(ftTorsion)

            child = FTAtomsNode()
            if _intelContest_:
                child._num = self._num
                self._num+=1
            
            self.nodes.append(child)
            ftTorsion.addChild(child)
            atIndices = [torTree.serialToIndex[x] for x in ttchild.atoms]
            #print 'for', atIndices,
            #print 'relative to', [torTree.serialToIndex[x] for x in ttchild.parent().atoms]
            child.setAtomIndexes(atIndices)
            #import pdb; pdb.set_trace()
            self.atomsAffectedByTorsion[ftTorsion] = atIndices
            #print '        ', atIndices
            if _intelContest_:
                self._f.write('NODE %d PARENT %d TORSION %d %d %d %d MOVING '%(
                    child._num, ftParent._num, a,b,c,d))
                for x in atIndices:
                    self._f.write('%d '%x)
                self._f.write('\n')
                

            for c in ttchild.children:
                affectedAtoms = _handleChild(torTree, c, child)
                self.atomsAffectedByTorsion[ftTorsion].extend(affectedAtoms)
                
            return self.atomsAffectedByTorsion[ftTorsion]

        #resnum = self.mol._ag.select('protein').getResnums()
        #self._fstResnum = min(resnum)
        #self._lstResnum = max(resnum)
        
        for ttchild in torTree.children:
            _handleChild(torTree, ttchild, ftCore)

        ramaTors = []
        if self.useRama:
            geneInd = 0
            phi = []
            psi = []
            for m in self.motions:
                if isinstance(m, PyFTTorsion):
                    #print m, self.atomsAffectedByTorsion[m]
                    m.setAffectedAtoms(self.atomsAffectedByTorsion[m]+[-1], atomSet)
                    # check is phi or psi angle
                    inds = m.getDihedralAtomIndexesPy()
                    a1, a2, a3, a4 = [self.mol._ag[x] for x in inds.tolist()]
                    if a2.getResname()!='PRO' and self.isPhi(a1, a2, a3, a4)[0]:
                        phi.append( (geneInd, (a1, a2, a3, a4) ))
                    if a2.getResname()!='PRO' and self.isPsi(a1, a2, a3, a4)[0]:
                        psi.append( (geneInd, (a1, a2, a3, a4) ))
                geneInd += m.getNumVariables()

            #import pdb; pdb.set_trace()
            if len(phi) and len(psi):
                residues = {}
                for n,atoms in phi:
                    res = '%s:%s%d'%(atoms[1].getChid(), atoms[1].getResname(), atoms[1].getResnum())
                    assert not residues.has_key(res)
                    residues[res] = [(n, atoms, atoms[1].getResname()), None]
                for n,atoms in psi:
                    res = '%s:%s%d'%(atoms[1].getChid(), atoms[1].getResname(), atoms[1].getResnum())
                    if residues.get(res, None):
                        residues[res][1] = (n, atoms, atoms[1].getResname())
                    else:
                        residues[res] = [None, (n, atoms, atoms[1].getResname())]
                for k, v in residues.items():
                    if len([x for x in v if x is not None])==2: # we have both phi and psi
                        ramaTors.append(v[0][0])
                        ramaTors.append(v[1][0])
                        if v[0][2]=='ALA': ramaTors.append(1)
                        elif v[0][2]=='GLY': ramaTors.append(2)
                        else: ramaTors.append(0)
                    else:
                        if v[0] is None:
                            print 'WARNING: PHI angle not found for residue %s '%(k,)
                        else:
                            print 'WARNING: PSI angle not found for residue %s '%(k,)
                ramaTors.append(-1)
	        nbvar = 0
        	for m in self.motions:
            	    nbvar += m.getNumVariables()
        	assert max(ramaTors)<nbvar
                print 'RAMATors:', ramaTors
        self._ramaTors = ramaTors

        #import pdb; pdb.set_trace()
        ftRoot.update() # force evaluation of torsion so that
                        # ftTorsion.getOrigAngle() returns the value
        #for m in self.motions:
        #    if isinstance(m, PyFTTorsion):
        #        m.idGenes = m.getGenesForValues([m.getOrigAngle()])
        if _intelContest_:
            f.close()
        #del self._fstResnum
        #del self._lstResnum
    
    def getValues(self):
        off = 0
        values = []
        for m in self.motions:
            values.extend( m.getVariablesPy())
        return values
        
    def configureMotions(self, values):
        off = 0
        for m in self.motions:
            nb = m.getNumVariables()
            m.setVariables( values[off:off+nb])
            off += nb
        
    def setBox(self, mini, maxi):
        trans = self.motions[0]
        assert isinstance(trans, FTDiscreteTranslation), "LigandFT:setBox: trying to set bounds for translation box on wrong motion object"
        trans.setBounds(mini, maxi)
        if _intelContest_:
            f = open('ligandFT.txt', 'a')
            f.write('TRANS_BOUNDS %.3f %.3f %.3f %.3f %.3f %.3f\n'%(mini[0],mini[1],mini[2],
                                                                    maxi[0],maxi[1],maxi[2]))
            f.close()
        trans.idGenes = trans.getGenesForValuesPy(self.rotationCenter)
        trans.setIdentityGenes(trans.getGenesForValuesPy(self.rotationCenter))
        trans.mini = mini
        trans.maxi = maxi

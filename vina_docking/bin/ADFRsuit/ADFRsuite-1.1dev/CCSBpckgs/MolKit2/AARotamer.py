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

########################################################################
#
# Date: Jan 2014 Authors: Michel Sanner
#
#   sanner@scripps.edu
#       
#   The Scripps Research Institute (TSRI)
#   Molecular Graphics Lab
#   La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

#
# $Header: /mnt/raid/services/cvs/MolKit2/AARotamer.py,v 1.1.4.2 2018/02/20 19:09:12 sanner Exp $
#
# $Id: AARotamer.py,v 1.1.4.2 2018/02/20 19:09:12 sanner Exp $
#
import numpy, os

from math import pi, sqrt
degtorad = pi/180.

from mglutil.math.rotax import rotax
from mglutil.math.torsion import torsion
from prody.atomic import Atom, AtomGroup
from rotamerLib import RotamerLib

class AARotamer:
    """
Objet representing a amino acid side chain and allowing to get
coordinates of side chain atoms for rotameric side chain conformations
    """

    def setOrigCoords(self):
        self._mol._ag.setACSIndex(0)
        
    def __init__(self, atoms, angleDef, angleList, mol):
        """constructor

        object <- Rotamer(atoms, angleDef)
        
atoms   : amino acid side atoms including N and C-alpha
angleDef: the list of torsion definitions. For each torsion we get
          a list of 4 atoms defining the angle and a list of atoms
          moved only by this angle [ [[a1,a2,a3,a4], [a5,a6,..]], ...]
angleList: list of Chi angles from rotamer library
           NOTE: angleList is not used if we use getCoordsFromAngles
"""       
        self.angleDef = angleDef # for each angle [a,b,c,d], [e,f,..]]
        self.angleList = angleList # for each rotamer CH1, CH2, CH3 ...
        self.nbChi = len(self.angleDef)
        self.atoms = atoms
        self.originalAngles = []
        self._mol = mol
        self.setOrigCoords()

        self.atomOrder = range(len(atoms))
        self.alignBBmatrix = numpy.identity(4)
        
        # compute Chi angles in incomming residue
        for ang in self.angleDef:
            coords = []
            for atName in ang[0]:
                at = atoms.select('name %s'%atName)
                if at is not None and len(at)==1:
                    coords.append(at.getCoords()[0])
            if len(coords)==4:
                t = torsion(coords[0], coords[1], coords[2], coords[3])
                self.originalAngles.append(t)
            else:
                self.originalAngles.append(None)
        # build index into atoms so so that lu[122] == 0 if atmo 122 in atoms
        # if atom 122 in the molecule containing res (atoms) if the first atom
        # in the residue
        self._atomsIndexLookup = {} # key actual atom index in atoms
                                    # value sequential atom number in atoms
        for i, a in enumerate(atoms):
            self._atomsIndexLookup[a.getIndex()] = i
        lu = self._atomsIndexLookup

        ## WHERE IS THAT USED ??
        self.rotBondIndices = [] # (i,j) into self.atoms for each CHI
        self.rotatedAtomsIndices = [] # (i,j,..) into self.atoms for moved atoms
        bdAts = []
        for ang in self.angleDef:
            if len(ang[0])==0:
                continue
            a1, a2, a3, a4 = ang[0]
            atom1 = atoms.select('name %s'%a2)
            atom2 = atoms.select('name %s'%a3)
            if atom2:
                i2 = atom2.getIndices()[0]
                self.rotBondIndices.append( (lu[atom1.getIndices()[0]], lu[i2]) )

                movedAtomIndices = []
                # added atms bounded to atom2
                for atom3 in Atom(atoms._ag, i2, 0).iterBonded():
                    if atom3.getIndex() in atoms._indices:
                        if atom3.getElement()=='H':
                            movedAtomIndices.append(lu[atom3.getIndex()])

                # add atoms from listed in rotamer library as moved by this  torsion
                # and they H atoms is any
                for atName in ang[1]:
                    atom = atoms.select('name %s'%atName)
                    ind = atom.getIndices()[0]
                    movedAtomIndices.append( lu[ind] )
                    for atom3 in Atom(atoms._ag, ind, 0).iterBonded():
                        if atom3.getIndex() in atoms._indices:
                            if atom3.getElement()=='H':
                                movedAtomIndices.append(lu[atom3.getIndex()])

                        ## if Atm1 not in self.atoms and Atm1.name =='SG':
                        ##     movedAtomIndices = []
                        ##     break
                        ## elif (Atm1.parent.name[:3] not in rotlib.residueNames):
                        ##     #print atom, atom.name,Atm1.parent.name
                        ##     movedAtomIndices = []
                        ##     break
                        ## else:
                        ##     if Atm1 in self.atoms:
                        ##         if Atm1.element == 'H':
                        ##             movedAtomIndices.append( self.atoms.index(Atm1))
                self.rotatedAtomsIndices.append( movedAtomIndices )

        # create array of original coordinates
        self.origCoords = atoms.getCoords()
        self.horigCoords = [ (x[0],x[1],x[2],1.0) for x in self.origCoords]

    def getCoordsFromAngles(self, angles, check=False):
        """
returns coordinates of side chain atoms for the specified CHI angles

        coords <- getCoords(angles)
        """
        self.setOrigCoords()
        # create array of transformed coordinates (initially orig)
        coords = numpy.array(self.horigCoords, 'd')
        oc = self.origCoords
        hoc = self.horigCoords
        # used to propagate tansformation along the chain
        #cmat = identity(4)
        cmat = None
        
        #print 'angles:', index, deviations
        for i, angle in enumerate(angles):
            #print 'origAngle', self.originalAngles[i]
            #print 'ANGLE1', angle
            if self.originalAngles[i]:
                angleDelta = angle - self.originalAngles[i]
            else:
                angleDelta = 0.0
            if angle < 0:
                angle += 360.
                angleDelta += 360.
                #angleDelta = self.originalAngles[i] - angle
            #print 'angleDeltaorig', angleDelta
            #if angleDelta>180.:
            #    angleDelta = -360.+angleDeltaangles
            #elif angleDelta<-180.:
            #    angleDelta = 360.+angleDelta
            #print 'ANGLEDelta', angleDelta
            a1Ind, a2Ind = self.rotBondIndices[i]

            # compute xform matrix for Chi x with this angle
            #mat = rotax(oc[a1Ind], oc[a2Ind], (angle-self.originalAngles[i])*degtorad,
            #            transpose=1)
            mat = rotax(oc[a1Ind], oc[a2Ind], angleDelta*degtorad, transpose=1)
            
            # add mat to cmat
            if cmat is None:
                cmat = mat
            else:
                cmat = numpy.dot(mat, cmat)

            # transform atoms effected by Chi x
            for j in self.rotatedAtomsIndices[i]:
                coords[j] = numpy.dot( [hoc[j]], cmat)
                #x,y,z = oc[j]
                #coords[j][0] = cmat[0][0]*x + cmat[1][0]*y + cmat[2][0]*z + cmat[3][0]
                #coords[j][1] = cmat[0][1]*x + cmat[1][1]*y + cmat[2][1]*z + cmat[3][1]
                #coords[j][2] = cmat[0][2]*x + cmat[1][2]*y + cmat[2][2]*z + cmat[3][2]

            #sanity check
            if check:
                atoms = []
                for atName in self.angleDef[i][0]:
                    atoms.append( [x for x in self.atoms if x.name==atName][0] )
                t = torsion(coords[self.atoms.index(atoms[0])][:3],
                            coords[self.atoms.index(atoms[1])][:3],
                            coords[self.atoms.index(atoms[2])][:3],
                            coords[self.atoms.index(atoms[3])][:3])
                if t< 0.:
                    t += 360.0
                #sa1 = min( abs(angle), abs(angle)-180, abs(angle)-360) 
                #sa2 = min( abs(t), abs(t)-180, abs(t)-360)
                #sa1 = min( angle, abs(angle)-180,abs(angle)-360.0) 
                #sa2 = min( t,360-abs(t))
                #print 'FUGU', angle, t
                #if abs(sa1-sa2) > 1.0:
                #print 'WANTED', angle, 'GOTTEN', t
                if abs(angle-t) > 1.0 and abs(angle-(360.0-t)) > 1.0:
                    print 'WANTED', angle, 'GOTTEN', t, 'DIFF', abs(angle-t)
                    raise RuntimeError
                    #import pdb
                    #pdb.set_trace()


        #print coords[3:]
        # move to target residue (if defined, else alignBBmatrix is identity)
        coords = numpy.dot( coords, self.alignBBmatrix )
        return coords[self.atomOrder,:3]

    def getCoordsForRotamer(self, rotNum):
        assert rotNum<len(self.angleList)
        self.setOrigCoords()
        return self.getCoordsFromAngles(self.angleList[rotNum])

    def scoreRotamers(self, colliderAtoms):
        from prody.kdtree import KDTree
        self.setOrigCoords()
        kdtree = KDTree(colliderAtoms.getCoords())
        return self._scoreRotamers(kdtree, colliderAtoms.getRadii())

    def _scoreRotamers(self, collidersKDTree, colliderRadii):
        # scan rotamer
        self.setOrigCoords()
        scores = []
        bestScore = None
        sc = self.atoms.select('sc and not name CB')
        scRad = sc.getRadii()
        dcut = max(colliderRadii) + max(sc.getRadii()) + 0.25
        # loop over rotamers
        scores = []
        clashes = []
        favorable = []
        for i in range(len(self.angleList)):
            # set coordinates for rotamer
            scCoords = self.getCoordsForRotamer(i)[5:]
            score = 0
            # loop over ligand atoms to compute score
            lind = 0
            lf = [] # favorable
            lc = [] # clash
            for c, r in zip(scCoords, scRad):
                result = collidersKDTree(dcut, c)
                if result[0] is not None:
                    for atInd, dist in zip(result[0],result[1]):
                        tdist = r + colliderRadii[atInd]
                        if dist < tdist-0.1:
                            score += int(10*(dist-tdist))
                            lc.append((lind, atInd, dist, dist-tdist))
                        elif dist < tdist+0.25:
                            score += 1
                            lf.append((lind, atInd, dist, dist-tdist))
                lind += 1
            favorable.append(lf)
            clashes.append(lc)
            scores.append(score)
            if bestScore is None or score > bestScore:
                bestScore = score
                bestRotIndex = i
        return bestRotIndex, scores, favorable, clashes
    
    ## def build4Points(self, atoms):
    ##     """
    ##     build the list of coordinates for N, CA, CB, C in that order from atoms
    ##     """
    ##     points = []
    ##     for name in ['N', 'CA', 'CB', 'C']:
    ##         atom = atoms.select('name %s'%name)
    ##         assert len(atom)==1
    ##         points.append(atom.getCoords()[0])
    ##     return points

    def build4Points(self, atoms):
        """
        build the list of coordinates for N, CA, C in that order from atoms
        """
        self.setOrigCoords()
        points = []
        #for name in ['N', 'CA', 'C']:
        res = atoms.getHierView().iterResidues().next()
        n_ = atoms.select('name N')
        if n_ is None: raise RuntimeError('atom N not found in %s:%s`%d'%(
            res.getChid(), res.getResname(), res.getResnum()))
        points.append(n_.getCoords()[0])

        ca = atoms.select('name CA')
        if ca is None: raise RuntimeError('atom CA not found in %s:%s`%d'%(
            res.getChid(), res.getResname(), res.getResnum()))
        points.append(ca.getCoords()[0])

        cb = atoms.select('name CB')
        if cb is None:
            c_ = atoms.select('name C')
            if c_ is None:
                raise RuntimeError('atom C not found in %s:%s`%d'%(
                    res.getChid(), res.getResname(), res.getResnum()))
            aamut = AARotamerMutator()
            ag = atoms.getAtomGroup()
            cbc = aamut.buildCB(res)
            newsc = AtomGroup('cb')
            newsc._setCoords(numpy.array([cbc]))
            newsc.setNames(['CB'])
            newsc.setElements(['C'])
            newsc.setResnames([res.getResname()])
            newsc.setResnums([res.getResnum()])
            newsc.setChids([res.getChid()])
            newsc.setFlags('hetatm', [False])
            newsc.setFlags('pdbter', [False])
            newsc.setAltlocs([''])
            newsc.setIcodes([''])
            newsc.setSegnames([res.getSegname()])
            newsc.setSerials([max(ag.getSerials())+1])
            ag.extend(newsc)
        else:
            cbc = cb.getCoords()[0]
        points.append(numpy.array(cbc))
            
        v1 = (points[0]-points[1])/numpy.linalg.norm(points[0]-points[1])
        v2 = (points[2]-points[1])/numpy.linalg.norm(points[2]-points[1])
        v3 = numpy.cross(v1, v2)/numpy.linalg.norm(numpy.cross(v1, v2))

        return [points[1]+v1, points[1], points[1]+v2, points[1]+v3]

    def alignRotToResBB(self, atoms):
        """
        build self.alignBBmatrix which is a 4x4 transformation matrix that will
        transform (N,CA,CB) of this rotamer to the equivalent atoms from the
        provided list.
        """
        # build self.atomOrder, i.e. atomOrder[i] = j where i is the index of
        # the atom in the list of atoms used to build the rotamer and j the list
        # of atoms in the list of atom used to align
        names = self.atoms.getNames().tolist()
        self.atomOrder = []
        for name in atoms.getNames():
            ar = atoms.select('name %s'%name)
            self.atomOrder.append(names.index(name))
        while len(self.atomOrder) < len(self.atoms):
            self.atomOrder.append(len(self.atomOrder))
        
        # build the transformation that transforms (N,CA,CB,C) of self.atoms
        # onto the corresponding atoms in the provided list
        movingPts = self.build4Points(self.atoms)
        targetPoints = self.build4Points(atoms)
        from mglutil.math.rigidFit import RigidfitBodyAligner
        aligner = RigidfitBodyAligner()
        aligner.setRefCoords( numpy.array( targetPoints))
        aligner.rigidFit( numpy.array( movingPts) )        
        self.alignBBmatrix = aligner.getMatrix()

class RotamerMol:
    """ provide a Pmv molecule for an amino acod for which we can flip
    between rotamers"""
    
    _residueNames = set([
        'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
        'THR', 'TRP', 'TYR', 'VAL'])
    _rotlib = RotamerLib()

    def __init__(self, restype):
        restype = restype.upper()
        assert restype in self._residueNames
        
        import MolKit2
        self._mol = MolKit2.Read(os.path.join(MolKit2.__path__[0],
                                      'Data', 'canRot', '%s.pdb'%restype))
        self._mol.defaultRadii()
        self.rotamer = AARotamer(self._mol._ag.all,
                                 self._rotlib.getAngleDef(restype),
                                 self._rotlib.getAngles(restype),
                                 self._mol)
    
class CanonicalAARotamers:
    """
    """
    residueNames = [
        'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
        'THR', 'TRP', 'TYR', 'VAL']
    
    def __init__(self):
        import MolKit2
        self._mol = MolKit2.Read(os.path.join(MolKit2.__path__[0],
                                      'Data', 'canonicalAA.pdb'))
        self._mol.defaultRadii()
        from MolKit2.rotamerLib import RotamerLib
        self._rotlib = RotamerLib()
        self._canRot = {}
        
    def findMissingAndExtraAtoms(self, atoms, rtype=None):
        """
        check is all atoms needed for the rotamer rtype are in atoms
        """
        rotlib = self._rotlib
        assert len(numpy.unique(atoms.getResnums()))==1
        assert len(numpy.unique(atoms.getChids()))==1
        if rtype is None:
            rtype = atoms.getResnames()[0]
        rotAtomNames = {}
        for angAtoms, moveAtoms in rotlib.getAngleDef(rtype):
            rotAtomNames.update({}.fromkeys(angAtoms))
            rotAtomNames.update({}.fromkeys(moveAtoms))

        atomNames = {}.fromkeys(atoms.getNames())
        missing = []
        for key in rotAtomNames:
            if not atomNames.has_key(key):
                missing.append(key)
            else:
                del atomNames[key]
        return missing, atomNames.keys()

    def __call__(self, residueName):
        residueName = residueName.upper()
        assert residueName in self.residueNames
        canR = self._canRot.get(residueName, None)
        if canR is None:
            canR = AARotamer(self._mol.select('resname %s'%residueName),
                             self._rotlib.getAngleDef(residueName),
                             self._rotlib.getAngles(residueName),
                             self._mol)
            self._canRot[residueName] = canR
        return canR
    
class AARotamerMutator:
    """
    """
    residueNames = [
        'ALA','ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL']

    def __init__(self, canRotFactory=None):
        if canRotFactory is None:
            canRotFactory = CanonicalAARotamers()
        else:
            assert isinstance(canRotFactory, CanonicalAARotamers)
        self.canRotamers = canRotFactory

        
    def mutate(self, mol, chid, resnum, newResType, origAngles=None,
               bbatoms=None, scatoms=None, icode=''):
        
        ## mutate side chaind residue resnum in chain chid of mol to a new residue type
        ## new residue type can be any standadrd 3-letters amino acid code except PRO
        

        from MolKit2.molecule import Molecule
        assert isinstance(mol, Molecule)
        assert resnum in mol._ag.getResnums()
        res = mol.select('chid %s resnum `%d`'%(chid, resnum))
        for res in res.getHierView().iterResidues():
            self.mutateRes(mol, res, chid, resnum, newResType,
                           origAngles=origAngles,
                           bbatoms=bbatoms, scatoms=scatoms, icode=icode)

    def buildCB(self, res, chirality='L'):
        cacb = 1.532 # ca-cb bond length

        # compute normalize CA-N and CA->C vectors
        ca_ = res['CA'].getCoords()
        can = res['N'].getCoords()-ca_
        can /= numpy.linalg.norm(can)
        cac = res['C'].getCoords()-ca_
        cac /= numpy.linalg.norm(cac)
        zz = can+cac

        xy = numpy.cross(can, cac)
        q = numpy.dot(can, cac)

        # -CA- angles are flexible synchronous scissors
        p = -.54835422911833237471
        q = .56844147803928796138 - .75614562856111715263 * q

        # levo-dextro isomerism alters one sign
        if chirality=='D':
            q = -q

        aa = p * cacb
        bb = q * cacb
        cb = [ca_[0] + aa * zz[0] + bb * xy[0],
              ca_[1] + aa * zz[1] + bb * xy[1],
              ca_[2] + aa * zz[2] + bb * xy[2]]
        return cb

    def mutateRes(self, mol, res, chid, resnum, newResType, origAngles=None,
                  bbatoms=None, scatoms=None, icode='', altloc=None):

        from time import time
        t0 = time()
        newResType = newResType.upper()
        assert newResType in self.residueNames
        self._deletedAtoms = mol.emptySelection()
        self._addedAtoms = mol.emptySelection()

        # set the residue name to the new type so that original CHI angles can be computed
        res.setResnames(newResType)
        # make all atoms not HETATM in the residues, useful to set backbone flag in atom for modified residues
        res.setFlags('hetatm', [False]*len(res))
        ##res.select('name N CA C O').setFlags('bb', [True]*4) # readonly flag :(
        #mol._ag._flags['bb'][res.select('name N CA C O').getIndices()] = [True]*4
        #scatoms = res.select('not name N CA C O')
        ##scatoms.setFlags('backbone', [True]*len(scatoms))
        #mol._ag._flags['sc'][sc.getIndices()] = [True]*len(scatoms)
        if origAngles is None:
            origAngles = mol.measureCHIsInRes(res, resname=newResType)[0]

        resind = res.getData('resindex')[0]
        angles = []

        if newResType=='ALA':
            # delete everything up to CB atoms
            self._deletedAtoms = res.select('sc and not name CB')
            if self._deletedAtoms:
                mol.deleteAtoms(self._deletedAtoms)
            # if CB was not there we nee to build it
            if res.select('sc and name CB') is None:
                coords = self.buildCB(res)
                newsc = AtomGroup('cb')
                newsc._setCoords(numpy.array([coords]))
                newsc.setNames(['CB'])
                newsc.setElements(['C'])
                newsc.setResnames(['ALA'])
                newsc.setResnums([resnum])
                newsc.setChids([chid])
                newsc.setFlags('hetatm', [False])
                newsc.setFlags('pdbter', [False])
                newsc.setAltlocs([''])
                newsc.setIcodes([''])
                newsc.setSegnames([res.getSegnames()[0]])
                newsc.setSerials([1])
            else:
                newsc = []
        else:
            # delete side chain atoms
            if scatoms is None:
                sc = res.select('sc')
            else:
                sc = scatoms

            if sc is not None:
                mol.deleteAtoms(sc)
                self._deletedAtoms = sc

            if newResType=='GLY':
                return []

            # create canonical rotamer of new residue type
            # this is a factory. the canR objects are create a cached
            self.rotamer = canR = self.canRotamers(newResType)

            # align it target residue in mol
            if bbatoms is None:
                bbatoms = res.select('bb or name CB')
            canR.alignRotToResBB(bbatoms)

            # create a new atom group with side chain atoms from canonical rotamer
            newsc = canR.atoms.select('sc').toAtomGroup()
            newsc.setResnums([resnum]*len(newsc))
            newsc.setChids([chid]*len(newsc))
            res.setResnames([newResType]*len(res))
            newsc.setData('segment', numpy.array([res.getSegname()]*len(newsc)))
            if icode not in ['_', ' ']:
                newsc.setData('icode', [icode]*len(newsc))

            if newResType=='PRO':
                newsc.setCoords(canR.getCoordsForRotamer(0)[4:])
                angles = canR.angleList[0]
            else:
                nbChiBefore = len(origAngles)
                nbChiAfter = len(canR.angleDef)
                if nbChiBefore >= nbChiAfter:
                    angles = origAngles[:nbChiAfter]
                else:
                    angles = origAngles + [180.0]*(nbChiAfter-nbChiBefore)
                newsc.setCoords(canR.getCoordsFromAngles(angles)[4:])

        if newsc:
            mol._ag.extend(newsc)

            # set protein flag for new side chain atoms in mol

            #print "ADDED ATOMS:", self._addedAtoms
            # _addedatons is None
            # fixme : the mutated residue atoms get different values for the 'resindex' data label. This causes wrong results in HierView.iterResidues().
            natoms = len(newsc)
            if icode=='':
                icode = '_'
            selstr = 'chid %s resnum `%d` icode %s'%(chid, resnum, icode)
            if altloc:
                selstr += ' altloc _ %s'%altloc

            res.select('sc') # needed for sc flags to exist
            mol._ag._flags['sc'][-natoms:] = [True]*natoms
            mol._ag._flags['protein'][-natoms:] = [True]*natoms
            mol._ag._data['resindex'][-natoms:] = [resind]*natoms
            # the initial selection was empty because sc was not set
            self._addedAtoms = mol.select('sc and ( %s )'%selstr)
            #import pdb; pdb.set_trace()

            # add CA CB bond
            cai = mol.select('%s name CA'%selstr).getIndices()[0]
            cbi = mol.select('%s name CB'%selstr).getIndices()[0]
            mol._ag._bonds = numpy.concatenate( (mol._ag._bonds, [[cai,cbi]]))
            if mol._ag._bondOrder is not None:
                mol._ag._bondOrder["%d %d"%(cai,cbi)] = 1
            for j in range(len(mol._ag._bmap[cai])):
                n  = mol._ag._bmap[cai][j]
                if n==-1 or mol._ag._flags['deleted'][n]:
                    mol._ag._bmap[cai][j] = cbi
            for j in range(len(mol._ag._bmap[cbi])):
                n = mol._ag._bmap[cbi][j]
                if n==-1 or mol._ag._flags['deleted'][n]:
                    mol._ag._bmap[cbi][j] = cai
            mol._ag._bondIndex["%d %d"%(cai,cbi)] = len(mol._ag._bonds)-1

        if newResType=='PRO':
            # add CD N bond
            cdi = mol.select('%s name CD'%selstr).getIndices()[0]
            ni = mol.select('%s name N'%selstr).getIndices()[0]
            mol._ag._bonds = numpy.concatenate( (mol._ag._bonds, [[ni,cdi]]))
            if mol._ag._bondOrder is not None:
                mol._ag._bondOrder["%d %d"%(ni,cdi)] = 1
            for j in range(len(mol._ag._bmap[ni])):
                n  = mol._ag._bmap[ni][j]
                if n==-1 or mol._ag._flags['deleted'][n]:
                    mol._ag._bmap[ni][j] = cdi
            for j in range(len(mol._ag._bmap[cdi])):
                n = mol._ag._bmap[cdi][j]
                if n==-1 or mol._ag._flags['deleted'][n]:
                    mol._ag._bmap[cdi][j] = ni
            mol._ag._bondIndex["%d %d"%(ni,cdi)] = len(mol._ag._bonds)-1
        #print 'BEST %s %3d/%3d %4d in %.4f (s)'%(newResType, bestRotIndex,
        #                                         len(canR.angleList), bestScore, time()-t0)
        #return the angles
        return angles
    
if __name__=='__main__':
    from MolKit2 import Read
    import prody
    
    # score TYR rotamers for TYR29 in 1crn
    from MolKit2.rotamerLib import RotamerLib
    rotlib = RotamerLib()
    prot = Read('../1crn.pdb')
    rot = AARotamer(prot.select('chid A resnum 29'),
                    rotlib.getAngleDef('TYR'),
                    rotlib.getAngles('TYR'))
    colliderAtoms = prot.select('not (chid A and resnum 29)')
    best, scores, fav, clash = rot.scoreRotamers(colliderAtoms)
    print 'best:', best
    print 'scores:', scores
    print 'favorables:',fav
    print 'clashes:', clash
    #import pdb; pdb.set_trace()
    #
    canRot = CanonicalAARotamers()

    res = prot.select('chid A resnum 29')
    missing, extra = canRot.findMissingAndExtraAtoms(res)
    missing1, extra1 = canRot.findMissingAndExtraAtoms(
        res.select('bb or name CB'))
    print 'MISSING', missing1

    mutator = AARotamerMutator()
    #mutator.mutate(prot, 'A', 29, 'ARG')
    #import pdb; pdb.set_trace()
    #print prot._ag._bonds[-1]
    #print res.getNames()
    #print res.getIndices()
    #print res.getBonds()[1]
    #import pdb; pdb.set_trace()
    ## prody.writePDB('1crnY_V29.pdb', prot.select())

    ## mutator.mutate(prot, 'A', 29, 'ARG')
    ## prody.writePDB('1crnY_R29.pdb', prot.select())
    
    prot = Read('../1jff.pdb')
    #mutator.mutate(prot, 'B', 283, 'ARG')
    #print prot.select('chid B resnum 283').getNames()
    #import pdb; pdb.set_trace()

    for restype in mutator.residueNames:
        mutator.mutate(prot, 'B', 283, restype)
        #if restype!='GLY':
        #    print restype, prot.select('chid B resnum 283 sc').getNames()
        #else:
        #    print restype
        prody.writePDB('1jffY_%s283.pdb'%restype, prot.select())
    
    #mutator.mutate(prot, 'B', 283, 'ARG')
    #prody.writePDB('1jffY_R283.pdb', prot.select())
    
    ## canR = rotFactory('TYR')
    ## c = canR.getCoordsFromAngles([0,]*len(rotlib.getAngleDef('TYR')))
    ## print c
    ## prot = Read('1crn.pdb')
    ## res = prot.select('resnum 29')
    ## mat = canR.alignRotToResBB(res)

    ## chiAngles, resType = prot.measureCHIs('A', 29, rotlib)
    
    ## # try to reproduce residue coordinates from rotamer object built
    ## # in canonical form
    ## nc = canR.getCoordsFromAngles(chiAngles)
    ## pc = res.getCoords()
    ## names = res.getNames()
    ## for i in range(len(nc)):
    ##     print '%3s %9.3f (%9.3f %9.3f %9.3f) (%9.3f %9.3f %9.3f)'%(
    ##         names[i], numpy.linalg.norm(pc[i]-nc[i]),
    ##         pc[i][0], pc[i][1], pc[i][2], nc[i][0], nc[i][1], nc[i][2])


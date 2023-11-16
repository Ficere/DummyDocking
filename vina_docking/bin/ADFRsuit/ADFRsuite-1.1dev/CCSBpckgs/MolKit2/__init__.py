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
# Copyright: M. Sanner TSRI 2015
#
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/MolKit2/__init__.py,v 1.2.2.2 2017/11/10 18:49:08 sanner Exp $
# 
# $Id: __init__.py,v 1.2.2.2 2017/11/10 18:49:08 sanner Exp $
#
from __future__ import absolute_import
import string, openbabel, gzip, mmap

import mmtf as MMTF
from .babelElements import babel_elements as elements
from .cif import CIFparser
from MolKit2.cifInterface import CIFtoMolecule

for k,v in elements.items():
   if len(k)==2:
       elements[k.upper()] = v

import os
import prody
# put error as default level in ProDyDist/prody/__init__.py
prody.confProDy(verbosity='error', auto_secondary=True, silent=True)
from .molecule import Molecule, MultiMolecule

_KNOW_FORMATS = ["pdb", "pdbqt", "ent", "ent.gz",
                 "mmtf", "mol2", "sdf", "cif",
                 "cif.gz",]#, "pqr", "cif"]

def isMultiMol(filename, fileFormat):
    """check if there ismore than 1 molecule in a mol2, sdf, cif or cif.gz file
    """
    from time import time
    t0 = time()
    if filename[-3:].lower() == ".gz":
        fp = gzip.open(cifFile)
    else:
        fp = open(filename)
    data = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
    if fileFormat=='mol2':
        start = data.find("@<TRIPOS>MOLECULE\n", 0)
        nextMolStart = data.find("@<TRIPOS>MOLECULE\n", start+18)

    elif fileFormat=='sdf':
        start = data.find("$$$$\n", 0)
        nextMolStart = data.find("$$$$\n", start+5)

    elif fileFormat in ['cif', 'cif,gz']:
        start = data.find("data_", 0)
        nextMolStart = data.find("data_", start+5)

    print 'isMultiMol', time()-t0

    if nextMolStart == -1:
       fp.close()
       return False, None, None
    else:
       return True, fp, data

def Read(filename, fileFormat='auto', **kw):
    if fileFormat == 'auto':
        fileFormat = os.path.splitext(filename)[1][1:]
        if fileFormat=='gz':
            ext = os.path.basename(filename).split('.')[1]
            assert ext in ['ent', 'pdb'], "ERROR: invalid extention %s.gz"%ext
            fileFormat='ent.gz'
          
    indexingMode = kw.get('indexingMode', 'foreground')
    buildBondsByDistance = kw.get('buildBondsByDistance', True)
    assert fileFormat in _KNOW_FORMATS, 'ERROR: bad fileFormat parameter, got ".%s" expected one of %s'%(
        fileFormat,_KNOW_FORMATS)
    if fileFormat in ["pdb", "pdbqt", "ent", "ent.gz"]:
        mol = readWithPrody(filename, buildBondsByDistance=buildBondsByDistance, **kw)
    elif fileFormat == "mol2":
        mol = readWithOB(filename, fileFormat, indexingMode=indexingMode)
    elif fileFormat == "sdf":
        mol = readWithOB(filename, fileFormat, indexingMode=indexingMode)
    elif fileFormat in ["cif", "cif.gz"]:
        mol = readCIF(filename, fileFormat, indexingMode=indexingMode)
    elif fileFormat == "mmtf":
        decoder = kw.get('decoder', None)
        mol = readMMTF(filename, decoder)
    ## elif fileExt == "pqr":
    ##     mols = readPQR(filename)
    ## elif fileExt == "cif":
    ##     mols = readMMCIF(filename)
    return mol

def readMMTF(filename, decoder=None):
    if decoder is None:
        mmtfDecoder = MMTF.parse(filename)
    else:
        mmtfDecoder = decoder
    from .mglmmtf import MMTFtoPrody
    name = os.path.splitext(os.path.basename(filename))[0]
    atGroup = MMTFtoPrody(mmtfDecoder, name=name)
    mol = Molecule(name, atGroup, filename)
    mol.name = name
    return mol

def readWithPrody(filename, buildBondsByDistance=True,**kw):
    #prody.confProDy(verbosity='error')
    header = kw.pop('header', False)
    model = kw.get('model', False)
    if not kw.has_key('altloc'):
        kw['altloc'] = string.printable[:-6]
    if not os.path.exists(filename):
        raise AssertionError , "%s doesn't exist" %filename
    ext = os.path.splitext(filename)[1]
    #name = os.path.splitext(os.path.split(filename)[1])[0]
    name = os.path.split(filename)[1].split(os.path.extsep)[0]
    if ext.lower()=='.pdbqt':
        ag = prody.parsePDB(filename, format='PDBQT', **kw)
        mol = Molecule(name, ag, filename=filename)
        if buildBondsByDistance:
           mol.buildBondsByDistance()
        mol.defaultRadii()
        if header is True:
            try:
               mol.pdbHeader = prody.parsePDB(filename, model=0, header=True)
            except ValueError:
               mol.pdbHeader = {}
        findHbAcceptorsDonors(mol)
    else: #if ext.lower() in ['.pdb', '.ent', 'ent.gz']:
        ag = prody.parsePDB(filename, **kw)
        mol = Molecule(name, ag, filename=filename)
        if buildBondsByDistance:
           mol.buildBondsByDistance()
        mol.defaultRadii()
        if header is True:
            try:
                mol.pdbHeader = prody.parsePDB(filename, model=0, header=True)
            except ValueError:
                mol.pdbHeader = {}
    return mol

def readWithOB(filename, _format, group=None, indexingMode='foreground'):
    assert _format in ["mol2", "sdf"], "ERROR: bad file format, expected mol2 or sdf got %s"%_format
    obconv = openbabel.OBConversion()
    obconv.SetInFormat(_format)
    obmol = openbabel.OBMol() 
    thereIsMore = obconv.ReadFile(obmol, filename)
    # check if there are more molecules in the file
    obmolDum = openbabel.OBMol()
    thereIsMore = obconv.Read(obmolDum)
    if thereIsMore: # MultiMolecule
        mol = MultiMolecule(filename, indexingMode=indexingMode)
        mol._multi = 'molecules'

    else: # single molecule in Mol2 file
        from MolKit2.openBabelInterface import OBMolToPrody
        mol = OBMolToPrody(obmol)
        mol._multi = 'False'
        mol.name = obmol.GetTitle()
        mol._basename = os.path.splitext(os.path.basename(filename))[0]
    return mol
   
def readCIF(filename, _format, group=None, indexingMode='foreground'):
    assert _format in ["cif", "cif.gz"], "ERROR: bad file format, expected cif or cif.gz got %s"%_format
    parser = CIFparser()
    multi, fp, data = isMultiMol(filename, _format)
    if multi:
        mol = MultiMolecule(filename, indexingMode=indexingMode)
        mol._multi = 'molecules'
    else:
        if filename[-3:].lower() == ".gz":
            fp = gzip.open(cifFile)
        else:
            fp = open(filename)

        parser.parseString(fp.read())
        # we should only have one entry in the data dict
        for name, data in parser.data.items():
            mol = CIFtoMolecule(name[5:], data)
            mol._multi = 'False'
            mol.name = name
            mol._basename = os.path.splitext(os.path.basename(filename))[0]
    return mol

## def readPQR(filename, group=None):
##    from MolKit2.pdbParser import PQRParser
##    newparser = PQRParser(filename, modelsAs=modelsAs)
##    mols = newparser.parse()
##    if mols is None :
##       del newparser
##       return
##    newmol = []

##    for m in mols:
##       mol = self.addMolecule(m, group=group, filename=filename)
##       if mol is None:
##          del newparser
##          return mols.__class__([])
##       newmol.append(mol)
##    return mols.__class__(newmol)


## def readMMCIF(filename, group=None):
##    from MolKit2.mmcifParser import MMCIFParser
##    newparser = MMCIFParser(filename)
##    mols = newparser.parse()
##    if mols is None: return
##    newmol = []
##    for m in mols:
##       mol = self.addMolecule(m, group=group, filename=filename)
##       if mol is None:
##          del newparser
##          return mols.__class__([])
##       newmol.append(mol)
##    return mols.__class__(newmol)

from MolKit2.PDBresidueNames import AAnames
def getSequence(atoms):
   """return sequence string with 1 character amino acids names"""
   try:
       return "".join([AAnames[x] for x in atoms.select('ca').getResnames()]) 
   except KeyError:
       return None

def _crystalMatesMatrices(ctof, ftoc, symMats=None, n=1, contactOnly=None):
    """
    build ans returns instance matrices givem symmetry matrices (typically
    obtained from mol.pdbHeader['space_group']['symMats']), and matrices
    for going between carteesian and fractional space (typically obtained
    from mol.pdbHeader['SCALE']).

    if contactOnly is specified it has to be a tuple (coords, cutoff), where
    coords are the set of coordinates from which distances are measured and
    cutoff is the upper bound of distance for which crystal mates matrices
    will be returned.

    The first returned matrix is identity, to facilitate removing it if desired
    """
    cmMats = [] # crystal mat matrices

    # build a list of displacements in fractional space fpr the required shells
    offset = range(-n, n + 1) # e.g. -1, 0, 1
    vectors = [(x, y, z, 0) for x in offset for y in offset for z in offset]
    vectors.remove((0, 0, 0, 0)) # remove identify
    vectors.insert(0, (0, 0, 0, 0)) # insert it as first matrix

    if contactOnly:
        # get the reference coordinates and the distance cutoff
        _coords, cutoff = contactOnly
        cutoff2 = cutoff*cutoff
        # make uniform cordinates
        ucoords = np.ones((len(_coords),4), 'f')
        ucoords[:, :3] = _coords
        _fcoords = _coords.astype('f')
        # build a BHTree to checking distances
        bht = bhtreelib.BHtree( _fcoords, None, 10)
        # allocate arrays for coordinate indices and distances returned by
        # _bht.closePointsDist
        ind = np.array((len(_coords)+1,), 'i')
        dist2 = np.array((len(_coords)+1,), 'f')
    else:
        cutoff= None

    if symMats is None:
        symMats = [np.identity(4, 'f')]

    # loop over symmetry related copies matrices
    for n1, symMat in enumerate(symMats):
        # apply symmetry matrix followed by to fractional matrix
        mat0 = np.dot(ctof, symMat)
        for vec in vectors:
            vecMat = np.identity(4)
            vecMat[:3, 3] = vec[:3]
            mat1 = np.dot(vecMat, mat0) # apply translation in fractional space
            mat2 = np.dot(ftoc, mat1)   # then go back to cartesian
            if cutoff:
                # check if this crystal mate has atoms within cutoff of ref
                XMcoords = np.dot(ucoords, mat2.transpose())[:, :3]
                contact = False
                for pt in XMcoords:
                    nb = bht.closePointsDist(tuple(pt), cutoff, ind, dist2)
                    if nb:
                        if min(dist2)<=cutoff2:
                            contact = True
                if contact:
                    cmMats.append(mat2)
            else:
                cmMats.append(mat2)
        
    return cmMats

def getCrystalMatesMatrices(mol, numShells=1, cutoff=None):
    """
    return instanceMatrices for a specified number of shells around the
    asymmetric unit.
    if cutoff is specified, only matrices of copies with atoms within the
    specified cutoff from the asymmetric unit (or symmetry related
    molecules) will be returned.
    The molecule is expected to have a pdbHeader attribute which 'SCALE' and
    'space_group' keys
    """
    if not hasattr(mol, 'pdbHeader'):
        print 'ERROR: molecule %s is missing pdbHeader record. Either the molecule was not read from a PDB file, or it was parsed without the header=True'%self.name
        return []

    # get matrices for symmetry related molecules
    if mol.pdbHeader.has_key('space_group'):
        symMats = mol.pdbHeader['space_group']['symMats']
    else:
        print 'WARNING: molecule %s is missing REMARK 290. Assuming no symmetry related molecules in crystal asymetric unit. Crystal lattice might be missing copies'%self.name
        symMats = [np.identity(4)]

    # get matrices for switch between cartesian and fractional space
    if mol.pdbHeader.has_key('SCALE'):
        ctof = mol.pdbHeader['SCALE']['ctof']
        ftoc = mol.pdbHeader['SCALE']['ftoc']
    else:
        print 'ERROR: molecule %s is missing SCALEx record'%self.name
        return []

    if cutoff:
        contact = (mol._ag.getCoords(), cutoff)
    else:
        contact = None
    mats = _crystalMatesMatrices(ctof, ftoc, n=numShells,
                                 symMats=symMats, contactOnly=contact)
    return mats


## HBOND functions for PDBQT files
from .selection import Selection
import numpy as np

from bhtree import bhtreelib

def findHbAcceptorsDonors(mol):
   
    acceptor_types = {}.fromkeys(['OA', 'NA', 'SA'])
    donor_types = {}.fromkeys(['N', 'O', 'OA', 'NA', 'SA'])
    donors = [False]*len(mol._ag)
    acceptors = [False]*len(mol._ag)
    for n, atom in enumerate(mol._ag):
        atype = atom.getType().strip()
        index = atom.getIndex()
        if acceptor_types.has_key(atype):
            acceptors[n] = True
        if donor_types.has_key(atype):
           hasH = False
           for neighbor in atom.iterBonded():
               if neighbor.getType()=='HD':
                   donors[n] = True

    mol._ag.setFlags('donor', donors)
    mol._ag.setFlags('acceptor', acceptors)

def hbGeom(atoms1, atoms2, cut=3.21, foff=0):
    cut2 = cut*cut
    hbCoords = []
    faces = []
    n = 0
    if len(atoms2) > len(atoms1):
       tmp = atoms1
       atoms1 = atoms2
       atoms2 = tmp
    ag1 = atoms1.getAtomGroup()
    ag2 = atoms2.getAtomGroup()
    d1 = atoms1.select('donor')
    if d1:
        #rcoords = d1.getCoords() this give the coord set fro when atom was created
        rcoords = ag1.getCoords()[d1.getIndices()]
        a2 = atoms2.select('acceptor')
        if a2:
            for lc in ag2.getCoords()[a2.getIndices()]:
                #lc = a.getCoords()
                delta = lc-rcoords
                dist2 = np.sum(delta*delta, 1)
                for i in np.where(dist2 < cut2)[0]:
                    hbCoords.append( lc )
                    hbCoords.append( rcoords[i] )
                    faces.append( (n+foff, n+1+foff) )
                    n += 2

    a1 = atoms1.select('acceptor')
    if a1:
        #rcoords = a1.getCoords()
        rcoords = ag1.getCoords()[a1.getIndices()]
        d2 = atoms2.select('donor')
        if d2:
            #for a in d2:
                #lc = a.getCoords()
            for lc in ag2.getCoords()[d2.getIndices()]:
                delta = lc-rcoords
                dist2 = np.sum(delta*delta, 1)
                for i in np.where(dist2 < cut2)[0]:
                    hbCoords.append( lc )
                    hbCoords.append( rcoords[i] )
                    faces.append( (n+foff, n+1+foff) )
                    n += 2
    return hbCoords, faces

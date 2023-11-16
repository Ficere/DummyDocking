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
# $Header: /mnt/raid/services/cvs/MolKit2/openBabelInterface.py,v 1.2.4.1 2017/07/26 22:03:40 annao Exp $
# 
# $Id: openBabelInterface.py,v 1.2.4.1 2017/07/26 22:03:40 annao Exp $
#
import openbabel as ob

import prody
import numpy as np

from prody.atomic.atomgroup import AtomGroup
from prody.atomic import ATOMIC_FIELDS
from prody.atomic.selection import Selection as ProDySel 

from MolKit2.selection import Selection as MolKitSel
from MolKit2.molecule import Molecule

from mglutil.util.io import StringIO

#from MolKit2.selection import Selection
#from prody.atomic.selection import Selection as ProdySelection

# list of elements the OBGasteiger method can handle
OBGastKnownElements = ['H','C','N','O','F','P','S','Cl','Br','I','Al']

_pdbToOBMmol = ob.OBConversion()
_pdbToOBMmol.SetInFormat('pdb')

def ProdyToOBMol(atoms, title = "ProdyToOBMOL"):
    """
    lossy conversion from prody AtomGroup or Selection to OBMol
    """
    if isinstance(atoms, ProDySel):
        cset = atoms.getACSIndex()
    elif isinstance(atoms, AtomGroup):
        sel = atoms.all
        cset = atoms.getACSIndex()
    elif isinstance(atoms, MolKitSel):
        sel = Molecule._ag[atoms.indices]
        cset = Molecule._ag.getACSIndex()
    elif isinstance(atoms, Molecule):
        sel = Molecule._ag.all
        cset = Molecule._ag.getACSIndex()
        
    prodyToObIndex = {}
    obToProdyIndex = np.zeros( len(atoms), 'int' )
    for n, atom in enumerate(atoms):
        prodyToObIndex[atom.getIndex()] = n
        obToProdyIndex[n] = atom.getIndex()

    #import pdb; pdb.set_trace()
    # at this point atoms is a prody selection
    stream = StringIO()
    # create the OBMOL
    obmol = ob.OBMol()
    # create PDB records
    stream = StringIO()
    prody.writePDBStream(stream, atoms, csets=[cset])
    # fill OBMol based on PDB lines
    _pdbToOBMmol.ReadString(obmol, ''.join(stream.readlines()))

    obmol.prodyToObIndex = prodyToObIndex
    obmol.obToProdyIndex = obToProdyIndex
    obmol.prodyAtomList = atoms
    
    return obmol

## this version tries to build the OBMol from scratch allowing to push
## bond types etc into the obmol but it produces OBMols that fail smarts patterns matching 
def _ProdyToOBMol(sel, title = "ProdyToOBMOL"):
    """sel needs to have a .getHierViewer() and .getAtomGroup() method.
    MolKit.selection.Selection and prody.atomic.selection.Selection are OK"""
    try:
        getattr(sel, 'getHierView')
    except AttributeError:
        raise RuntimeError, "ERROR: bad argument"
    if isinstance(sel, ProDySel):
        sel = MolKitSel(sel.getAtomGroup(), sel.getIndices(), '')
    obmol = ob.OBMol()
    obmol.BeginModify()
    eTable = ob.OBElementTable()
    atoms = []
    hv = sel.getHierView()

    # bidirectional atom indices lookup between OB and ProDy mol
    prodyToObIndex = np.zeros( len(sel.getAtomGroup()), 'int' )
    obToProdyIndex = np.zeros( len(sel), 'int' )
    nbOBAtoms = 0
    for chain in hv.iterChains():
        cId = chain.getChid()
        for res in chain.iterResidues():
            obres = obmol.NewResidue()
            #print res.getResname(), res.getResnum()
            obres.SetName(res.getResname() )
            #obres.SetNum( int(res.getResnum()) ) # faile wit rednum -1 4k3q
            obres.SetNum( str(res.getResnum()) )
            obres.SetChain(cId[0])
            #print "RESNAME", obres.GetName(), res.getResname()
            #print "RESNUM", obres.GetNum(), res.getResnum()
            for atom in res.iterAtoms():
                obatom = obmol.NewAtom()
                #pind = ob.OBPairData()
                #pind.SetAttribute('prodyIndex')
                #pind.SetValue(atom.getIndex())
                #obatom.SetData(pind)
                #obatom.SetIdx(atom.getIndex()+1)
                obatom.SetAtomicNum( int(atom.getData('atomicNumber') ))
                prodyToObIndex[atom.getIndex()] = nbOBAtoms
                obToProdyIndex[nbOBAtoms] = atom.getIndex()
                nbOBAtoms += 1
                obres.AddAtom(obatom)
                obatom.SetResidue(obres)
                obres.SetAtomID(obatom, atom.getName())
                atoms.append(atom)
                #obatom.SetAtomicNum( int(atom.getData('atomicNumber') ))
                coords = atom.getCoords()
                vec = ob.vector3()
                vec.SetX(coords[0])
                vec.SetY(coords[1])
                vec.SetZ(coords[2])
                obatom.SetVector(vec)
                obres.SetHetAtom(obatom, bool(atom.getFlag('hetatm')))
                #obres.SetSerialNum(obatom, int(atom.getSerial()))
    
    ag = sel.getAtomGroup()
    if ag._bonds is not None:
        bonds = sel.getBonds()
        for bo in [1,2,3]:
            for i,j in bonds[bo]:
                if ag._bondOrder:
                    bondOrder  = ag._bondOrder['%d %d'%(i,j)]
                else:
                    bondOrder = 1
                #print i, j, int(prodyToObIndex[i]+1), int(prodyToObIndex[j]+1), bondOrder
                obmol.AddBond(int(prodyToObIndex[i]+1),
                              int(prodyToObIndex[j]+1), bondOrder)
            #print i, j,bondOrder
            #obmol.AddBond(i+1, j+1, bondOrder)
            #obmol.AddBond(int(i+1), int(j+1), bondOrder)

    ## if pmol._bonds is not None:
    ##     for n in range(len(pmol._bonds)):
    ##         i,j = pmol._bonds[n]
    ##         if pmol._bondOrder:
    ##             bondOrder  = pmol._bondOrder['%d %d'%(i,j)]
    ##         else:
    ##             bondOrder = 1
    ##         obmol.AddBond(int(i+1), int(j+1), bondOrder)

    #print 'FAGA', obmol.NumResidues(), obmol.NumBonds()

    resdat = ob.OBResidueData()
    resdat.AssignBonds(obmol,ob.OBBitVec())
    obmol.EndModify(True)
    # MS: adds a bond in chain A of 1jff
    #obmol.ConnectTheDots() # build bonds based on covalent radii
    obmol.PerceiveBondOrders()

    obmol.prodyToObIndex = prodyToObIndex
    obmol.obToProdyIndex = obToProdyIndex
    obmol.prodyAtomList = atoms

    #typer = ob.OBAtomTyper()
    #typer.Init()
    #typer.AssignHyb(obmol)

    #for res in ob.OBResidueIter(obmol):
    #    print res.GetName(), type(res.GetName()), len(res.GetName())

    #obmolNew = ob.OBMol(obmol)
    #obmolNew.prodyToObIndex = prodyToObIndex
    #obmolNew.obToProdyIndex = obToProdyIndex
    #import pybel
    #pybel.Molecule(obmol).write('pdb', 'debug.pdb', overwrite=1)
    #return obmolNew
    return obmol


def OBMolToProdyAG(obmol):
    # return Prody AtomGroup-Based MolKit.molecule.Moelcule object
    assert isinstance(obmol, ob.OBMol), "ERROR: bad argument"
    name = obmol.GetTitle()
    atomgroup = AtomGroup(name)

    asize = obmol.NumAtoms()
    prodyToObIndex = np.zeros( asize, 'int' )
    obToProdyIndex = np.zeros( asize, 'int' )

    # create numpy arrays for
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    atomtypes = np.zeros(asize, dtype=ATOMIC_FIELDS['type'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)

    #termini = np.zeros(asize, dtype=bool)
    #altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)

    charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
    radii = np.zeros(asize, dtype=ATOMIC_FIELDS['radius'].dtype)

    #bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
    #occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)

    #print obmol.NumAtoms()
    
    #maxAtypLen = 0
    elementTable = ob.OBElementTable()
    atnum = 0

    for atnum, atom in enumerate(ob.OBMolAtomIter(obmol)):
        res = atom.GetResidue()
        if res is None: # e.g. https://files.rcsb.org/ligands/view/00C_ideal.sdf
            resName = name
            resNum = 1
            atomID = atnum
        else:
            resName = res.GetName()
            resNum = res.GetNum()
            atomID = res.GetAtomID(atom)
        aIndex = atom.GetIdx()-1
        #print 'ATMO', atnum, atom.GetIdx(), resName, resNum, residue.GetAtomID(atom)
        prodyToObIndex[atnum] = aIndex
        obToProdyIndex[aIndex] = atnum

        coordinates[aIndex] = atom.GetX(),atom.GetY(),atom.GetZ()
        atomnames[aIndex] = atomID.strip()
        resnames[aIndex] = resName[:3]
        resnums[aIndex] = resNum
        chainids[aIndex] = 'A'
        hetero[aIndex] = atom.IsHeteroatom()
        #termini[aIndex] = False
        serials[aIndex] = aIndex
        #segnames[aIndex] = 'noSegm'
        atype = atom.GetType()
        #if len(atype)> maxAtypLen:
        #    maxAtypLen = len(atype)
        atomtypes[aIndex] = atype
        #print atype # different from file :( "." is missing
        elements[aIndex] = elementTable.GetSymbol(atom.GetAtomicNum())

        charges[aIndex] = atom.GetPartialCharge()
        radii[aIndex] = 1.0

        #altlocs[aIndex] = ' '
        #bfactors[aIndex] = 0.0
        #occupancies[aIndex] = 0.0

        #atId = atom.GetIdx()
        #name = res.GetAtomID(atom)
        #coords = atom.GetCoordinate()
        #atomicNumber = atom.GetAtomicNum()
        #valence = atom.GetValence()
        #hyb = atom.GetHyb()
        #atype = atom.GetType()
        #pcharge = atom.GetPartialCharge()
        #print resName, resNum, name, atId, coords, atomicNumber, valence, hyb, atype, pcharge
        atnum += 1

    atomgroup._setCoords(coordinates)
    atomgroup.setNames(atomnames)
    atomgroup.setTypes(atomtypes)
    atomgroup.setData('atomType', atomtypes)
    atomgroup.setResnames(resnames)
    atomgroup.setResnums(resnums)
    atomgroup.setChids(chainids)
    atomgroup.setFlags('hetatm', hetero)
    #atomgroup.setFlags('pdbter', termini)
    #atomgroup.setAltlocs(altlocs)
    #atomgroup.setIcodes(np.char.strip(icodes))
    atomgroup.setSerials(serials)
    #atomgroup.setBetas(bfactors)
    #atomgroup.setOccupancies(occupancies)
    #atomgroup.setSegnames(np.char.strip(segnames))
    atomgroup.setElements(np.char.strip(elements))
    atomgroup.setCharges(charges)
    atomgroup.setRadii(radii)

    # add bonds and bondOrders
    pairs = []
    bo = []
    #print 'AAAA', obmol.NumAtoms(), obmol.NumBonds()
    #print 'FAGA', obmol.NumAtoms(), obmol.NumBonds(), id(obmol)
    for bond in ob.OBMolBondIter(obmol):
        #i = bond.GetBeginAtomIdx()-1
        #j = bond.GetEndAtomIdx()-1
        i = bond.GetBeginAtomIdx()-1
        j = bond.GetEndAtomIdx()-1
        #print 'BO', i, j, bond.GetBondOrder()
        if i>j:
            a=i
            i=j
            j=a
        pairs.append( (i,j) )
        if bond.IsAromatic():
            bondOrder = 4
        elif bond.IsAmide():
            bondOrder = 5
        else:
            bondOrder = bond.GetBondOrder()
        bo.append(bondOrder)
        #print 'BOND', bond.GetBeginAtomIdx()-1, bond.GetEndAtomIdx()-1, bondOrder
        #atom1 = obmol.GetAtom(bond.GetBeginAtomIdx())
        #atomname1 = atom1.GetResidue().GetAtomID(atom1)
        #atom2 = obmol.GetAtom(bond.GetEndAtomIdx())
        #atomname2 = atom2.GetResidue().GetAtomID(atom2)
        #print atomname1, atomname2, bond.GetBondOrder()
    #print 'FAGA END'
    #print pairs
    if len(pairs):
        atomgroup.setBonds(pairs, bo)
    obmol.prodyToObIndex = prodyToObIndex
    obmol.obToProdyIndex = obToProdyIndex
    obmol.prodyAtomList = atomgroup
    return atomgroup

def OBMolToPrody(obmol):
    return Molecule(obmol.GetTitle(), OBMolToProdyAG(obmol))

  
def assignPartialChargesOB(atoms, chargeModel='gasteiger'):
    ag = atoms.getAtomGroup()
    if ag.getCharges() is None:
        ag.setCharges([0.0]*len(ag))
    
    if chargeModel=='gasteiger':
        okAtoms = atoms.select('element %s'%' '.join(OBGastKnownElements))
        noCharge = atoms.select('charge 0.0')
        if noCharge:
            from prody.atomic.selection import Selection
            notokAtoms = Selection(ag, list(set(noCharge.getIndices())-
                                            set(okAtoms.getIndices())), '')
        else:
            notokAtoms = []
            
        if len(notokAtoms):
            print 'WARNING: %s atoms have no charge. Elements %s'%(len(notokAtoms), set(notokAtoms.getElements()))
        obmol = ProdyToOBMol(okAtoms)

        charger = ob.OBChargeModel.FindType('gasteiger')
        success = charger.ComputeCharges(obmol)
        if not success:
            raise IOError("ERROR: Failed to assign gasteiger charges")
        totalCharge = obmol.GetTotalCharge ()

        # copy charges from obmol to prody molecule
        inds = []
        charges = []
        for obatom in ob.OBMolAtomIter(obmol):
            obIndex = obatom.GetIndex()
            inds.append(obmol.obToProdyIndex[obIndex])
            charges.append(obatom.GetPartialCharge())
            patom = ag[obmol.obToProdyIndex[obIndex]]
            res = obatom.GetResidue()
            assert patom.getIndex() == inds[-1]
            assert res.GetAtomID(obatom).strip()==patom.getName()

        # set the Gasteiger charges for these atoms
        ag._data['charge'][inds] = charges

        # print nonIntegral charge
        for res in atoms.getHierView().iterResidues():
            resCharge = np.sum(res.getData('charge'))
            if abs(resCharge-int(resCharge)) > 0.1:
                print 'Non integral charge' , res, resCharge
            if abs(resCharge) > 0.1:
                print 'Charged residue' , res, resCharge

    else:
        raise RuntimeError("ERROR: chargeModel %s not yet implemented, use 'gasteiger'")

import numpy
import openbabel as ob

from prody.atomic.atom import Atom

from MolKit2.molecule import Molecule
from MolKit2.openBabelInterface import (ProdyToOBMol, OBMolToProdyAG,
                                        OBMolToPrody, OBGastKnownElements)
from MolKit2.openBabelInterface import assignPartialChargesOB

from ADFR.recFromPDB import protonate

from ADFR import ADelements, PDBQTlinesForAtoms

AD4ElemsSet = set(ADelements)

from ADFR.utils.ligTT import LigTT

def addHydrogenAtomsOB(obmol, pH=7.4):
    obmol.DeleteHydrogens()
    obmol.UnsetFlag(ob.OB_PH_CORRECTED_MOL)
    for a in ob.OBMolAtomIter(obmol):
        a.SetFormalCharge(0)
    obmol.SetAutomaticFormalCharge(True)
    obmol.AddHydrogens(False, True, pH)

def assignChargesOB(obmol, chargeModel='gasteiger'):
    chargeAssigner = ob.OBChargeModel.FindType(chargeModel)
    obmol.UnsetPartialChargesPerceived()
    obmol.SetAutomaticPartialCharge(False)   
    success = chargeAssigner.ComputeCharges(obmol)
    if not success:
        msg = ('%s [calculateCharges] WARNING: charge model is '
               'missing parameters for some atoms'%mol.name)
        return None, msg
    else:
        return sum([a.GetPartialCharge() for a in ob.OBMolAtomIter(obmol)]), 'OKAY'

##
## prody-based helpers
##
def autoroot(atoms):
    """
    for a ProDy atomset find the the root of the flexibility tree
    """
    mol = atoms._ag.getMolecule()
    bestBranch = len(atoms)
    bestList = []
    for atom1 in atoms:
        if atom1.getElement()=='H': continue
        maxbranch = 0
        for atom2 in atom1.iterBonded():
            if atom2.getElement()=='H': continue
            thistree = mol.subTree(atom1, atom2).select('not hydrogen')
            thisbranch = len(thistree)
            if thisbranch>maxbranch:
                maxbranch = thisbranch
        if maxbranch<bestBranch:
            bestList = []
            bestList.append(atom1)
            bestBranch = maxbranch
        elif maxbranch==bestBranch:
            bestList.append(atom1)
    if len(bestList)==0:
        return None
    else:
        return bestList[0]

def getNonPolarHydrogens(atoms):
    """
    return the list of non polar hydrogen atoms
    """
    nph = []
    hatoms = atoms.select('hydrogen')
    if hatoms is None:
        return nph
    for h in hatoms:
        for atom in h.iterBonded():
            if atom.getElement()=='C':
                nph.append(h.getIndex())
    return nph

def findNonPolarHydrogenCharges(atoms):
    """
    return the list of indices for non polar hydrogen atoms and
    add their charges to the carbon they are attached to
    """
    nph = []
    hatoms = atoms.select('hydrogen')
    if hatoms is None:
        return nph
    for h in hatoms:
        for atom in h.iterBonded():
            if atom.getElement()=='C':
                atom.setCharge(atom.getCharge()+h.getCharge())
                h.setCharge(0.0)
                nph.append(h.getIndex())
    return nph

def unsupportedADelements(atoms):
    """
    returns list elements that are not supported by AutoDock4
    """
    return list(set(atoms.getElements()) - AD4ElemsSet)


##
## helper using both pmol and obmol
##
def setADElements(atoms, obmol=None):
    if obmol is None:
        obmol = ProdyToOBMol(atoms)

    ag = atoms.getAtomGroup()
    if 'AD_element' not in ag._data.keys():
        ag.setData('AD_element', ag.getElements())
    else:
        atoms.setData('AD_element', atoms.getElements())

    for atom in ob.OBMolAtomIter(obmol):
        index = obmol.obToProdyIndex[atom.GetIndex()]
        #print atom.GetIndex(), index
        #import pdb; pdb.set_trace()
        if atom.IsHydrogen(): ag[index].setData('AD_element', 'HD')
        elif atom.IsOxygen(): ag[index].setData('AD_element', 'OA')
        else:
            aromatic = atom.IsAromatic()
            acceptor = atom.IsHbondAcceptor()
            if atom.IsCarbon() and aromatic: ag[index].setData('AD_element', 'A')
            elif atom.IsNitrogen() and acceptor: ag[index].setData('AD_element', 'NA')
            elif atom.IsSulfur() and acceptor: ag[index].setData('AD_element', 'SA')

class ADLigand:
    """
    class to process a set of atoms describing a ligand for docking with the
    AutoDock Suite docking engines
    """

    def __init__(self, atoms, frozenBonds=None, addH='gasteiger', pH=7.4,
                 chargeModel='gasteiger', root=None):
        self.molPrep = None # will be prody molecule ready to be saved as PDBQT
        self.processAtoms(atoms, frozenBonds=frozenBonds, addH=addH, pH=pH,
                          chargeModel=chargeModel, root=root)

    def processAtoms(self, atoms, frozenBonds=None, addH='gasteiger', pH=7.4,
                     chargeModel='gasteiger', root=None):

        ## save current parameters so we can writ them in the output file
        if addH is not None:
            self._atoms = atoms.select('not hydrogen')
        else:
            self._atoms = atoms
        self._frozenBonds = frozenBonds  # [[indAt1, indAt2]]
        self._root = root
        self._addH = addH
        self._pH = pH
        self._chargeModel = chargeModel
        self.obmol = ProdyToOBMol(self._atoms)

        # check that we know all atom types
        unsup = unsupportedADelements(self._atoms)
        if len(unsup):
            raise RuntimeError('ERROR: ligand contains atoms of type %s for which ther are not parameters in the AD4.1 forcefield'%".=,".join(unsup))

        # add hydrogen atoms
        if addH=='gasteiger':
            #import pdb; pdb.set_trace()
            addHydrogenAtomsOB(self.obmol, pH=7.4)
            self._atomsH = OBMolToProdyAG(self.obmol)
        elif addH=='reduce':
            #import pdb; pdb.set_trace()
            recH = protonate(self._atoms, 'lig_H')
            recH.buildBondsByDistance()
            self._atomsH = recH._ag
            self.obmol = ProdyToOBMol(recH._ag.all)

        # partial charges
        assignPartialChargesOB(self._atomsH.all)

        # create prody molecule of charged protonated ligand
        self.molPrep = Molecule('molPrep', self._atomsH)

        # freezebonds
        bonds = [] # [[at1Ind, at2ind]] in self.molPrep
        if frozenBonds:
            d = {}
            n = 0
            for x,y,z in self.molPrep._ag.getCoords():
                d['%9.3f %9.3f %9.3f'%(x,y,z)] = n
                n += 1
            ag1 = self._atomsH
            blist = [list(b) for b in self.molPrep._ag._bonds]
            for i,j in frozenBonds:
                olda1 = ag1[i]
                at1ind = d['%9.3f %9.3f %9.3f'%tuple(olda1.getCoords())] 
                olda2 = ag1[j]
                at2ind = d['%9.3f %9.3f %9.3f'%tuple(olda2.getCoords())] 
                bonds.append([at1ind, at2ind])

        # merge non-polar hydrogens
        self.nphIndices = findNonPolarHydrogenCharges(self.molPrep._ag)
        self.molPrep._ag.setFlags("nph", [False]*len(self.molPrep._ag))
        self.molPrep._ag._flags["nph"][self.nphIndices] = True
        ligAtoms = self.molPrep._ag.select("not nph")

        from MolKit2.selection import Selection
        sel = Selection(self.molPrep._ag, ligAtoms.getIndices(), '')
        # create the receptor atom group
        self.molPrep._ag._bondData = {} # if this dict exists sel.toAtomGroup fails:(
        ag = sel.toAtomGroup('ligand')
        self.molPrep = Molecule('%s_adlig'%atoms.getAtomGroup().getTitle(), ag)

        # build torsion tree
        self.ttbuilder = LigTT(self.molPrep, frozenBonds=bonds, #noScore=self.nphIndices,
                               flexRings=False)

        # FIXME finding autoRoot should be in LigTT.buildTree is root is None 
        #get the root atom for TorTree
        if root is None:
            self._root = autoroot(self.molPrep.select('not hydrogen and not deleted'))
            if self._root is None:
                print("WARNING: automatic root atom detection failed, using first atom as root for ligand torsion tree")
                self._root = self.molPrep._ag[0]
        else:
            assert isinstance(root, Atom)
            self._root = root

        self.TTroot = self.ttbuilder.buildTree(self._root, flexRings=False)

        setADElements(self.molPrep._ag.all)
        # get acceptors
        ## acceptors = []
        ## for atom in ob.OBMolAtomIter(self.obmol):
        ##     acceptors.append(atom.IsHbondAcceptor())
        ## adElem = []
        ## import pdb; pdb.set_trace()
        ## for a, acceptor in zip(self.molPrep._ag, acceptors):
        ##     if a.ishydrogen: adElem.append('HD')
        ##     elif a.iscarbon and a.isaromatic: adElem.append('A')
        ##     elif a.isoxygen: adElem.append('OA')
        ##     elif a.isnitrogen and acceptor: adElem.append('NA')
        ##     elif a.issulfur and acceptor: adElem.append('SA')
        ##     else: adElem.append(a.getElement())

        #setADElements(atoms, obmol=None):
        #self.molPrep._ag.setData("AD_element", adElem)

    def getPDBQTlines(self):

        def _depthFirst(node):
            self._PDBQTlines.append("BRANCH %3d %3d"%tuple(node.bond))
            atomIndices = node.atoms
            atomIndices.sort()
            self._PDBQTlines.extend(
                PDBQTlinesForAtoms(self.molPrep._ag[atomIndices]))
            for i in atomIndices:
                self.renumberLU[i] = self._n
                self._n += 1

            for child in node.children:
                _depthFirst(child)

            self._PDBQTlines.append("ENDBRANCH %3d %3d"%tuple(node.bond))

        self.renumberLU = {} # used to renumber the atoms
        self._n = 1 # global "new atom number" count
        
        self._PDBQTlines = []
        self._PDBQTlines.append("ROOT")
        atomIndices = self.TTroot.atoms
        atomIndices.sort()
        self._PDBQTlines.extend(PDBQTlinesForAtoms(self.molPrep._ag[atomIndices]))
        self._PDBQTlines.append("ENDROOT")
        for i in atomIndices:
            self.renumberLU[i] = self._n
            self._n += 1

        for child in self.TTroot.children:
            _depthFirst(child)

        # now renumber the atoms and Branches
        from Support import version
        newLines = [
            "REMARK 850 file prepared by AGFR version %s"%(version.__version__),
            "REMARK 850 source %s.pdb"% (self._atomsH.getTitle(),),
            "REMARK %d active torsions:"%self.ttbuilder._nbTorsions,
            "REMARK  status: ('A' for Active; 'I' for Inactive)"]
        n = 0
        for i1, i2 in self.ttbuilder._allRotatableBonds:
            a1 = self.molPrep._ag[i1]
            a2 = self.molPrep._ag[i2]
            newLines.append("REMARK %4d  A    between atoms: %s and %s"%(
                n, '%s_%d'%(a1.getName(), self.renumberLU[i1]),
                '%s_%d'%(a2.getName(), self.renumberLU[i2])))
            n += 1

        for line in self._PDBQTlines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                serial = int(line[6:11])
                line = 'HETATM%5d%s'%(self.renumberLU[serial], line[11:])
            elif line.startswith('BRANCH'):
                s1 = int(line[6:11])
                s2 = int(line[11:15])
                line = "BRANCH %3d %3d"%(self.renumberLU[s1],self.renumberLU[s2])
            elif line.startswith('ENDBRANCH'):
                s1 = int(line[9:14])
                s2 = int(line[14:18])
                line = "ENDBRANCH %3d %3d"%(self.renumberLU[s1],self.renumberLU[s2])
            newLines.append(line)

        newLines.append("TORSDOF %d"%self.ttbuilder._torsdof)

        del self._PDBQTlines
        del self.renumberLU
        del self._n
        return newLines


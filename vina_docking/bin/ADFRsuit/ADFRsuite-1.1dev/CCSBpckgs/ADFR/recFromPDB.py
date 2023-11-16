import os, sys, subprocess, platform, urllib2, httplib, numpy
from collections import defaultdict, OrderedDict

import openbabel as ob
import prody
from prody.measure.contacts import findNeighbors
from prody.proteins.header import buildBiomolecules
from prody.atomic.selection import Selection

from Bio import pairwise2

from mglutil.util.io import StringIO

import MolKit2
from MolKit2 import Read, Molecule
from MolKit2.AARotamer import AARotamer, AARotamerMutator
from MolKit2.openBabelInterface import OBGastKnownElements, ProdyToOBMol
from MolKit2.PDBresidueNames import AAnames, stdAAset
from MolKit2.PDBresData import resData
from MolKit2.selection import PSelToMKSel

from ADFR import ADelements, PDBQTlinesForAtoms
from ADFR.cofactors import cofactorDict
from ADFR.PDBLigDicts import pdb2lig
from ADFR.PDBChemComp import chemComp
#from ADFR.utils.ADprep import mergeNonPolarHydrogenCharges, setADElements

def have_internet():
    conn = httplib.HTTPConnection("www.google.com", timeout=5)
    try:
        conn.request("HEAD", "/")
        conn.close()
        return True
    except:
        conn.close()
        return False

def getResidueKeyAlt(res, alt=None):
    # generate "chid:resname:resnum:icode:altloc" where icode and altloc
    # are '_' when blank
    if not hasattr(res, 'getChid'): # Selection of Atoms, not Residue
        hv = res.getHierView()
        assert hv.numResidues()==1
        res = hv.iterResidues().next()
    chid = res.getChid()
    resnum = res.getResnum()
    resname = res.getResname()          
    icodeSel = icode = res.getIcode()
    if icodeSel == '': icodeSel='_'
    alts = set(res.getAltlocs())-set([' '])
    if len(alts) == 0: # if there are alternate locations
        alt = '_'
    else:
        if alt is None:
            alt = list(alts)[0]
        else:
            assert alt in alts
    return '%s:%s:%s:%s:%c'%(chid, resname, resnum, icodeSel, alt)


def getResidueKey(res):
    # generate "chid:resname:resnum:icode" where icode is '__' when blank
    if not hasattr(res, 'getChid'): # Selection of Atoms, not Residue
        hv = res.getHierView()
        assert hv.numResidues()==1
        res = hv.iterResidues().next()
    chid = res.getChid()
    resnum = res.getResnum()
    resname = res.getResname()          
    icodeSel = icode = res.getIcode()
    if icodeSel == '': icodeSel='_'
    return '%s:%s:%s:%s'%(chid, resname, resnum, icodeSel)

def atomSelStr(atom):
    # generate "chid:resname:resnum:icode" where icode is '__' when blank
    chid = atom.getChid()
    resnum = atom.getResnum()
    resname = atom.getResname()          
    icodeSel = icode = atom.getIcode()
    if icodeSel == '': icodeSel='_'
    return 'chid %s resname %s resnum `%s` icode %s'%(chid, resname, resnum, icodeSel)

def parseResMutationstring(ag, resStr):
    # "A:LDH_110_#LYS; B:MLZ_505_#LYS" --> [['A', 'LDH', 110, None, 'LYS'], ['B', 'MLZ', 505 , None, 'LYS']
    #[('A:LDH:110', ['LDH', 'LYS'])]
    if not resStr: return None
    reslist = {}#[]
    for expr in resStr.split(';'):
        chid, resname, resnum, icode, newres
        _reslist = []
        for res in residues.replace(',', ' ').split():
            print res
            try:
                res, newres = res.split("#")
            except:
                print "ERROR: unsupported format for residues string %s. \nFormat: Chid:Resname_Resnum_Icode#ModResname,...;Chid:Resname_Resnum_#ModResname\nNote: Icode is optional.\n Example: %s" % (resStr, "A:LDH_110_A#LYS;B:MLZ_505_#LYS")
                raise RuntimeError
            icode = None
            items = res.split("_")
            if len(items) == 2:
                resname, resnum = items
                _rstr = "%s:%s:%s:_"%(chid, resname, resnum)
            elif len(items) == 3:
                resname, resnum, icode = items
                _rstr = "%s:%s:%s:%s"%(chid, resname, resnum, icode)
                if not len(icode):
                    icode = None
            else:
                print "ERROR: unsupported format for residues string %s. \nFormat: Chid:Resname_Resnum_Icode#ModResname,...;Chid:Resname_Resnum_#ModResname\nNote: Icode is optional.\n Example: %s" % (resStr, "A:LDH_110_A#LYS;B:MLZ_505_#LYS")
                raise RuntimeError
                return []
            #_rstr = "%s:%s:%s:_:*"%(chid, resname, resnum)
            #reslist.append([chid, resname, int(resnum), icode, newres])
            reslist[_rstr] = [resname, newres]
        #
    return reslist

class PDBProfiler:
    """
    classifies residues into the following categories:
       protein, cofactors, charge, ligands, additives, containing unsupported atoms, and water.
       These are represented by AtomGRoups or lists of residues.
       Binding sites are also identified for ligands, and cofactors
    cofactors, ligands and binding sites are identified using PDBe rest API.
    additives are identified using a list of known residues names.
    charge are identified as atoms of supported by AD4.1 but not gasteriger
    and therefore need a user-specified charge
    """
    # set of elements that are not metals
    _nonMetal = set([
        'H', 'H', 'C', 'N', 'O', 'F', 'Ne', 'P', 'S', 'Cl', 'Ar', 'Se', 'Br',
        'Kr', 'I', 'Xe', 'Rn' 'B', 'Si', 'Ge', 'As', 'Sb', 'Te', 'At'])

    # set of residue names known to be additives
    _additives = set([
        'P4G', 'EOH', '1PE', 'DMS', 'ACT', 'ACY', 'PEG', 'PGE', 'BME', 'FMT',
        'BUD', 'GOL', 'MRY', 'BU3', 'PO4', 'MLI', 'PDO', 'PG4', 'SO4', '15P',
        'CO3', 'PE8', 'PE4', 'EDO'])

    def __init__(self):
        self._receptor = None # will be the ProDy molecule from the PDB file
        self._pdbid = None # pdb id of the file, used to find cofactors and ligands
        self._assemblies = [] # will be [['ASU', ['A', 'B']], ['Biomat1', ['A']], ...]
        self._sites = {}
        self._resToMutate = {} # will contain residues that need to be mutated {getResidueKey(res): (curtype, newtype)}
        self._missingAtoms = {} # {getResidueName(res): list of missing sidechain atom names}
        self._altlocsRes = OrderedDict() # { getResidueKey(): [0, (alt1,alt2,..), (occup1, occup2,...) ]}

        self._addH = True
        self._currentBM = None
        self._customCharges  = {}
        self.resetDefaultsReceptorIncludes()
        self.setDefaultsReceptor(only=False)
        self._manualCharges = {}
        self.resetManualCharges()
        
    def resetManualCharges(self):
        # define charges for AD4.1 atoms for which OpenBabel cannot assign
        # charges using Gasteiger
        ## AD4.1  :  H  C  N  O F  MG P  S  CL CA MN FE ZN BR I  
        ## OBGast :  H  C  N  O F     P  S  CL             BR I AL]
        
        self._manualCharges['*:MN:*:*'] =  4.0 # Manganese
        self._manualCharges['*:CA:*:*'] =  2.0 # Calcium
        self._manualCharges['*:MG:*:*'] =  2.0 # Magnesium
        self._manualCharges['*:FE:*:*'] =  2.0 # Iron
        self._manualCharges['*:ZN:*:*'] =  2.0 # Zinc
        self._manualCharges['*:BR:*:*'] =  1.0 # Bromine
        #self._manualCharges['*:I:*:*']  = -1.0  # Iodine
        #self._manualCharges['*:CL:*:*'] = -1.0 # Chlorine

    def updateManualChargesTable(self, **kw):
        for k,v in kw.items():
            assert k.count(':')==3, "ERROR: incorrect ion specification string chid:resname:resnum:icode expected"
            assert isinstance(v, float)
            self._manualCharges[k] = v
        
    def getChargeforResidueKey(self, key):
        if self._manualCharges.has_key(key):  # exact match
            return self._manualCharges[key]
        # try relaxing match progressively
        chid, resname, resnum, icode = key.split(':')

        # using 1 star
        key1 = "%s:%s:%s:*"%(chid, resname, resnum)
        if self._manualCharges.has_key(key1):  # exact match
            return self._manualCharges[key1]
        key1 = "%s:%s:*:%s"%(chid, resname, icode)
        if self._manualCharges.has_key(key1):  # exact match
            return self._manualCharges[key1]
        key1 = "*:%s:%s:%s"%(resname, resnum, icode)
        if self._manualCharges.has_key(key1):  # exact match
            return self._manualCharges[key1]

        # using 2 stars
        key2 = "%s:%s:*:*"%(chid, resname)
        if self._manualCharges.has_key(key2):  # exact match
            return self._manualCharges[key2]
        key2 = "*:%s:%s:*"%(resname, resnum)
        if self._manualCharges.has_key(key2):  # exact match
            return self._manualCharges[key2]
        key2 = "*:%s:*:%s"%(resname, icode)
        if self._manualCharges.has_key(key2):  # exact match
            return self._manualCharges[key2]

        # using 3 stars
        key3 = "*:%s:*:*"%(resname)
        if self._manualCharges.has_key(key3):  # exact match
            return self._manualCharges[key3]

        raise RuntimeError, 'no ion charge found for residue %s'%key

    def assignDefaultManualCharges(self):
        manual = self._receptor.select('category manch')
        if manual:
            for res in manual.getHierView().iterResidues():
                key = getResidueKey(res)
                charge = self.getChargeforResidueKey(key)
                res.setCharges([charge])
    
    def resetDefaultsReceptorIncludes(self):
        self.includeInReceptor = {'prote' : True,
                                  'nucle' : True,
                                  'modif' : True,
                                  'cmpsc' : True,
                                  'cofac' : True,
                                  'missa' : True,
                                  'manch' : False, # manual charge
                                  'addit' : False,
                                  'ligan' : False,
                                  '_wate' : False,
                                  'other' : False,

                                  #'missS' : True, # classify2, include if category is true as well
                                  #'missA' : True, # classify2, include if category is true as well
                                  }
        
    def setDefaultsReceptor(self, only=False, **kw):
        """keys can be 'prote', 'nucle', 'modif', 'cmpsc', 'cofac', 'missa',
'manch', 'addit, 'ligan', 'water' and 'other'. Values are True of False
if only is specified, all keys not present to kw are set to False"""
        if only:
            for k in self.includeInReceptor.keys():
                self.includeInReceptor[k] = False
                
        for k,v in kw.items():
            assert k in self.includeInReceptor.keys(), "ERROR: invalud key %s for includeInReceptor\n"
            assert isisntance(v, bool)
            self.includeInReceptor[k] = v

    def readFromFile(self, filename):
        mol = Read(filename, header=True)
        self.setReceptor(mol)

    def readFromLines(self, pdbid, lines):
        ag = prody.parsePDBStream(lines)
        mol = Molecule(pdbid, ag)
        self.setReceptor(mol)

    def _makeAssemblyList(self):
        self._assemblies = OrderedDict()

        self._receptor._ag.all.setSegnames('A')
        self._assemblies['ASU'] = [set(self._receptor._ag.getChids()), self._receptor]
        self._currentBM = self._receptor
        bmTrans = self._receptor.pdbHeader.get('biomoltrans', None)
        nbAtASU = len(self._receptor._ag)
        if bmTrans:
            bms = buildBiomolecules(self._receptor.pdbHeader, self._receptor._ag)
            # bms could be a list of atom groups or just one atom group.
            if not isinstance (bms, (list, tuple)):
                bms = [bms]
            keys = bmTrans.keys()
            keys.sort()
            for k in keys:
                v = bmTrans[k]
                n = int(k)
                chidlist = []
                for data in v:
                    nbMat = (len(data)-1)/3
                    chids = data[0]*nbMat
                    chidlist.extend(chids)
                chids = set(chidlist)
                firstname = [x for x in self._assemblies.keys() if x.startswith('ASU')][0]
                chset, mol = self._assemblies[firstname]
                ## using only chids does not work (eg 2dpr.pdb) as same chaisn can be used
                ## with different transformations
                ## if chset==chids and len(chidlist)==len(chset):
                ##     #del self._assemblies[firstname]
                ##     self._assemblies.pop(firstname)
                ##     self._assemblies['%s and Biomol %s'%(firstname, n)] = [chset, mol]
                ## else:
                ##     self._assemblies['Biomol %d'%n] = [chids, Molecule('Biomol %d'%n, bms[n-1])]
                if len(bms[n-1])!=nbAtASU:
                    self._assemblies['Biomol %d'%n] = [chids, Molecule('Biomol %d'%n, bms[n-1])]
                else:
                    # match up atoms. if all ,atch up it is hte same construct
                    match = findNeighbors(bms[n-1], 0.1, self._receptor._ag)
                    if len(match)==nbAtASU:
                        self._assemblies.pop(firstname)
                        self._assemblies['%s and Biomol %s'%(firstname, n)] = [chset, mol]
                    else:
                        self._assemblies['Biomol %d'%n] = [chids, Molecule('Biomol %d'%n, bms[n-1])]

        #print "mols in _assemblies:", self._assemblies.keys()
            
    def setReceptor(self, mol, defaultReceptor=None, defaultLigand=None):
        assert isinstance(mol, MolKit2.Molecule)
        self._receptor = mol
        self._pdbid = self._receptor.name
        self._classify(self._pdbid, defaultReceptor=defaultReceptor,
                       defaultLigand=defaultLigand)
        self._receptor._ag.setData('charge', [0.0]*len(self._receptor._ag))
        self.assignDefaultManualCharges()
        self._makeAssemblyList()
        
    def isManualCharge(self, res):
        return res.numAtoms()==1 and self._manualCharges.has_key('*:%s:*:*'%res.getElements()[0].upper())

    def getCofactors(self, pdbid):
        """use PDBe rest API to get cofactors for a PDB entry"""
        cofactors = []
        for res in self._receptor._ag.getHierView().iterResidues():
            if cofactorDict.has_key(res.getResname()):
                cofactors.append(res)
        return cofactors
      
        ## if not have_internet():
        ##     return cofactors

        ## url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/cofactor/%s'%pdbid
        ## try:
        ##     site = urllib2.urlopen(url)
        ##     d = eval(site.read())
        ## except urllib2.HTTPError:
        ##     return cofactors

        ## for d1 in d[pdbid]:
        ##     chid = d1['chain_id']
        ##     resname = d1['chem_comp_id']
        ##     resnum = str(d1['author_residue_number'])+d1['author_insertion_code']
        ##     res = self._receptor._ag.select('resname %s resnum %s'%(resname, resnum))
        ##     res0 = res.getHierView().iterResidues().next()
        ##     cofactors.append(res0)
        ##     #cofactors.append('%s:%s:%s'%(chid, resname, resnum))
        ## return cofactors

    def getLigands(self, pdbid):
        """use PDBe rest API to get ligands for a PDB entry
        only ligands that are not additives, cofactors or manual charges are reported
        """
        ligands = []
        ligandNames = {}
        for lig in pdb2lig.get(pdbid, []):
            atoms = self._receptor._ag.select('resname %s'%lig)
            if atoms is None:
                print 'Warning: ligand %s not found in structure. Consider getting the latest version of the pdb file'%lig
                continue
            # if a residues is already tagged as missing something for instance it should not be a ligand
            for frag in prody.findFragments(atoms):
                if len(frag)<5: continue
                if list(set(atoms.getData('category')))[0]=='other':
                    for res in atoms.getHierView().iterResidues():
                        if res.select('category cofac addit manch') is None:
                            ligands.append(res)
                            ligandNames[lig] = '%s (%s)'%(chemComp[lig]['name'],chemComp[lig]['formula']) 
        return ligands, ligandNames
        #if not have_internet():
        #    return ligands, ligandNames
        
        #url='https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/%s'%pdbid
        #try:
        #    site = urllib2.urlopen(url)
        #    d = eval(site.read())
        #except urllib2.HTTPError:
        #    return ligands, ligandNames

        #for d1 in d[pdbid]:
        #    chid = d1['chain_id']
        #    resname = d1['chem_comp_id']
        #    #import pdb; pdb.set_trace()
        #    ligandNames[d1['chem_comp_id']] = d1['chem_comp_name']
        #    resnum = str(d1['author_residue_number'])
        #    icode = d1['author_insertion_code'] if d1['author_insertion_code'] else '_'
        #    if not resnum: continue
        #    res = self._receptor._ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
        #    if res is None:
        #        continue # file might be aclled 1b2d but SO4 in chain B was removed
        #    ## what is ligand has multiple residues ??
        #    ## should we iterate here ??
        #    if res.getHierView().numResidues()>1:
        #        raise RuntimeError('pdb with ligand containing more than 1 residue')
        #    res0 = res.getHierView().iterResidues().next()
        #    if res0.select('category cofac addit metal') is None:
        #        ligands.append(res0)
        #import pdb; pdb.set_trace()
        #return ligands, ligandNames

    ## def getSites(self, pdbid, ligandNames):
    ##     """use PDBe rest API to get binding site for a PDB entry
    ##     only sites for ligands that are not additives, cofactors or metals are
    ##     reported
    ##     """
    ##     sites = {}
    ##     if not have_internet():
    ##         return sites

    ##     url ='https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/%s'%pdbid
    ##     try:
    ##         site = urllib2.urlopen(url)
    ##         d = eval(site.read().replace('null', '""'))
    ##     except urllib2.HTTPError:
    ##         return sites
        
    ##     for d1 in d[pdbid]:
    ##         d2 = d1['ligand_residues'][0]
    ##         #import pdb; pdb.set_trace()
    ##         chid = d2['chain_id'] if d2['chain_id'] else '_'
    ##         resname = d2['chem_comp_id']
    ##         resnum = str(d2['author_residue_number'])
    ##         icode = d2['author_insertion_code'] if d2['author_insertion_code'] else '_'
    ##         if not resnum:continue
    ##         res = self._receptor._ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
    ##         if not res: continue
    ##         #try: # (1cwk.pdb)
    ##         #    res = self._receptor._ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
    ##         #except: continue
    ##         res0 = res.getHierView().iterResidues().next()
    ##         if res0.select('category cofac addit metal'):
    ##             continue
    ##         siteRes = []
    ##         siteWaters = []
    ##         for resDict in d1['site_residues']:
    ##             chid = resDict['chain_id']
    ##             resname = resDict['chem_comp_id']
    ##             resnum = str(resDict['author_residue_number'])
    ##             icode = resDict['author_insertion_code'] if resDict['author_insertion_code'] else '_'
    ##             res = self._receptor._ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
    ##             res0 = res.getHierView().iterResidues().next()
    ##             if res0.water:
    ##                 siteWaters.append(res0)
    ##             else:
    ##                 siteRes.append(res0)
    ##         #import pdb; pdb.set_trace()
    ##         chem_comp = d1['ligand_residues'][0]['chem_comp_id']
    ##         if ligandNames.has_key(chem_comp):
    ##             sites['%s %s'%(chem_comp,ligandNames[chem_comp])] = (siteRes, siteWaters)
    ##     return sites

    def getSites(self, pdbid, ligands):
        """use distance cutoff from ligand to define sites
        """
        sites = {}
        rads = self._receptor._ag.getRadii()
        ag = self._receptor._ag
        for lig in ligands:
            others = ag[list(set(range(len(ag)))-set(lig.getIndices()))]
            neigh = findNeighbors(lig, 4.0, others)
            atomInds = {}
            for a1,a2,d in neigh:
                a1i = a1.getIndex()
                a2i = a2.getIndex()
                if d < 1.2*(rads[a1i]+rads[a2i]):
                    atomInds[a2i] = True
            if len(atomInds):
                siteRes = []
                siteWaters = []
                for res in self._receptor._ag[atomInds.keys()].getHierView().iterResidues():
                    if res.water:
                        siteWaters.append(res)
                    else:
                        siteRes.append(res)
                rname = '%s:%s:%d'%(lig.getChid(), lig.getResname(),lig.getResnum())
                sites['%s %s'%(rname,chemComp[lig.getResname()]['name'])] = (siteRes, siteWaters)
        return sites

    def _gaps(self, polymers):
        # find gaps: for each gap we will have the residue before and after the gap
        # and the list of missing residues inbetween.
        # missing N- and C-termini are dealt with separately 
        #from time import time
        #t0 = time()
        ag = self._receptor._ag
        self._gapsPerChain = {}
        if polymers:
            for chain in ag.getHierView().iterChains():
                chid = chain.getChid()
                fasta = ''
                reslist = []
                sel = chain.select('protein')
                if sel is None: continue # chain can be only water
                for res in chain.select('protein').getHierView().iterResidues():
                    gaps = [] # list of [[pdbRes1, restyp, ... restyp, pdbRes2], ...]
                    rname = res.getResname()
                    fasta += AAnames.get(rname, '?')
                    reslist.append(res)

                aln = pairwise2.align.globalxx(polymers[chid].sequence, fasta, one_alignment_only=True)[0]
                #print "ALIGNMENT:"
                #print aln[0]
                #print aln[1]
                if aln[1].count('-'): # there are missing residues
                    i = 0
                    ri = i # index in pdb sequence
                    while aln[1][i]=='-': i+=1
                    missingNterm = i
                    j = 0
                    while aln[1][-j-1]=='-': j+=1
                    missingCterm = j
                    #ri = i # index in pdb sequence
                    gap = []
                    for i in range(missingNterm, len(aln[1])-missingCterm):
                        if aln[1][i]=='-':
                            if len(gap)==0:
                                gap = [reslist[ri-1]]
                            gap.append(aln[0][i])
                        else:
                            if len(gap):
                                gap.append(reslist[ri])
                                # check that there is a physical gap (i.e. 1jff GLY34-HIS61 covalent bond)
                                c = gap[0].select('name C')
                                n = gap[-1].select('name N')
                                cc = None
                                if n and c:
                                    for b in ag[c.getIndices()[0]].iterBonded():
                                        if b.getIndex()==n.getIndices()[0]: # C and N are bonded, i.e. no gap in structure
                                            cc = c.getCoords()[0]
                                            break
                                    cc = c.getCoords()[0]
                                    cn = n.getCoords()[0]
                                    cutoff = 12.25 # 3.5**2. i.e. vdw(C) + vdw(N)
                                else: ## either C or N was not found, so we look at CA-CA distance
                                    ca1 = gap[0].select('name CA')
                                    ca2 = gap[-1].select('name CA')
                                    if ca1 and ca2:
                                        cc = ca1.getCoords()[0]
                                        cn = ca2.getCoords()[0]
                                        cutoff = 16 # 4**2
                                    else: # jsut pick first atom in each of the 2 residues
                                        ca1 = gap[0].all[0]
                                        ca2 = gap[-1].all[0]
                                        cc = ca1.getCoords()[0]

                                if cc is not None:
                                    d2 = numpy.sum((cn-cc)*(cn-cc))
                                    if d2>cutoff:
                                        # get the gap residues name, num and icode
                                        gapData = [[gap[0].getResname(), gap[0].getResnum(), gap[0].getIcode()]]
                                        gapData.extend(gap[1:-1])
                                        gapData.append([gap[-1].getResname(), gap[-1].getResnum(), gap[-1].getIcode()])
                                        gaps.append(gapData)
                                gap = []
                            ri += 1
                    mrn = aln[0][:missingNterm] # sequence of missing Nterminus
                    mrc = aln[0][-missingCterm:] # sequence of missing Cterminus
                    if len(mrn):
                        rr = reslist[0]
                        mrn = [[rr.getResname(), rr.getResnum(), rr.getIcode()], mrn]
                    if len(mrc):
                        rr = reslist[-1]
                        mrc = [[rr.getResname(), rr.getResnum(), rr.getIcode()], mrc]
                else:
                    mrn = mrc = []
                    gaps = []
                if len(mrn) or len(gaps) or len(mrc):
                    self._gapsPerChain[chid] = [mrn, gaps, mrc]
        
    ## def _classify2(self, pdbid, defaultReceptor=None, defaultLigand=None, completeSidechains=True):
    ##     # pdbid is isued to lookup cofactors and ligands
    ##     # defaultReceptor has to be a dictionary specifying for each category flags ? operation ?
    ##     # defaultLigand prody selection string specifying a ligand
        
    ##     # This function creates a "category: data field for each atom
    ##     # The category label is a 5-character string that can be
    ##     #     'water' : for water molecules
    ##     #     'modif' : for atoms in modified residues
    ##     #     'cofac' : for atom sin cofactors
    ##     #     'metal' : for metal atoms
    ##     #     'addit' : for atoms in additives
    ##     #     'ligan' : for atoms in residues declared as ligand by PDB
    ##     #     'nucle' : for atom sin nucleic acids
    ##     #     'prote' : for atoms in amino acids
    ##     #     'other' : for anything else
    ##     #
    ##     # flags:
    ##     #     'receptor': include this atom in the receptor, mutually exlusive with ligand
    ##     #     'ligand'  : include this atom in the ligand, mutally exclusive with receptor
    ##     #     'missS' : for atoms in standard amino acids missing sidechain atoms
    ##     #     'missA' : for atoms in residues with missing atoms other thans missS
    ##     #     'hasAlt'  : for residues alternate locations
    ##     #     'unsupType':for residues with atoms not known in the AD4 forcefield

    ##     if defaultReceptor is None:
    ##         defaultReceptor = self.includeInReceptor
        
    ##     ag = self._receptor._ag
    ##     # set category to 'other for everythingt
    ##     ag.setData('category', ['other']*len(ag)) # create array of string len 5

    ##     # create the flags and set them to False by default
    ##     for flag in ['receptor', 'ligand', 'modified', 'altlocs', 'unsupType',
    ##                   'missS', 'missA']:
    ##         ag.setFlags(flag, [False]*len(ag))

    ##     # set category data
    ##     # NOTE: use '_wate' as category for water as water is a reserved keyword in prody
    ##     for cat in ['protein', 'nucleic', 'water']:
    ##         if cat=='water': cat5 = '_wate'
    ##         else: cat5 = cat[:5]
    ##         atoms = ag.select(cat)
    ##         if atoms:
    ##             atoms.setData('category', cat5)
    ##             if defaultReceptor[cat5]:
    ##                 atoms.setFlags('receptor', True)

    ##     # set cofactor category
    ##     cofactors = self.getCofactors(pdbid)
    ##     if cofactors:
    ##         for res in cofactors:
    ##             res.setData('category', 'cofac')
    ##             if defaultReceptor['cofac']:
    ##                 res.setFlags('receptor', True)

    ##     # set ligand category
    ##     ligands, ligandNames = self.getLigands(pdbid)
    ##     if ligands:
    ##         for l in ligands:
    ##             l.setData('category', 'ligan')
    ##             if defaultReceptor['ligan']:
    ##                 l.setFlags('receptor', False)

    ##     self._sites = self.getSites(pdbid, ligands)

    ##     # annotate metals and additives so that we will not consider them
    ##     # as ligands or building their sites
    ##     #atoms = ag.select('not protein and not water')
    ##     atoms = ag.select('category other')
    ##     if atoms:
    ##         hv1 = atoms.getHierView()
    ##         for res in hv1.iterResidues():
    ##             #resnum = str(res.getResnum())+res.getIcode() 
    ##             resname = res.getResname()
    ##             if self.isMetal(res):
    ##                 res.setData('category', 'metal')
    ##                 if defaultReceptor['metal']:
    ##                     res.setFlags('receptor', True)                    
    ##             else:
    ##                 if resname in self._additives:
    ##                     res.setData('category', 'addit')
    ##                     if defaultReceptor['addit']:
    ##                         res.setFlags('receptor', True)

    ##     # set unsupType flag for atoms not known in the AD4 forcefield
    ##     recElems = set(ag.getElements())
    ##     unknownElems = recElems - set(ADelements.keys())
    ##     if len(unknownElems):
    ##         atoms = ag.select('element %s'%' '.join(list(unknownElems)))
    ##         atoms.setFlag('unsupType', lem(atoms))

    ##     # handle modified residues
    ##     self._modifiedRes = {}
    ##     _polymers = {}
    ##     if len(self._receptor.pdbHeader)==0:
    ##         print 'WARNING: this PDB file has not header records'

    ##     if self._receptor.pdbHeader.has_key('polymers'):
    ##         for poly in self._receptor.pdbHeader['polymers']:
    ##             _polymers[poly.chid] = poly
    ##             if poly.modified is None: continue
    ##             for resname, resnum, baseResName, name in poly.modified:
    ##                 if resnum[-1].isalpha():
    ##                     icode = resnum[-1]
    ##                     resnum = resnum[:-1]
    ##                 else:
    ##                     icode = '_'

    ##                 atoms = ag.select('chid %s resnum `%s` icode %s'%(
    ##                     poly.chid, resnum, icode))
    ##                 if atoms: # residue in chain not present in this assembly
    ##                     atoms.setData('category', 'modif')
    ##                     supported = len(atoms.select('unsupType'))==0
    ##                     if defaultReceptor['modif'] and supported:
    ##                         atoms.setFlags('receptor', True)
    ##                         mutResName = None
    ##                         for r in atoms.getHierView().iterResidues():
    ##                             allKnown = len(set(atoms.getElements())- set(ADelements.keys()))==0
    ##                             self._modifiedRes[getResidueKeyAlt(r)] = [baseResName,allKnown,mutResName]
    ##                     elif not supported:
    ##                         pass
                        
    ##     ## at this points atoms have a category that is 'protein', 'nucleic', '_wate',
    ##     ## 'cofac', 'ligan', 'metal', 'addit' or 'modif'
    ##     ## now set dome of the flags

    ##     # set the ligand flag
    ##     if defaultLigand:
    ##         ligAtoms = mol._ag.select(defaultLigand)
    ##         if ligAtoms:
    ##             ligAtoms.setFlags('ligand', [True]*len(ligAtoms))

    ##     # set alternate location flag and chose deault alternate location
    ##     # get all atoms that have altloc other than ' '
    ##     alt = ag.select('not altloc _')
    ##     if alt:
    ##         altres = self._altlocsRes
    ##         for res in alt.getHierView().iterResidues():

    ##             key = getResidueKey(res)
    ##             # get all altloc characters for this residue
    ##             alts = list(set(res.getAltlocs()))

    ##             # we assume that the altloc characters are consistent within a residue
    ##             # get set of atoms for this residue for each altloc
    ##             cat = res.getData('category')[0]
    ##             resind = res.getResindices()[0]
    ##             resAtomsForAlt = []
    ##             nbat = []
    ##             occup = []
    ##             resAllAtoms = ag.select('resindex %d'%resind)
    ##             resAllAtoms.setFlags('altlocs', True)
    ##             resAllAtoms.setFlags('receptor', False)
    ##             # loop over residues for each altloc to find atoms, 
    ##             # and occupancies of atoms with altloc
    ##             for alt in alts:
    ##                 # complete the residue using atoms with and without altlocs
    ##                 atoms = resAllAtoms.select('altloc _ or altloc %s'%alt)
    ##                 resAtomsForAlt.append(atoms)
    ##                 nbat.append(len(atoms))
    ##                 _occup = []
    ##                 for a in atoms:
    ##                     if a.getAltloc()!=' ':
    ##                         _occup.append(a.getOccupancy())
    ##                 occup.append(max(_occup))

    ##             # find max number of atoms in res for each altloc
    ##             maxn = max(nbat)
    ##             # favor residues with more atoms (i.e more complete)
    ##             # break ties using occupancy values
    ##             # check how many have maxn atoms
    ##             if nbat.count(maxn)==1: # if only one select this one
    ##                 ind = nbat.index(maxn)
    ##                 if not res.water: # water is not included by default
    ##                     resAtomsForAlt[ind].setFlags('receptor', True)
    ##             else:
    ##                 # get all with maxn atoms
    ##                 inds = [i for i in range(len(nbat)) if nbat[i] == maxn] 
    ##                 # pick the one with maxn atoms and highest occupancy
    ##                 occmax = max([occup[j] for j in inds])
    ##                 defal = None
    ##                 for n, al in enumerate(alts):
    ##                     if occup[n]==occmax:
    ##                         defal = al
    ##                 ind = alts.index(defal)

    ##             # save the list of alternate location and occupancies and
    ##             # the index the once selected by default
    ##             altres[key] = [ind, alts, occup]
    ##             catego = res.getData('category')[0]
    ##             if defaultReceptor[catego]:
    ##                 resAtomsForAlt[ind].setFlags('receptor', True)

    ##     # handle missing atoms declared in the PDB header
    ##     mutator = AARotamerMutator()
    ##     if self._receptor.pdbHeader.has_key('missing_atoms'):
    ##         # loop over list os MISSING ATOMS records
    ##         for chid, resname, resnum, icode, missAtNames in self._receptor.pdbHeader['missing_atoms']:
    ##             if icode=='': icode = '_'
    ##             allres = ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
    ##             cat = allres.getData('category')[0]
    ##             if cat=='prote' and len(set(missAtNames) & set(['N', 'CA', 'C', 'O']))==0:
    ##                 # this is a standard amino acid and missing atoms are NOT in the backbone
    ##                 for res in allres.getHierView().iterResidues():
    ##                     if res.getResname() in mutator.residueNames:
    ##                         res.setFlags('missS', ['True']*len(res))
    ##                         if defaultReceptor['missS'] and defaultReceptor[cat]:
    ##                             res.setFlags('receptor', True)
    ##                             key = getResidueKey(res)
    ##                             self._missingAtoms[key] = missAtNames
    ##                             if completeSidechains:
    ##                                 self._resToMutate[key] = [resname,  resname]
    ##                     else:
    ##                         raise RuntimeError("ERROR: %s unmutatable amino acids with missing side chain atoms", res.getResname())
    ##                 else: # e.g. 1qj6 I:TYS63 is missing O3
    ##                     # missing side chain atoms that we canot rebuild (yet)
    ##                     res.setFlags('missA', ['True']*len(res))
    ##                     if defaultReceptor['missA'] and defaultReceptor[cat]:
    ##                         res.setFlags('receptor', True)
    ##                         self._missingAtoms[getResidueKey(res)] = missAtNames

    ##     # handle gaps
    ##     self._gaps(_polymers)

    ##     expType = self._receptor.pdbHeader.get('experiment', None)
    ##     hydrogens = ag.select('hydrogen')
    ##     if not expType:
    ##         self._addH = False
    ##     else:
    ##         if expType == 'SOLUTION NMR' and hydrogens:
    ##             self._addH = False
        
    def _classify(self, pdbid, defaultReceptor=None, defaultLigand=None, completeSidechains=True):
        # This function creates a "category: data field for each atom
        # The category label is a 5-character string that can be
        #     'water' : for water molecules
        #     'modif' : for atoms in modified residues
        #     'cmpsc' : for atoms in residues with missing sidechain atoms
        #               that can be reconstructed. Currentyl only amino acids
        #     'missa' : for atoms in residues with other missing atoms
        #     'cofac' : for atoms in cofactors
        #     'manch' : for atoms needing a charge sprecified by the user
        #     'addit' : for atoms in additives
        #     'ligan' : for atoms in ligands
        #     'nucle' : for atoms in nucleic acids
        #     'prote' : for atoms in amino acids
        #     'other' : for anything else
        #
        # in what will be the receptor molecule. The categories included by default are:
        # prote, nucle, modif, missa, cmpsc, cofac, unsre, unsat
        #
        # if alternate position are present an additional flag is created 
        # The function also creates a structure describing gaps (i.e sequences of missing residues) in the chains
        # unsat is set for atoms not parametrized int he AD4.1 forcefield
        # unsre is set for residues with atoms not parametrized int he AD4.1 forcefield

        if defaultReceptor is None:
            defaultReceptor = self.includeInReceptor

        ag = self._receptor._ag
        # set everythingto be category 'other'
        ag.setData('category', ['other']*len(ag)) # create array of string len 5
        # create a receptor and ligand flag both set to False for all atoms
        ag.setFlags('receptor', [False]*len(ag)) # will be true for atoms to keep in receptor
        ag.setFlags('ligand', [False]*len(ag)) # true for atoms to use as ligand molecule
        ag.setFlags('altlocs', [False]*len(ag)) # true for residues with atoms with alternate locations
        ag.setFlags('unsat', [False]*len(ag)) # for atoms with unknown altom types for DA4.1
        ag.setFlags('unsre', [False]*len(ag)) # for residues with unknown altom types for DA4.1
        
         # set unsupType flag for atoms not known in the AD4 forcefield
        recElems = set(ag.getElements())
        unknownElems = recElems - set(ADelements.keys())
        if len(unknownElems):
            atoms = ag.select('element %s'%' '.join(list(unknownElems)))
            atoms.setFlags('unsat', [True]*len(atoms))
            for res in atoms.getHierView().iterResidues():
                #ectend selection to entire residue
                resAtoms = ag.select("resindex %d"%res.getResindices()[0])
                resAtoms.setFlags('unsre', [True]*len(atoms))

        atoms = ag.select('water')
        if atoms:
            atoms.setData('category', '_wate') # because water is a reserved keyword in prody
            if defaultReceptor['_wate']:
                atoms.setFlags('receptor', True)

        # tag protein
        atoms = ag.select('protein')
        if atoms:
            atoms.setData('category', 'prote')
            if defaultReceptor['prote']:
                atoms.setFlags('receptor', True)

        # tag nucleic acids
        atoms = ag.select('nucleic')
        if atoms:
            atoms.setData('category', 'nucle')
            if defaultReceptor['nucle']:
                atoms.setFlags('receptor', True)

        # set the ligand flag
        if defaultLigand:
            ligAtoms = ag.select(defaultLigand)
            if ligAtoms:
                #ligAtoms.setFlags('ligand', [True]*len(ligAtoms))
                ligAtoms.setFlags('ligand', True)

        # get list of cofactors
        #if have_internet():
        cofactors = self.getCofactors(pdbid)
        if cofactors:
            for res in cofactors:
                res.setData('category', 'cofac')
                if defaultReceptor['cofac']:
                    res.setFlags('receptor', True)

        # annotate manual charges and additives so that we will not consider them
        # as ligands or building their sites
        #atoms = ag.select('not protein and not water')
        atoms = ag.select('category other')
        if atoms:
            hv1 = atoms.getHierView()
            for res in hv1.iterResidues():
                #resnum = str(res.getResnum())+res.getIcode() 
                resname = res.getResname()
                if self.isManualCharge(res):
                    res.setData('category', 'manch')
                    if defaultReceptor['manch']:
                        res.setFlags('receptor', True)                    
                else:
                    if resname in self._additives:
                        res.setData('category', 'addit')
                        if defaultReceptor['addit']:
                            res.setFlags('receptor', True)

        mutator = AARotamerMutator()
        # tag modified residues
        _polymers = {}
        if len(self._receptor.pdbHeader)==0:
            print 'WARNING: this PDB file has not header records'    
        if self._receptor.pdbHeader.has_key('polymers'):
            for poly in self._receptor.pdbHeader['polymers']:
                _polymers[poly.chid] = poly
                if poly.modified is None: continue
                for resname, resnum, baseResName, name in poly.modified:
                    if resnum[-1].isalpha():
                        icode = resnum[-1]
                        resnum = resnum[:-1]
                    else:
                        icode = '_'

                    atoms = ag.select('chid %s resnum `%s` icode %s'%(
                        poly.chid, resnum, icode))
                    #print "SETTING category modif to res", poly.chid, resname, resnum
                    if atoms:
                        atoms.setData('category', 'modif')
                        for res in atoms.getHierView().iterResidues():
                            # need at least N CA C or N CA CB in order to perform mutation
                            # e.g. 1gbb P:B2A or 1p02 P:B2A:1
                            #import pdb; pdb.set_trace()
                            set1 = res.select('name N CA C and protein')
                            set2 = res.select('name N CA CB and protein')
                            if (set1 and len(set1)==3) or (set2 and len(set2)==3):
                                # if all atoms are known we can keep it as is
                                key = getResidueKey(res)
                                if not res.getFlags('unsre')[0]:
                                    if defaultReceptor['modif']:
                                        #res.setFlags('receptor', [True]*len(res))
                                        res.setFlags('receptor', True)
                                elif baseResName in mutator.residueNames:
                                    # contains unsupported atoms but can be mutated to something we know
                                    if defaultReceptor['modif']:
                                        res.setFlags('receptor', [True]*len(res)) # keep it
                                        if completeSidechains:
                                            self._resToMutate[key] = [resname,  baseResName] # set it to be mutated
                                else: # unsupported atoms and can't be mutated, so it can't be included
                                    res.setFlags('receptor', [False]*len(res))

        # tag res with missing atoms
        if self._receptor.pdbHeader.has_key('missing_atoms'):
            # loop over list os MISSING ATOMS records
            for chid, resname, resnum, icode, missAtNames in self._receptor.pdbHeader['missing_atoms']:
                if icode=='': icode = '_'
                allres = ag.select('chid %s resnum `%s` icode %s'%(chid, resnum, icode))
                if len(set(missAtNames) & set(['N', 'CA', 'C', 'O'])): # missing backbone atoms
                    # we can include as is (default) or remove
                    for res in allres.getHierView().iterResidues():
                        res.setData('category', ['missa']*len(res))
                        if defaultReceptor['missa']:
                            res.setFlags('receptor', True)
                            self._missingAtoms[getResidueKey(res)]=missAtNames

                else: # NOT missing backbone atoms
                    for res in allres.getHierView().iterResidues():
                        #print "SETTING category cmpsc to res", chid, resname, resnum
                        if res.getResname() in mutator.residueNames:
                            res.setData('category', ['cmpsc']*len(res))
                            res.setData('category', ['cmpsc']*len(res))
                            if defaultReceptor['cmpsc']:
                                res.setFlags('receptor', True)
                                key = getResidueKey(res)
                                self._missingAtoms[key] = missAtNames
                                if completeSidechains:
                                    self._resToMutate[key] = [resname,  resname]
                                
                        else: # e.g. 1qj6 I:TYS63 is missing O3
                            # missing side chain atoms that we canot rebuild (yet)
                            res.setData('category', ['missa']*len(res))
                            if defaultReceptor['missa']:
                                res.setFlags('receptor', True)
                                self._missingAtoms[getResidueKey(res)] = missAtNames

        # additional rough check for missing atoms
        # e.g. cryoEM structures that have only CA but no MISSING ATOMS records
        for res in self._receptor._ag.iterResidues():
            resname = res.getResname()
            if resData.has_key(resname):
                nbAtoms, nbHeavy, pepLink, smilessrc, smiles, name, fch, arom, leav, chiral = resData[resname]
                heavy = res.select('not hydrogen')
                if res.getResname() not in mutator.residueNames: continue
                if len(heavy) < nbHeavy:
                    if res.getData('category')[0] in ['missa', 'cmpsc']:
                        continue
                    # do not include by default residues with no bonds
                    # and missing atoms (e.g. 1eg0.pdb)
                    if len(PSelToMKSel(heavy).getBonds()[1])==0: # no neighbors
                        res.setData('category', ['missa']*len(res))
                        res.setFlags('receptor', False)
                        self._missingAtoms[getResidueKey(res)] = ['?']
                        continue
                    if pepLink: #res.protein:
                        bbat = res.select('name N CA C O')
                        if bbat and len(bbat)==4:
                            if bbat.bb is None:
                                # force backbone flag
                                self._receptor._ag._flags['bb'][bbat.getIndices()] = True
                            if len(heavy)+1==nbHeavy: # missing OXT
                                continue
                            res.setData('category', ['cmpsc']*len(res))
                            if defaultReceptor['cmpsc']:
                                res.setFlags('receptor', True)
                                key = getResidueKey(res)
                                self._missingAtoms[key] = ['?']
                                if completeSidechains:
                                    self._resToMutate[key] = [resname,  resname]
                        else:
                            res.setData('category', ['missa']*len(res))
                            if defaultReceptor['missa']:
                                res.setFlags('receptor', True)
                                self._missingAtoms[getResidueKey(res)] = ['?']
                    else:
                        res.setData('category', ['missa']*len(res))
                        if defaultReceptor['missa']:
                            res.setFlags('receptor', True)
                            self._missingAtoms[getResidueKey(res)] = ['?']
        
        ligands, ligandNames = self.getLigands(pdbid)
        if ligands:
            for l in ligands:
                l.setData('category', 'ligan')
                if defaultReceptor['ligan']:
                    l.setFlags('receptor', True)
            
        # handle gaps
        self._gaps(_polymers)
        
        #self._sites = self.getSites(pdbid, ligandNames)
        self._sites = self.getSites(pdbid, ligands)
        
        expType = self._receptor.pdbHeader.get('experiment', None)
        hydrogens = ag.select('hydrogen')
        if not expType:
            self._addH = False
        else:
            if expType == 'SOLUTION NMR' and hydrogens:
                self._addH = False

        # by default we only select atoms wit alternate location ' ' or 'A'
        #alt = ag.select('not altloc' _' and not altloc A)
        #if alt:
        #    alt.setFlags('receptor', False)
        # code above does not work because we canot assume first altloc is
        # always A e.g. 4a1x:B:CYS52 altlocs B C
        # we nee to loop over all residues in each chain and select an alternate
        # location tag for this residue
        # get all atoms that have altloc other than ' '
        alt = ag.select('not altloc _')
        if alt:
            altres = self._altlocsRes
            for chain in alt.getHierView().iterChains():
                for res in chain.iterResidues():
                    key = getResidueKey(res)

                    # get all altloc characters for this residue
                    alts = list(set(res.getAltlocs()))
                    
                    # here we assume that the altloc character is consistent within a residue
                    # now get set of atoms for this residue for each altloc
                    resind = res.getResindices()[0]
                    resAtomsForAlt = []
                    nbat = []
                    occup = []
                    resAllAtoms = ag.select('resindex %d'%resind)
                    resAllAtoms.setFlags('altlocs', True)
                    resAllAtoms.setFlags('receptor', False)
                    # loop over residues for each altloc to find atoms, and occupancies
                    # of atoms with altloc
                    for alt in alts:
                        atoms = resAllAtoms.select('altloc _ or altloc %s'%alt)
                        resAtomsForAlt.append(atoms)
                        nbat.append(len(atoms))
                        _occup = []
                        for a in atoms:
                            if a.getAltloc()!=' ':
                                _occup.append(a.getOccupancy())
                        occup.append(max(_occup))

                    # find max number of atoms in res for each altloc
                    maxn = max(nbat)

                    # check how many have maxn atoms
                    if nbat.count(maxn)==1: # if only one select this one
                        ind = nbat.index(maxn)
                        if not res.water: # water is not included by default
                            resAtomsForAlt[ind].setFlags('receptor', True)
                    else:
                        # get all with maxn atoms
                        inds = [i for i in range(len(nbat)) if nbat[i] == maxn] 
                        # pick the one with maxn atoms and highest occupancy
                        occmax = max([occup[j] for j in inds])
                        defal = None
                        for n, al in enumerate(alts):
                            if occup[n]==occmax:
                                defal = al
                        ind = alts.index(defal)

                    altres[key] = [ind, alts, occup]
                    catego = res.getData('category')[0]
                    if defaultReceptor[catego]:
                        resAtomsForAlt[ind].setFlags('receptor', True)
                        
    def getSummary(self):
        ag = self._receptor._ag
        recProt = ag.select('category prote and receptor')
        if not recProt:
            recProt = []
        notrecProt = ag.select('category prote and not receptor')
        if not notrecProt:
            notrecProt = []

        def resname(res):
            return '%s:%s%d%s'%(res.getChid(), res.getResname(), res.getResnum(), res.getIcode())

        cofactors = ag.select('category cofac and receptor')
        cofactorshv = cofactors.getHierView().iterResidues() if cofactors else []

        ligands = ag.select('category ligan and receptor')
        ligandshv = ligands.getHierView().iterResidues() if ligands else []

        additives = ag.select('category addit and receptor')
        additiveshv = additives.getHierView().iterResidues() if additives else []
        manuals = ag.select('category manch and receptor')
        manualshv = metals.getHierView().iterResidues() if manuals else []

        modified = ag.select('category modif and receptor')
        modifiedhv = modified.getHierView().iterResidues() if modified else []

        missingAtoms = ag.select('category missi and receptor')
        missingAtomshv = missingAtoms.getHierView().iterResidues() if missingAtoms else []

        #water = ag.select('water and receptor')
        #waterhv = water.getHierView().iterResidues() if water else []
        
        other = ag.select('category other and receptor')
        otherhv = other.getHierView().iterResidues() if other else []

        unsup = ag.select('category unsre and receptor')
        unsuphv = unsup.getHierView().iterResidues() if unsup else []
        
	altloc = ag.select('not altloc _ and receptor')
        altlochv = altloc.getHierView().iterResidues() if altloc else []

        lines = ["PDB entry: %s\n"%self._pdbid,
                 " Included in receptor:\n"
                 "  protein  : %d atoms in chains %s\n"%(len(recProt), list(set(recProt.getChids()) if recProt else "")),
                 "  cofactors: %s\n"%[resname(x) for x in cofactorshv if x.select('receptor')],
                 "  ligands  : %s\n"%[resname(x) for x in ligandshv if x.select('receptor')], 
                 "  manCharge: %s\n"%[resname(x) for x in manualshv if x.select('receptor')], 
                 "  modified : %s\n"%[resname(x) for x in modifiedhv if x.select('receptor')], 
                 "  additives: %s\n"%[resname(x) for x in additiveshv if x.select('receptor')], 
#                 "  water    : %s\n"%[resname(x) for x in waterhv if x.select('receptor')], 
                 "  unsupport: %s\n"%[resname(x) for x in unsuphv if x.select('receptor')], 
#                 "  altlocs  : %s\n"%[resname(x) for x in altlochv if x.select('receptor')], 
                 "  other    : %s\n"%[resname(x) for x in otherhv if x.select('receptor')], 
                 ]
        lines.append("  defAlt:\n")
        n = 0
        s = ''
        for k,v in self._altlocsRes.items():
            ind, alts, occ = v
            s += '%s:(%c %.2f), '%(k, alts[ind], occ[ind])
            n += 1
            if n==5:
                s+='\n'
                lines.append(s)
                s = ''
                n = 0
        if len(s):
            s+='\n'
            lines.append(s)
            
        cofactors = ag.select('category cofac and not receptor')
        cofactorshv = cofactors.getHierView().iterResidues() if cofactors else []

        ligands = ag.select('category ligan and not receptor')
        ligandshv = ligands.getHierView().iterResidues() if ligands else []

        additives = ag.select('category addit and not receptor')
        additiveshv = additives.getHierView().iterResidues() if additives else []

        manuals = ag.select('category mandch and not receptor')
        manualshv = metals.getHierView().iterResidues() if manuals else []

        modified = ag.select('category modif and not receptor')
        modifiedhv = modified.getHierView().iterResidues() if modified else []

        missingAtoms = ag.select('category missi and not receptor')
        missingAtomshv = missingAtoms.getHierView().iterResidues() if missingAtoms else []

#        water = ag.select('water and not recegotoptor')
#        waterhv = water.getHierView().iterResidues() if water else []
        
        other = ag.select('category other and not receptor')
        otherhv = other.getHierView().iterResidues() if other else []

        unsup = ag.select('category unsre and not receptor')
        unsuphv = unsup.getHierView().iterResidues() if unsup else []
        lines.extend([
                 " Excluded from receptor:\n"
                 "  protein  : %d atoms in chains %s\n"%(len(notrecProt), list(set(notrecProt.getChids()) if notrecProt else "")),
                 "  cofactors: %s\n"%[resname(x) for x in cofactorshv if x.select('not receptor')],
                 "  ligands  : %s\n"%[resname(x) for x in ligandshv if x.select('not receptor')], 
                 "  manCharge: %s\n"%[resname(x) for x in manualshv if x.select('not receptor')], 
                 "  modified : %s\n"%[resname(x) for x in modifiedhv if x.select('not receptor')], 
                 "  additives: %s\n"%[resname(x) for x in additiveshv if x.select('not receptor')], 
#                 "  water    : %s\n"%[resname(x) for x in waterhv if x.select('not receptor')], 
                 "  unsupport: %s\n"%[resname(x) for x in unsuphv if x.select('not receptor')], 
                 "  other    : %s\n"%[resname(x) for x in otherhv if x.select('not receptor')], 
            ])
        lines.append('Sites:\n')
        for k,v in self._sites.items():
            lines.append("  %s : %d residues %s\n"%(k, len(v[0]), [resname(x) for x in v[0]]))
            lines.append("       %d water(s) %s\n"%(len(v[1]), [resname(x) for x in v[1]]))

        def resname(val):
            return '%s:%s%d%s'%(val[:4])

        lines.append('missing residues:\n    ')
        if self._receptor.pdbHeader.has_key('missing_residues'):
            lines.append(' %s\n'%' '.join([resname(x) for x in self._receptor.pdbHeader['missing_residues']]))
        else:
            lines.append('\n')
        lines.append('missing atoms:\n')
        for k,v in self._missingAtoms.items():
            lines.append('  %s: %s'%(k, ' '.join(v)))
        lines.append('\n')

        return lines

    def buildSelStr(self, atoms, exclude=[]):
        # build a selection string for atoms exluding the ones specified in exclude
        # exclude is a list of keys chid:resnum[Icodealtloc] where Icode and altloc are optional
        if atoms is None:
            return None
        
        excludeKeys = []
        for chid, rnum, icode in exclude:
            excludeKeys.append('%s:%d%s'%(chid, rnum, icode))
        skip = set(excludeKeys)

        # get all alternate location
        altlocs = set(atoms.getAltlocs())
        altlocres = defaultdict(list)
        selstr = ""
        for res in atoms.getHierView().iterResidues():
            chid = res.getChid()
            key = '%s:%d%s'%(chid, res.getResnum(), res.getIcode())
            icode = res.getIcode()
            if icode=='': icode = '_'
            key += icode
            alt = list(set(res.getAltlocs()))
            # replace ' ' with '_' in alt and skip the ones in exclude
            _alt = []
            for a in alt:
                if key+a in skip: continue
                if a==' ': _alt.append('_')
                else: _alt.append(a)

            selstr += '(chid %s resnum `%d` icode %c altloc %s) or'%(
                chid, res.getResnum(), icode, ' '.join(_alt))
 
        assert len(self._currentBM.select(selstr[:-3]))%len(atoms)==0
        return selstr[:-3] # remove last 'or '
        
    def getProtToExcludeString(self):
        """create the prody selection string for residues that are protein
        but are not selected to be part of the receptor
        """
        return self.buildSelStr(self._currentBM.select('(category prote cmpsc missa) and not receptor'))
    
    def getHeteroToKeepString(self, exclude=[]):
        """create the prody selection string for residues that are hetero
        but are not selected to be part of the receptor
        """
        return self.buildSelStr(self._currentBM.select('(not category prote cmpsc missa) and receptor'), exclude)

    ## def handleMSE(self):
    ##     # handle MSE residues by replacing Se by S
    ##     mse = self._currentBM.select('resname MSE')
    ##     if mse:
    ##         sele = mse.select('element SE')
    ##         sele.setElements(['S']*len(sele))
    ##         sele.setNames(['S']*len(sele))
    ##         mse.setResnames(['MET']*len(mse))
    ##     self._MSEtoMET = mse

    def getPDBQTHeaderLines(self):
        from Support import version
        lines = ["REMARK 850 AGFR version %s"%(version.__version__),
                 "REMARK 850 AGFR source %s.pdb"% self._receptor.name]

        # biomolecule
        #lines.append("REMARK 850 AGFR biomol %s:%s%s %s->%s"%(chid, resnum, icode, resname, mresname))

        # mutation
        if len(self._resToMutate):
            for key, resnames in self._resToMutate.items():
                resname, mresname = resnames
                if resname != mresname:
                    lines.append("REMARK 850 AGFR mutation %s -> %s"%(key, mresname))

        # completed sidechains
        if len(self._resToMutate):
            for key, resnames in self._resToMutate.items():
                resname, mresname = resnames
                if resname == mresname:
                    lines.append("REMARK 850 AGFR side chain completed %s"%(key))

        # selected alternate location
        bm = self._currentBM
        for key,values in self._altlocsRes.items():
            # check that this residue is in rececptor
            chid, resname, resnum, icode = key.split(':')
            if bm.select('chid %s resnum `%s` icode %s and receptor'%(chid, resnum, icode)) is None: continue
            ind, alts, occ = values
            lines.append('REMARK 850 AGFR altlocs %s:%s:%s:%s @%s (%.2f)'%(
                chid, resname, resnum, icode, alts[ind], occ[ind]))

        # protonation
        lines.append("REMARK 850 AGFR protonated using %s"%'reduce')
        lines.append("REMARK 850 AGFR partial charge assigned using %s"%'OpenBabel Gasteiger model')
        # charges
        manuals = self._receptor.select('category manch and receptor')
        if manuals:
            for res in manuals.getHierView().iterResidues():
                lines.append("REMARK 850 AGFR custom charge for %s:%s:%s:%s = %.4f"%(
                    res.getChid(), res.getResname(), res.getResnum(), res.getIcode(), res.getCharges()[0]))
        lines.append("REMARK 850 AGFR end")
        return lines


def protonate(atoms, molName):
    """execute reduce to protonate a set of Prody Atoms
    handle segments by temporarily renaming chains to have unique chain names
    The original chain bnames are restored after reduce added hydrogen atoms
    """
    heavy = atoms.select("not hydrogen")

    # since reduce does not like chains with the same name in different segments
    # and creates unbound hydrogen atoms we make sure each chain has its own name
    #allsegs = numpy.unique(added[comp]['_contactingSegments']+[v['pepSeg']])
    allsegs = numpy.unique(heavy.getSegnames())
    taken = {}
    if len(allsegs)>1:
        chids = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
                 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
                 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

        for chid in numpy.unique(atoms.getChids()):
            chids.remove(chid)

        for chain in heavy.getHierView().iterChains():
            chid = chain.getChid()
            seg = chain.getSegname()
            if taken.has_key(chid):
                #we need to rename this chain
                if len(chids)==0:
                    # we ran out of chain ids to rename e.e. 5vf1.pdb 6 chains and 24 BIOMAT for 144-mer
                    print 'ERROR: protonate, too many chains'
                    return None
                newchid = chids.pop(0)
                chain.setChids(newchid)
                taken[newchid] = "%c %c"%(seg, chid)
            else:
                taken[chid] = "%c %c"%(seg, chid)

    ## run reduce as a sbprocess to avoid having to write a file
    #reducePath = os.path.join(os.path.split(sys.executable)[0], 'reduce')
    from ADFR.utils.MakeGrids import findBinary
    reducePath = findBinary("reduce")
    #print "REDUCE PATH:", reducePath
    system_info = platform.uname()
    if system_info[0] == 'Windows':
        shell=False
    else:
        shell=True

    try:
        #import pdb; pdb.set_trace()
        proc = subprocess.Popen('%s -build -DB %s_wwPDB_het_dict.txt -'%(reducePath,reducePath),
                                stdin=subprocess.PIPE , 
                                stdout=subprocess.PIPE , 
                                stderr=subprocess.PIPE, 
                                bufsize=1, shell=shell)
        prody.writePDBStream(proc.stdin, atoms.select('not hydrogen'))
    except IOError:
        print 'ERROR: cmd %s failed'%(reducePath+' -build -')
        return None

    (stdout, stderr) = proc.communicate()

    # reduce command failed
    if stdout[:16]!="USER  MOD reduce":
        print "ERROR: output of reduce is wrong"
        return None

    #out = stdout.split('\n')
    #for l in out: print l
    agH = prody.parsePDBStream(StringIO(stdout.split('\n')))
    molh = Molecule(molName, agH)
    heavyH = molh._ag.select("not hydrogen")
    assert len(heavy)==len(heavyH)
    hat = molh._ag.select("hydrogen")
    if hat is None: # reduce failed to protonate e.g. 6ann.pdb
        print "ERROR: no hydrogen atoms after executing reduce"
        return None
    numBonds = [x.numBonds() for x in hat]
    if min(numBonds)==0:
        raise RuntimeError, '%d H atoms with no bond after reduce'%syst

    # restore CHIDs
    ## v = taken
    ## atomsets = {}
    ## if len(v):
    ##     # build atoms sets for each chain before resetting any chid to
    ##     # avoid affect selection of atoms later no
    ##     for newchid, old in v.items():
    ##         atomsets[newchid] = molh._ag.select('chid %s'%newchid)
    ##     for newchid, old in v.items():
    ##         oldSeg, oldChid = old.split()
    ##         _atoms = atomsets[newchid]
    ##         _atoms.setSegnames([oldSeg]*len(atoms))
    ##         _atoms.setChids([oldChid]*len(atoms))

    return molh

## class pdbToBioMolPDBQT:
##     """
##     class to generate PDBQT for the BioMolecule of selected parts
##     The bio molecule will be built using the biomoltrans entry in mol.pdbHeader
##     if the header is missing the symetric unit will be used
##     """

##     def __init__(self):
##         self._fromPDB = None # willbe a ProDy molecuel read from PDB with header
##         self._bm = None # will be biomolecule
##         self._proteinToRemove = None
##         self._heteroToKeep = None
##         self._addH = 'auto'
##         self._receptor = None # will be the BM molecules that can be converted
##                               # to PDBQT
##         self._recPDBQT = None

##         # atom sets for various categories
##         self._protein = None
##         self._cofactors = None
##         self._ligands = None
##         self._additives = None
##         self._siteWaters = None
##         self._metals = None
##         self._waters = None
##         self._other = None
##         # a flag indicating whether we need to build the biomolecule or not
        
##     def setPDBfromFile(self, filename, proteinToRemove=[], heteroToKeep=[],
##                        addH='auto'):
##         self._proteinToRemove = proteinToRemove
##         self._heteroToKeep = heteroToKeep
##         self._addH = addH
##         self._fromPDB = Read(filename, header=True)
##         self.buildReceptor()


##     def setPDBfromMolecule(self, mol, bioMol, proteinToRemove=None, heteroToKeep=None,
##                            addH='auto', pdbqtFileName=None, resToMutate=None, charges={},
##                            pdbqtHeaderLines=[]):
##         assert isinstance(mol, Molecule)
##         self._proteinToRemove = proteinToRemove
##         self._heteroToKeep = heteroToKeep
##         #print "In recFromPdb.setPDBfromMolecule():"
##         #print "proteinToRemove", self._proteinToRemove
##         #print "heteroToKeep:", self._heteroToKeep
##         #print "mutateRes:", mutateRes
##         self._addH = addH
##         self._fromPDB = mol
##         self._resToMutate = resToMutate
##         self._pdbqtFileName = pdbqtFileName
##         self._bm = bioMol._ag
##         self._charges = charges
##         self._pdbqtHeaderLines = pdbqtHeaderLines
##         self.processReceptor()

##     def _getPDBQTatoms(self):
##         mol = self._receptor 
##         mol._ag.setCharges([0.]*len(mol._ag))
##         for resStr, ch in self._charges.items():
##             chid, resname, resnum, icode = resStr.split(":")
##             if chid == ' ': chid = '_'
##             if icode:
##                 mtl = mol._ag.select('chid %s resnum `%s` icode %s'%(chid,resnum,icode))
##             else:
##                 mtl = mol._ag.select('chid %s resnum `%s`'%(chid,resnum))
##             if mtl:
##                 mtl.setCharges([ch])
 
##         okAtoms = mol.select('element %s'%' '.join(OBGastKnownElements))
##         obmol = ProdyToOBMol(okAtoms)

##         charger = ob.OBChargeModel.FindType('gasteiger')
##         success = charger.ComputeCharges(obmol)
##         if not success:
##             raise IOError("Failed to assign gasteiger charges")

##         totalCharge = obmol.GetTotalCharge ()
##         # copy charges from obmol to prody molecule
##         inds = []
##         charges = []
##         for obatom in ob.OBMolAtomIter(obmol):
##             obIndex = obatom.GetIndex()
##             inds.append(obmol.obToProdyIndex[obIndex])
##             charges.append(obatom.GetPartialCharge())
##             patom = mol._ag[obmol.obToProdyIndex[obIndex]]
##             res = obatom.GetResidue()
##             #aIndex = obatom.GetIdx()-1
##             #print '%s%d:%s %s%d:%s %.2f'%(res.GetName(), res.GetNum(), res.GetAtomID(obatom),
##             #                             patom.getResname(), patom.getResnum(), patom.getName(),
##             #                             charges[-1])
##             assert patom.getIndex() == inds[-1]
##             assert res.GetAtomID(obatom)==patom.getName()

##         # set the Gasteiger charges for these atoms
##         mol._ag._data['charge'][inds] = charges

##         # print nonIntegral charge
##         for res in mol._ag.iterResidues():
##             resCharge = numpy.sum(res.getData('charge'))
##             if abs(resCharge-int(resCharge)) > 0.1:
##                 print 'Non integral charge' , res, resCharge
##             if abs(resCharge) > 0.1:
##                 print 'Charged residue' , res, resCharge

##         #print 'totalCharge:', numpy.sum(mol._ag._data['charge'])
##         setADElements(mol._ag.all)
##         mol.buildBondsByDistance()
##         nphind = mergeNonPolarHydrogenCharges(mol._ag.all)
##         mol._ag.setFlags("nph", [False]*len(mol._ag))
##         mol._ag._flags["nph"][nphind] = True

##         recAtoms =  mol._ag.select("not nph")
##         # write PDBQT file
##         # put remarks for provenance data
##         lines = self._pdbqtHeaderLines
        
##         #import pdb; pdb.set_trace()
##         lines.extend(PDBQTlinesForAtoms(recAtoms))
##         if not self._pdbqtFileName:
##             self._pdbqtFileName = '%s_rec.pdbqt'%self._fromPDB.name
##         f = open(self._pdbqtFileName,  'w')
##         [f.write('%s\n'%l) for l in lines]
##         f.close()

##     def addAG(self, ag1, ag2):
##         # returns ag1 + ag2
##         if not ag1: return ag2
##         if not ag2: return ag1
##         return Selection(self._bm, list(set(ag1._getIndices()) |
##                                             set(ag2._getIndices())), '')
    
##     def subAG(self, ag1, ag2):
##         # returns ag1 - ag2
##         if not ag1:
##             return Selection(self._bm, [], "")
##         if not ag2: return ag1
##         return Selection(self._bm, list(set(ag1._getIndices()) -
##                                             set(ag2._getIndices())), '')
        
##     def processReceptor(self):

##         # mutate residues
##         if self._resToMutate:
##             mol = self.mutateResidues(self._bm)

##         # select atoms
##         ag1 = self._bm.select('protein')
##         if self._proteinToRemove:
##             ag2 = self._bm.select(self._proteinToRemove)
##             ag = self.subAG(ag1, ag2)
##         else:
##             ag = ag1
##         if self._heteroToKeep:
##             ag3 = self._bm.select(self._heteroToKeep)
##             atoms = self.addAG(ag, ag3)
##         else:
##             atoms = ag
        
##         atoms = atoms.select('not deleted')
##         # reset all alternate locations (else reduce will drop the one that are not ' ' or 'A')
##         atoms.setData('altloc', [' ']*len(atoms))

##         if atoms:
##             from MolKit2.selection import Selection
##             sel = Selection(self._bm, atoms.getIndices(), '')
##             ag = sel.toAtomGroup('receptor')

##         # hydrogen atoms
##         hydrogens = ag.select('hydrogen')
##         if hydrogens is not None and len(hydrogens) == 0: hydrogens = None
        
##         ## metals = ag.select('category metal and receptor')
##         ## if metals:
##         ##      self._metals = [(res.getChid(), res.getResname(), res.getResnum(), res.getIcode()) for res in metals.getHierView().iterResidues() ]
##         ## else:
##         ##     self._metals = []
##         ## print "METALS:", self._metals
##         #import pdb; pdb.set_trace()

##         # as input to reduce and reads PDB records back from reduce
##         if self._addH is True or (self._addH is 'auto' and hydrogens is None):
##             recH = protonate(ag.select('not hydrogen'), 'rec_H')
##             if recH is None:
##                 print 'WARNING: reduced failed to protonate receptor'
##                 if hydrogens is not None and len(hydrogens) > 0:
##                     # the receptor has hydrogen atoms, we proceed whit what is there
##                     recH = Molecule('rec_H', ag)
##                 else:
##                     raise RuntimeError, 'ERROR: receptor has no hydrogen atoms, consider turning off adding hydrogen atoms'                    
##         else:
##             recH = Molecule('rec_H', ag)

##         self._receptor = recH
##         self._getPDBQTatoms()
##         #stream = pdbSrcStream(self._recPDBQTstr.split('\n'))
##         #self.rec_pdbqt = Molecule(self._fromPDB.name, prody.parsePDBStream(stream, format='PDBQT'))
##         #return self.rec_pdbqt
               
##     def _printModResAtoms(self, obj, msgstr=""):
##         if self._resToMutate:
##             for key, resnames in self._resToMutate.items():
##                 chid, resname, resnum, icode = key.split(':')
##                 resname, mutrtype = resnames
##                 if icode=='': icode = '_'
##                 res = obj.select("chid %s resnum `%d` icode %s altloc %c" % (chid, resnum, icode, alt))
##                 if res:
##                     print msgstr, chid, resnum
##                     print res.getResnames(), len(res.getResnames())
##                     print res.getNames()
##                     print res.getCoords()

##     def mutateResidues(self, ag):
##         #from MolKit2.AARotamer import AARotamer, AARotamerMutator
##         mutator = AARotamerMutator()
##         if hasattr (ag, "_molecule"):
##             mol = ag._molecule()
##         else:
##             mol = Molecule('tmp', ag)
##         mol._ag.delData('anisou')
##         # loop over residues tages to be mutated
##         chids = set(ag.getChids())
##         for key, resnames in self._resToMutate.items():
##             chid, resname, resnum, icode,  = key.split(':')
##             oldres, mutrtype = resnames
##             # multiple residues might be selected when BM has segments
##             if chid not in chids: continue # skip mutations outside this BM
##             atoms = mol.select("chid %s resnum `%s` icode %s and receptor" % (chid, resnum, icode))
##             if atoms is None or len(atoms) == 0:
##                 print "ERROR: In mutateResidues: residue %s not found" % (key)
##                 raise RuntimeError
##                 continue
            
##             ## loops to mutate the residue (possibly in multiple segments)
##             ## to a new residue type or to same residue type to fix missing atoms
##             for res in atoms.getHierView().iterResidues():
##                 segment = res.getSegname()
##                 ## print 
##                 if res.getResname()==mutrtype:
##                     print "COMPLETE sidechain res /%s/%s"%(segment, key)
##                 else:
##                     print "MUTATE residue /%s/%s from %s to %s"%(
##                         segment, key, res.getResname(), mutrtype)

##                 resbb = res.select("name N CA C O")
##                 ressc = res.select("not (name N CA C O)")
##                 # REMARK .. original CHI angles should be used to bias rotamer selection
##                 # but currently are only use to build side chains and will be over written by best rotamer
##                 # which can have very different angles
##                 origChi = mol.measureCHIsInRes(res, resname=mutrtype)[0]
##                 angles = mutator.mutateRes(mol, res, chid, int(resnum), mutrtype, origChi,
##                                            bbatoms=resbb, scatoms=ressc,icode=icode)
##                 mutator._addedAtoms.setFlags('receptor', [True]*len(mutator._addedAtoms))

##                 ## find best rotamer
##                 ## fisrt select the rebuilt residue
##                 if mutrtype=='ALA':
##                     continue

##                 # re-select to include rebuilt atoms
##                 allres = mol.select("segment %s chid %s resnum `%s` icode %s and receptor" % (
##                     segment, chid, resnum, icode))

##                 for res in allres.getHierView().iterResidues():
##                     seg = res.getSegname()

##                     ## build a rotamer object for this residue
##                     rotamer = AARotamer(res, mutator.rotamer.angleDef,
##                                         mutator.rotamer.angleList, mol)

##                     ## build collider atoms to identify best rotamer
##                     colliders = mol.select('not deleted and not water and (not (segment %s chid %s resnum `%s` icode %s) or (segment %s chid %s and resnum `%s` icode %s name N O))'%(segment, chid, resnum, icode, segment, chid, resnum, icode))

##                     ## score all rotamers
##                     result = mutator.rotamer.scoreRotamers(colliders)
##                     bestRotIndex, scores, favorable, clashes = result

##                     #TO DO: check for negative score( scores[bestRotIndex])
##                     res.setCoords(rotamer.getCoordsForRotamer(bestRotIndex))
##                     inds = res.getIndices()

##                     # the 'protein" flag for the atoms that are not bb is set to False after mutation.
##                     #set 'protein' flag to "True" (prody does not let use setFlags() method for flag "protein")
##                     mol._ag._flags['protein'][inds] = numpy.ones(len(inds), dtype='bool')
##                     mol._ag._flags['receptor'][inds] = numpy.ones(len(inds), dtype='bool')
##         return mol
        

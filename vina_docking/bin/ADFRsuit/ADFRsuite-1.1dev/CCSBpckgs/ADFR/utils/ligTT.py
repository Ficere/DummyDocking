from ADFR.torTreeFromPDBQT import TreeNode
from MolKit2.openBabelInterface import ProdyToOBMol

from rings import RingConformers, FlexibleRings
#from puckerCoords import HillReillyShort as HillReilly
from puckerCoords import HillReilly


class LigTT:

    def numRings(self):
        return len(self._6rings)+len(self._5rings)
    
    def numTorsions(self):
        return len(self._torsions)
    
    def freezeBonds(self, frozenBonds, noScore=None, flexRings=True):
        self._rotatableBonds = {}

        # make single bonds rotatable
        for k,v in self.mol._ag._bondOrder.items():
            if v==1:
                self._rotatableBonds[k] = True
            else:
                self._rotatableBonds[k] = False

        # make ring bonds not rotatable
        self._obrings = self.obmol.GetSSSR() # get the Smallest Set of Smallest Rinfor ring in torTreeBuilder.rings:
        for ring in self._obrings:
            atIndices = self.obmol.obToProdyIndex[[x-1 for x in list(ring._path)]]

            #print atIndices, self.mol._ag._data['name'][list(atIndices)]
            l = len(atIndices)
            for i in range(l):
                if i==l-1:
                    j = 0
                else:
                    j = i+1
                a1 = atIndices[i] 
                a2 = atIndices[j] 
                if a1>a2:
                    tmp=a1
                    a1=a2
                    a2=tmp
                #print '   ', a1, a2, self.mol._ag._data['name'][[a1,a2]]
                self._rotatableBonds['%d %d'%(a1,a2)] = False

        # freeze single bonds that only rotate non scorable atoms
        if noScore:
            noScore = set(noScore)
            for bond in mol._ag._bondOrder.keys():
                if self._rotatableBonds[bond]:
                    i1, i2 = [int(i) for i in bond.split()]
                    indices = mol.subTree(mol._ag[i1],mol._ag[i2]).getIndices()[2:]
                    if len(noScore & set(indices))==len(indices):
                        self._rotatableBonds[bond] = False
                    indices = mol.subTree(mol._ag[i2],mol._ag[i1]).getIndices()[2:]
                    if len(noScore & set(indices))==len(indices):
                        self._rotatableBonds[bond] = False
        if flexRings:
            self._flexRings = FlexibleRings(self.obmol, self._obrings)
        else:
            self._flexRings = None

        # freeze user specified bonds
        self._frozenBonds = frozenBonds
        frozenKeys = {}
        if frozenBonds:
            for i,j in frozenBonds:
                if j<i:
                    tmp = i
                    i = j
                    j = tmp
                frozenKeys['%d %d'%(i,j)] = True

        for k,v in self.mol._ag._bondOrder.items():
            if frozenKeys.has_key(k):
                self._rotatableBonds[k] = False
                i, j = [int(x) for x in k.split()]
                print 'freezing bond', self.mol._ag[i], self.mol._ag[j]
            
    def __init__(self, mol, frozenBonds=None, noScore=None, flexRings=True):
        self.mol = mol
        self._nbTorsions = 0
        self._torsdof = 0
        self._nb6rings = 0
        self._nb5rings = 0
        self._rotatableBonds = {}
        self._frozenBonds = None

        self.obmol = ProdyToOBMol(mol.select())
        self.freezeBonds(frozenBonds, noScore=noScore, flexRings=flexRings)

            
    def _rigidBody(self, startIndex):
        """recursively find all neighbor atoms not connect by rotatable bond"""

        self._seen[startIndex] = True
        self._rbAtoms.append(startIndex)
        for at2Ind in self.mol._ag._bmap[startIndex]:
            if at2Ind==-1:
                break
            if self._seen[at2Ind]:
                continue
            if self.mol._ag._data['numbonds'][at2Ind]==1: # leaf belongs to this rigid body
                self._seen[at2Ind] = True
                self._rbAtoms.append(at2Ind)
                continue
            i1 = startIndex
            i2 = at2Ind
            if i1>i2:
                i1 = at2Ind
                i2 = startIndex
            #print i1, i2, self._rotatableBonds['%d %d'%(i1,i2)], self.mol._ag._data['name'][[i1,i2]]
            if not self._rotatableBonds['%d %d'%(i1,i2)]:
                self._rigidBody( at2Ind)
            else:
                #print i1, i2, self._rotatableBonds['%d %d'%(i1,i2)], self.mol._ag._data['name'][[i1,i2]], self.mol._ag._bondOrder['%d %d'%(i1,i2)]
                # parent is -1 because there might be another ftParent if there is a ring in this rigid body
                self._rbBonds.append( [i1, i2, -1] )
        return 

    def _walk(self, atom):
        if self._ringAtom is not None: return # we found it already
        self._seen1[atom] = True # we saw this atom
        if self._ring1[atom]:
            self._ringAtom = atom
            return
        for at2 in self.mol._ag._bmap[atom]:
            seen = self._seen1.get(at2)
            if seen is not None: # at2 is in "all"
                if seen is False:
                    self._walk(at2)

    def firstAtomInRing(self, start, all, ring):
        # return the first atom in "ring" encountered when
        # traversing the molecule starting from "start"
        # while only considering atoms in "all"
        # start is an int and all and ring are lists of int
        # corresponding to atom indices in self.mol
        self._seen1 = {}.fromkeys(all, False)
        self._ring1 = {}.fromkeys(all, False)
        self._ringAtom = None
        for i in ring:
            self._ring1[i] = True
        self._walk(start)
        del self._ring1
        del self._seen1
        value = self._ringAtom
        del self._ringAtom
        return value
    
    def rigidBodyRings(self, rigidBodyIndices, bonds, fromAtom, ftparent):
        if self._flexRings is None:
            return 0
        obmol = self._flexRings.obmol
        rings, attached = self._flexRings.findFlexibleRing(
            #[obmol.prodyToObIndex[i] for i in rigidBodyIndices])
            [i+1 for i in rigidBodyIndices])
        if len(rings)==0:
            self._nodeDescr.append([ftparent, ('ADFRcc.adfr','FTAtomsNode'),
                                    ('setAtomIndexes', [rigidBodyIndices], {})])
            for n in range(len(bonds)-1,-1,-1):
                if bonds[n][2]==-1:
                    bonds[n][2] = len(self._nodeDescr)-1
                    
            #print 'rigid'
            return 0
        elif len(rings)==1:
            #ringInd =  [obmol.obToProdyIndex[i] for i in rings[0]]
            ringInd = rings[0]
            # loop over children in tor tree to find atom indices that produce children in
            # torsion tree
            withChild = []
            attachedInd = attached[0]

            # shift ring so that first atom is point we enter from
            # which will become first verex of triangle in 6 ring
            # when flapping corners
            if fromAtom !=-1:
                shift = None
                if fromAtom in ringInd:
                    shift = ringInd.index(fromAtom)
                else:
                    # look which ring atom is where we start from
                    fromAt = self.firstAtomInRing(fromAtom, rigidBodyIndices, ringInd)
                    shift = ringInd.index(fromAt)

                if shift is None:
                    #raise RuntimeError,"%s, Ring entering atom not found for ring: "%self.mol.name+str(ringInd)
                    print "%s, Ring entering atom not found for ring: "%self.mol.name+str(ringInd)
                if shift > 0:
                    ringInd = ringInd[shift:]+ringInd[:shift]
                    attachedInd = attachedInd[shift:] + attachedInd[:shift]
                    
            n = 0
            bondsWichChildInd = []
            for i1, i2, p in bonds:
                if i1 in rigidBodyIndices:
                    withChild.append(i1)
                    bondsWichChildInd.append(n)
                elif i2 in rigidBodyIndices:
                    withChild.append(i2)
                    bondsWichChildInd.append(n)
                n+=1
            withChild = set(withChild)
            # convert attached indices from OB to prody and build subtree list
            subTreeInd = []
            for b in bonds:
                if b[2]==-1:
                    if b[0] in ringInd:
                        subTreeInd.append(ringInd.index(b[0]))
                    elif b[1] in ringInd:
                        subTreeInd.append(ringInd.index(b[1]))

            r = RingConformers(self.mol._ag.getCoords(), rigidBodyIndices,
                               ringInd, attachedInd, subTreeInd)
            #print self.mol.name, r._origAngles, len([len(x) for x in attachedInd if len(x)>0])
            allCoords, matrices = r.conformerCoords(HillReilly)
            coords = self.mol._ag.getCoords()
            # DEBUG
            ## import prody
            ## for n, c in enumerate(allCoords):
            ##     coords[ rigidBodyIndices] = c
            ##     self.mol._ag.addCoordset(coords)
            ##     prody.writePDB('rconf/%02f.pdb'%n, self.mol._ag, [n])
            self._nodeDescr.append([ftparent, ('ADFR.LigandFT','PyFTDiscreteConformation'),
                                    ('setAtomIndexes', [rigidBodyIndices], {}),
                                    ('setCoords', [allCoords.tolist()], {})] )
            if len(ringInd)==6:
                self._6rings.append( [len(self._nodeDescr)-1, r._origAngles, ringInd,
                                      len([len(x) for x in attachedInd if len(x)>0])] )
            elif len(ringInd)==5:
                self._5rings.append( [len(self._nodeDescr)-1, r._origAngles, ringInd,
                                      len([len(x) for x in attachedInd if len(x)>0])] )
            else:
                raise RuntimeError, "ring of length %d"%len(ringInd)
            
            parent = len(self._nodeDescr)-1
            for n in range(len(subTreeInd)):
                self._nodeDescr.append([parent, ('ADFRcc.adfr','FTDiscreteTransformation'),
                                        ('setMatricesData', [matrices[n]], {})])
                for ind in bondsWichChildInd:
                    i1,i2,p = bonds[ind]
                    if ringInd[subTreeInd[n]]==i1:
                        bonds[ind][2] = len(self._nodeDescr)-1
                    elif ringInd[subTreeInd[n]]==i2:
                        bonds[ind][2] = len(self._nodeDescr)-1
                    
            #print '%d conformers'%len(allCoords)
            return 1
        else:
            #print "multiple rings in root, Not handled. Using rigid root"
            return 0
        #return ftparents

    def _buildTree(self, parentNode, atInParent, at2, ftparent):
        self._rbAtoms = []
        self._rbBonds = []
        self._rigidBody(at2)
        atInd = self._rbAtoms
        #print 'RB', atInParent, at2, atInd
        #print bonds
        self._nodeNum += 1
        node = TreeNode(self._nodeNum, parent=parentNode)
        self._allRigidBodies.append(node)
        parentNode.children.append(node)
        node.atoms = atInd
        node.bond = [atInParent, at2]
        node._rbBonds = self._rbBonds
        self._nbTorsions += 1

        # check if this bond only rotates a polar hydrogen
        moved = self.mol._ag[atInd].select('hydrogen')
        if moved and len(moved)==2:
            movH = moved.select('hydrogen')
            movNotH = moved.select('not hydrogen')
            # FIXME not sure what the rules are for torsdof
            # if more than 1 H
            if len(movH)==1 and len(moved)==1:
                hat = moved.iterAtoms().next()
                heavy = hat.iterBonded().next()
                if heavy.getElement() not in ['N', 'O']:
                    self._torsdof += 1
#            if not len(moved.select('hydrogen')==1) or \
#                   moved.select('not hydrogen').getElement()=='C':
                self._torsdof += 1
        else:
            self._torsdof += 1
            
        # find 2 more atoms to define torsion
        for a0 in self.mol._ag._bmap[atInParent]:
            if a0!=at2:
                break
        for a3 in self.mol._ag._bmap[at2]:
            if a3!=atInParent:
                break

        #print 'node', self._nodeNum, 'child of node', parentNode.num, self.mol._ag._data['name'][[atInParent,at2]], self.mol._ag._data['name'][atInd], self.mol._ag.getNames()[[a0, atInParent, at2, a3]]
        self._nodeDescr.append([ftparent, ('ADFR.LigandFT', 'PyFTTorsion'),
                                ('setDihedralAtomIndexes', [a0, atInParent, at2, a3], {})])
        self._torsions.append(len(self._nodeDescr)-1)
        parent = len(self._nodeDescr)-1
        nbfr = self.rigidBodyRings(atInd, self._rbBonds, at2, parent)

        #print self.mol._ag.getNames()[atInd]
        for i1, i2, p in node._rbBonds:
            if self._seen[i1] and self._seen[i2]:
                continue
            if i1 in atInd:
                self._allRotatableBonds.append([i1,i2])
                self._buildTree(node, i1, i2, p)
            elif not self._seen[i1]:
                self._allRotatableBonds.append([i2,i1])
                self._buildTree(node, i2, i1, p)
        
    def buildTree(self, rootAtom, frozenBonds=None, noScore=None,
                  flexRings=True):
        self._seen = [False]*len(self.mol._ag)
        self.freezeBonds(frozenBonds, noScore=noScore, flexRings=flexRings)
        root = TreeNode(0)
        root.mol = self.mol
        self._nodeNum = 0
        self._nbTorsions = 0
        self._torsions = []
        self._6rings = []
        self._5rings = []

        # these are used to store inkdices of atoms and rotatable bond
        # for each rigid body foudn during traversal
        self._rbAtoms = []
        self._rbBonds = []
        self._allRotatableBonds = []
        self._allRigidBodies = [root]
        
        self._rigidBody(rootAtom.getIndex())
        root.atoms = self._rbAtoms
        root._rbBonds = self._rbBonds
        
        g = root.mol._ag.getCoords()[root.atoms[0]].tolist()
        self._nodeDescr = [[-1, ('ADFR.LigandFT','PyFTDiscreteTranslation'),
                            ('setPreferredPoints', [[[0,0,0]]], {}) ],
                           [0, ('ADFR.LigandFT','PyFTRotationAboutPointQuat'),
                            ('setRotPoint', g, {})]
                           ]
        #print 'Root', self.mol._ag.getNames()[atoms]
        #print 'Root', atoms
        self.rigidBodyRings(root.atoms, root._rbBonds, -1, 1)

        for i1, i2, p in root._rbBonds:
            if self._seen[i1] and self._seen[i2]:
                continue
            if i1 in root.atoms:
                self._allRotatableBonds.append([i1,i2])
                self._buildTree(root, i1, i2, p)
            else:
                self._allRotatableBonds.append([i2,i1])
                self._buildTree(root, i2, i1, p)

        serialToIndex = {}
        for a in self.mol.select():
            serialToIndex[a.getSerial()] = a.getIndex()
        root.serialToIndex = serialToIndex
        root.nodeDescr = self._nodeDescr
        del self._nodeDescr
        del self._rbAtoms
        del self._rbBonds

        return root
        #print [mol._ag.getNames()[list(b)] for b in bonds]

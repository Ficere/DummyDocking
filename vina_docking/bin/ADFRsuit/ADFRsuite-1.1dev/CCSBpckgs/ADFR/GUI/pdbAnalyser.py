import sys, os, weakref, prody
from collections import OrderedDict
from PySide import QtCore, QtGui

from ADFR.GUI import ICONPATH, waiting_effects
from ADFR.recFromPDB import getResidueKey, getResidueKeyAlt

from PmvApp.Pmv import RefreshDisplayEvent
from PmvApp.pmvPalettes import AtomElements
from prody.atomic.selection import Selection
from MolKit2.AARotamer import RotamerMol

residueNames = [
        'ALA','ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL']
    
class PDBAnalyserWidget(QtGui.QWidget):
    # Widget containing tools for displaying and selecting fragments of a PDB file
    # to be included in the receptor
    def __init__(self, app, pdbprofile, parent=None):
        super(PDBAnalyserWidget, self).__init__(parent)
        self.PDBprofile = pdbprofile
        self.app = weakref.ref(app) # GridGUI instance
        self.pmv = weakref.ref(app.pmv)
        self.fileName = None
        self.buildUI(parent=parent)
        self.setMinimumWidth(300)

    def highlightObj(self, atoms, mol):
        #print 'higlighting', atoms
        from MolKit2.selection import Selection
        from PmvApp.pmvPalettes import AtomElements
        self.pmv().undisplaySB(mol)
        self.pmv().undisplayCPK(mol)
        if atoms is None:
            self.pmv().customColor(mol.select('element C'), [(AtomElements['C'])])
            return
        sel = Selection(mol._ag, atoms.getIndices(), "")
        if len(atoms)==1: # e.g. ions
            self.pmv().displayCPK(sel)
        else:
            self.pmv().customColor(mol.select('element C'), [(AtomElements['C'])])
            self.pmv().displaySB(sel)
            self.pmv().customColor(sel.select('element C'), [(0.,1.,0.)])

    def focusObj(self, atoms, mol):
        # if the focus button below the BM tree is checked, create an animation
        # to translate and zoom in or out to bring atoms in focus
        from MolKit2.selection import Selection
        if self.app().autofocus():
            self.pmv().focusScene(obj=Selection(mol._ag, atoms.getIndices(), ''))

    def buildUI(self, parent=None):
        # create master layout
        self.pdbAnaLayout = QtGui.QVBoxLayout()

        # title text widget
        self.text = QtGui.QTextEdit()
        self.text.setReadOnly(True)
        html = """<title>PDB profiler</title>
</head>
<body>
<h3>PDB analyzer</h3>
<hr/>
Select biological unit, fragments to include in receptor, alternate locations. Mutate/complete residue sidechains.
</body>
</html>"""
        self.text.insertHtml(html)
        self.text.setFixedHeight(120)
        self.pdbAnaLayout.addWidget(self.text)

        expType = self.PDBprofile._receptor.pdbHeader.get('experiment', None)
        bmTrans = self.PDBprofile._receptor.pdbHeader.get('biomoltrans', None)
        hydrogens = self.PDBprofile._receptor._ag.select('hydrogen')
        if hydrogens is not None and len(hydrogens) == 0: hydrogens = None

        # descripition label
        labelTxt = "%s.pdb " % self.PDBprofile._pdbid
        if expType:
            if expType.find("NMR") >= 0:
                # some NMR files have no MODEL/ENDMDL e.g. 103d.pdb
                modelNum = self.PDBprofile._receptor.pdbHeader.get('n_models', 1)
                labelTxt +=  "%s, n models: %d" % (expType, modelNum)
            elif  expType.find("CRYSTAL") >= 0 or expType.find("X-RAY") >= 0 :
                labelTxt += "%s, resolution: %0.2f" % (expType, self.PDBprofile._receptor.pdbHeader['resolution'])
        label = QtGui.QLabel(labelTxt)
        self.pdbAnaLayout.addWidget(label)

        # create layout for tree and tool bar under it
        layout = QtGui.QVBoxLayout()
        self.pdbAnaLayout.addLayout(layout)

        # create notebook for trees for various BMs
        self.treesNoteBook = QtGui.QTabWidget()
        self.treeWidgets = OrderedDict()

        for name in self.PDBprofile._assemblies.keys():
            self.recItemsTreeW = widget = ReceptorItemsTree(self.PDBprofile, name, self)
            self.treeWidgets[name] = widget
            self.treesNoteBook.addTab(widget, name)
            
        self.treesNoteBook.currentChanged.connect(self.switchBM)
        layout.addWidget(self.treesNoteBook)

        # tool bar for BM trees
        self.toolbar = QtGui.QToolBar(self)
        # populate tool bar
        #b = self._focusToggleButton = QtGui.QPushButton(
        #    QtGui.QIcon(os.path.join(ICONPATH, 'focus.png')),'', self)
        #b.setCheckable(True)
        #b.setChecked(True)
        #b.setStatusTip('toggle focussing 3D view on selected tree element')
        #self.toolbar.addWidget(b)
        self.toolbar.addSeparator()
        
        b = self._showReceptorButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'magnifying-glass.png')),
            'rec.', self)
        b.setStatusTip('Show atoms included in the receptor')
        self.toolbar.addWidget(b)

        b = self._showLigandButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'magnifying-glass.png')),
            'lig.', self)
        b.setStatusTip('Show atoms included in ligand(s)')
        self.toolbar.addWidget(b)

        b = self._showOtherButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'magnifying-glass.png')),
            'other', self)
        b.setStatusTip('Show atoms not included in receptor or ligand')
        self.toolbar.addWidget(b)

        self._showReceptorButton.clicked.connect(self.showCurrentReceptor)
        self._showLigandButton.clicked.connect(self.showCurrentLigand)
        self._showOtherButton.clicked.connect(self.showOther)
        layout.addWidget(self.toolbar)

        # # build group for saving receptor and/or ligand and move to AGFRgui
        b = self.buildBMGroup = QtGui.QGroupBox(self)
        #b.setStyleSheet("QGroupBox { padding: 0px;}")
        self.saveRecAsBt =  QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'folder.png')),'rec. ...')
        self.saveRecAsBt.clicked.connect(self.makeAndSaveRecFile)

        self.saveLigAsBt =  QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'folder.png')), 'lig. ...')
        self.saveLigAsBt.clicked.connect(self.makeAndSaveLigFile)

        # create the create to proceed to AGFRgui to specify the docking box
        self.proceedToMapsBt = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH,'iso.png')), 'maps')
        self.proceedToMapsBt.clicked.connect(self.app().finalizeMolecules)
        
        gblayout = QtGui.QHBoxLayout()
        self.buildBMGroup.setLayout(gblayout)
        gblayout.addWidget(self.saveRecAsBt)
        gblayout.addWidget(self.saveLigAsBt)
        gblayout.addWidget(self.proceedToMapsBt)
        self.pdbAnaLayout.addWidget(self.buildBMGroup)

        # maybe we put this in settings panel since we can specify different
        # protonation methods, at least for the ligand 
        #self.addHcheckB = cb = QtGui.QCheckBox("add Hydrogens")
        #cb.stateChanged.connect(self._toggleAddH)

        ## FIXME WE HAVE TO SEE HOW TO HANDLE THIS
        ## addH = True
        ## if not expType:
        ##     addH = False
        ## else:
        ##     if expType == 'SOLUTION NMR' and hydrogens:
        ##         addH = False
        ## cb.setChecked(addH)
        ## self.PDBprofile._addH= addH

        self.setLayout(self.pdbAnaLayout)

        # enable/disable buttons based on receptor/ligand flags
        self._configureButtons()
        
    ## def _toggleAddH(self, val=None):
    ##     # callback of the "add hydrogens" check button
    ##     self.PDBprofile._addH = {0:False, 2:True}.get(val)

    def _configureButtons(self):
        # enable/disable buttons based on receptor/ligand flags
        rec = self.PDBprofile._currentBM
        recAtoms = rec.select('receptor')
        ligAtoms = rec.select('ligand')
        otherAtoms = rec.select('not receptor and not ligand')
        if not recAtoms:
            self._showReceptorButton.setDisabled(True)
            self.proceedToMapsBt.setDisabled(True)
            self.saveRecAsBt.setDisabled(True)
        else:
            self._showReceptorButton.setDisabled(False)
            self.proceedToMapsBt.setDisabled(False)
            self.saveRecAsBt.setDisabled(False)

        if not ligAtoms:
            self._showLigandButton.setDisabled(True)
            self.saveLigAsBt.setDisabled(True)
        else:
            self._showLigandButton.setDisabled(False)
            self.saveLigAsBt.setDisabled(False)
            
        if not otherAtoms:
            self._showOtherButton.setDisabled(True)
        else:
            self._showOtherButton.setDisabled(False)
        
    # FIXME try using transparency of what we do not want to see
    def showCurrentReceptor(self):
        rec = self.PDBprofile._currentBM
        self.pmv().undisplaySB(rec)
        self.pmv().customColor(rec.select('element C'), [(AtomElements['C'])])
        sel = rec.select('receptor')
        if sel:
            self.pmv().displaySB(sel)
            self.pmv().customColor(sel.select('element C'), [(0.27,0.62,0.11)])
            self.focusObj(sel, rec)

    # FIXME try using transparency of what we do not want to see
    def showCurrentLigand(self):
        rec = self.PDBprofile._currentBM
        self.pmv().undisplaySB(rec)
        self.pmv().customColor(rec.select('element C'), [(AtomElements['C'])])
        sel = rec.select('ligand')
        if sel:
            self.pmv().displaySB(sel)
            self.pmv().customColor(sel.select('element C'), [(0.27,0.62,0.11)])
            self.focusObj(sel, rec)
        
    # FIXME try using transparency of what we do not want to see
    def showOther(self):
        rec = self.PDBprofile._currentBM
        self.pmv().undisplaySB(rec)
        self.pmv().customColor(rec.select('element C'), [(AtomElements['C'])])
        sel = rec.select('not ligand and not receptor')
        if sel:
            self.pmv().displaySB(sel)
            self.pmv().customColor(sel.select('element C'), [(0.8,0.62,0.3)])
            self.focusObj(sel, rec)

    @waiting_effects
    def _processReceptor(self):
        mol = self.app().agfr._BMProcessed
        if mol is None:
            mol = self.app().agfr.processBM()

        atoms = mol.select('receptor')
        if atoms is None:
            print 'WARNING: no atoms selected for the receptor'
            return
        return self.app().agfr.processReceptor(atoms, charges=self.getCharges())
    
    def makeAndSaveRecFile(self):
        # process the receptor ro have inof needed for pdbqt
        # and ask for filename and save remember processed receptor
        agfr = self.app().agfr
        processedReceptor = self._processReceptor()
        
        fdialog = QtGui.QFileDialog()
        fdialog.setDefaultSuffix("pdbqt")
        name = self.PDBprofile._pdbid+'_rec.pdbqt'
        filepath = fdialog.getSaveFileName(
            self, "receptor filename", name, 
            filter="PDBQT files (*.pdbqt);; Any files (*)")[0]

        if filepath is None: # cancelled dialog
            return
        
        name, _format = os.path.splitext(filepath)
        if _format=='':
            _format = '.pdbqt'
            filepath += _format

        agfr.saveReceptorAsPDBQT(processedReceptor, filepath)
        agfr._processedReceptor = processedReceptor

    @waiting_effects
    def _processLigands(self, ligatoms):
        return self.app().agfr.processLigands(ligatoms)

    def _writeLines(self, lines, filename):
        f = open(filename, 'w')
        [f.write('%s\n'%l) for l in lines]
        f.close()
        
    def makeAndSaveLigFile(self, filename=None):
        # process the ligand to have info needed for pdbqt
        # and ask for filename and save remember processed receptor
        agfr = self.app().agfr
        mol = agfr._BMProcessed
        if mol is None:
            mol = agfr.processBM()
 
        ligAtoms = mol.select('ligand')
        frags = prody.findFragments(ligAtoms)
        if len(frags)>1:
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("INFO: ligand atoms form %d connected fragments.\nWe will create a file for each fragment"%len(frags))
            msgBox.exec_()

        agfr._processedLigands = self._processLigands(ligAtoms)
        if filename is None:
            fdialog = QtGui.QFileDialog()
            fdialog.setDefaultSuffix("pdbqt")
            name = self.PDBprofile._pdbid+'_lig.pdbqt'
            filepath = fdialog.getSaveFileName(self, "ligand filename", name, 
                                               filter="PDBQT files (*.pdbqt);; Any files (*)")[0]

            if filepath is None: # cancelled dialog
                return
        else:
            filepath = filename

        name, _format = os.path.splitext(filepath)
        if _format=='':
            _format = '.pdbqt'
            filepath += _format

        if len(agfr._adLigs)>1:
            basename = os.path.splitext(filepath)[0]
            for num, adlig in enumerate(agfr._adLigs): 
                filename = '%s_%03d.pdbqt'%(basename, num+1)
                self._writeLines(adlig.getPDBQTlines(), filename)
                agfr._processedLigands[num].filename = filename
        else:
            self._writeLines(agfr._adLigs[0].getPDBQTlines(), filepath)
            agfr._processedLigands[0].filename = filepath

    @waiting_effects
    def switchBM(self, *args): # added *args for waiting effects to work
        treeW = self.treesNoteBook.currentWidget()
        name = treeW._name
        if not self.treeWidgets[name]._hasUI:
            self.treeWidgets[name].buildUI()
            updateAltlocs = False
        else:
            updateAltlocs = True

        #print "switching to ", name
        for nn in self.PDBprofile._assemblies.keys():
            mol = self.treeWidgets[nn]._mol
            if nn == name:
                # update alternate location selections in this tree to reflect
                # self.PDBprofile._altlocsRes which might have been altered in another tree
                if updateAltlocs:
                    for chid in self.treeWidgets[nn]._chaindInBM:
                        if self.treeWidgets[nn]._altlocItems.has_key(chid): # 
                            for key, widget in self.treeWidgets[nn]._altlocItems[chid].items():
                                ind, alts, occ = self.PDBprofile._altlocsRes[key]
                                widget.altWidget.buttons[ind].setChecked(True)

                self.pmv().showMolecules([x._mol for x in self.treeWidgets[nn]._rotamols], True)
                self.pmv().showMolecules(mol, True)
                focusMol = mol
                self.PDBprofile._currentBM = mol
                self.app().addGeomsForGaps(self.PDBprofile._gapsPerChain)
            else:
                self.pmv().showMolecules([x._mol for x in self.treeWidgets[nn]._rotamols], False)
                self.pmv().showMolecules(mol, False)
        self.pmv().focusScene(obj=focusMol, duration=0.1, padding=0)

    def getCharges(self):
        #import pdb; pdb.set_trace()
        charges = self.treesNoteBook.currentWidget()._getCharges()
        self.PDBprofile._customCharges = charges
        return charges

    def getAltlocResidues(self):
        return self.treesNoteBook.currentWidget()._getAltlocResidues()
    

class CustomTristateItem(QtGui.QTreeWidgetItem):
    def __init__(self,parent=None):
        super(CustomTristateItem, self).__init__(parent)
        self.partiallyCheckedItms= []
        self._checkBColumn = 1
        self._childCheckBColumn = 1
        
    def setData(self, a, b, c):
        state = self.checkState(self._checkBColumn )
        #print "set data:", a, b, c, "state",state
        if state == QtCore.Qt.Unchecked:
            # make this click for partial selection
            #import pdb; pdb.set_trace()
            if self.childCount() > 0:
                for item, st in self.partiallyCheckedItms:
                    #print "setting state", item, st
                    item.setCheckState(self._childCheckBColumn, st)
            else:
                super(CustomTristateItem, self).setData(a,b,c)
        else:
            super(CustomTristateItem, self).setData(a, b, c) 

    def setpartiallyCheckedItms(self, items):
        #items is a list of lists [TreeItemObject, QtCore.Qt.Checked(Unchecked) ]
        self.partiallyCheckedItms = items


class AltlocsChooser(QtGui.QWidget):
    """
    Widget used to display alternate locations wih occupancies and alow to chose
    one of them exclusively (i.e. radio button)
    """
    
    #def __init__(self, ind, alts, occ, altAtoms, resItem, profile, focus_cb, parent=None):
    def __init__(self, resItem, alts, occ, current, tree, parent=None):
        # alts: list of alternate location characters ['A', 'B']
        # occ : list of alternate location occupancies [.65, .25]
        # current: the 0-based index of the currently selected alternate location
        assert len(alts)==len(occ)
        assert current >= 0 and current<len(alts)
        assert isinstance(tree, ReceptorItemsTree)
        
        super(AltlocsChooser, self). __init__(parent)
        self.resItem = resItem
        self._altlocs = alts       # hold list of alternate location character
        self._occupancies = occ    # hold list of alternate location occupancies
        self._tree = tree
        self._current = current # store index of currently selected alternate location
        self._altloc = alts[current] # hold the current alternate location character
        self._active = False # used to prevent onclick from doing anything upon construction
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        self.buttons = []
        for i in range(len(alts)):
            button = QtGui.QRadioButton('%s %.2f'%(alts[i], occ[i]))
            button.toggled.connect(self.onclick)
            self.buttons.append(button)
            if i==current: button.setChecked(True)
            layout.addWidget(button)

    def onclick(self):
        if not self._active: return
        button = self.sender()
        if button.isChecked():
            alt = button.text().split()[0]
            self._current = self._altlocs.index(alt)
            self._altloc = self._altlocs[self._current]
            self._tree._setItemFlags(self.resItem, 'receptor', self.resItem.checkState(0))
            self._tree._setItemFlags(self.resItem, 'ligand', self.resItem.checkState(2))
            # update display
            atoms = self.resItem._pmvObj.select('altloc _ %s'%self._altloc)
            self._tree.visuallyIdentify(atoms)
            key = getResidueKey(self.resItem._pmvObj)
            # update PDBprofile._altlocsRes
            self._tree.app().PDBprofile._altlocsRes[key][0] =  self._current

    def getAltlocChar(self):
        return self._altloc

    def getAltlocIndex(self):
        return self._current
    
from ADFR.GUI import ICONPATH
from mglutil.util.callback import CallbackFunction

class ReceptorItemsTree(QtGui.QWidget):
    """
    Tree widget for picking parts of the PDB file to include in the receptor
    or ligand, define charge for atom snot handles by Gasteiger, mutate
    residues etc...
    """
    def __init__(self, profile, name, analyserGUI, parent=None):
        super(ReceptorItemsTree, self).__init__(parent)
        self._chaindInBM, self._mol = profile._assemblies[name]
        self._ag = self._mol._ag
        self._suspend = False
        self._selectionToItem = {}
        self.profile = profile
        self.showgapscb = None
        self._typItems = {}
        self.fileName = None
        self.app = weakref.ref(analyserGUI) # instance of PDBAnalyserWidget
        self.pmv = weakref.ref(analyserGUI.app().pmv)
        self._name = name
        self._rotamols = [] # list of RotamerMol object used to display rotamer for completed side chains
        self._hasUI = False # we will lazy buikd when it gets selected
        self._altlocItems = {} # {chid: {getResidueKey(res): resItem}}

    def buildUI(self):
        # lazy build triggered by switchBM()
        if self._hasUI: return
        self._hasUI = True

        layout =  QtGui.QVBoxLayout()
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        #self.toolbar = QtGui.QToolBar(self)
        layout.addWidget(self.treeWidget)
        #layout.addWidget(self.toolbar)
        self.setLayout(layout)
        
        treeWidget.setColumnCount(3)
        treeWidget.headerItem().setText(0, self._mol.name+" receptor")
        #treeWidget.headerItem().setText(1, "rec.")
        treeWidget.headerItem().setText(1, "params")
        treeWidget.headerItem().setText(2, "lig.")
        treeWidget.setColumnWidth(0, 200)
        treeWidget.setColumnWidth(1, 70)
        treeWidget.setColumnWidth(2, 30)
        treeWidget.itemClicked.connect(self.onClick)
        self.treeWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeWidget.customContextMenuRequested.connect(self.openResMenu)

        chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
        resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
        atIcon = QtGui.QIcon(os.path.join(ICONPATH, 'atom.png'))

        for chain in self._ag.select('segment A').getHierView().iterChains():

            chidSel = chid = chain.getChid()
            if chidSel == ' ': chidSel = '_' # used for prody selections

            # only consider chains from ASU that appear in this BM
            if chid not in self._chaindInBM:
                continue
                        
            # create the tree item for the chain
            root = chItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), chid)
            chItem.setExpanded(True)
            chItem.setIcon(0,chIcon)
            
            # dict to save tree item by typ
            self._typItems[chid] = typItems = {}

            altlocNO = chain.select('not altlocs and not unsre')
            altloc = chain.select('altlocs and not unsre')
            unsup = chain.select('unsre')

            for _iter, chainAtoms in enumerate([altlocNO, altloc, unsup]):
                # on second iteration creat the altloc tree item that will be the parent
                # item for all categories
                if _iter==1: # we iterate or set with altlocs
                    self._altlocItems[chid] = altlocItems = {}
                    if chainAtoms is None or len(chainAtoms)==0:
                        #self._altlocItems[chid] = {}
                        continue
                    chItem = QtGui.QTreeWidgetItem(chItem)
                    chItem.setText(0, 'altloc')
                    chItem.setExpanded(True)
                    typItems['altloc'] = chItem
                    chItem.setFlags(chItem.flags() | QtCore.Qt.ItemIsTristate)
                    chItem._pmvObj = chainAtoms

                if _iter==2: # we iterate or set with unsupported atom types
                    if chainAtoms is None or len(chainAtoms)==0:
                        #self._altlocItems[chid] = {}
                        continue
                    chItem = QtGui.QTreeWidgetItem(root)
                    chItem.setText(0, 'unknown atom types')
                    chItem.setExpanded(True)
                    typItems['unsre'] = chItem
                    chItem.setFlags(chItem.flags() | QtCore.Qt.ItemIsTristate)
                    
                # loop over all categories
                for typ in ['protein', 'nucleic', 'sidechains to complete',
                            'missing atoms', 'cofactors', 'modified', 'manual charges',
                            'ligands', 'additives', 'water', 'other']:
                    if typ=='sidechains to complete':
                        cat = 'cmpsc'
                    elif typ=='missing atoms':
                        cat = 'missa'
                    elif typ=='manual charges':
                        cat = 'manch'
                    elif typ=='water':
                        cat = '_wate'
                    else:
                        cat = typ[:5]

                    # check if we have residues in this category
                    atoms = chainAtoms.select('category %s'%cat)
                    if atoms:
                        # create the top node
                        parentItem = typItems[typ] = QtGui.QTreeWidgetItem(chItem )
                        parentItem.setText(0, '%s'%(typ))
                        parentItem._pmvObj = atoms
                        parentItem.setFlags(typItems[typ].flags() | QtCore.Qt.ItemIsTristate)
                        if typ in ['protein', 'nucleic', 'sidechains to complete']:
                            parentItem.setCheckState(0, QtCore.Qt.Checked)
                        elif not _iter==2:
                            parentItem.setCheckState(0, QtCore.Qt.Unchecked)

                        if _iter==0 or _iter==2: # no alternate location
                            alt = '_' # these atoms have no alternate locations
                            for res in atoms.getHierView().iterResidues():
                                # make sure the residue has no alternate locations
                                allres = self._ag.select('segment A chid %s resindex %d'%(chidSel, res.getResindices()[0]))
                                if len(allres)>len(res):
                                    continue # some atoms in this residue have altlocs
                                resnum = res.getResnum()
                                resname = res.getResname()
                                icode = icodeSel = res.getIcode()
                                if not icode:
                                    icodeSel = "_"
                                resItem = QtGui.QTreeWidgetItem(parentItem)
                                resItem.setIcon(0, resIcon)
                                resItem._pmvObj = res
                                key = getResidueKey(res)
                                if res.select('category modif'):
                                    old, new = self.profile._resToMutate.get(key, (None, None))
                                    if old is not None:
                                        resItemTxt =  "%s'%s -> %s" %(old, resnum, new)
                                    else:
                                        resItemTxt =  "%s'%s" %(resname, resnum)
                                        
                                elif res.select('category cmpsc'):
                                    missingAtomNames = self.profile._missingAtoms[key]
                                    resItemTxt = "%s'%s%s build:%s"%(resname, resnum, icode, ' '.join(missingAtomNames))
                                    # Add molecule to PMV for to show completed sidechain
                                    self.addRotamerMol(resItem, res, chidSel, resname, resnum, icodeSel)
                                elif res.select('category missa'):
                                    # amino acid with missing atoms
                                    atoms = self.profile._missingAtoms[key]
                                    resItemTxt = "%s'%s%s missing:%s"%(resname, resnum, icode, ' '.join(atoms))
                                else:
                                    resItemTxt = "%s'%s%s"%(resname, resnum, icode)

                                if res.select('receptor'):
                                    resItem.setCheckState(0, QtCore.Qt.Checked)
                                elif not _iter==2:
                                    resItem.setCheckState(0, QtCore.Qt.Unchecked)

                                if res.select('category manch'):
                                    resItem.setIcon(0, atIcon)
                                    resItem.spinBox = QtGui.QDoubleSpinBox()
                                    charge = res.getCharges()[0]
                                    resItem.spinBox.setMaximum( 10.0 )
                                    resItem.spinBox.setMinimum( -10.0 )
                                    resItem.spinBox.setSingleStep( .05 )
                                    resItem.spinBox.setValue( charge )
                                    treeWidget.setItemWidget(resItem, 1, resItem.spinBox )

                                resItem.setText(0, resItemTxt)
                                if res.select('ligand'):
                                    resItem.setCheckState(2, QtCore.Qt.Checked) # lig check box
                                    resItem.setCheckState(0, QtCore.Qt.Unchecked) # rec check box
                                    print "ligand res:", res.getResname(), res.getResnum()
                                elif _iter==2: # unsupported
                                    # remove user checkable flag
                                    resItem.setFlags(resItem.flags() ^ QtCore.Qt.ItemIsUserCheckable)
                                    #resItem.setCheckState(2, QtCore.Qt.Unchecked) # lig check box
                                    #resItem.setCheckState(0, QtCore.Qt.Unchecked) # rec check box
                                    #import pdb; pdb.set_trace()
                                else:
                                    resItem.setCheckState(2, QtCore.Qt.Unchecked) # lig check box
                                    
                                if _iter==2: # for residues with unsupported atom types
                                    # add a child to the residue for the unsupported atom
                                    resItem.setExpanded(True)
                                    atoms = res.select('unsat')
                                    if atoms:
                                        for atom in atoms:
                                            atomItem = QtGui.QTreeWidgetItem(resItem)
                                            atomItem.setText(0, '%s (%s)'%(atom.getName(), atom.getElement()))
                                            atomItem._pmvObj = Selection(self._ag, [atom.getIndex()], '')
                                            # remove user checkable flag
                                            atomItem.setFlags(atomItem.flags() ^ QtCore.Qt.ItemIsUserCheckable)


                        else: # with alternate location
                            # get dict of altres descriptions for this chain
                            altres = self.profile._altlocsRes
                            # loop over residues
                            for res in atoms.getHierView().iterResidues():
                                # atoms in this residue that have altloc == ' ' are not included
                                # we need to expand atoms to contains all atoms in this residue
                                #resAtoms = chain.select('resindex %d'%res.getResindices()[0])
                                #res = resAtoms.getHierView().iterResidues().next()
                                key = getResidueKey(res)
                                resItem = QtGui.QTreeWidgetItem(parentItem)
                                altlocItems[key] = resItem
                                resItem._pmvObj = res
                                name = "%s'%s%s"%(res.getResname(), res.getResnum(),
                                                 res.getIcode())
                                resItem.setText(0 , name)
                                resItem.setCheckState(2, QtCore.Qt.Unchecked) # lig checkbox
                                if res.select('receptor'):
                                    resItem.setCheckState(0, QtCore.Qt.Checked)
                                else:
                                    resItem.setCheckState(0, QtCore.Qt.Unchecked)
                                # get index of currently selected altloc, list of alts and occups
                                ind, alts, occ = altres[key]
                                altAtoms = [res.select('altloc _ or altloc %s'%alt) for alt in alts]
                                #if key=="A:SER:166:_":
                                #    import pdb; pdb.set_trace()
                                #altWidget = AltlocsChooser(ind, alts, occ, altAtoms, resItem,
                                #                           self.profile, self.focusObj)
                                altWidget = AltlocsChooser(resItem, alts, occ, ind, self)
                                resItem.altWidget = altWidget # so the object does not get garbage coll
                                treeWidget.setItemWidget(resItem, 1, altWidget)
                                altWidget._active = True

                                if _iter==2: # for residues with unsupported atom types
                                    # add a child to the residue for the unsupported atom
                                    resItem.setExpanded(True)
                                    atoms = res.select('unsat')
                                    if atoms:
                                        for atom in atoms:
                                            atomItem = QtGui.QTreeWidgetItem(resItem)
                                            atomItem.setText(0, '%s (%s)'%(atom.getName(), atom.getElement()))
                                            atomItem._pmvObj = atom
                                        # remove user checkable flag
                                        atomItem.setFlags(atomItem.flags() ^ QtCore.Qt.ItemIsUserCheckable)

                # add entry for gaps
                chgaps = self.profile._gapsPerChain.get(chid, None)
                if _iter==0 and chgaps:
                    misResItem = QtGui.QTreeWidgetItem(chItem)
                    misResItem.setText(0, 'missing residues')
                    if len(chgaps[0]):
                        ntermItem = QtGui.QTreeWidgetItem(misResItem)
                        ntermItem.setText(0, "N term: %s" % chgaps[0][1])
                        rName, rNum, rIcode = chgaps[0][0]
                        if rIcode:
                            # do not want to use the _pmbObj attribute , so that it is not used to set
                            # 'receptor' flag to the residues .
                            ntermItem._pmvGapObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, rNum, rIcode))
                        else:
                            ntermItem._pmvGapObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, rNum))
                    for gap in chgaps[1]:
                        gapItem = QtGui.QTreeWidgetItem(misResItem)
                        #r1 = gap[0]
                        #r2 = gap[-1]
                        r1Name, r1Num, r1Icode = gap[0]
                        r2Name, r2Num, r2Icode = gap[-1]
                        rr = tuple(gap[1:-1])
                        gapItem.setText(0, 'gap: %s%d--'%(r1Name, r1Num) + "%s"*len(rr)%rr + "--%s%d"%(r2Name, r2Num))
                        if r1Icode:
                            r1 = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, r1Num, r1Icode))
                        else:
                            r1 = self._ag.select("segment A chid %s resnum `%d`" %(chid, r1Num))
                        if r2Icode:
                            r2 = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, r2Num, r2Icode))
                        else:
                            r2 = self._ag.select("segment A chid %s resnum `%d`" %(chid, r2Num))
                        gapItem._pmvGapObj = r1+r2                    
                    if len(chgaps[2]):
                        ctermItem = QtGui.QTreeWidgetItem(misResItem)
                        ctermItem.setText(0, "C term: %s" % chgaps[2][1])
                        cName, cNum, cIcode = chgaps[2][0]
                        if cIcode:
                            ctermItem._pmvGapObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, cNum, cIcode))
                        else:
                            ctermItem._pmvGapObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, cNum))
            # end loop over altlocNO,altloc
            
    def setRecFlag(self, item, on):

        def _setRecFlag(item, on):
            if hasattr(item, 'altWidget'):
                # this node represents a residue with alternate location
                # we set the flag only for the atoms with altloc _ or currenly selected
                if on:
                     # none of the atoms of this residue can be ligand
                    item._pmvObj.setFlags('ligand', False)
                    if item.flags() & QtCore.Qt.ItemIsUserCheckable:
                        item.setCheckState(2, QtCore.Qt.Unchecked)

                    # set receptor=False for all atoms in the residue
                    item._pmvObj.setFlags('receptor', False) # set all
                    # get the atoms for the current altloc and set their receptor flag to True
                    atoms = item._pmvObj.select('altloc _ %s'%item.altWidget.getAltlocChar())
                    atoms.setFlags('receptor', True)
                else:
                    # set receptor=False for all atoms in the residue
                    item._pmvObj.setFlags('receptor', False) # set all
            else:
                if on:
                    item._pmvObj.setFlags('ligand', False)
                    if item.flags() & QtCore.Qt.ItemIsUserCheckable:
                        item.setCheckState(2, QtCore.Qt.Unchecked)
                    item._pmvObj.setFlags('receptor', True) # set all
                else:
                    item._pmvObj.setFlags('receptor', False) # set all
                    
        nbChildren = item.childCount()
        if nbChildren:
            #iterate over children
            for n in range(item.childCount()):
                self.setRecFlag(item.child(n), on)
        else:
            _setRecFlag(item, on)

    def setLigFlag(self, item, on):

        def _setLigFlag(item, on):
            if hasattr(item, 'altWidget'):
                # this node represents a residue with alternate location
                # we set the flag only for the atoms with altloc _ or currenly selected
                if on:
                     # none of the atoms of this residue can be receptor
                    item._pmvObj.setFlags('receptor', False)
                    if item.flags() & QtCore.Qt.ItemIsUserCheckable:
                        item.setCheckState(0, QtCore.Qt.Unchecked)

                    # set ligand=False for all atoms in the residue
                    item._pmvObj.setFlags('ligand', False) # set all
                    # get the atoms for the current altloc and set their receptor flag to True
                    atoms = item._pmvObj.select('altloc _ %s'%item.altWidget.getAltlocChar())
                    atoms.setFlags('ligand', True)
                else:
                    # set receptor=False for all atoms in the residue
                    item._pmvObj.setFlags('ligand', False) # set all False
            else:
                if on:
                    item._pmvObj.setFlags('receptor', False)
                    if item.flags() & QtCore.Qt.ItemIsUserCheckable:
                        item.setCheckState(0, QtCore.Qt.Unchecked)
                    item._pmvObj.setFlags('ligand', True) # set all
                else:
                    item._pmvObj.setFlags('ligand', False) # set all
                    
        nbChildren = item.childCount()
        if nbChildren:
            #iterate over children
            for n in range(item.childCount()):
                self.setLigFlag(item.child(n), on)
        else:
            _setLigFlag(item, on)
        
    def _setItemFlags(self, item, flag_name, on):
        # item is a node in the tree that has a _pmvObj attribute
        # toggle receptor and ligand flags in the currentBM

        assert flag_name in ['receptor', 'ligand']
        #assert on in [True, False, 1, 0] # can aso be PySide.QtCore.Qt.CheckState.Checked

        if flag_name=='receptor':
            self.setRecFlag(item, on)
            # reset processed receptor
            self.app().app().agfr._processedReceptor = None
        
        if flag_name=='ligand':
            self.setLigFlag(item, on)
            # reset list of processed ligands
            self.app().app().agfr._adLigs = []
            self.app().app().agfr._processedLigands = []
            
        self.app().app().agfr._BMProcessed = None
        
        # enable/disable buttons based on receptor/ligand flags
        self.app()._configureButtons()

    def onClick(self, item, column):
        if column==0: # clicked on name or receptor check button
            if hasattr(item, '_pmvObj'):
                obj = item._pmvObj
                # manage the receptor flag
                if item.flags() & QtCore.Qt.ItemIsUserCheckable:
                    self._setItemFlags(item, 'receptor', item.checkState(column))

                if hasattr(item, 'altWidget'):
                    obj = obj.select('receptor')

            elif hasattr(item, '_pmvGapObj'):
                obj = item._pmvGapObj
            else:
                obj = None

            if obj:
                self.visuallyIdentify(obj)
                
        elif column==2: # click on ligand check box
            # handle ligand flag
            if hasattr(item, '_pmvObj') and item.flags() & QtCore.Qt.ItemIsUserCheckable:
                self._setItemFlags(item, 'ligand', item.checkState(column))

    def visuallyIdentify(self, obj):
        self.app().highlightObj(obj, self._mol)
        self.app().focusObj(obj, self._mol)
            
    def _setSCcompletion(self, act, item):
        def setName(item, act):
            res = item._pmvObj
            resnum = res.getResnum()
            resname = res.getResname()
            icode = res.getIcode()
            key = getResidueKey(res)
            atoms = self.profile._missingAtoms[key]
            if act=='complete':
                resItemTxt = "%s%s%s build:%s"%(resname, resnum, icode, ' '.join(atoms))
                self.profile._resToMutate[key] = [resname,  resname]
            else:
                resItemTxt = "%s%s%s use as is missing:%s"%(resname, resnum, icode, ' '.join(atoms))
                if self.profile._resToMutate.has_key(key):
                    self.profile._resToMutate.pop(key)
            item.setText(0, resItemTxt)

        if str(item.text(0)) in ['sidechains to complete', 'partial backbone']:
            for n in range(item.childCount()):
                setName(item.child(n), act)
        else:
            setName(item, act)
        

    def openResMenu(self, position):
        itms = self.treeWidget.selectedItems()
        if len(itms) > 0: # len(items) is always 1 currently because we
                          # cannot select more than one item at a time
            level = 0
            selectedItem = item = itms[0]
            while item.parent() and item.text(0)!='altloc':
                print level, item.text(0),
                if item.parent():
                    item.parent().text(0)
                else:
                    print
                item = item.parent()
                level += 1
        # now level is 0:chain, 1:category, 2:item in category

        def sideChainCompleteMenu(menu, selectedItem, forCategory):
            cb = CallbackFunction(self._setSCcompletion, 'complete', selectedItem)
            act1 = QtGui.QAction('Complete sidechain', self, checkable=True, triggered=cb)       
            cb = CallbackFunction(self._setSCcompletion, 'include as is', selectedItem)
            act2 = QtGui.QAction('Include as is', self, checkable=True, triggered=cb)
            ag = QtGui.QActionGroup(self)#, exclusive=True)
            ag.addAction(act1)
            ag.addAction(act2)
            if forCategory:
                act1.setChecked(True)
            else:
                if self.profile._resToMutate.has_key(getResidueKey(selectedItem._pmvObj)):
                    act1.setChecked(True)
                else:
                    act2.setChecked(True)

            menu.addAction(act1)
            menu.addAction(act2)
            menu.exec_(self.treeWidget.viewport().mapToGlobal(position))

        menu = QtGui.QMenu()
        if level == 1: # we click on a category
            prtTxt = selectedItem.text(0)
            if prtTxt in ["sidechains to complete", "partial backbone"]:
                # for the side chain to complete category we have 2 choices: complete or keepAsIs
                sideChainCompleteMenu(menu, selectedItem, True)

        elif level == 2:
            prtTxt = selectedItem.parent().text(0)
            if prtTxt in ["sidechains to complete"]:
                sideChainCompleteMenu(menu, selectedItem, False)
                
            elif prtTxt in [ "modified", "protein"]:
                itemTxt = selectedItem.text(0).split(" -> ")
                res = selectedItem._pmvObj
                resname = res.getResname()

                if prtTxt == "modified":
                    baseResName,allKnown, mutResName = self.profile._modifiedRes[getResidueKeyAlt(res)]
                    resToExclude = [baseResName]
                    if allKnown:
                        if len(itemTxt) > 1:
                            names = [resname, "Base Res: %s" % baseResName]
                        else:
                            names = ["Base Res: %s" % baseResName]
                    else:
                        names = ["Base Res: %s" % baseResName]
                    resNames = list(set(residueNames) - set(resToExclude))
                    resNames.sort()
                    names.extend(resNames)

                elif  selectedItem.parent().text(0) =="protein":
                    if len(itemTxt) > 1:
                        names = ["Base Res: %s" % resname]
                        resNames = list(set(residueNames) - set([resname]))
                        resNames.sort()
                        names.extend(resNames)
                    else:
                        names = list(set(residueNames) - set([resname]))
                        names.sort()

                for _name in names:
                    cb = CallbackFunction(self._selected,  _name, selectedItem)
                    act = QtGui.QAction(_name, self, triggered = cb)
                    menu.addAction(act)
                menu.exec_(self.treeWidget.viewport().mapToGlobal(position))

    def addRotamerMol(self, resItem, res, chid, resname, resnum, icode):
        # create RotamerMol object and add rotamol._mol to PMV
        # add combo box allowing to browse rotamers
        pmv = self.pmv()
        if resname=='ALA':
            # AlA with CB missing e.g. 2p4n A:ALA265
            return
        natoms = len(res._ag)
        rotmol = RotamerMol(resname)
        rotmol._mol.name = 'rotamer for %s:%s:%s'%(chid, resname, resnum)
        self._rotamols.append(rotmol)
        rotmol.rotamer.alignRotToResBB(res.select('bb or name CB'))
        resItem._rotmol = rotmol
        natoms1 = len(res._ag)
        if natoms1 != natoms: # alignRotToResBB() calls build4Points() that can add a CB atom to the _ag
            # update allCoords array in the geometry container of the molecule
            mol = res._ag.getMolecule()
            gc = mol.geomContainer
            gc.allCoords = mol._ag._coords[mol._ag._acsi].astype('f')
            for atomset in gc.atoms.values():
                atomset._ag = mol._ag
            coords = mol._ag.getCoords().astype('f')
            # fix geometry vertices
            for geom in gc.geoms.values():
                if hasattr(geom, '_hasAllVertices') and geom._hasAllVertices:
                    geom.vertexSet.vertices.array = coords
        ## build a set of collider atoms to identify best rotamer and save in resItem, just in case
        segment = res.getSegnames()[0]
        resItem.colliders = self._ag.select('not deleted and not water and (not (segment %s chid %s resnum `%s` icode %s) or (segment %s chid %s and resnum `%s` icode %s name N O))'%(segment, chid, resnum, icode, segment, chid, resnum, icode))

        ## score all rotamers and save results in resItem
        result = rotmol.rotamer.scoreRotamers(resItem.colliders)
        bestRotIndex, scores, favorable, clashes = result
        # save rotamer scores
        resItem._autoRotIndex = bestRotIndex
        resItem._rotaScores = [scores, favorable, clashes]
        # set side chain coordinates
        indices = rotmol._mol._ag.select('sc').getIndices()
        rotmol._mol._ag[indices].setCoords(rotmol.rotamer.getCoordsForRotamer(bestRotIndex)[indices])
        # set the coordinate of the CA atom
        rotmol._mol._ag.select('ca').setCoords(res.select('ca').getCoords())

        ## add molecule to PmvApp
        pmv.addMolecule(rotmol._mol)
        pmv.customColor(rotmol._mol.select('element C'), [(.2,8.,6.)])
        pmv.undisplayLines(rotmol._mol.select('name N C O'))

        ## add the rotamer combobox
        resItem.comboBox = QtGui.QComboBox()
        resItem.comboBox.addItems(['auto %d'%(bestRotIndex+1)]+
                                  [str(x) for x in range(1, len(rotmol.rotamer.angleList)+1)])
        self.treeWidget.setItemWidget(resItem, 1, resItem.comboBox)
        cb = CallbackFunction(self._setRotamer, resItem)
        resItem.comboBox.currentIndexChanged.connect(cb)
        
    def _setRotamer(self, resItem, index):
        ## callback function for residue with rotamer object on side chain
        rotmol = resItem._rotmol
        pmv = self.pmv()
        if index==0:
            # in dex 0 is the bstRotamer identified automatically
            index = resItem._autoRotIndex
        else:
            # else we go from 1-based rotamer to 0-based indexing
            index = index - 1
        # get the coordinates for this rotamer and set them for side chain atoms
        newcoords = rotmol.rotamer.getCoordsForRotamer(index)
        rotmol._mol._ag.setCoords(newcoords)
        indices = rotmol._mol._ag.select('sc').getIndices()
        # update the allCoords array un the geometry container
        rotmol._mol.geomContainer.allCoords[indices] = newcoords[indices]
        # display message about the current rotamer
        # FIXME would better handled with signal
        self.app().app().statusbar.showMessage('Score: %.1f, contacts: %d favorable %d clashes'%(
            resItem._rotaScores[0][index], len(resItem._rotaScores[1][index]),
            len(resItem._rotaScores[2][index])))
        #pmv.undisplayLines(rotmol._mol)
        #pmv.displayLines(rotmol._mol.select('not name N C O'))

        # update the Viewer to reflex the new coordinates, whatever the geometry is
        pmv = self.pmv()
        pmv.gui().viewer.suspendRedraw = True
        event = RefreshDisplayEvent(molecule=rotmol._mol)
        pmv.eventHandler.dispatchEvent(event)
        pmv.gui().viewer.selectionEventHandler()
        pmv.gui().viewer.suspendRedraw = False
        pmv.gui().viewer.OneRedraw()
       
    def _selected(self, newresname, item):
        # callback function of the residue item's drop down menu.
        res = item._pmvObj
        resnum = res.getResnum()
        resname = res.getResname()
        
        #print "menu item selected: ", newresname, "RESNAME:", resname, "RESNUM:", resnum
        key = getResidueKey(res)
        newresname = newresname.split("Base Res:")
        if len(newresname) > 1:
            mutResName = newresname[1].strip()
        else:
            mutResName = newresname[0].strip()
        if mutResName == resname:
            mutResName = None
            item.setText(0, "%s%s"%(resname, resnum))
            if self.profile._resToMutate.has_key(key):
                self.profile._resToMutate.pop(key)
        else:
            item.setText(0, "%s%s -> %s"%(resname,resnum, mutResName))
            self.profile._resToMutate[key] = [resname,  mutResName]

    def _getCharges(self):
        charges = {}
        for chid, chItems in self._typItems.items():
            item = chItems.get('manual charges', None)
            if item:
                for nn in range(item.childCount()):
                    mItem = item.child(nn)
                    if  mItem.checkState(1) == QtCore.Qt.Checked:
                        charge = mItem.spinBox.value()
                        res = mItem._pmvObj
                        key = getResidueKey(res)
                        charges[key] = charge
        return charges

    def _getAltlocResidues(self):
        # return a dictionary of residue data (name, num, icode, alt) per chain included in the receptor and that have alternate locations.  
        altresDict = {}
        for chid in self._chaindInBM:
            reslist = []
            for key, resItem in self._altlocItems[chid].items():
                if resItem._pmvObj.select('receptor') is None:
                    continue # ignore things not in the receptor
                chid, resname, resnum, icode = key.split(':')
                if icode=='_': icode = ''
                reslist.append((resname, resnum, icode))
            altresDict[chid] = reslist
        #import pdb; pdb.set_trace()
        return altresDict

        ## for chid, resdict in self._altlocResPerChain.items():
        ##     allres = []
        ##     for resStr, reslist in resdict.items():
        ##         for alt, res in reslist:
        ##             sel = res.select('receptor')
        ##             #print chid, res.getResname(), res.getResnum(), alt, len(res), res.getFlags('receptor'),
        ##             #if sel : print len(sel)
        ##             #else: print 0
                    
        ##             if  sel and len(sel)==len(res): 
        ##                 allres.append([res.getResname(), res.getResnum(), res.getIcode(), alt])
        ##     altresDict[chid] = allres
        ## #print "altloc residues:", altresDict
        ## return altresDict


    ## def _subtree(self, atomsets, parentItem, altloc):
    ##     if altloc:
    ##         item = QtGui.QTreeWidgetItem(parentItem)
    ##         item._pmvObj = atoms
    ##         item.setText(0, 'alternate locations')
    ##         parentItem = item
            
    ##     for name, atoms in atomsets.items():
    ##         if atoms and len(atoms):
    ##             item = QtGui.QTreeWidgetItem(parentItem)
    ##             item._pmvObj = atoms
    ##             item.setText(0, name)
    ##             item.setFlags(parentItem.flags() | QtCore.Qt.ItemIsTristate)
    ##             item.setCheckState(0, QtCore.Qt.Checked)

    ##             for res in atoms.getHierView().iterResidues():
    ##                 resnum = res.getResnum()
    ##                 resname = res.getResname()
    ##                 icode = icodeSel = res.getIcode()
    ##                 if not icode:
    ##                     icodeSel = "_"
    ##                 resItem = QtGui.QTreeWidgetItem(item)
    ##                 #resItem.setIcon(0, resIcon)
    ##                 resItem._pmvObj = res
    ##                 key = getResidueKey(res)
    ##                 resItemTxt = "%s'%s%s"%(resname, resnum, icode)
    ##                 resItem.setText(0, resItemTxt)
    ##                 if res.select('receptor'):
    ##                     resItem.setCheckState(0, QtCore.Qt.Checked)
    ##                 else:
    ##                     resItem.setCheckState(0, QtCore.Qt.Unchecked)
                    
    ## def _subtreeForCategory(self, chid, catatoms, name, root):

    ##     # create the top node
    ##     parentItem = self._typItems[chid][name] = QtGui.QTreeWidgetItem(root)
    ##     parentItem.setText(0, name)
    ##     parentItem._pmvObj = catatoms
    ##     parentItem.setFlags(parentItem.flags() | QtCore.Qt.ItemIsTristate)
    ##     parentItem.setCheckState(0, QtCore.Qt.Checked)
    ##     parentItem.setExpanded(True)
        
    ##     noalt = catatoms.select('not altlocs')
    ##     alt = catatoms.select('altlocs')
    ##     for altFlag, atoms in enumerate([noalt, alt]):
    ##         if not atoms: continue
    ##         if name in ["metal", "water"]:
    ##             if altFlag==0:
    ##                 d = {'without alternate location': atoms}
    ##             else:
    ##                 d = {'with alternate locations': atoms}
    ##         elif name=='protein':
    ##             d = {'complete': atoms.select('not missS and not missA'),
    ##                  'missing sidechain atoms': atoms.select('missS'),
    ##                  'other missing atoms': atoms.select('missA')}
    ##         else:
    ##             d = {'complete': atoms.select('missS and not missA'),
    ##                  'missing atoms': atoms.select('missA')}

    ##         self._subtree( d, parentItem, altFlag)
        
    ## def buildUI2(self):
    ##     # relies on classify2
    ##     if self._hasUI: return
    ##     self._hasUI = True

    ##     layout =  QtGui.QVBoxLayout()
    ##     self.treeWidget = treeWidget = QtGui.QTreeWidget()
    ##     self.toolbar = QtGui.QToolBar(self)
    ##     layout.addWidget(self.treeWidget)
    ##     layout.addWidget(self.toolbar)

    ##     # build tree
    ##     self.treeWidget = treeWidget = QtGui.QTreeWidget()
    ##     treeWidget.setColumnCount(3)
    ##     treeWidget.headerItem().setText(0, self._mol.name)
    ##     #treeWidget.headerItem().setText(1, "rec.")
    ##     treeWidget.headerItem().setText(1, "params")
    ##     treeWidget.headerItem().setText(2, "lig.")
    ##     treeWidget.setColumnWidth(0, 200)
    ##     treeWidget.setColumnWidth(1, 70)
    ##     #treeWidget.setColumnWidth(2, 40)
    ##     treeWidget.setColumnWidth(2, 30)
    ##     treeWidget.itemClicked.connect(self.onClick)
    ##     self.treeWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    ##     self.treeWidget.customContextMenuRequested.connect(self.openResMenu)
        
    ##     #chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
    ##     #resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
    ##     #atIcon = QtGui.QIcon(os.path.join(ICONPATH, 'atom.png'))

    ##     # loop over chains
    ##     for chain in self._ag.select('segment A').getHierView().iterChains():
    ##         chid = chain.getChid()
    ##         self._altlocItems[chid] = {}

    ##         # only consider chains from ASU that appear in this BM
    ##         if chid not in self._chaindInBM:
    ##             continue
                        
    ##         # create the tree item for the chain
    ##         chItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), chid)
    ##         chItem.setExpanded(True)
    ##         #chItem.setIcon(0,chIcon)

    ##         # dict to save tree item by category
    ##         self._typItems[chid] = {}

    ##         protein = chain.protein
    ##         if protein:
    ##             self._subtreeForCategory(chid, protein, 'protein', chItem)

    ##         nucleic = chain.nucleic
    ##         if nucleic:
    ##             self._subtreeForCategory(chid, nucleic, 'nucleic', chItem)

    ##         other = chain.select('not protein and not nucleic')
    ##         supAtoms = other.select('not unsupType')
    ##         unsupAtoms = other.select('unsupType')
    ##         for unsupFlag, atoms in enumerate([supAtoms, unsupAtoms]):
    ##             if not atoms: continue
    ##             item = QtGui.QTreeWidgetItem(chItem)
    ##             if unsupFlag:
    ##                 item.setText(0, 'Unsupported atom type')
    ##             else:
    ##                 item.setText(0, 'hetero')
    ##             item.setFlags(item.flags() | QtCore.Qt.ItemIsTristate)
    ##             item.setCheckState(0, QtCore.Qt.Unchecked)
    ##             item.setExpanded(True)
                
    ##             for cat in ['cofactor', 'metal', 'ligand', 'additives', 'water', 'other']:
    ##                 if cat=='water':
    ##                     atoms = chain.select('category _wate')
    ##                 else:
    ##                     atoms = chain.select('category %s'%cat)
    ##                 if atoms:
    ##                     self._subtreeForCategory(chid, atoms, cat, item)

    ##         # add entry for gaps
    ##         chgaps = self.profile._gapsPerChain.get(chid, None)
    ##         if chgaps:
    ##             misResItem = QtGui.QTreeWidgetItem(chItem)
    ##             misResItem.setText(0, 'missing residues')
    ##             if len(chgaps[0]):
    ##                 ntermItem = QtGui.QTreeWidgetItem(misResItem)
    ##                 ntermItem.setText(0, "N term: %s" % chgaps[0][1])
    ##                 rName, rNum, rIcode = chgaps[0][0]
    ##                 if rIcode:
    ##                     ntermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, rNum, rIcode))
    ##                 else:
    ##                     ntermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, rNum))
    ##             for gap in chgaps[1]:
    ##                 gapItem = QtGui.QTreeWidgetItem(misResItem)
    ##                 #r1 = gap[0]
    ##                 #r2 = gap[-1]
    ##                 r1Name, r1Num, r1Icode = gap[0]
    ##                 r2Name, r2Num, r2Icode = gap[-1]
    ##                 rr = tuple(gap[1:-1])
    ##                 gapItem.setText(0, 'gap: %s%d--'%(r1Name, r1Num) + "%s"*len(rr)%rr + "--%s%d"%(r2Name, r2Num))
    ##                 if r1Icode:
    ##                     r1 = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, r1Num, r1Icode))
    ##                 else:
    ##                     r1 = self._ag.select("segment A chid %s resnum `%d`" %(chid, r1Num))
    ##                 if r2Icode:
    ##                     r2 = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, r2Num, r2Icode))
    ##                 else:
    ##                     r2 = self._ag.select("segment A chid %s resnum `%d`" %(chid, r2Num))
    ##                 gapItem._pmvObj = r1+r2                    
    ##                 #gapItem._gapData = gap
    ##             if len(chgaps[2]):
    ##                 ctermItem = QtGui.QTreeWidgetItem(misResItem)
    ##                 ctermItem.setText(0, "C term: %s" % chgaps[2][1])
    ##                 cName, cNum, cIcode = chgaps[2][0]
    ##                 if cIcode:
    ##                     ctermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, cNum, cIcode))
    ##                 else:
    ##                     ctermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, cNum))
    ##         # end loop over altlocNO,altloc

    ##     layout = QtGui.QVBoxLayout()
    ##     layout.addWidget(self.treeWidget)
    ##     self.setLayout(layout)

## class ReceptorItemsDialog(QtGui.QDialog):
##     def __init__(self, profile, title='receptor items selector', parent=None):
##         super(ReceptorItemsDialog, self).__init__(parent)
##         self.setWindowTitle(title)
##         self.recItemsWidget = ReceptorItemsTree(profile, parent)
##         #self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
##         self.buttonBox = QtGui.QDialogButtonBox()
##         self.buttonBox.accepted.connect(self.accept)

##         layout = QtGui.QVBoxLayout()
##         layout.addWidget(self.recItemsWidget)
##         layout.addWidget(self.buttonBox)
##         # Set dialog layout
##         self.setLayout(layout)

##     def reject(self):
##         QtGui.QDialog.reject(self)
##         self.closeEvent()
        
##     def closeEvent(self, evnt=None):
##         print 'closed'

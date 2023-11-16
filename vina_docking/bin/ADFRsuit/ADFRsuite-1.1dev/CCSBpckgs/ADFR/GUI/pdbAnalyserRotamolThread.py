import sys, os
from collections import OrderedDict
from PySide import QtCore, QtGui
from threading import Thread

from ADFR.GUI import ICONPATH
from ADFR.recFromPDB import getResidueKey, getResidueKeyAlt
from ADFR.GUI.gridGUI import waiting_effects

from PmvApp.Pmv import RefreshDisplayEvent

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
        self.app = app
        self.fileName = None
        self.buildUI(parent=parent)

    def buildUI(self, parent=None):
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
        self.pdbAnaLayout = QtGui.QVBoxLayout()
        self.treesNoteBook = QtGui.QTabWidget()
        self.treeWidgets = OrderedDict()

        for name in self.PDBprofile._assemblies.keys():
            self.recItemsTreeW = widget = ReceptorItemsTree(self.PDBprofile, name, self)
            widget.highlightcb = self.app._highlightPDBParts
            widget.focuscb = self.app._focusPDBParts
            self.treeWidgets[name] = widget
            self.treesNoteBook.addTab(widget, name)
        
        self.treesNoteBook.currentChanged.connect(self.switchBM)
        
        self.pdbAnaLayout.addWidget(self.text)
        self.pdbAnaLayout.addWidget(self.treesNoteBook)

        # create the create receptor button
        self.createReceptorBt = QtGui.QPushButton('create receptor PDBQT, proceed to AGFR')
        self.createReceptorBt.clicked.connect(self.createReceptor)

        # tool bar widget
        self.toolbar = QtGui.QToolBar(self)
        b = self._focusToggleButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'crosshair.png')),
            'auto focus', self)
        b.setCheckable(True)
        b.setChecked(True)
        b.setStatusTip('toggle focussing 3D view on selected tree element')
        self.toolbar.addWidget(b)
        self.toolbar.addSeparator()
        
        b = self._showReceptorButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'receptor2_NO.png')),
            'show rec.', self)
        b.setStatusTip('Show atoms included in receptor')
        self.toolbar.addWidget(b)

        b = self._showNotReceptorButton = QtGui.QPushButton(
            QtGui.QIcon(os.path.join(ICONPATH, 'ligand64_NO.png')),
            'show not rec.', self)
        b.setStatusTip('Show atoms not included in receptor')
        self.toolbar.addWidget(b)
        self.pdbAnaLayout.addWidget(self.toolbar)
        self._showReceptorButton.clicked.connect(self.app._showCurrentReceptor)
        self._showNotReceptorButton.clicked.connect(self.app._showCurrentNotReceptor)
        self._completionCB = self.app.completeReceptorCreationFromPDB

        self.buildBMGroup = QtGui.QGroupBox(self)
        # widget for specifying name under which to save PDBQT file
        self.saveAsBt =  QtGui.QPushButton('save receptor as ...')
        self.saveAsBt.clicked.connect(self.selectPdbqtFile)
        self.saveAsLabel = QtGui.QLabel("%s.pdbqt" % self.PDBprofile._pdbid)

        gblayout = QtGui.QGridLayout()#QtGui.QHBoxLayout()
        gblayout.setContentsMargins(0,0,0,0)
        self.buildBMGroup.setLayout(gblayout)

        self.addHcheckB = cb = QtGui.QCheckBox("add Hydrogens")
        cb.stateChanged.connect(self._toggleAddH)
        expType = self.PDBprofile._receptor.pdbHeader.get('experiment', None)
        bmTrans = self.PDBprofile._receptor.pdbHeader.get('biomoltrans', None)
        hydrogens = self.PDBprofile._receptor._ag.select('hydrogen')
        if hydrogens is not None and len(hydrogens) == 0: hydrogens = None
        #print "header:", self.PDBprofile._receptor.pdbHeader.keys()
        #print "EXPERIMENT:", expType
        buildBM = False
        addH = True
        if not expType:
            addH = False
        else:
            if expType == 'SOLUTION NMR' and hydrogens:
                addH = False
        cb.setChecked(addH)
        self.PDBprofile._buildBM=buildBM
        self.PDBprofile._addH= addH
        
        gblayout.addWidget(self.saveAsBt, 0,0,1,1)
        gblayout.addWidget(self.saveAsLabel, 0,1,1,1)
        gblayout.addWidget(cb, 1, 0, 1,1)
        gblayout.addWidget(self.createReceptorBt, 2, 0,1,2)

        layout = QtGui.QVBoxLayout()
        labelTxt = "%s.pdb " % self.PDBprofile._pdbid
        if expType:
            if expType.find("NMR") >= 0:
                # some NMR files have no MODEL/ENDMDL e.g. 103d.pdb
                modelNum = self.PDBprofile._receptor.pdbHeader.get('n_models', 1)
                labelTxt +=  "%s, n models: %d" % (expType, modelNum)
            elif  expType.find("CRYSTAL") >= 0 or expType.find("X-RAY") >= 0 :
                labelTxt += "%s, resolution: %0.2f" % (expType, self.PDBprofile._receptor.pdbHeader['resolution'])
        label = QtGui.QLabel(labelTxt)
        layout.addWidget(label)

        self.pdbAnaLayout.addLayout(layout)       
        self.pdbAnaLayout.addWidget(self.buildBMGroup)

        self.setLayout(self.pdbAnaLayout)

    def _toggleAddH(self, val=None):
        # callback of the "add hydrogens" check button
        self.PDBprofile._addH = {0:False, 2:True}.get(val)

    def createReceptor(self):
        self.createReceptorBt.setDisabled(True)
        if self._completionCB:
            self._completionCB()

    def selectPdbqtFile(self):
        file_dialog = QtGui.QFileDialog(parent=self,
            caption="Select PDBQT file name",
            #directory=,
            filter="PDBQT file (*.pdbqt)")
        file_dialog.setDefaultSuffix("pdbqt")
        file_dialog.setLabelText(QtGui.QFileDialog.Accept, "OK")
        file_dialog.setLabelText(QtGui.QFileDialog.Reject, "Cancel")
        #file_dialog.setOption(QtGui.QFileDialog.DontResolveSymlinks)
        file_dialog.setOption(QtGui.QFileDialog.DontUseNativeDialog)
        if not file_dialog.exec_():
            return
        fileName = file_dialog.selectedFiles()[0]
        if fileName:
            self.fileName = fileName.encode('ascii', 'replace')
            import os
            self.saveAsLabel.setText("%s" % os.path.basename(self.fileName))

    @waiting_effects
    def switchBM(self, *args): # added *args for waiting effects to work
        treeW = self.treesNoteBook.currentWidget()
        name = treeW._name
        if not self.treeWidgets[name]._hasUI:
            #self.treeWidgets[name].buildUI2()
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
                        for key, widget in self.treeWidgets[nn]._altlocItems[chid].items():
                            ind, alts, occ = self.PDBprofile._altlocsRes[key]
                            widget.altWidget.buttons[ind].setChecked(True)

                self.app.pmv.showMolecules([x._mol for x in self.treeWidgets[nn]._rotamols], True)
                self.app.pmv.showMolecules(mol, True)
                focusMol = mol
                self.PDBprofile._currentBM = mol
                self.app.addGeomsForGaps(self.PDBprofile._gapsPerChain)
            else:
                self.app.pmv.showMolecules([x._mol for x in self.treeWidgets[nn]._rotamols], False)
                self.app.pmv.showMolecules(mol, False)
        self.app.pmv.focusScene(obj=focusMol)

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

    def __init__(self, ind, alt, occ, altAtoms, resItem, profile, focus_cb, parent=None):
        super(AltlocsChooser, self). __init__(parent)
        self.altAtoms = altAtoms
        self.resItem = resItem
        self.focus_cb = focus_cb
        self.profile = profile
        self._active = False # used to prevent onclick from doing anything upon construction
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,4,0,4)
        self.buttons = []
        for i in range(len(alt)):
            button = QtGui.QRadioButton('%s %.2f'%(alt[i], occ[i]))
            self.buttons.append(button)
            button.ind = i
            button.toggled.connect(self.onclick)
            if i==ind:
                button.setChecked(True)
                resItem._pmvObj = altAtoms[i]
            layout.addWidget(button)

    def onclick(self):
        if not self._active: return
        button = self.sender()
        if button.isChecked():
            # reset flag for current altloc residue
            self.resItem._pmvObj.setFlags('receptor', False)
            res = self.altAtoms[button.ind]
            # set new altloc res to be the on for this residue
            self.resItem._pmvObj = res
            # set the receptor flag for the new residue
            res.setFlags('receptor', True)
            res.setFlags('ligand', False)
            # set the selected altloc for this residue in profile
            key = getResidueKey(res)
            self.profile._altlocsRes[key][0] = button.ind
            # update display
            self.focus_cb(self.resItem._pmvObj)
            
from ADFR.GUI import ICONPATH
from mglutil.util.callback import CallbackFunction

class ReceptorItemsTree(QtGui.QWidget):
    """
    Tree widget for picking parts of the PDB file to include in the receptor
    and set charges of kept metals
    """
    def __init__(self, profile, name, gui, parent=None):
        super(ReceptorItemsTree, self).__init__(parent)
        self._chaindInBM, self._mol = profile._assemblies[name]
        self._ag = self._mol._ag
        self._suspend = False
        self._selectionToItem = {}
        self.profile = profile
        self.highlightcb = None
        self.focuscb = None
        self.showgapscb = None
        self._completionCB = None
        self._typItems = {}
        self.fileName = None
        self.gui = gui
        self._name = name
        self._rotamols = [] # list of RotamerMol object used to display rotamer for completed side chains
        self._hasUI = False # we will lazy buikd when it gets selected
        self._altlocItems = {} # {chid: {getResidueKey(res): resItem}}
        #self.buildUI()

    def _subtree(self, atomsets, parentItem, altloc):
        if altloc:
            item = QtGui.QTreeWidgetItem(parentItem)
            item._pmvObj = atoms
            item.setText(0, 'alternate locations')
            parentItem = item
            
        for name, atoms in atomsets.items():
            if atoms and len(atoms):
                item = QtGui.QTreeWidgetItem(parentItem)
                item._pmvObj = atoms
                item.setText(0, name)
                item.setFlags(parentItem.flags() | QtCore.Qt.ItemIsTristate)
                item.setCheckState(0, QtCore.Qt.Checked)

                for res in atoms.getHierView().iterResidues():
                    resnum = res.getResnum()
                    resname = res.getResname()
                    icode = icodeSel = res.getIcode()
                    if not icode:
                        icodeSel = "_"
                    resItem = QtGui.QTreeWidgetItem(item)
                    #resItem.setIcon(0, resIcon)
                    resItem._pmvObj = res
                    key = getResidueKey(res)
                    resItemTxt = "%s'%s%s"%(resname, resnum, icode)
                    resItem.setText(0, resItemTxt)
                    if res.select('receptor'):
                        resItem.setCheckState(0, QtCore.Qt.Checked)
                    else:
                        resItem.setCheckState(0, QtCore.Qt.Unchecked)
                    
    def _subtreeForCategory(self, chid, catatoms, name, root):

        # create the top node
        parentItem = self._typItems[chid][name] = QtGui.QTreeWidgetItem(root)
        parentItem.setText(0, name)
        parentItem._pmvObj = catatoms
        parentItem.setFlags(parentItem.flags() | QtCore.Qt.ItemIsTristate)
        parentItem.setCheckState(0, QtCore.Qt.Checked)
        parentItem.setExpanded(True)
        
        noalt = catatoms.select('not altlocs')
        alt = catatoms.select('altlocs')
        for altFlag, atoms in enumerate([noalt, alt]):
            if not atoms: continue
            if name in ["metal", "water"]:
                if altFlag==0:
                    d = {'without alternate location': atoms}
                else:
                    d = {'with alternate locations': atoms}
            elif name=='protein':
                d = {'complete': atoms.select('not missS and not missA'),
                     'missing sidechain atoms': atoms.select('missS'),
                     'other missing atoms': atoms.select('missA')}
            else:
                d = {'complete': atoms.select('missS and not missA'),
                     'missing atoms': atoms.select('missA')}

            self._subtree( d, parentItem, altFlag)
        
    def buildUI2(self):
        # relies on classify2
        if self._hasUI: return
        self._hasUI = True
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        treeWidget.setColumnCount(3)
        treeWidget.headerItem().setText(0, self._mol.name)
        #treeWidget.headerItem().setText(1, "rec.")
        treeWidget.headerItem().setText(1, "params")
        treeWidget.headerItem().setText(2, "lig.")
        treeWidget.setColumnWidth(0, 200)
        treeWidget.setColumnWidth(1, 70)
        #treeWidget.setColumnWidth(2, 40)
        treeWidget.setColumnWidth(2, 30)
        treeWidget.itemClicked.connect(self.itemClicked)
        self.treeWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeWidget.customContextMenuRequested.connect(self.openResMenu)
        
        #chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
        #resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
        #atIcon = QtGui.QIcon(os.path.join(ICONPATH, 'atom.png'))

        # loop over chains
        for chain in self._ag.select('segment A').getHierView().iterChains():
            chid = chain.getChid()
            self._altlocItems[chid] = {}

            # only consider chains from ASU that appear in this BM
            if chid not in self._chaindInBM:
                continue
                        
            # create the tree item for the chain
            chItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), chid)
            chItem.setExpanded(True)
            #chItem.setIcon(0,chIcon)

            # dict to save tree item by category
            self._typItems[chid] = {}

            protein = chain.protein
            if protein:
                self._subtreeForCategory(chid, protein, 'protein', chItem)

            nucleic = chain.nucleic
            if nucleic:
                self._subtreeForCategory(chid, nucleic, 'nucleic', chItem)

            other = chain.select('not protein and not nucleic')
            supAtoms = other.select('not unsupType')
            unsupAtoms = other.select('unsupType')
            for unsupFlag, atoms in enumerate([supAtoms, unsupAtoms]):
                if not atoms: continue
                item = QtGui.QTreeWidgetItem(chItem)
                if unsupFlag:
                    item.setText(0, 'Unsupported atom type')
                else:
                    item.setText(0, 'hetero')
                item.setFlags(item.flags() | QtCore.Qt.ItemIsTristate)
                item.setCheckState(0, QtCore.Qt.Unchecked)
                item.setExpanded(True)
                
                for cat in ['cofactor', 'metal', 'ligand', 'additives', 'water', 'other']:
                    if cat=='water':
                        atoms = chain.select('category _wate')
                    else:
                        atoms = chain.select('category %s'%cat)
                    if atoms:
                        self._subtreeForCategory(chid, atoms, cat, item)

            # add entry for gaps
            chgaps = self.profile._gapsPerChain.get(chid, None)
            if chgaps:
                misResItem = QtGui.QTreeWidgetItem(chItem)
                misResItem.setText(0, 'missing residues')
                if len(chgaps[0]):
                    ntermItem = QtGui.QTreeWidgetItem(misResItem)
                    ntermItem.setText(0, "N term: %s" % chgaps[0][1])
                    rName, rNum, rIcode = chgaps[0][0]
                    if rIcode:
                        ntermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, rNum, rIcode))
                    else:
                        ntermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, rNum))
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
                    gapItem._pmvObj = r1+r2                    
                    #gapItem._gapData = gap
                if len(chgaps[2]):
                    ctermItem = QtGui.QTreeWidgetItem(misResItem)
                    ctermItem.setText(0, "C term: %s" % chgaps[2][1])
                    cName, cNum, cIcode = chgaps[2][0]
                    if cIcode:
                        ctermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d` icode %s" %(chid, cNum, cIcode))
                    else:
                        ctermItem._pmvObj = self._ag.select("segment A chid %s resnum `%d`" %(chid, cNum))
            # end loop over altlocNO,altloc
            
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.treeWidget)
        self.setLayout(layout)

    def buildUI(self):
        if self._hasUI: return
        self._hasUI = True
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        treeWidget.setColumnCount(3)
        treeWidget.headerItem().setText(0, self._mol.name+" receptor")
        #treeWidget.headerItem().setText(1, "rec.")
        treeWidget.headerItem().setText(1, "params")
        treeWidget.headerItem().setText(2, "lig.")
        treeWidget.setColumnWidth(0, 200)
        treeWidget.setColumnWidth(1, 70)
        #treeWidget.setColumnWidth(2, 40)
        treeWidget.setColumnWidth(2, 30)
        treeWidget.itemClicked.connect(self.itemClicked)
        self.treeWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeWidget.customContextMenuRequested.connect(self.openResMenu)
        
        chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
        resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
        atIcon = QtGui.QIcon(os.path.join(ICONPATH, 'atom.png'))
        getBestRotaList = [] # list used to call getBestRotamer in the background

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

            altlocNO = chain.select('not altlocs and not unsupType')
            altloc = chain.select('altlocs and not unsupType')
            unsup = chain.select('unsupType')

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

                if _iter==2: # we iterate or set with unsupported atom types
                    if chainAtoms is None or len(chainAtoms)==0:
                        #self._altlocItems[chid] = {}
                        continue
                    chItem = QtGui.QTreeWidgetItem(root)
                    chItem.setText(0, 'unknown atom types')
                    chItem.setExpanded(True)
                    typItems['unsup'] = chItem
                    chItem.setFlags(chItem.flags() | QtCore.Qt.ItemIsTristate)
                    
                # loop over all categories
                for typ in ['protein', 'nucleic', 'sidechains to complete',
                            'missing atoms', 'cofactors', 'modified', 'metals',
                            'ligands', 'additives', 'water', 'other']:
                    if typ=='sidechains to complete':
                        cat = 'cmpsc'
                    elif typ=='missing atoms':
                        cat = 'missa'
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
                                    # create combox but disaebl it
                                    resItem.comboBox = QtGui.QComboBox()
                                    resItem.comboBox.setDisabled(True)
                                    # Add molecule to PMV for to show completed sidechain
                                    self.addRotamerMol(resItem, res, chidSel, resname, resnum, icodeSel)
                                    # this is a little slow so we will do in in the background
                                    getBestRotaList.append([resItem, res, chidSel, resname, resnum, icodeSel])
                                    
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

                                if res.select('category metal'):
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
                                elif not _iter==2:
                                    resItem.setCheckState(2, QtCore.Qt.Unchecked) # lig check box
                                    
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
                                if res.receptor:
                                    resItem.setCheckState(0, QtCore.Qt.Checked)
                                else:
                                    resItem.setCheckState(0, QtCore.Qt.Unchecked)
                                # get index of currently selected altloc, list of alts and occups
                                ind, alts, occ = altres[key]
                                altAtoms = [res.select('altloc _ or altloc %s'%alt) for alt in alts]
                                #if key=="A:SER:166:_":
                                #    import pdb; pdb.set_trace()
                                altWidget = AltlocsChooser(ind, alts, occ, altAtoms, resItem,
                                                           self.profile, self.focusObj)
                                resItem.altWidget = altWidget # so the object does not get garbage coll
                                treeWidget.setItemWidget(resItem, 1, altWidget)
                                altWidget._active = True

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
            
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.treeWidget)
        self.setLayout(layout)
        # start a thread to find best rotamers
        if getBestRotaList:
            self.addRotamerMol_bg(getBestRotaList)
            #thread = Thread(target=self.addRotamerMol_bg, args = (getBestRotaList,))
            #thread.start()
            #thread.join()
            print 'thread started'
        
    def itemClicked(self, item, column):
        #print "itemClicked:", item, column, item.text(0)
        from MolKit2.selection import Selection
        self.highlightcb(None, self._mol)

        if column==0: # clicked on name or receptor check button
            if not hasattr(item, '_pmvObj'):
                obj= None
            else:
                obj = item._pmvObj
            if hasattr(item, '_pmvGapObj'):
                obj = item._pmvGapObj
            #print 'Highlight ', obj

            if self.highlightcb:
                self.highlightcb(obj, self._mol)
                
            if obj and self.focuscb and self.gui._focusToggleButton.isChecked():
                self.focuscb(obj, self._mol)

            # manage the receptor flag
            if hasattr(item, '_pmvObj'):
                if item.checkState(column):
                    item._pmvObj.setFlags('receptor', True)
                    item._pmvObj.setFlags('ligand', False)
                    item.setCheckState(2, QtCore.Qt.Unchecked)
                else:
                    item._pmvObj.setFlags('receptor', False)

        elif column==2: # click on ligand check box
            # handle ligand flag
            if hasattr(item, '_pmvObj'):
                if item.checkState(column):
                    item._pmvObj.setFlags('ligand', True)
                    item._pmvObj.setFlags('receptor', False)
                    item.setCheckState(0, QtCore.Qt.Unchecked)
                else:
                    item._pmvObj.setFlags('ligand', False)

    def focusObj(self, focusObj):
        if self.highlightcb:
            self.highlightcb(focusObj, self._mol)
        if self.focuscb and self.gui._focusToggleButton.isChecked():
            self.focuscb(focusObj, self._mol)
            
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

    def addRotamerMol_bg(self, rotalist):
        # function to add the RotamerMol objects in the background
        for resItem, res, chid, resname, resnum, icode in rotalist:
            ## build a set of collider atoms to identify best rotamer and save in resItem, just in case
            segment = res.getSegnames()[0]
            #import pdb; pdb.set_trace()
            resItem.colliders = self._ag.select('not deleted and not water and (not (segment %s chid %s resnum `%s` icode %s) or (segment %s chid %s and resnum `%s` icode %s name N O))'%(segment, chid, resnum, icode, segment, chid, resnum, icode))

            rotmol = resItem._rotmol
            result = rotmol.rotamer.scoreRotamers(resItem.colliders)
            bestRotIndex, scores, favorable, clashes = result
            # save rotamer scores
            resItem._autoRotIndex = bestRotIndex
            resItem._rotaScores = [scores, favorable, clashes]
            # set side chain coordinates
            indices = rotmol._mol._ag.select('sc').getIndices()
            rotmol._mol._ag[indices].setCoords(rotmol.rotamer.getCoordsForRotamer(bestRotIndex)[indices])
            # configure combobox
            resItem.comboBox.addItems(['auto %d'%(bestRotIndex+1)]+
                                      [str(x) for x in range(1, len(rotmol.rotamer.angleList)+1)])
            cb = CallbackFunction(self._setRotamer, resItem)
            resItem.comboBox.currentIndexChanged.connect(cb)
            resItem.comboBox.setDisabled(False)
            print 'ROTAMOL THREAD', res, bestRotIndex, scores
           
    def addRotamerMol(self, resItem, res, chid, resname, resnum, icode):
        # create RotamerMol object and add rotamol._mol to PMV
        # add combo box allowing to browse rotamers
        pmv = self.gui.app.pmv
        rotmol = RotamerMol(resname)
        rotmol._mol.name = 'rotamer for %s:%s:%s'%(chid, resname, resnum)
        self._rotamols.append(rotmol)
        rotmol.rotamer.alignRotToResBB(res.select('bb'))
        resItem._rotmol = rotmol

        # set the coordinate of the CA atom
        rotmol._mol._ag.select('ca').setCoords(res.select('ca').getCoords())

        ## add molecule to PmvApp
        pmv.addMolecule(rotmol._mol)
        pmv.customColor(rotmol._mol.select('element C'), [(.2,8.,6.)])
        pmv.undisplayLines(rotmol._mol.select('name N C O'))
        self._setRotamer(resItem,1)

        ## add the rotamer combobox
        resItem.comboBox = QtGui.QComboBox()
        self.treeWidget.setItemWidget(resItem, 1, resItem.comboBox)
        
    def _setRotamer(self, resItem, index):
        ## callback function for residue with rotamer object on side chain
        rotmol = resItem._rotmol
        pmv = self.gui.app.pmv
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
        # display message about the curent rotamer
        if hasattr(resItem, '_rotaScores'):
            self.gui.app.statusbar.showMessage('Score: %.1f, contacts: %d favorable %d clashes'%(
                resItem._rotaScores[0][index], len(resItem._rotaScores[1][index]),
                len(resItem._rotaScores[2][index])))
        #pmv.undisplayLines(rotmol._mol)
        #pmv.displayLines(rotmol._mol.select('not name N C O'))

        # update the Viewer to reflex the new coordinates, whatever the geometry is
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
            item = chItems.get('metals', None)
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


class ReceptorItemsDialog(QtGui.QDialog):
    def __init__(self, profile, title='receptor items selector', parent=None):
        super(ReceptorItemsDialog, self).__init__(parent)
        self.setWindowTitle(title)
        self.recItemsWidget = ReceptorItemsTree(profile, parent)
        #self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.recItemsWidget)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)

    def reject(self):
        QtGui.QDialog.reject(self)
        self.closeEvent()
        
    def closeEvent(self, evnt=None):
        print 'closed'

            

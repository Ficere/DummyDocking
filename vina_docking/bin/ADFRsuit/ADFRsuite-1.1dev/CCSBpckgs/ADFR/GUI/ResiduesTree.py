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
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/ResiduesTree.py,v 1.16.2.2 2017/08/02 00:27:34 annao Exp $
#
# $Id: ResiduesTree.py,v 1.16.2.2 2017/08/02 00:27:34 annao Exp $
#
from PySide import QtCore, QtGui
from ADFR.GUI import ICONPATH
from prody.measure.contacts import findNeighbors
import os, numpy

from MolKit2.selection import Selection
from mglutil.util.callback import CallbackFunction

class ResiduesTree(QtGui.QWidget):

    resChanged = QtCore.Signal()

    def __init__(self, app, rec, lig=None, checkedItems=None, parent=None, sites=None,altloc=None, mutated=None, mse=None, flexResiduesNames=None):
        # checkedItems is a a dict where the key is 'chid:resnumIcode' and the
        # value is True or False, only residues with a key in checked are shown

        # flexResiduesNames is a list of resnames that can be made flexible
        super(ResiduesTree, self).__init__(parent)
        self._suspend = False
        self._app = app
        self._resToItem = {}
        self._resCA = []
        self._CAindices = set([])
        self.receptor = rec
        self.ligand = None
        self.checkedItems = checkedItems
        self._sites = sites
        self._mutated = mutated
        self._altloc = altloc
        self._mse = mse
        self.resCategories = {}
        self._flexResiduesNames = flexResiduesNames
        self.buildUI()
        
        if lig:
            self.setLigand(lig)

    def onCheckGroup(self):
        if self.groupBox1.isChecked():
            self.selectClose()
        
    def buildUI(self):
        self.groupBox1 = QtGui.QGroupBox(self)
        self.groupBox1.setTitle(self.tr("residues within"))
        self.groupBox1.setCheckable(True)
        self.groupBox1.setChecked(False)
        self.groupBox1.clicked.connect(self.onCheckGroup)
        self.distanceWidget = w = QtGui.QDoubleSpinBox()
        w.setToolTip("select residues with atoms within this distance of ligand atoms")
        w.setMinimum(0.1)
        w.setSingleStep(0.2)
        w.setDecimals(2)
        w.setValue(3.0)
        self.groupBox1.setDisabled(True)
        w.valueChanged.connect(self.selectClose)

        self.buildTree()

        layout = QtGui.QVBoxLayout()
        layout1 = QtGui.QHBoxLayout()
        self.groupBox1.setLayout(layout1)
        layout.addWidget(self.groupBox1)
        layout1.addWidget(self.distanceWidget)
        layout1.addWidget(QtGui.QLabel(' of ligand'))
        if self._mse is not None:
            # mse is a selection
            mse = []
            for res in self._mse.getHierView().iterResidues():
                resKey = '%s:%s%s%s'%(res.getChid(), res.getResname(), res.getResnum(), res.getIcode())
                if self._resToItem.has_key(resKey):
                    mse.append(resKey)
            if len(mse):
                self.resCategories['MSE'] = mse
                #print "RES TREE MSE:", mse
        if self._mutated is not None:
            mut = []
            constr = []
            for key, rnames  in self._mutated.items():
                chid, resname, resnum, icode = key.split(":")
                if icode == "_": icode = ""
                if rnames[0] == rnames[1]:
                    resKey = '%s:%s%s%s'%(chid,rnames[0],resnum,icode)
                    if self._resToItem.has_key(resKey): constr.append(resKey)
                else:
                    resKey = '%s:%s%s%s'%(chid,rnames[1],resnum,icode)
                    if self._resToItem.has_key(resKey):
                        if rnames[0] != "MSE":
                            mut.append(resKey)
            if len(mut): self.resCategories['mutated'] = mut
            if len(constr): self.resCategories['constructed'] = constr
            #print 'RES TREE MUT:', mut
            #print 'RES TREE const:', constr
                    
        if self._altloc is not None:
            altloc = []
            for chid, reslist in self._altloc.items():
                for resname, resnum, icode in reslist:
                    resKey = '%s:%s%s%s'%(chid,resname,resnum,icode)
                    if self._resToItem.has_key(resKey): altloc.append(resKey)
            if len(altloc): self.resCategories['alternate locations'] = altloc
            #print 'RES TREE ALTLOC', altloc
        
            
        if len(self.resCategories):
            self.groupBox2 = QtGui.QGroupBox(self)
            self.groupBox2.setTitle("select/deselect residues with:")
            self.groupBox2.setCheckable(False)
            sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
            self.groupBox2.setSizePolicy(sizePolicy)
        
            layout2 = QtGui.QGridLayout()
            layout2.setContentsMargins(0,0,0,0)
            #layout2.setHorizontalSpacing(0)
            self.groupBox2.setLayout(layout2)
            buttons = []
            count = 0
            for name in ['alternate locations', 'mutated', 'constructed', "MSE"] :
                resdict = self.resCategories.get(name, None)
                if not resdict: continue
                lab = QtGui.QLabel("%s: "% name)
                txt = "+"
                b1 = QtGui.QPushButton("+")
                width = b1.fontMetrics().boundingRect(txt).width() + 12
                b1.setMaximumWidth(width)
                b1.setMaximumHeight(width)
                b2 = QtGui.QPushButton("-")
                b2.setMaximumWidth(width)
                b2.setMaximumHeight(width)
                b1.clicked.connect(CallbackFunction(self._toggleSelectResCategory, name, True))
                b2.clicked.connect(CallbackFunction(self._toggleSelectResCategory, name, False))
                layout2.addWidget(lab, count, 0, 1, 1)
                layout2.addWidget(b1, count, 1, 1, 1, QtCore.Qt.AlignLeft)
                layout2.addWidget(b2, count, 2, 1, 1, QtCore.Qt.AlignLeft)
                spacer = QtGui.QWidget()
                spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
                layout2.addWidget(spacer, count, 3, 1, 1, QtCore.Qt.AlignLeft)
                layout2.setRowStretch(count, 1)
                buttons.extend([b1, b2])
                count += 1
            layout.addWidget(self.groupBox2)
            for i, b in enumerate(buttons):
                # one of the buttons is always selected (highlighted). To prevent this:
                b.setAutoDefault(False)
                b.setDefault(False)
        layout.addWidget(self.treeWidget)
        self.setLayout(layout)
        #treeWidget.itemClicked is called when we click on box or label
        # treeWidget.itemClicked.connect(self.onClick)

        # treeWidget.itemChanged is only called when box toggles
        self.treeWidget.itemChanged.connect(self.onClick)



    def buildTree(self):
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        hv = self.receptor.getHierView()
        chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
        resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
        _CAindices = []
        #print "Checked Items:", self.checkedItems
        self._sitesResToItem ={}
        if self._sites:
            siteRootItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), ["binding sites for:"])
            siteRootItem.setExpanded(True)
            #we will add the items under this tree branch later
                    
        for chain in hv.iterChains():
            chid = chain.getChid()
            chItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), chid)
            chItem.setExpanded(True)
            chItem.setIcon(0,chIcon)
            
            for res in chain.iterResidues():
                resnum = res.getResnum()
                icode = res.getIcode()
                key1 = '%s:%d%s'%(chid,resnum,icode)
                if icode is None: icode = ""
                if self.checkedItems and key1 not in self.checkedItems:
                    continue
                resname = res.getResname()
                if self._flexResiduesNames and resname not in self._flexResiduesNames:
                    continue
                resItem = QtGui.QTreeWidgetItem(chItem )
                resItem.setText(0, '%s%d%s'%(resname,resnum,icode ))
                resItem.setIcon(0, resIcon)
                resItem._pmvObj = res
                if self.checkedItems:
                    if self.checkedItems['%s:%d%s'%(chid,resnum,icode)]:
                        #sel = res.select('ca')
                        #if sel is not None:
                        if len(icode):
                            sel = res.select("chain %s resnum `%d` icode %s"%(chid,resnum,icode))
                        else:
                            sel = res.select("chain %s resnum `%d`"%(chid,resnum))
                        #self._CAindices.append(sel.getIndices()[0])
                        _CAindices.extend(sel.getIndices())
                        resItem.setCheckState(0, QtCore.Qt.Checked)
                    else:
                        resItem.setCheckState(0, QtCore.Qt.Unchecked)
                else:
                    resItem.setCheckState(0, QtCore.Qt.Unchecked)
                key = '%s:%s%d%s'%(chid,resname,resnum,icode)
                self._CAindices = set(_CAindices)
                self._resToItem[key] = resItem
        if self._sites:
            # now add the site items
            for sitename, residues in self._sites.items():
                siteItem = QtGui.QTreeWidgetItem(siteRootItem)
                siteItem.setText(0, sitename)
                siteItem.setCheckState(0, QtCore.Qt.Unchecked)
                chainItems = {}
                for res in residues:
                    chid, resname, resnum, icode = res.getChid(), res.getResname(), res.getResnum(), res.getIcode()
                    key = '%s:%s%d%s'%(chid,resname,resnum,icode)
                    resItem = self._resToItem.get(key)
                    if resItem:
                        if not chainItems.has_key(chid):
                            chItem = QtGui.QTreeWidgetItem(siteItem, chid)
                            chainItems[chid] = chItem
                            chItem.setExpanded(False)
                            chItem.setIcon(0,chIcon)
                        else:
                            chItem = chainItems[chid]
                        siteResItem = QtGui.QTreeWidgetItem(chItem )
                        siteResItem.setText(0, '%s%d%s'%(resname,resnum,icode ))
                        siteResItem.setIcon(0, resIcon)
                        siteResItem._pmvObj = resItem._pmvObj
                        self._sitesResToItem[key] = siteResItem
                        siteResItem.setCheckState(0, resItem.checkState(0))

    def selectClose(self):
        #if self.ligand is None:
        #    return

        atoms = findNeighbors(self.ligand._ag,
                              self.distanceWidget.value(),
                              self.receptor)
        if len(atoms):
            self._suspend = True
            indices = numpy.unique([x[1].getIndices() for x in atoms])
            resnums = self.receptor._ag._data['resnum']
            icodes = self.receptor._ag._data['icode']
            chids = self.receptor._ag._data['chain']
            if self.checkedItems:
                indicesInBox = []
                for ind in indices:
                    if self.checkedItems.has_key('%s:%d%s'%(chids[ind],
                                                            resnums[ind],icodes[ind])):
                        indicesInBox.append(ind)
                indices = indicesInBox
            selStr = '(same residue as index %s) and ca'%' '.join([str(x) for x in indices])
            self._resCA = resCA = self.receptor.select(selStr)
            #for atom in resCA:
            #    print '%s%d '%(atom.getResname(), atom.getResnum())
            # uncheck all
            for item in self._resToItem.values():
                item.setCheckState(0, QtCore.Qt.Unchecked)
            # check close
            for i, ca in enumerate(resCA):
                key = '%s:%s%d%s'%(ca.getChid(), ca.getResname(), ca.getResnum(), ca.getIcode())
                self._resToItem[key].setCheckState(0, QtCore.Qt.Checked)
                if i==0:
                    # make first check residue visible
                    self.treeWidget.scrollToItem(self._resToItem[key])
            #print len(atoms), len(resCA), time()-t0
            #self._CAindices = resCA.getIndices().tolist()
            self._CAindices = set(resCA.getIndices())
            self._suspend = False
        else:
            self._resCA = []
            self._CAindices = set([])
        self.resChanged.emit()

    def _toggleSelectResCategory(self, category, select):
        #import pdb; pdb.set_trace()
        self._suspend = True
        _CAindices = []
        #print "TOGGLEselect in box:", self.resCategories[category]
        count = 1
        for key in self.resCategories[category]:
            resItem = self._resToItem.get(key)
            if resItem:
                if hasattr(resItem, "_pmvObj"):
                    res = resItem._pmvObj
                    if select:
                        #check if the residue inside the box:
                        #inside, outside, inds = self._app.pointsInBox(res.getCoords())
                        #if not len(outside):
                        _CAindices.extend(res.getIndices())
                        resItem.setCheckState(0, QtCore.Qt.Checked)
                        if count == 1:
                            self.treeWidget.scrollToItem(resItem)
                        count =+ 1
                    else:
                        _CAindices.extend(res.getIndices())
                        resItem.setCheckState(0, QtCore.Qt.Unchecked)
        if len(_CAindices):
            if select:
                self._CAindices.update(_CAindices)
            else:
                self._CAindices = self._CAindices - set(_CAindices)
            self.resChanged.emit()
        self._suspend = False
                    
    def getAtoms(self):
        if len(self._CAindices)==0:
            return None
        else:
            return self.receptor.select(
                'same residue as index %s'%' '.join(
                    [str(x) for x in self._CAindices]))

    def setLigand(self, ligand):
        self.ligand = ligand
        self.groupBox1.setDisabled(False)

    def onClick(self, item, column):
        if self._suspend: return
        checked = item.checkState(column)
        #print "ResiduesTree.onClick(), item text:", item.text(0), column
        level = 0
        _item = item
        while _item.parent():
            _item = _item.parent()
            level += 1
        if level == 1:
            if str(_item.text(0)).find("binding sites") >= 0:
                self._suspend = True
                inds = []
                for i in range(item.childCount()):
                    chItem = item.child(i)
                    for j  in range(chItem.childCount()):
                        siteResItem = chItem.child(j)
                        res = siteResItem._pmvObj
                        key = '%s:%s%d%s'%(res.getChid(), res.getResname(), res.getResnum(), res.getIcode())
                        resItem = self._resToItem.get(key, None)
                        if resItem:
                            resItem.setCheckState(0, checked)
                            inds.extend(resItem._pmvObj.getIndices())
                        siteResItem.setCheckState(0, checked)
                self._suspend = False
            else: # 
               res = item._pmvObj # this item is a residue item under a chain tree item(not in sites)
               inds = res.getIndices()
               # find out if we have sites and the sites have the same res item:
               key = '%s:%s%d%s'%(res.getChid(), res.getResname(), res.getResnum(), res.getIcode())
               sitesResItem = self._sitesResToItem.get(key, None)
               if sitesResItem:
                   self._suspend = True
                   sitesResItem.setCheckState(0, checked)
                   self._suspend = False
        elif level == 3: # this is a res item under "sites" tree branch
            res = item._pmvObj
            inds = res.getIndices()
            key = '%s:%s%d%s'%(res.getChid(), res.getResname(), res.getResnum(), res.getIcode())
            resItem = self._resToItem.get(key, None)
            if resItem: # check/uncheck the corresponding res item under the chains
                self._suspend = True
                resItem.setCheckState(0, checked)
                self._suspend = False 
        if checked == QtCore.Qt.Checked:
            #self._CAindices.append(ca.getIndices()[0])
            self._CAindices.update(inds)
        else:
            #self._CAindices.remove(ca.getIndices()[0])
            self._CAindices = self._CAindices - set(inds)
        self.resChanged.emit()

class ResiduesTreeDialog(QtGui.QDialog):
    closedSignal = QtCore.Signal()

    def __init__(self, app, rec, lig=None, checkedItems=None, title='No Name',
                 sites=None, altloc=None, mutated=None, mse=None, parent=None, flexResiduesNames=None):
        # sites - binding sites for residues (a dictionary that comes from PDBAnalyser) 
        
        super(ResiduesTreeDialog, self).__init__(parent)
        self.setWindowTitle(title)

        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)

        self.resTreeWidget = ResiduesTree(app, rec, lig, checkedItems, parent, sites, altloc=altloc,
                                          mutated=mutated, mse=mse, flexResiduesNames=flexResiduesNames)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.resTreeWidget)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)
        #self.setWindowFlags( self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)

    def reject(self):
        QtGui.QDialog.reject(self)
        self.closeEvent()
        
    def closeEvent(self, evnt=None):
        self.closedSignal.emit()
        
if __name__=='__main__':
    import sys

    from MolKit2 import Read
    rec = Read('4EK3_rec.pdbqt')
    lig = Read('4EK4_lig.pdbqt')
    app = QtGui.QApplication(sys.argv)
    from time import time
    t0 = time()
    widget = ResiduesTree(rec, lig)
    #widget.setLigand(lig)
    #def p(chid, resname, checked):
    #    print chid, resname, checked
    #widget.resClickedSignal.connect(p)
    
    #print time()-t0
    widget.resize(150,300)
    widget.show()
    sys.exit(app.exec_())

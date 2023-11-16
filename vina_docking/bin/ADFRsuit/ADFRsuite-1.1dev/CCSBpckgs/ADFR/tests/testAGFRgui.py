import os, sys
import unittest
import numpy as np

from PySide import QtCore, QtGui
from PySide.QtTest import QTest
from PySide.QtCore import Qt

from ADFR.GUI import makeGridGUI, _AGFRGUI_debug
_AGFRGUI_debug = True  # avoid @waiting hiding exception

app = QtGui.QApplication(sys.argv)

class GridGUITest(unittest.TestCase):
    '''Test the margarita mixer GUI'''
    def setUp(self):
        '''Create the GUI'''
        self.app = makeGridGUI(app)
        self.app.closeEvent = self.app.closeEventNoAsk # to avoid asking are you sure
        
    def tearDown(self):
        '''Create the GUI'''
        self.app.close()

    def test_loadReceptorPDBQT(self):
        self.app.loadReceptor('ADFR/tests/Data/4EK3_rec.pdbqt')
        assert self.app.receptor 
        assert len(self.app.pmv.Mols)==1
        assert self.app.pmv.Mols[0].name=='4EK3_rec'
        
    def test_loadReceptorPDB(self):
        self.app.loadReceptor('ADFR/tests/Data/6bhg.pdb')
        assert self.app.PDBprofile is not None

        # check receptor in PDBprofile
        rec = self.app.PDBprofile._receptor._ag
        assert len(rec.select('receptor'))==1750
        assert rec.select('ligand') is None
        assert len(rec.select('unsat'))==40
        assert len(rec.select('unsre')) == 40
        assert len(rec.select('altlocs')) == 175
        assert set(rec.getData('category')) == set(['_wate', 'cmpsc', 'missa', 'other', 'prote', 'modif'])
        assert len(rec.select('category prote')) == 1665
        assert len(rec.select('category _wate')) == 223
        assert len(rec.select('category cmpsc')) == 127
        assert len(rec.select('category missa')) == 2
        assert len(rec.select('category modif')) == 12
        assert len(rec.select('category other')) == 40

        assert len(self.app.PDBprofile._assemblies)==1
        misatkeys = self.app.PDBprofile._missingAtoms.keys()
        assert set(misatkeys) == set(
            ['B:SER:10:_', 'A:LYS:257:_', 'A:ALA:406:_', 'A:ARG:384:_',
             'A:LYS:290:_', 'A:LYS:194:_', 'A:ASP:256:_', 'A:LYS:288:_',
             'A:ASN:236:_', 'A:LYS:237:_', 'A:LYS:349:_', 'A:TRP:275:_',
             'A:LYS:233:_', 'A:ARG:209:_', 'A:LYS:364:_', 'A:SER:405:_',
             'B:LYS:18:_', 'A:LYS:239:_', 'A:ASP:270:_', 'B:ARG:17:_',
             'A:LYS:402:_'])
        mutateKeys = self.app.PDBprofile._resToMutate.keys()
        assert set(mutateKeys) == set([
            'A:ARG:209:_', 'A:ARG:384:_', 'A:ASN:236:_', 'A:ASP:256:_',
            'A:ASP:270:_', 'A:LYS:194:_', 'A:LYS:233:_', 'A:LYS:237:_',
            'A:LYS:239:_', 'A:LYS:257:_', 'A:LYS:288:_', 'A:LYS:290:_',
            'A:LYS:349:_', 'A:LYS:364:_', 'A:LYS:402:_', 'A:SER:405:_',
            'A:TRP:275:_', 'B:ARG:17:_' , 'B:LYS:18:_' , 'B:SER:10:_'])

        # check PDBAnalyzer widget
        assert self.app.analyserWidget is not None
        assert len(self.app.analyserWidget.treeWidgets) == 1
        tree = self.app.analyserWidget.treeWidgets.values()[0].treeWidget
        root = tree.invisibleRootItem()
        assert root.childCount()==2
        assert set([str(root.child(n).text(0)) for n in range(root.childCount())]) == set(['A', 'B'])
        AItem = root.child(0)
        catA = [str(AItem.child(n).text(0)) for n in range(AItem.childCount())]
        assert set(catA) == set(['protein', 'sidechains to complete', 'water',
                                 'missing residues', 'altloc',
                                 'unknown atom types'])
        altCatItem = AItem.child(catA.index('altloc'))
        altcats = [str(altCatItem.child(n).text(0)) for n in range(altCatItem.childCount())]
        assert set(altcats) == set(['protein', 'sidechains to complete',
                                    'missing atoms', 'water'])

        # check alternate location state of ILE 204 in chain A
        item = altCatItem.child(altcats.index('protein'))
        ile204 = [item.child(n) for n in range(item.childCount()) if str(item.child(n).text(0))=="ILE'204"][0]
        assert hasattr(ile204, 'altWidget')
        assert ile204.altWidget.getAltlocChar()=='B'
        assert len(ile204._pmvObj) == 12
        assert len(ile204._pmvObj.select('receptor'))==8
        assert set(ile204._pmvObj.select('receptor').getAltlocs()) == set([' ', 'B'])
        
        # check alternate location state of SER 10 in chain B
        BItem = root.child(1)
        catB = [str(BItem.child(n).text(0)) for n in range(BItem.childCount())]
        assert set(catB) == set(['protein', 'sidechains to complete', 'water',
                                 'missing residues', 'altloc',
                                 'modified'])
        altCatItem = BItem.child(catB.index('altloc'))
        altcats = [str(altCatItem.child(n).text(0)) for n in range(altCatItem.childCount())]
        assert set(altcats) == set(['sidechains to complete'])
        ser10 = altCatItem.child(0).child(0)
        assert hasattr(ser10, 'altWidget')
        assert ser10.altWidget.getAltlocChar()=='B'
        assert len(ser10._pmvObj) == 10
        assert len(ser10._pmvObj.select('receptor')) == 5
        assert set(ser10._pmvObj.select('receptor').getAltlocs()) == set(['B'])

        ## make chain B ligand (except for water)
        # get handle to ReceptorItemsTree
        recItemTree = self.app.analyserWidget.treeWidgets.values()[0]
        # get tree item for protein category
        proteinItem = BItem.child(0)

        # check that ChainB: protein is currently in receptor
        assert proteinItem.checkState(0)==QtCore.Qt.CheckState.Checked
        # and not in ligand
        assert proteinItem.checkState(2)==QtCore.Qt.CheckState.Unchecked
        
        ## make chain B ligand (except water)
        for child in [BItem.child(n) for n in range(BItem.childCount())]:
            ## simulate clicking on ligand check button
            # set ligand check button checkSate to checked and call associated cb
            if str(child.text(0))=='water': continue
            child.setCheckState(2, QtCore.Qt.Checked)
            recItemTree.onClick(child, 2)
            # no longer receptor
            assert child.checkState(0)==QtCore.Qt.CheckState.Unchecked
        
        assert len(rec.select('ligand')) == 54
        # makes sure we have not duplicates in B:SER10 which has altlocs
        assert set(rec.select('resnum 10 and ligand').getNames()) == set(
            ['N', 'CA', 'C', 'O', 'CB'])

        # delete file we expect to create (if they exist)
        if os.path.exists('6bhg_lig.pdbqt'):
            os.remove('6bhg_lig.pdbqt')
        if os.path.exists('6bhg_rec.pdbqt'):
            os.remove('6bhg_rec.pdbqt')

        ## proceed to maps
        QTest.mouseClick(self.app.analyserWidget.proceedToMapsBt, Qt.LeftButton)
        
        # ensure files were created
        assert os.path.exists('6bhg_lig.pdbqt')
        assert os.path.exists('6bhg_rec.pdbqt')

        # check receptor molecule
        assert self.app.receptor.name == '6bhg_rec'
        assert len(self.app.receptor._ag) == 2146

        # check ligand
        assert self.app.ligand.name == '6bhg_adlig'
        assert len(self.app.ligand._ag) == 85

        # make sure side chain of SER 10 was built
        ser10 = self.app.ligand._ag.select('resname SER')
        # FIXEM .. why are the names not stripped ?
        assert set(ser10.getNames())==set(['N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'H', 'H', 'H'])
        # make sure we got alternate location B for SER 10
        assert abs(np.sum(ser10.getCoords()[0]-[48.598, -3.4, 62.474])) < 0.0001

    def test_loadPDBSaveLigand(self):
        self.app.loadReceptor('ADFR/tests/Data/6bhg.pdb')
        assert self.app.PDBprofile is not None
        tree = self.app.analyserWidget.treeWidgets.values()[0].treeWidget
        root = tree.invisibleRootItem()
        BItem = root.child(1)
        recItemTree = self.app.analyserWidget.treeWidgets.values()[0]
        proteinItem = BItem.child(0)
        # make chain B ligand (except water)
        for child in [BItem.child(n) for n in range(BItem.childCount())]:
            if str(child.text(0))=='water': continue
            child.setCheckState(2, QtCore.Qt.Checked)
            recItemTree.onClick(child, 2)

        rec = self.app.PDBprofile._receptor._ag
        import prody
        frags = prody.findFragments(rec.select('ligand and not deleted'))

        # FIXME we should also test cancel
        #
        if os.path.exists('6bhg_lig.pdbqt'):
            os.remove('6bhg_lig.pdbqt')
        # click save ligand
        QTest.mouseClick(self.app.analyserWidget.saveLigAsBt, Qt.LeftButton)

        
        #import pdb; pdb.set_trace()
        
if __name__ == "__main__":
    unittest.main()

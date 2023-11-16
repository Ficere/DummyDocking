import os
from ADFR.utils.runAGFR import runAGFR
from MolKit2 import Read
from ADFR.tests.defaultkw import kw_base


def test1():
    # test that icons get default charges
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data/1jff.pdb', 'toPdbqt': True,
               "heteroToKeep":"chid A resnum 501 900"})
    if os.path.exists('1jff_rec.pdbqt'):
        os.remove('1jff_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)
    mol = Read('1jff_rec.pdbqt')
    zn = mol._ag.select('resname ZN')
    assert zn.getCharges()[0] == 2.0
    mg = mol._ag.select('resname MG')
    assert mg.getCharges()[0] == 2.0

def test2():
    # test that we can set the chat from the command line
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data/1jff.pdb', 'toPdbqt': True,
               "heteroToKeep":"chid A resnum 501 900",
               "ionCharges":":ZN:::1.5;:MG:501::-1.2"})
    if os.path.exists('1jff_rec.pdbqt'):
        os.remove('1jff_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)
    mol = Read('1jff_rec.pdbqt')
    zn = mol._ag.select('resname ZN')
    assert zn.getCharges()[0] == 1.5
    mg = mol._ag.select('resname MG')
    assert mg.getCharges()[0] == -1.2

    
if __name__=='__main__':
    test1()
    test2()

import os
from ADFR.utils.runAGFR import runAGFR
from MolKit2 import Read
from ADFR.tests.defaultkw import kw_base

def test1():
    # test that we can process a ligand
    kw = kw_base.copy()
    #kw.update({"receptorFile":'ADFR/tests/Data/1jff.pdb', 'toPdbqt': True,
    #           #"heteroToKeep":"chid B resname TA resnum 1601",
    #           "pdbLigand":"chid B resname TA1 resnum 601"})
    kw.update({"receptorFile":'ADFR/tests/Data/6ggd.pdb', 'toPdbqt': True,
               #"heteroToKeep":"chid A resname TA resnum 1601",
               "pdbLigand":"chid A resname EYB resnum 402"})
    if os.path.exists('6ggd_lig.pdbqt'):
        os.remove('6ggd_lig.pdbqt')
    runner = runAGFR()
    runner(**kw)
    assert os.path.exists('6ggd_lig.pdbqt')
    mol = Read('6ggd_lig.pdbqt')
    
if __name__=='__main__':
    test1()

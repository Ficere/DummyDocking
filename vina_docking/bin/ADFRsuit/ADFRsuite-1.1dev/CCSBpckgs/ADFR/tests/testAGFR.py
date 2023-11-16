import os
import numpy as np

from ADFR.utils.runAGFR import runAGFR
from MolKit2 import Read
from ADFR.tests.defaultkw import kw_base

def test1():
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data/4a1x.pdb', 'toPdbqt': True})
    if os.path.exists('4a1x_rec.pdbqt'):
        os.remove('4a1x_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)
    mol = Read('4a1x_rec.pdbqt')
    assert len(mol._ag)==1897
    # make sure arg 180 has been rebuilt
    arg180 = mol._ag.select('chid A resnum 180')
    assert len(arg180.getNames())==17

def test1_1():
    # test was we can exclude chain C
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data/4a1x.pdb', 'toPdbqt': True,
               'protToRemove': 'chid C'})
    if os.path.exists('4a1x_rec.pdbqt'):
        os.remove('4a1x_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)
    mol = Read('4a1x_rec.pdbqt')
    assert len(mol._ag)==1742
    chids = set(mol._ag.getChids())
    assert len(chids)==1 and list(chids)[0]=='A'
    
    
def test2():
    # test that we can select Biomol 2 rather than the default 1
    # keep hetero ZN 401 and EYB402
    # mutate residues with and without alternate locations
    kw = kw_base.copy()

    # keep ligand A:402 and ZN401
    heteroToKeep="chid B resnum 402 401"

    # mutate SER96 to ILE, PRO98 to ALA, GLN100 to PHE
    # mutate B:SER166 which has altlocs A and B (B selected by default)
    # mutate B:SER182 which has altlocs A and B (A selected by default)
    mutations = "B:SER:96::ILE; B:PRO:98::ALA; B:GLN:100::PHE; B:SER:166::ALA; B:CYS:182::ALA"
    #mutations = None

    kw.update({"receptorFile":'ADFR/tests/Data/6ggd.pdb', 'toPdbqt': True,
               'bioMol': '2', 'heteroToKeep':heteroToKeep, 'mutateRes':mutations})
    if os.path.exists('6ggd_rec.pdbqt'):
        os.remove('6ggd_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)

    mol = Read('6ggd_rec.pdbqt')
    # check that we selected Biomol2 whch is chain B
    set(mol._ag.getChids())==set(['B'])
    
    # check that B: ZN401 adn EYB402 have been kept
    zn = mol._ag.select('resnum 401')
    assert set(zn.getNames())==set(['ZN'])

    # fails because we have no charges from cmdline
    #assert zn.getCharges()[0] = 1.0 # fails 

    lig = mol._ag.select('resnum 402')
    assert lig is not None
    assert len(lig)==23
    assert set(lig.getResnames())==set(['EYB'])
         
    # make sure B:LYS120,ARG209,GLU221,ARG290,LYS291 has been completed
    assert set(mol._ag.select('resnum 120 sc').getNames())==set(['CB', 'CG', 'CE', 'CD', 'NZ', 'HZ1', 'HZ3', 'HZ2'])
    assert set(mol._ag.select('resnum 209 sc').getNames())==set(['CB', 'CD', 'NE', 'CG', 'NH1', 'NH2', 'CZ', 'HE', '1HH1', '2HH1', '1HH2', '2HH2'])
    assert set(mol._ag.select('resnum 221 sc').getNames())==set(['CB', 'CD', 'OE1', 'OE2', 'CG'])

    # check that the mutations occurred
    res96 = mol._ag.select('resnum 96 sc')
    assert set(res96.getNames())==set(['CB', 'CD1', 'CG1', 'CG2'])
    assert set(res96.getResnames())==set(['ILE'])

    res98 = mol._ag.select('resnum 98 sc')
    assert set(res98.getNames())==set(['CB'])
    assert set(res98.getResnames())==set(['ALA'])

    res100 = mol._ag.select('resnum 100 sc')
    assert set(res100.getNames())==set(['CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ'])
    assert set(res100.getResnames())==set(['PHE'])

    res166 = mol._ag.select('resnum 166 sc')
    assert set(res166.getNames())==set(['CB'])
    assert set(res166.getResnames())==set(['ALA'])

    res182 = mol._ag.select('resnum 182 sc')
    assert set(res182.getNames())==set(['CB'])
    assert set(res182.getResnames())==set(['ALA'])

def test3():
    # test that we can select Biomol 'ASU'
    # mutate residue with Icode
    # 1a3r chain L has SER27A LEU27B and LEU27C
    kw = kw_base.copy()

    # mutate L:LEU27C and make sure L:LEU27B is not affected
    mutations = "L:LEU:27:C:ALA"

    kw.update({"receptorFile":'ADFR/tests/Data/1a3r.pdb', 'toPdbqt': True,
               'bioMol': 'ASU', 'mutateRes':mutations})
    
    if os.path.exists('1a3r_rec.pdbqt'):
        os.remove('1a3r_rec.pdbqt')
    runner = runAGFR()
    runner(**kw)
    mol = Read('1a3r_rec.pdbqt')

    res27A = mol._ag.select('chid L resnum 27 icode A sc')
    assert set(res27A.getNames())==set(['CB', 'OG', 'HG'])
    assert set(res27A.getResnames())==set(['SER'])
    assert set(res27A.getIcodes())==set(['A'])

    res27B = mol._ag.select('chid L resnum 27 icode B sc')
    assert set(res27B.getNames())==set(['CB', 'CG', 'CD1', 'CD2'])
    assert set(res27B.getResnames())==set(['LEU'])
    assert set(res27B.getIcodes())==set(['B'])

    res27C = mol._ag.select('chid L resnum 27 icode C sc')
    assert set(res27C.getNames())==set(['CB'])
    assert set(res27C.getResnames())==set(['ALA'])
    assert set(res27C.getIcodes())==set(['C'])


def test4():
    # test that we compute the maps
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data/1jff.pdb', 'bioMol': '1', 'receptorGradient': False,
               "boxMode": ["user", 0.908, -16.044, 14.073, 15.0, 15.0, 15.0]})

    if os.path.exists('1jff_rec.trg'):
        os.remove('1jff_rec.trg')
    runner = runAGFR()
    runner(**kw)
    assert os.path.exists('1jff_rec.trg')
    from ADFR.utils.maps import MapsFile
    mf = MapsFile('1jff_rec.trg')
    data = mf.getData()
    assert np.sum(data['boxCenter']-np.array([0.908, -16.044, 14.073]))< 0.001 
    assert np.sum(data['boxLengths']-np.array([15.0, 15.0, 15.0]))< 0.001 
    assert data['mapGradients'] == False
    assert data['covalentBond'] is None
    assert data['boxMode'][0] == 'user'
    assert data['pocketmode'][0] == 'all'
    assert data['inputReceptor'] == '1jff_rec.pdbqt'
    assert data['flexResStr'] == None
    assert data['spacing'] == 0.375
    assert data['boxPadding'] == 0.0
    assert data['covalentLigandFile'] == ''
    assert data['covalentRes'] is None
    assert data['gradCutOff'] is None
    assert data['wMapEntropy'] == -0.2
    assert np.sum(data['boxSize']-np.array([40, 40, 40]))==0
    assert data['flexRecFile'] == ''
    assert data['covalentLigandAtomIndices'] is None
    assert data['wMapWeight'] == 0.6

def test5():
    # test flexible residues option:
    # --flexRes selectionString 
    kw = kw_base.copy()
    kw.update({"receptorFile":"ADFR/tests/Data/4EK3_rec.pdbqt",
               "ligandFile": "ADFR/tests/Data/4EK4_lig.pdbqt",
               "flexres": "A:ILE10,PHE82", 'receptorGradient': False, "outputFile": "4EK3_flx_str", }
            )
    if os.path.exists('4EK3_flx_str.trg'):
        os.remove('4EK3_flx_str.trg')
    runner = runAGFR()
    runner(**kw)
    assert os.path.exists('4EK3_flx_str.trg')
    from ADFR.utils.maps import MapsFile
    mf = MapsFile('4EK3_flx_str.trg')
    data = mf.getData()
    assert np.sum(data['boxCenter'] -np.array([ 23.332 ,  28.9215,  29.5975]))< 0.001 
    assert np.sum(data['boxLengths']-np.array([ 17.25,  16.5 ,  12.  ]))< 0.001
    assert data['flexResStr'] == 'A:ILE10,PHE82'
    assert data['mapGradients'] == False
    assert data['covalentBond'] is None
    assert data['boxMode'][0] == 'ligand'
    assert data['pocketmode'][0] == 'all'
    assert data['inputReceptor'] == '4EK3_rec.pdbqt'
    assert data['spacing'] == 0.375
    assert data['boxPadding'] == 4.0
    assert data['covalentLigandFile'] == ''
    assert data['covalentRes'] is None
    assert data['gradCutOff'] is None
    assert data['wMapEntropy'] == -0.2
    assert np.sum(data['boxSize']-np.array([46, 44, 32]))==0
    assert data['flexRecFile'] == 'flexRec.pdbqt'
    assert data['covalentLigandAtomIndices'] is None
    assert data['wMapWeight'] == 0.6

def test6():
    # test --boxMode user ligand 10 13 15 (center on ligand and cpecify box lengths)
    kw = kw_base.copy()
    kw.update({"receptorFile":"ADFR/tests/Data/4EK3_rec.pdbqt", "ligandFile": "ADFR/tests/Data/4EK4_lig.pdbqt",
               'receptorGradient': False, "outputFile": "4EK3_box_lig", "padding":5, "pocketMode":["best"],
               "boxMode": ["user", "ligand", 10, 13, 15]}
            )
    if os.path.exists('4EK3_box_lig.trg'):
        os.remove('4EK3_box_lig.trg')
    runner = runAGFR()
    runner(**kw)
    assert os.path.exists('4EK3_box_lig.trg')
    from ADFR.utils.maps import MapsFile
    mf = MapsFile('4EK3_box_lig.trg')
    data = mf.getData()
    assert np.sum(data['boxCenter'] -np.array([ 23.332 ,  28.9215,  29.5975]))< 0.001 
    assert np.sum(data['boxLengths']-np.array([ 10.5,  13.5,  15.]))< 0.001
    assert data['flexResStr'] == None
    assert data['mapGradients'] == False
    assert data['covalentBond'] is None
    assert data['boxMode'][0] == 'ligand'
    assert data['pocketmode'][0] == 'best'
    assert data['inputReceptor'] == '4EK3_rec.pdbqt'
    assert data['spacing'] == 0.375
    assert data['boxPadding'] == 0.
    assert data['covalentLigandFile'] == ''
    assert data['covalentRes'] is None
    assert data['gradCutOff'] is None
    assert data['wMapEntropy'] == -0.2
    assert np.sum(data['boxSize']-np.array([28, 36, 40]))==0
    assert data['flexRecFile'] == ''
    assert data['covalentLigandAtomIndices'] is None
    assert data['wMapWeight'] == 0.6


def test7():
    # test --boxMode residues B:VAL23,ASP26,GLU27,LEU217,HIS229 (box around specified residues)
    kw = kw_base.copy()
    kw.update({"receptorFile":"ADFR/tests/Data/1jff.pdb", 
               'receptorGradient': False, "outputFile": "1jff_box_res",  "pocketMode":["best"],
               "boxMode": ["residues", "B:VAL23,ASP26,GLU27,LEU217,HIS229"]}
            )
    if os.path.exists('1jff_box_res.trg'):
        os.remove('1jff_box_res.trg')
    runner = runAGFR()
    runner(**kw)
    #import pdb; pdb.set_trace()
    assert os.path.exists('1jff_box_res.trg')
    from ADFR.utils.maps import MapsFile
    mf = MapsFile('1jff_box_res.trg')
    data = mf.getData()
    assert np.sum(data['boxCenter'] -np.array([ -0.035 ,  -14.581  ,  13.591])) < 0.001 
    assert np.sum(data['boxLengths']-np.array([21.75, 25.5, 25.5]))< 0.001
    assert data['flexResStr'] == None
    assert data['mapGradients'] == False
    assert data['covalentBond'] is None
    assert data['boxMode'][0] == 'residues'
    assert data['pocketmode'][0] == 'best'
    assert data['inputReceptor'] == '1jff_rec.pdbqt'
    assert data['spacing'] == 0.375
    assert data['boxPadding'] == 4.0
    assert data['covalentLigandFile'] == ''
    assert data['covalentRes'] is None
    assert data['gradCutOff'] is None
    assert data['wMapEntropy'] == -0.2
    assert np.sum(data['boxSize']-np.array([58, 68,68]))==0
    assert data['flexRecFile'] == ''
    assert data['covalentLigandAtomIndices'] is None
    assert data['wMapWeight'] == 0.6
    #print "MAPS:", set(data["mapTypes"])
    assert set(data["mapTypes"]) == set(['C', 'GA', 'Fe', 'A', 'e', 'd', 'S', 'Q', 'F', 'SA', 'G', 'P', 'Mg', 'Cl', 'Zn', 'OA', 'NS', 'Ca', 'Br', 'Mn', 'OS', 'NA', 'J', 'I', 'HS', 'HD', 'H', 'Z', 'N', 'W', 'constr'])

def test8():
    #create  maps for N  top rancking pockets , and maptypes "ligand"
    kw = kw_base.copy()
    kw.update({"receptorFile":"ADFR/tests/Data/4EK3_rec.pdbqt", "ligandFile": "ADFR/tests/Data/4EK4_lig.pdbqt",
               'receptorGradient': False, "outputFile": "4EK3_ats", "padding":5, "pocketMode":["forTop"],
               "pocketCutoff":3,  "mapTypes": "ligand", #autoSiteVersion": 1.1
               }
            )
    from glob import glob
    for f in glob('4EK3_ats_pocket00?.trg'):
        os.remove(f)
    runner = runAGFR()
    runner(**kw)
    from ADFR.utils.maps import MapsFile
    numFPoints = [300, 24, 44]
    for i in range(3):
        assert os.path.exists('4EK3_ats_pocket00%d.trg'%i)

        mf = MapsFile('4EK3_ats_pocket00%d.trg'%i)
        data = mf.getData()
        assert data['flexResStr'] == None
        assert data['mapGradients'] == False
        assert data['covalentBond'] is None
        assert data['boxMode'][0] == 'ligand'
        assert data['pocketmode'][0] == "forTop"
        assert data['inputReceptor'] == '4EK3_rec.pdbqt'
        assert data['spacing'] == 0.375
        assert data['boxPadding'] == 5.
        assert data['covalentLigandFile'] == ''
        assert data['covalentRes'] is None
        assert data['gradCutOff'] is None
        assert data['wMapEntropy'] == -0.2
        assert np.sum(data['boxSize']-np.array([52, 50, 38]))==0
        assert np.sum(data['boxCenter'] -np.array([ 23.332 ,  28.9215,  29.5975])) < 0.001 
        assert np.sum(data['boxLengths']-np.array([19.5 ,  18.75,  14.25]))< 0.001
        assert data['flexRecFile'] == ''
        assert data['covalentLigandAtomIndices'] is None
        assert data['wMapWeight'] == 0.6
        assert data['AutoSiteVersion'] == '1.0'
        assert set(data["mapTypes"]) == set(['C', 'A', 'e', 'd', 'OA', 'Br', 'NA', 'HD', 'N', 'W','constr'])
        #print "FPOINTS:", i, data['nbFillPoints']
        assert data['nbFillPoints'] == numFPoints[i]

def test9():
    # test that the receptorFile (-r) option supports the name of the data folder
    # containing pdb and pdbqt files
    kw = kw_base.copy()
    kw.update({"receptorFile":'ADFR/tests/Data', 'toPdbqt': True,})
    files = []
    pdbf =[]
    pdbqt=[]
    pwd = os.getcwd()
    os.chdir('ADFR/tests/Data')
    for f in os.listdir('.'):
        if os.path.isfile(f):
            file, ext = os.path.splitext(f)
            if ext == ".pdb":
                pdbf.append(file)
            elif ext == ".pdbqt" and file.find("lig") < 0:
                pdbqt.append(file)
    os.chdir(pwd)
    #print pdbf
    #print pdbqt
    for f in pdbf:
        if os.path.exists('%s_rec.pdbqt'% (f,)):
            print "removing:", '%s_rec.pdbqt' % (f,)
            os.remove('%s_rec.pdbqt' % (f,))
    for f in pdbqt:
        if os.path.exists('%s.trg'% (f,)):
            print "removing:", '%s.trg' % (f,)
            os.remove('%s.trg' % (f,))
    #import pdb; pdb.set_trace()
    f = open("ADFR/bin/runAGFR.py")
    sc= f.read()
    import sys 
    sys.argv = ["ADFR/bin/runAGFR.py", "-r", 'ADFR/tests/Data', '--toPdbqt', "-ng"]
    exec(sc)
    #runner = runAGFR()
    for f in pdbf:
        assert os.path.exists('%s_rec.pdbqt'% (f,))
    for f in pdbqt:
        assert os.path.exists('%s.trg'% (f,))
    
    
if __name__=='__main__':
    test1()
    test1_1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()


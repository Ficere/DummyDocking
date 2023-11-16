import prody
import numpy as np
from prody.atomic.atomgroup import AtomGroup
from prody.atomic import ATOMIC_FIELDS
from prody.atomic.selection import Selection as ProDySel 

from MolKit2.selection import Selection as MolKitSel

def CIFToProdyAG(name, data):
    # cif is assumed to be  the data dict for a molecule as
    # returned by 
    atomgroup = AtomGroup(name)
    atdata = data.get('chem_comp_atom', None)
    if atdata is None:
        return atomgroup, {}
    asize = len(atdata['comp_id'])

    # check for missing coords
    coordinates = np.zeros((asize, 3), dtype=float)
    # DUM says no on flags but has ? for coords
    if data['chem_comp']['pdbx_ideal_coordinates_missing_flag'][0]=='N':
        coordinates[:,0] = ['0.' if x=='?' else x for x in atdata['pdbx_model_Cartn_x_ideal']]
        coordinates[:,1] = ['0.' if x=='?' else x for x in atdata['pdbx_model_Cartn_y_ideal']]
        coordinates[:,2] = ['0.' if x=='?' else x for x in atdata['pdbx_model_Cartn_z_ideal']]
    elif data['chem_comp']['pdbx_model_coordinates_missing_flag'][0]=='N':
        coordinates[:,0] = ['0.' if x=='?' else x for x in atdata['model_Cartn_x']]
        coordinates[:,1] = ['0.' if x=='?' else x for x in atdata['model_Cartn_y']]
        coordinates[:,2] = ['0.' if x=='?' else x for x in atdata['model_Cartn_z']]
    else:
        coordsLabel =  'pdbx_model_Cartn_%c_ideal'
        coordinates[:,0] = ['0.' if x=='?' else x for x in atdata[coordsLabel%'x']]
        coordinates[:,1] = ['0.' if x=='?' else x for x in atdata[coordsLabel%'y']]
        coordinates[:,2] = ['0.' if x=='?' else x for x in atdata[coordsLabel%'z']]

    atomgroup._setCoords(coordinates)
    atomgroup.setNames(atdata['atom_id'])
    atomgroup.setResnames(atdata['comp_id'])
    atomgroup.setResnums([1]*asize)
    atomgroup.setChids(['A']*asize)
    atomgroup.setSerials(atdata['pdbx_ordinal'])
    atomgroup.setElements(atdata['type_symbol'])
    atomgroup.setData('atomType', atdata['type_symbol'])
    atomgroup.setData('formalCharge', atdata['charge'])
    atomgroup.setFlags('hetatm', [True]*asize)
    atomgroup.setFlags('aromatic_', [x=='Y' for x in atdata['pdbx_aromatic_flag']])
    atomgroup.setFlags('leaving', [x=='Y' for x in atdata['pdbx_leaving_atom_flag']])
    atomgroup.setFlags('chiral', [x=='Y' for x in atdata['pdbx_stereo_config']])

    # add bonds and bondOrders
    nameToIndex = {}
    for n, name in enumerate(atdata['atom_id']):
        nameToIndex[name] = n

    bdata = data.get('chem_comp_bond', None)
    if bdata:
        bonds = []
        bondOrder = []
        bondOrderNum = {'SING':1, 'DOUB':2, 'TRIP':3}
        for a1, a2, order, arom in zip(bdata['atom_id_1'], bdata['atom_id_2'],
                                     bdata['value_order'], bdata['pdbx_aromatic_flag']):
            bonds.append( (nameToIndex[a1], nameToIndex[a2]) )
            if arom=='Y':
                bondOrder.append(4)
            else:
                bondOrder.append(bondOrderNum[order])

        if len(bonds):
            atomgroup.setBonds(bonds, bondOrder)

    header = {}
    header['name'] = data['chem_comp']['name'][0]
    header['comp_chem_type'] = data['chem_comp']['type'][0]
    header['comp_chem_formula'] = data['chem_comp']['formula'][0]
    header['comp_chem_weight'] = data['chem_comp']['formula_weight'][0]
    d = data.get('pdbx_chem_comp_descriptor', None)
    if d:
        for n in range(len(d['comp_id'])):
            key = '%s %s'%(d['type'][n], d['program'][n])
            value = d['descriptor'][n]
            header[key] = value
        
    return atomgroup, header

def CIFtoMolecule(name, cifdata):
    from MolKit2.molecule import Molecule
    ag, header = CIFToProdyAG(name, cifdata)
    mol = Molecule(name, ag)
    mol.cifHeader = header
    return mol

if __name__=='__main__':
    from cif import CIFparser
    parser = CIFparser()
    f = open('mse.cif')
    lines = f.readlines()
    f.close()
    from time import time
    t0 = time()
    parser.parseString(''.join(lines))
    for name, data in parser.data.items():
        mol  = CIFtoMolecule(name[5:], data)
    #if cifFile[-3:].lower() == ".gz": parser.parse(gzip.open(cifFile))
    print time()-t0

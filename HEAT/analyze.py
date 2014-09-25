"""
Script for calculating the atomization energies and comparing them with the
CCSDTQ results in the 2004 HEAT paper.

Usage:
(1) Make sure the .csv file containing the CCSDTQ results is in the ./ref
    directory.
(2) Modify the script to change the entry names, etc.
(3) Run the script from a specific test directory, e.g. ./def2-QZVP.pbe .
"""
import os
import re
import numpy  as np
import pandas as pd


KJ_MOL_IN_AU = 2625.499630
KJ_MOL_IN_KCAL_MOL = 627.5095/2625.499630

funs = ['neo','rpa','axk','sosex','pbe','pbe0','b3-lyp']

def get_energy(m, fun):
    """
    Extract the (pbe, non-scf-acpbe, non-scf-tpss) total energies for a given
    molecule m and given functional fun.
    """
    if fun == 'neo':
        out = open(m + '/ri/ridft.out', 'rU').read()
        cor = open(m + '/neo.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'RPA correlation energy\s*\:\s*(-\d+\.\d+)D([-+]\d+)',cor)
        assert match 
        assert cormatch 
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'rpa':
        out = open(m + '/ri/ridft.out', 'rU').read()
        cor = open(m + '/axk.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'RPA correlation energy\s*\:\s*(-\d+\.\d+)E([-+]\d+)',cor)
        assert match 
        assert cormatch
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'axk':
        out = open(m + '/ri/ridft.out', 'rU').read()
        cor = open(m + '/axk.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'AXK\s*correlation energy\s*\:\s*(\-\d\.\d+)E([-+]\d+)',cor)
        assert match 
        assert cormatch 
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'sosex':
        out = open(m + '/ri/ridft.out', 'rU').read()
        cor = open(m + '/axk.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'ACSOSEX\s*correlation energy\s*\:\s*([-\s]\d\.\d+)E([-+]\d+)',cor)
        assert match 
        assert cormatch 
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    else:
        out = open(m + '/dscf.' + fun, 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        assert match
        energy = float(match.group(1)) 

    return energy


def get_atoms(m):
    """
    Return the atoms in a molecule m.
    """
    coord = open(m + '/coord', 'rU')
    atoms = [line.split()[-1] for line in coord if not line.startswith('$')]

    return atoms


# CCSDTQ calc. from Tajti et al. JCP 121, 11599 (2004)
# Take the CCSDTQ nonrelativistic calc. without vibration correction as ref.
# Energies in kJ/mol
mol_data = pd.read_csv('../ref/HEAT_CCSDTQ_rel.csv')
mol_data['Ea CCSDTQ'] = mol_data['Total'] - mol_data['ZPE'] - mol_data['DBOC'] \
                      - mol_data['scalar rel'] - mol_data['SO rel']
mols = [m.lower() for m in mol_data['molecule']]

mae = {}
for fun in funs:
    ae_mols = np.zeros(len(mols))
    for i, mol in enumerate(mols):
        en_mol = get_energy(mol, fun)
        atoms_in_mol = get_atoms(mol)
        sum_en_atoms = sum([get_energy(a, fun) for a in atoms_in_mol])
        ae_mols[i] = sum_en_atoms - en_mol
    # convert from Hartree to kJ/mol
    mol_data[fun] = ae_mols * KJ_MOL_IN_AU
    mol_data[fun + ' err.'] = mol_data[fun] - mol_data['Ea CCSDTQ']
    mae[fun] = np.mean(np.abs(mol_data[fun + ' err.']))

pd.set_option('display.max_columns', 100)
print mol_data
print '\n'
for fun in funs:
    print "MAE of", fun, "=", mae[fun]*KJ_MOL_IN_KCAL_MOL, 'kcal/mol' 

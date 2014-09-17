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


KJ_MOL_IN_AU = 2625.49962

funs = ['pbe', 'pbe0', 'b-lyp', 'b3-lyp', 'acgga', 'acgga0', 'b-acgga', 'b-acgga0']
for f in ['acgga0', 'b-acgga0']:
    for exx in np.arange(0.12, 0.23, 0.01):
        suffix = "{:.3f}".format(exx).split(".")[1]
        funs.append(f + "." + suffix)

# funs = ['pbe',
#         'b-acgga0.100', 'b-acgga0.150', 'b-acgga0.200', 'b-acgga0.250',
#         'acgga0.100', 'acgga0.150', 'acgga0.200', 'acgga0.250']

def get_energy(m, fun):
    """
    Extract the (pbe, non-scf-acpbe, non-scf-tpss) total energies for a given
    molecule m and given functional fun.
    """
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
    print "fun =", fun

pd.set_option('display.max_columns', 100)
print mol_data
print '\n'
for fun in funs:
    print "MAE of", fun, "=", mae[fun], 'kJ/mol'





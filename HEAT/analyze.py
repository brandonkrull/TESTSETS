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


def get_energies(m):
    """
    Extract the (pbe, non-scf-acpbe, non-scf-tpss) total energies for a given
    molecule m.
    """
    energies = np.zeros(3)
    for i, calc in enumerate(['pbe', 'acpbe', 'tpss']):
        out = open(m + '/dscf.' + calc, 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        assert match
        energies[i] = float(match.group(1))

    return energies


def get_atoms(m):
    """
    Return the atoms in a molecule m.
    """
    coord = open(m + '/coord', 'rU')
    atoms = [line.split()[-1] for line in coord if not line.startswith('$')]

    return atoms


atoms = ['h', 'c', 'n', 'o', 'f']

en_atoms = {}
for atom in atoms:
    en_atoms[atom] = get_energies(atom)

# CCSDTQ SC+SO relativistic calc. from Tajti et al. JCP 121, 11599 (2004)
# Energies in kJ/mol
mol_data = pd.read_csv('../ref/HEAT_CCSDTQ_rel.csv')
mols = [m.lower() for m in mol_data['molecule']]

ae_mols = np.zeros([len(mols), 3])
for i, mol in enumerate(mols):
    en_mol = get_energies(mol)
    atoms_in_mol = get_atoms(mol)
    sum_en_atoms = np.sum([en_atoms[a] for a in atoms_in_mol], axis=0)
    ae_mol = sum_en_atoms - en_mol
    ae_mols[i] = ae_mol

# convert from Hartree to kJ/mol
ae_mols = ae_mols * 2625.49962
mol_data['Ea CCSTDQ'] = mol_data['Total'] - mol_data['ZPE'] - mol_data['DBOC'] \
                      - mol_data['scalar rel'] - mol_data['SO rel']
mol_data['Ea PBE'] = ae_mols[:, 0]
mol_data['Ea acPBE@PBE'] = ae_mols[:, 1]
mol_data['Ea TPSS@PBE'] = ae_mols[:, 2]
mol_data['err. PBE'] = mol_data['Ea PBE'] - mol_data['Ea CCSTDQ']
mol_data['err. acPBE@PBE'] = mol_data['Ea acPBE@PBE'] - mol_data['Ea CCSTDQ']
mol_data['err. TPSS@PBE'] = mol_data['Ea TPSS@PBE'] - mol_data['Ea CCSTDQ']

print mol_data
print '\n'
print 'MAE of PBE', np.mean(np.abs(mol_data['err. PBE'])), ' kJ/mol'
print 'MAE of acPBE@PBE', np.mean(np.abs(mol_data['err. acPBE@PBE'])), ' kJ/mol'
print 'MAE of TPSS@PBE', np.mean(np.abs(mol_data['err. TPSS@PBE'])), ' kJ/mol'




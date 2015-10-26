"""
Script for calculating the electron affinity and comparing them with the
CCSDTQ results in the 2004 HEAT paper.

Usage:
(1) Make sure the .csv file containing the CCSDTQ results is in the ./ref
    directory.
(2) Modify the script to change the entry names, etc.
(3) Run the script from a specific test directory, e.g. ./def2-QZVP.pbe .
"""
import os
import re
import glob 
import numpy  as np
import pandas as pd 

KCAL_MOL_IN_AU = 627.5095

fexst = os.getcwd()+'/c/dscf.'
funs =  [f.replace(fexst,'') for f in glob.glob(fexst+'*')]

def get_energy(m, fun):
    """
    Extract the (pbe, non-scf-acpbe, non-scf-tpss) total energies for a given
    molecule m and given functional fun.
    """
    
    if fun == 'neo':
        out = open(m + '/ridft.hf', 'rU').read()
        cor = open(m + '/neo.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'RPA correlation energy\s*\:\s*(-\d+\.\d+)D([-+]\d+)',cor)
        print 'mol=',m
        assert match 
        assert cormatch 
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'rpa':
        out = open(m + '/ridft.hf', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        assert match 
        cor = open(m + '/rpacf.out', 'rU').read()
        cormatch = re.search(r'RPA correlation energy\s*\:\s*(-\d+\.\d+)D([-+]\d+)',cor)
        print 'mol=',m
        assert cormatch
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'axk':
        out = open(m + '/ridft.out', 'rU').read()
        cor = open(m + '/axk.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'AXK\s*correlation energy\s*\:\s*(\-\d\.\d+)D([-+]\d+)',cor)
        assert match 
        assert cormatch 
        energy = float(match.group(1))
        energy = energy+float(cormatch.group(1))*10**(float(cormatch.group(2)))

    elif fun == 'sosex':
        out = open(m + '/ridft.out', 'rU').read()
        cor = open(m + '/axk.out', 'rU').read()
        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
        cormatch = re.search(r'ACSOSEX\s*correlation energy\s*\:\s*([-\s]\d\.\d+)D([-+]\d+)',cor)
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

#reference data is in kcal/mol
mol_data = pd.read_csv('../ref/ref.csv',comment='#',skip_blank_lines=True)
mols = [m.lower() for m in mol_data['molecule']]

mse = {}
mae = {}
max_p = {}
max_m = {}
for fun in funs:
    ea_mols = np.zeros(len(mols))   # zero eas
    for i, mol in enumerate(mols): 
        en_mol = get_energy(mol, fun)   # get neutral total egy
        en_mol_neg = get_energy(mol+'-', fun) # get anion total egy
        ea_mols[i] = (en_mol - en_mol_neg)* KCAL_MOL_IN_AU # difference in kcal/mol
    mol_data[fun] = ea_mols 
    mol_data[fun + ' err.'] = (mol_data[fun] - mol_data['total'])
    mse[fun] = np.mean(mol_data[fun + ' err.'])
    mae[fun] = np.mean(np.abs(mol_data[fun + ' err.']))
    max_p[fun] = np.max(mol_data[fun + ' err.'])
    max_m[fun] = np.min(mol_data[fun + ' err.'])

pd.set_option('display.max_columns', 1500)
print mol_data
print '\n'
for fun in funs:
    print "MAE of", fun, "=", mae[fun], 'kcal/mol' 
    print "MSE of", fun, "=", mse[fun], 'kcal/mol' 
    print "MAX- of", fun, "=", max_m[fun], 'kcal/mol' 
    print "MAX+ of", fun, "=", max_p[fun], 'kcal/mol\n' 

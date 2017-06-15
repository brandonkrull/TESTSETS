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
import glob 
import numpy  as np
import pandas as pd 

KCAL_MOL_IN_AU = 627.5095
KJ_MOL_IN_AU = 2625.499630
KJ_MOL_IN_KCAL_MOL = 627.5095/2625.499630


fexst = os.getcwd()+'/h/dscf.'
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

#    elif fun.find('ksps') or fun.find('lmf') or fun.find('prop'): 
#        out = open(m + '/dscf.' + fun, 'rU').read()
#        match = re.search(r'total energy\s+=\s*(-\d+\.\d+)', out)
#        assert match
#        energy = float(match.group(1))
      
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


def get_method():
  wd = os.getcwd()
  test = wd.split('/')
  num = len(test)
  test = test[num-2]
  return test

#reference data is in kcal/mol
mol_data = pd.read_csv('../../allref/ref.'+get_method()+'.csv',comment='#')
mols = [m.lower() for m in mol_data['molecule']]

mse = {}
mae = {}
max_p = {}
max_m = {}
for fun in funs:
    prop_mols = np.zeros(len(mols))   # zero ips
    for i, mol in enumerate(mols): 
        en_mol = get_energy(mol, fun)   # get neutral total egy
        get_method()
        compute_prop(test)
        prop_mols[i] = prop
    mol_data[fun] = prop_mols 
    mol_data[fun + ' err.'] = (mol_data[fun] - mol_data['total'])
    mol_data[fun + ' err.'] = (mol_data[fun] - mol_data[fun+'ref'])
    mse[fun] = np.mean(mol_data[fun + ' err.'])
    mae[fun] = np.mean(np.abs(mol_data[fun + ' err.']))
    max_p[fun] = np.max(mol_data[fun + ' err.'])
    max_m[fun] = np.min(mol_data[fun + ' err.'])

def compute_prop(test):

    if test == "g21ea":
        en_mol = get_energy(mol, fun)   # get neutral total egy
        en_mol_neg = get_energy(mol+'-', fun) # get anion total egy
        prop = (en_mol - en_mol_neg)* KCAL_MOL_IN_AU # difference in kcal/mol
    elif test == "g21ip":
        en_mol = get_energy(mol, fun)   # get neutral total egy
        if mol == 'h':
            ip_mols[i] = -en_mol * KCAL_MOL_IN_AU # if hydrogen, there is no cation energy
        else: 
            en_mol_plus = get_energy(mol+'+', fun) # get cation total egy
    elif test == "bh76":
        en_sys1 = mol_data.stoich1[i] * get_energy(mol_data.sys1[i], fun)  
        en_sys2 = mol_data.stoich2[i] * get_energy(mol_data.sys2[i], fun) 
        if not pd.isnull(mol_data.sys3[i]):
           en_sys3 = mol_data.stoich3[i] * get_energy(mol_data.sys3[i], fun)
        else:
           en_sys3 = 0.
        prop = (en_sys1 + en_sys2 + en_sys3) * KCAL_MOL_IN_AU # convert au to kcal/mol
    elif test == "bh76rc":
        en_sys1 = mol_data.stoich1[i] * get_energy(mol_data.sys1[i], fun)  
        en_sys2 = mol_data.stoich2[i] * get_energy(mol_data.sys2[i], fun) 
        if not pd.isnull(mol_data.sys3[i]):
            en_sys3 = mol_data.stoich3[i] * get_energy(mol_data.sys3[i], fun)
        else:
            en_sys3 = 0.
        if not pd.isnull(mol_data.sys4[i]):
            en_sys4 = mol_data.stoich4[i] * get_energy(mol_data.sys4[i], fun)
        else:
            en_sys4 = 0.
        prop = (en_sys1 + en_sys2 + en_sys3 + en_sys4) * KCAL_MOL_IN_AU
    elif test == "AE6":
        en_mol = get_energy(mol, fun)
        atoms_in_mol = get_atoms(mol)
        sum_en_atoms = sum([get_energy(a, fun) for a in atoms_in_mol])
        prop = (sum_en_atoms - en_mol)* KJ_MOL_IN_AU
    elif test == "HEAT":
        en_mol = get_energy(mol, fun)
        atoms_in_mol = get_atoms(mol)
        sum_en_atoms = sum([get_energy(a, fun) for a in atoms_in_mol])
        ae_mols[i] = (sum_en_atoms - en_mol)* KJ_MOL_IN_AU*KJ_MOL_IN_KCAL_MOL

    return prop

pd.set_option('display.max_columns', 1500)
print mol_data
print '\n'
for fun in funs:
    print "MAE of", fun, "=", mae[fun], 'kcal/mol' 
    print "MSE of", fun, "=", mse[fun], 'kcal/mol' 
    print "MAX- of", fun, "=", max_m[fun], 'kcal/mol' 
    print "MAX+ of", fun, "=", max_p[fun], 'kcal/mol\n' 

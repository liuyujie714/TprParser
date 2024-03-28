""" @brief A test script for common tpr file
"""

from TprParser.TprReader import TprReader, SimSettings
from glob import glob
import sys
import numpy as np

tprlist = {
    '1EBZ.tpr' : [3218], 
    '2020.4_gra.tpr' : [4536], 
    '2022.tpr' : [165706], 
    '2023demo.tpr' : [165766],
    '2lyz_gmx_2021.tpr' : [2263], 
    '2lyz_gmx_2021_double.tpr' : [2263], 
    '2lyz_gmx_4.0.tpr' : [2263], 
    'ab42_gmx_4.6.1.tpr' : [44052], 
    'annealing.tpr' : [347443], 
    'benchMEM.tpr' : [81743], 
    'double_2023.tpr' : [16844], 
    'em.tpr' : [252], 
    'Inter-2019.6.tpr' : [157488], 
    'inter-md.tpr' : [13749], 
    'large_2021_aa_posres.tpr' : [34466], 
    'md.tpr' : [2520], 
    'md2024.tpr' : [58385], 
    'md_cg.tpr' : [8], 
    'pull.tpr' : [94560], 
    'semiP.tpr' : [4608],
    'CO2_LineAngle.tpr' : [3000]
}


def test_get_xvf(handle, ftype):
    try:
        ret = handle.get_xvf(ftype)
    except:
        sys.exit(f'Can not execute test_get_xvf("{ftype}") function')

def test_get_bonded(handle, ftype):
    try:
        ret = handle.get_bonded(ftype)
    except:
        sys.exit(f'Can not execute get_bonded("{ftype}") function')

def test_get_mq(handle, ftype):
    try:
        ret = handle.get_mq(ftype)
    except:
        sys.exit(f'Can not execute get_mq("{ftype}") function')

def test_get_name(handle, ftype):
    try:
        ret = handle.get_name(ftype)
    except:
        sys.exit(f'Can not execute get_name("{ftype}") function')

def test_tot_atoms(handle, natoms, fname):
    assert natoms == len(handle.get_name('res')), f"The number of atoms is wrong in file {fname}"

def do_test():
    for index, name in enumerate(tprlist.keys()):
        print(f'do test {index+1}', flush=True)
        fname = 'test/' + name
        try:
            reader = TprReader(fname)
        except:
            sys.exit(f'Can not init tpr handle for file: {fname}')

        # total atoms
        test_tot_atoms(reader, tprlist[name][0], fname)

        test_get_xvf(reader, 'x')
        test_get_xvf(reader, 'v')
        
        test_get_bonded(reader, 'bonds')
        # pure water use settle, no angle
        if 'semiP.tpr' not in fname:
            test_get_bonded(reader, 'angles')
        # these tpr has not dihedrals
        if ('md.tpr' not in fname) and ('md_cg.tpr' not in fname) and \
            ('semiP.tpr' not in fname):
            test_get_bonded(reader, 'dihedrals')
            test_get_bonded(reader, 'impropers')

        test_get_mq(reader, 'm')
        test_get_mq(reader, 'q')

        test_get_name(reader, 'res')
        test_get_name(reader, 'atom')
        
        # need delete obj
        del reader

def do_test2():
    fout = 'output.tpr'
    for index, name in enumerate(tprlist.keys()):
        print(f'do test {index+1}', flush=True)
        fname = 'test/' + name

        # get precision of tpr
        reader = TprReader(fname)
        prec = reader.get_prec()
        del reader
        
        with SimSettings(fname, fout) as writer:
            # change 
            writer.set_dt(0.002)
            writer.set_nsteps(100)
            # unsupport set_mdp_integer for gmx < 4.6
            if '4.0' not in name:
                writer.set_mdp_integer('nstxout', 100)
                writer.set_mdp_integer('nstenergy', 100)
                writer.set_mdp_integer('nsttcouple', 1)
                writer.set_mdp_integer('nstpcouple', 1)
                writer.set_mdp_integer('nstxout_compressed', 1032)

            if prec==4:
                writer.set_pressure('CRescale', 'Isotropic', 3.0, 
                                    [
                                        100,0, 0,
                                        0, 100,0,
                                        0, 0, 100
                                    ],
                                    [
                                        1,0,0,
                                        0,1,0,
                                        0,0,1,
                                    ]
                                    )
                newX = 152*np.ones(shape=(tprlist[name][0], 3))
                newV = 110*np.ones(shape=(tprlist[name][0], 3))
                writer.set_xvf('x', newX)
                writer.set_xvf('v', newV)
            
        # assert modify parameters
        reader = TprReader(fout)
        x = reader.get_xvf('x')
        v = reader.get_xvf('v')

        if prec==4:
            assert np.all(newX==x)
            assert np.all(newV==v)  
            
        if '4.0' not in name:
            assert reader.get_mdp_integer('nstxout') == 100
            assert reader.get_mdp_integer('nstenergy') == 100
            assert reader.get_mdp_integer('nsttcouple') == 1
            assert reader.get_mdp_integer('nstpcouple') == 1
            assert reader.get_mdp_integer('nstxout_compressed') == 1032
        
        del reader

if __name__ == '__main__':
    do_test()
    print('<'*10+'Passed All TprParser Tests'+'>'*10, flush=True)

    do_test2()
    print('<'*10+'Passed All SimSettings Tests'+'>'*10, flush=True)

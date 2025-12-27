""" @brief A test script for common tpr file
"""

from TprParser.TprReader import TprReader, SimSettings
import sys
import numpy as np

# All test tpr file: [natoms, prec]
tprlist = {
    # No dihedrals
    'md.tpr' :              [2520, 4], 
    'md_cg.tpr' :           [8, 4], 
    'semiP_water.tpr' :     [4608, 4],
    'CO2_LineAngle.tpr' :   [3000, 4],
    # gmx2025-beta
    'npt2025-beta_water.tpr':[2652,4],
    # gmx 2026-beta
    'md2026-beta_water.tpr' :     [2652,4],
    # enforced rotation
    'enforced_rotation_water.tpr' : [7306, 4],
    # test [ exclusion ] 
    'one_water_tip3p_excls.tpr' : [3, 4], 
    'one_water_tip4p_excls.tpr' : [4, 4], 
    'two_water_tip4p_excls.tpr' : [8, 4], 
    'CO2_vsites.tpr'            : [5000, 4], 

    '1EBZ.tpr' :                [3218, 4], 
    '2020.4_gra.tpr' :          [4536, 4], 
    '2022.tpr' :                [165706, 4], 
    '2023demo.tpr' :            [165766, 4],
    "2lyz_gmx_4.0.2.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.3.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.4.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.5.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.6.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.7.tpr"      : [2263, 4],
    "2lyz_gmx_4.0.tpr"        : [2263, 4],
    "2lyz_gmx_4.5.1.tpr"      : [2263, 4],
    "2lyz_gmx_4.5.2.tpr"      : [2263, 4],
    "2lyz_gmx_4.5.3.tpr"      : [2263, 4],
    "2lyz_gmx_4.5.4.tpr"      : [2263, 4],
    "2lyz_gmx_4.5.5.tpr"      : [2263, 4],
    "2lyz_gmx_4.5.tpr"        : [2263, 4],
    "2lyz_gmx_5.0.2.tpr"      : [2263, 4],
    "2lyz_gmx_5.0.4.tpr"      : [2263, 4],
    "2lyz_gmx_5.0.5.tpr"      : [2263, 4],
    "2lyz_gmx_5.1.tpr"        : [2263, 4],
    "2lyz_gmx_2016.tpr"       : [2263, 4],
    "2lyz_gmx_2018.tpr"       : [2263, 4],
    "2lyz_gmx_2019-beta3.tpr" : [2263, 4],
    "2lyz_gmx_2020.tpr"       : [2263, 4],
    "2lyz_gmx_2020_double.tpr": [2263, 8],
    # "2lyz_gmx_2020-beta2.tpr" : [2263, 4], # unsupport this beta version
    "2lyz_gmx_2021.tpr"       : [2263, 4],
    "2lyz_gmx_2021_double.tpr": [2263, 8],
    "2lyz_gmx_2022-rc1.tpr"   : [2263, 4],
    "2lyz_gmx_2023.tpr"       : [2263, 4],
    'ab42_gmx_4.6.tpr'   :      [44052, 4], 
    'ab42_gmx_4.6.1.tpr' :      [44052, 4], 
    'annealing.tpr' :           [347443, 4], 
    'benchMEM.tpr' :            [81743, 4], 
    'double_2023_cg.tpr' :      [16844, 8], 
    'em.tpr' :                  [252, 4], 
    'Inter-2019.6.tpr' :        [157488, 4], 
    'inter-md.tpr' :            [13749, 4], 
    'large_2021_aa_posres.tpr' : [34466, 4], 
    'md2024.tpr' :          [58385, 4], 
    'pull.tpr' :            [94560, 4],
    'nobox.tpr' :           [13, 4],
    'cg_big.tpr':           [290482, 4],
    # electric-field
    'elec5.1.2.tpr':        [45, 4], # along x
    'elec2019.tpr':         [45, 4], # along z
    'elecxyz.tpr':          [45, 4], # along xyz
    'elecxyz_2024.tpr':     [45, 4], # along xyz
    # FEP
    'benchBFC_FEP.tpr' :    [43952,4],
    # No lj parameters
    'extra-interactions-2018.tpr' : [17, 4],
}
NoDihedrals = [k for k in list(tprlist.keys())[0:11]]


rand_int = lambda : np.random.randint(0, 100000)
# for test mdp set and get
mdp_integer_data = {
    'nstlog' : rand_int(),
    'nstxout': rand_int(),
    'nstvout': rand_int(),
    'nstfout': rand_int(),
    'nstcalcenergy': rand_int(),
    'nstenergy': rand_int(),
    'nsttcouple': rand_int(),
    'nstpcouple': rand_int(),
    'nstxout_compressed': rand_int(),
    'nstlist': rand_int(),
    'nstcomm': rand_int(),
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

def test_get_ivector(handle, ftype, fname:str=""):
    try:
        ret = handle.get_ivector(ftype)
    except:
        sys.exit(f'Can not execute get_ivector("{ftype}") function for {fname}')

def test_tot_atoms(handle, natoms, fname):
    assert natoms == len(handle.get_name('res')), f"The number of atoms is wrong in file {fname}"

def test_precision(handle:TprReader, fname, prec=4):
    assert prec == handle.get_prec(), f"The precision is not euqal to {prec} for file: {fname}"

def test_exclusions(handle:TprReader, fname):
    excls_map = {
        'tip3p' : [[0, 1, 2], [0, 1, 2], [0, 1, 2]],
        'tip4p' : [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3]],
        'two_tip4p' : [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3], [4, 5, 6, 7], [4, 5, 6, 7], [4, 5, 6, 7], [4, 5, 6, 7]],
    }
    if 'tip3p' in fname:
        assert excls_map['tip3p'] == handle.get_exclusions(), "The exclusions is not euqal for file: {fname}"
    elif 'tip4p' in fname:
        if 'one' in fname:    
            assert excls_map['tip4p'] == handle.get_exclusions(), "The exclusions is not euqal for file: {fname}"
        elif 'two' in fname:
            assert excls_map['two_tip4p'] == handle.get_exclusions(), "The exclusions is not euqal for file: {fname}"

def test_get_vsites(handle:TprReader, fname):
    try:
        ret = handle.get_vsites()
    except:
        sys.exit(f'Can not execute get_vsites() function for {fname}')

def test_make_top_from_tpr(tpr, top):
    from TprParser.TprMakeTop import make_top_from_tpr
    try:
        make_top_from_tpr(tpr, top)
    except:
        sys.exit(f'Can not execute test_make_top_from_tpr for file: {tpr}')

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

        # test tpr precision
        if 'double' in fname:
            test_precision(reader, fname, 8)
        else:
            test_precision(reader, fname, 4)

        # test coords/velocity
        test_get_xvf(reader, 'x')
        test_get_xvf(reader, 'v')

        # test electric field to get
        if 'elec' in fname:
            test_get_xvf(reader, 'ef')

        # test bonded
        test_get_bonded(reader, 'bonds')
        # pure water use settle, no angle
        # CO2 use constraints, no angle
        if 'water' not in fname and 'CO2' not in fname:
            test_get_bonded(reader, 'angles')
        # these tpr has not dihedrals
        if name not in NoDihedrals:
            test_get_bonded(reader, 'dihedrals')
            test_get_bonded(reader, 'impropers')

        # test atomic mass and charge
        test_get_mq(reader, 'm')
        test_get_mq(reader, 'q')

        # test resname, atomname, atomtype
        test_get_name(reader, 'res')
        test_get_name(reader, 'atom')
        test_get_name(reader, 'type')
        
        # test resid, atomtypenumber(filever>128 not do_atomtypes), atomic number
        test_get_ivector(reader, "resid")
#        test_get_ivector(reader, "atnum")
        
        # elec/cg/some low version tpr all atom number == -1 or 0
        if ('elec' not in fname) and ('cg' not in fname) and \
            ('benchMEM' not in fname) and ('nobox' not in fname) and \
            ('extra' not in fname):
            test_get_ivector(reader, "atomicnum", fname)

        # test exclusions
        if 'excls' in fname:
            test_exclusions(reader, fname)
        
        # test virtual sites
        test_get_vsites(reader, fname)

        # need delete obj
        del reader

        # test write gromacs top from tpr
        # not include atomnumber of atomtype ['double_2023.tpr', 'double_2023.tpr', 'md2024.tpr']:
        # Not LJ parameters ['extra-interactions-2018.tpr']
        test_make_top_from_tpr(fname, 'md.top')

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
                for key, val in mdp_integer_data.items():
                    writer.set_mdp_integer(key, val)

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
                # add deform
                writer.set_pressure('Berendsen', 'Anisotropic', 1.0, 
                                    [
                                        100,0, 0,
                                        0, 100,0,
                                        0, 0, 100
                                    ],
                                    [
                                        1,0,0,
                                        0,1,0,
                                        0,0,1,
                                    ],
                                    [
                                        0, 0, 0,
                                        0, 0, 0,
                                        0.01, 0, 0
                                    ]
                                    )
                newX = np.random.uniform(-999, 999, tprlist[name][0]*3).reshape(-1, 3)
                newV = np.random.uniform(-999, 999, tprlist[name][0]*3).reshape(-1, 3)
                writer.set_xvf('x', newX)
                writer.set_xvf('v', newV)

                # test modify electric field
                if 'elecxyz' in fname:
                    # E0, omega, t0, sigma for each dim
                    ef = [
                        10, 0, 0,   0,
                        10, 0, 1.5, 0,
                        10, 0, 0,   2.0
                    ]
                    ef = np.array(ef, dtype=np.float32)
                    writer.set_xvf('ef', ef)
            
        # assert modify parameters
        reader = TprReader(fout)
        x = reader.get_xvf('x')
        v = reader.get_xvf('v')

        if prec==4:
            assert np.allclose(newX, x, atol=1E-3)
            assert np.allclose(newV, v, atol=1E-3)
            # assert electric-field
            if 'elecxyz' in fname:
                assert np.all(ef==reader.get_xvf('ef').flatten())
            
        if '4.0' not in name:
            for key, val in mdp_integer_data.items():
                assert reader.get_mdp_integer(key) == val, f'get_mdp_integer {key} should be {val}'
        
        del reader

if __name__ == '__main__':
    do_test()
    print('<'*10+'Passed All TprParser Tests'+'>'*10, flush=True)

    do_test2()
    print('<'*10+'Passed All SimSettings Tests'+'>'*10, flush=True)

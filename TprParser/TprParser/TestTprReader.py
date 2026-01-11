""" @brief A test script for common tpr file
"""

from TprParser.TprReader import TprReader, SimSettings
import MDAnalysis as mda
from copy import deepcopy
import numpy as np
import sys

import warnings
warnings.filterwarnings("ignore")

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
    'gmx_3.3.1_water.tpr'       : [864, 4], 
    'CO2_vsites.tpr'            : [5000, 4], 

    'lzm_gmx_3.2.tpr' :         [23207, 4], 
    'lzm_gmx_3.3.3.tpr' :       [23207, 4], 
    'dppc_sol_gmx_3.2.tpr' :    [121856, 4], 
    'dppc_sol_gmx_3.3.3.tpr' :  [121856, 4], 
    '1EBZ_4.5.5.tpr' :          [3218, 4], 
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
    # BHAM & No lj parameters, for test all extra interactions
    'mda_extra-interactions-2018.tpr' : [17, 4],
    # nyself, from dummy_2025&6.top
    'extra-interactions-2021.tpr' : [18, 4],
    'extra-interactions-2022.tpr' : [18, 4],
    'extra-interactions-2023.tpr' : [18, 4],
    'extra-interactions-2024.tpr' : [18, 4],
    'extra-interactions-2025.tpr' : [18, 4],
    'extra-interactions-2026.tpr' : [18, 4],
    # 计算电生理学
    'swapcoords_gmx_2018.tpr'     : [32681, 4],
    'swapcoords_gmx_2019.tpr'     : [32681, 4],
    'swapcoords_gmx_2019.6_double.tpr' : [32681, 8],
    'swapcoords_gmx_2021.tpr'     : [32681, 4],
    'swapcoords_gmx_2024.tpr'     : [32681, 4],
    'swapcoords_gmx_2025.tpr'     : [32681, 4],
    'swapcoords_gmx_2026-rc.tpr'  : [32681, 4],
}
NoDihedrals = [k for k in list(tprlist.keys())[0:12]]


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

def test_set_mq(handle, ftype:str):
    ret = handle.set_mq(ftype)
    if ret is not True:
        sys.exit(f'Can not execute set_mq("{ftype}") function')

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

def read_top_ff(fname:str):
    labels = [
        'virtual_sites1',
        'virtual_sites2',
        'virtual_sites3',
        'virtual_sites4',
        'virtual_sitesn',
        'pairs',
        'bonds',
        'constraints',
        'angles',
        'dihedrals',
    ]
    ff_items = {key: [] for key in labels}
    with open(fname, 'r') as f:
        lines = f.readlines()
        for label in labels:
            found = False
            for line in lines:
                if line.startswith(';'):
                    continue
                if line.startswith('[') and label in line:
                    found = True
                    continue
                if found and len(line.strip()) < 2:
                    found = False
                if found:
                    data = list(map(float, line.split(';')[0].strip().split()))
                    ff_items[label].append(data)

    return ff_items
def test_mda_top(handle:TprReader, fname:str):
    u = mda.Universe(fname, topology_format='itp')
    # below items always exists
    # compare bonds
    _bonds = handle.get_bonded('bonds') # 1-based
    bonds = np.array([[t[0], t[1]] for t in _bonds], dtype=int)
    mda_bonds = u.atoms.bonds.to_indices()+1
    assert np.array_equal(bonds, mda_bonds)

    # compare angles
    _angles = handle.get_bonded('angles')
    # mda can not read func type=9 (linear angle)
    angles = np.array([[t[0], t[1], t[2]] for t in _angles if t[3]!=9], dtype=int)
    mda_angles = u.atoms.angles.to_indices()+1
    assert np.array_equal(angles, mda_angles)

    # compare propers dihedrals
    _dihedrals = handle.get_bonded('dihedrals')
    dihedrals = np.array([[t[0], t[1], t[2], t[3]] for t in _dihedrals], dtype=int)
    mda_dihedrals = u.atoms.dihedrals.to_indices()+1
    assert np.array_equal(dihedrals, mda_dihedrals)

    # compare improper dihedrals
    _imps = handle.get_bonded('impropers')
    impropers = np.array([[t[0], t[1], t[2], t[3]] for t in _imps], dtype=int)
    # sort impropers for compare, because tprparser can not sort for improper, keep it original order
    adjusted = np.array([(d,c,b,a) if a>d else (a,b,c,d) for a,b,c,d in impropers])
    impropers = adjusted[np.lexsort(adjusted.T[::-1])]
    mda_impropers = u.atoms.impropers.to_indices()+1
    assert np.array_equal(impropers, mda_impropers)

    # comare atoms
    atomnames, resnames, types = handle.get_name('atom'), handle.get_name('res'), handle.get_name('type')
    assert np.array_equal(atomnames, u.atoms.names), "The atomnames is not euqal for file: {fname}"
    assert np.array_equal(resnames, u.atoms.resnames), "The resnames is not euqal for file: {fname}"
    assert np.array_equal(types, u.atoms.types), "The types is not euqal for file: {fname}"
    charges, masses = handle.get_mq('q'), handle.get_mq('m')
    assert np.array_equal(charges, u.atoms.charges), "The charges is not euqal for file: {fname}"
    assert np.array_equal(masses, u.atoms.masses), "The masses is not euqal for file: {fname}"

    # test ff
    ff_items = read_top_ff(fname)
    vsites = handle.get_vsites()
    for vsite in vsites:
        name = vsite[0]
        if name=="virtual_sitesn":
            continue
        tpr_vsite = vsite[1:]
        top_vsites = ff_items[name]
        found = False
        for top_vsite in top_vsites:
            size = len(top_vsite)
            if len(tpr_vsite)>=size and np.allclose(tpr_vsite[:size], top_vsite):
                found = True
        if not found:
            sys.exit(f'Can not match virtual site {name} in file: {fname}')
    # bonds ff
    for tpr_bond in _bonds:
        top_bonds = deepcopy(ff_items['bonds'])
        # bonds from constraints
        top_bonds.extend(ff_items['constraints'])
        found = False
        for top_bond in top_bonds:
            size = min(len(top_bond), len(tpr_bond))
            if top_bond[0]>top_bond[1]:
                top_bond[0], top_bond[1] = top_bond[1], top_bond[0]
            if np.allclose(tpr_bond[:size], top_bond[:size]):
                found = True
        if not found:
            sys.exit(f'Can not match bond {tpr_bond} in file: {fname}')
    # angles ff
    for tpr_angle in _angles:
        top_angles = deepcopy(ff_items['angles'])
        found = False
        for top_angle in top_angles:
            size = min(len(top_angle), len(tpr_angle))
            # exchange atom order
            if top_angle[0]>top_angle[2]:
                top_angle[0], top_angle[2] = top_angle[2], top_angle[0]
            if np.allclose(tpr_angle[:size], top_angle[:size]):
                found = True
        if not found:
            sys.exit(f'Can not match angle {tpr_angle} in file: {fname}')
    # dihedrals ff
    for tpr_dihedral in _dihedrals:
        top_dihedrals = deepcopy(ff_items['dihedrals'])
        found = False
        for top_dihedral in top_dihedrals:
            size = min(len(top_dihedral), len(tpr_dihedral))
            # exchange atom order
            if top_dihedral[0]>top_dihedral[3]:
                top_dihedral[0], top_dihedral[3] = top_dihedral[3], top_dihedral[0]
                top_dihedral[1], top_dihedral[2] = top_dihedral[2], top_dihedral[1]
            # type=1 and 9 is same
            if top_dihedral[4] == 1:
                top_dihedral[4] = 9
            if np.allclose(tpr_dihedral[:size], top_dihedral[:size]):
                found = True
        if not found:
            sys.exit(f'Can not match dihedral {tpr_dihedral} in file: {fname}')
    # impropers ff
    for tpr_improper in _imps:
        top_impropers = deepcopy(ff_items['dihedrals'])
        found = False
        for top_improper in top_impropers:
            size = min(len(top_improper), len(tpr_improper))
            # keep original atom order
            if np.allclose(tpr_improper[:size], top_improper[:size]):
                found = True
        if not found:
            sys.exit(f'Can not match improper {tpr_improper} in file: {fname}')

def test_make_top_from_tpr(tpr, top):
    from TprParser.TprMakeTop import make_top_from_tpr
    try:
        make_top_from_tpr(tpr, top)
    except:
        sys.exit(f'Can not execute test_make_top_from_tpr for file: {tpr}')

def do_reader():
    for index, name in enumerate(tprlist.keys()):
        print(f'do_reader {name}', flush=True)
        fname = 'test/' + name
        try:
            reader = TprReader(fname)
        except:
            sys.exit(f'Can not init tpr handle for file: {fname}')

        # total atoms
        test_tot_atoms(reader, tprlist[name][0], fname)

        # test tpr precision
        test_precision(reader, fname, tprlist[name][1])

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
            ('extra' not in fname) and ('gmx_3' not in fname):
            test_get_ivector(reader, "atomicnum", fname)

        # test exclusions
        if 'excls' in fname:
            test_exclusions(reader, fname)
        
        # test virtual sites
        test_get_vsites(reader, fname)

        # compare tpr and top
        if name.startswith('extra-interactions'):
            test_mda_top(reader, 'test/dummy_2025&6.top')

        # need delete obj
        del reader

        # test write gromacs top from tpr
        # not include atomnumber of atomtype ['double_2023.tpr', 'double_2023.tpr', 'md2024.tpr']:
        # Not LJ parameters ['extra-interactions-2018.tpr']
        test_make_top_from_tpr(fname, 'md.top')

def do_writer():
    fout = 'output.tpr'
    Verlet, Group = range(0, 2)
    for index, name in enumerate(tprlist.keys()):
        # skip very old tpr
        if 'gmx_3' in name:
            continue
        print(f'do_writer {name}', flush=True)
        fname = 'test/' + name

        expected_q = -666
        expected_m = 9999
        with SimSettings(fname, fout) as writer:
            # change 
            writer.set_dt(0.002)
            writer.set_nsteps(100)
            # change atomic charge & mass
            writer.set_mq('q', [expected_q]*tprlist[name][0])
            writer.set_mq('m', [expected_m]*tprlist[name][0])

            # unsupport set_mdp_integer for gmx < 4.6
            if '4.0' not in name:
                for key, val in mdp_integer_data.items():
                    writer.set_mdp_integer(key, val)

                if '4.5' not in name:
                    #writer.set_mdp_integer('cutoff_scheme', 0) # verlet
                    writer.set_mdp_integer('cutoff_scheme', Group) # group

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
        q = reader.get_mq('q')
        m = reader.get_mq('m')

        assert np.allclose(newX, x, atol=1E-3)
        assert np.allclose(newV, v, atol=1E-3)
        # assert electric-field
        if 'elecxyz' in fname:
            assert np.all(ef==reader.get_xvf('ef').flatten())
        # test charge & mass
        assert np.allclose(q, expected_q, atol=1E-3), f'q should be {expected_q}'
        assert np.allclose(m, expected_m, atol=1E-3), f'm should be {expected_m}'
            
        if '4.0' not in name:
            for key, val in mdp_integer_data.items():
                assert reader.get_mdp_integer(key) == val, f'get_mdp_integer {key} should be {val}'
            if '4.5' not in name:
                assert reader.get_mdp_integer('cutoff_scheme')==Group, f'get_mdp_integer cutoff_scheme should be {Group}'
        
        del reader

if __name__ == '__main__':
    do_reader()
    print('<'*10+'Passed All TprParser Tests'+'>'*10, flush=True)

    do_writer()
    print('<'*10+'Passed All SimSettings Tests'+'>'*10, flush=True)

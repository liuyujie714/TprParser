# Description

`TprParser` is a convenient Python module for reading and setting simulation parameters of gromacs `tpr` file. It does **NOT** rely on the GROMACS library and `only` pure Python environment.

This module mainly aimed to modify **atom property** of `tpr` and create a new tpr file (named `new.tpr`) after use any one `set_` method. 

The module supports many functions to get topology properties, such as atoms coordinates, velocity and force by `module.get_xvf(...)` function. 

Many properties can be set up by this module, such as total simulation time `nsteps`, simulation integrator interval `dt`, output control parameters (`nstxout, nstvout, etc.`) and temperature/pressure coupling parameters.



# Compatibility

System Support: Linux & Windows

GROMACS tpr version should between `3.2` to `2026`, too old tpr to be read/write by this module.

| TPX format | TPX generation | Gromacs release                               | read/write |
| :--------: | :------------: | :-------------------------------------------- | :--------: |
| 31         | 4              | 3.2                                           | read-only  |
| 40         | 7              | 3.3.1, 3.3.3                                  | read-only  |
| 58         | 17             | 4.0, 4.0.2, 4.0.3, 4.0.4, 4.0.5, 4.0.6, 4.0.7 | read-only  |
| 73         | 23             | 4.5.0, 4.5.1, 4.5.2, 4.5.3, 4.5.4, 4.5.5      | read-write |
| 83         | 24             | 4.6, 4.6.1                                    | read-write |
| 100        | 26             | 5.0, 5.0.1, 5.0.2, 5.0.3, 5.0.4, 5.0.5        | read-write |
| 103        | 26             | 5.1                                           | read-write |
| 110        | 26             | 2016                                          | read-write |
| 112        | 26             | 2018                                          | read-write |
| 116        | 26             | 2019                                          | read-write |
| 119        | 27             | 2020                                          | read-write |
| 122        | 28             | 2021                                          | read-write |
| 127        | 28             | 2022                                          | read-write |
| 129        | 28             | 2023                                          | read-write |
| 133        | 28             | 2024                                          | read-write |
| 137        | 28             | 2025                                          | read-write |
| 138        | 29             | 2026                                          | read-write |


# Installation

* Requirements

  * `Python >= 3.8`

  * `Numpy`


* Install

  The module is installed by `pip` method:

  ```
  pip install TprParser -i https://pypi.org/simple
  ```

  Please **ALWAYS** install Latest version.

  Update this module if you have installed:
  ```
  pip install TprParser --upgrade -i https://pypi.org/simple
  ```



# Usage

Write your python program like this:

```python
# import this module
from TprParser.TprReader import TprReader	
```



## Get atom property


```python
# get atom coords or velocity
reader = TprReader("your.tpr")
coords = reader.get_xvf('x')
velocity = reader.get_xvf('v')

# get atom charge or mass
charge = reader.get_mq('q')
mass = reader.get_mq('m')

# get residue or atom name or atomtype of each atom
resnames  = reader.get_name('res')
atomnames = reader.get_name('atom')
atomtypes = reader.get_name('type')

# get resids
resids = reader.get_ivector('resid')
```

Note: The behavior of `get_ivector('resid')` is different when `TprParser >= 0.1.58`, the resids is consistent with output of `gmx editconf -o xxx.gro`. If `TprParser < 0.1.58`, the resids is unique for global atoms, same as `MDAnalysis`



## Get bonds/angles/dihedrals(proper and impropers) forcefield parameters

```python
# get all bond pairs (1-based index) and it's parameters
bonds = reader.get_bonded('bonds')

# get all angles pairs (1-based index) and it's parameters
angles = reader.get_bonded('angles')

# get all proper dihedrals pairs (1-based index) and it's parameters
propers = reader.get_bonded('dihedrals')

# get all improper dihedrals pairs (1-based index) and it's parameters
impropers = reader.get_bonded('impropers')

```



## Get non-bonded paramaters

```python
# get pairs (1-based index) and it's parameters
pairs = reader.get_nonbonded('pairs')

# get atomtypes lj parameters for each atoms [sigma, epsion], which crossbonding to reader.get_name('type')
pairs = reader.get_nonbonded('lj')[:, 1:]
```



## Get virtual sites parameters

NOTE: `TprParser` must be >= `0.1.60`

```python
# get virtual site type + atom index (1-based index) + it's parameters
vsites = reader.get_vsites()
```



# Get mdp parameters

```python
# integer value
nstxout = reader.get_mdp_integer('nstxout')
nsttcouple = reader.get_mdp_integer('nsttcouple')

# float value, TprParser must be >= 0.1.64
rvdw = reader.get_mdp_float('rvdw')
taup = reader.get_mdp_float('tau_p')
```



## Modify atom property

### Modify atomic coords/velocity/box
```python
newcoords = np.array([[1,2,3], [4,5,6], [...]], dtype=np.float32) # shape= N*3
# The step will create new.tpr that used newcoords
reader.set_xvf('x', newcoords)
# set velocity or box, same as above
reader.set_xvf('v', newvelocity)
```

### Modify atomic charges/masses
Note: `TprParser` must be >= `0.1.63`

Atomic charges and masses are stored on a per molecule type in `.tpr` file compared to coordinates. So you must set it for each molecule type.

For example, if your simulation system topolgy is:
```bash
[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL                1000
Protein_chain_A     1
```
And each molecule type `Protein_chain_A` has `x` atoms, `SOL` has `y` atoms, you can set new charges and masses:
```python
pro_charges = [1.0,2.0,3.0, ...]  # length= x
sol_charges = [2.0,3.0,4.0]       # length= y, such as y==3 for TIP3P
# sum the number of molecules, the order of molecule is important
# the length(newcharges)==total number of atoms in system
newcharges = pro_charges*1 + sol_charges*1000 + pro_charges*1
reader.set_mq('q', newcharges)
reader.set_mq('m', newmasses) # same as above
```
>Tips: You can use `gmx check` to compare old and new tpr differences and make sure the changes are correct:
```bash
gmx check -s1 old.tpr -s2 new.tpr
```

## Modify system pressure

Such as define a function to do this work:

```Python
def change_pressure(fname):
    reader = TprReader(fname)
    # 100 bar
    ref_p = [
        100, 0, 0,
        0, 100, 0,
        0, 0, 100
    ]
    # compressibility 4.5E-5
    compress = [
        4.5E-5, 0, 0,
        0, 4.5E-5, 0,
        0, 0, 4.5E-5
    ]
    assert len(ref_p) == 9
    assert len(compress) == 9
    # use ParrinelloRahman algorithm and Isotropic pressure coupling method
    reader.set_pressure('ParrinelloRahman', 'Isotropic', 1.0, ref_p, compress)
```

Also can change deform parameters for `TprParser` must be >= `0.1.52`
```python
# optional modify deform
deform = [
  0, 0, 0,
  0, 0, 0,
  0.01, 0, 0
]
reader.set_pressure('Berendsen', 'anisotropic', 1.0, ref_p, compress, deform)
```



## Modify system temperature

```python
# set Berendsen algorithm and tau_t=0.2, ref_t=400 K for one temperature coupling group
reader.set_temperature(etc='Berendsen', tau_t=[0.2], ref_t=[400])
```



## Modify electric field parameters

NOTE: `TprParser` must be >= `0.1.53`
```python
# The modify must be matched to old tpr electric-field dimension
# E0, omega, t0, sigma for each dimension, use gromacs unit
newEF = [
  10, 1, 0, 0.1,  # x direction
  0,  0, 0, 0,    # y direction
  0,  0, 0, 0     # z direction
]
reader.set_xvf('ef', newEF)
```



## Modify multiple parameters
Use `SimSettings` class to do this work
```python
from TprParser.TprReader import SimSettings

with SimSettings('input.tpr', 'output.tpr') as writer:
    writer.set_dt(0.001)
    writer.set_mdp_integer('nstxout', 100)
    writer.set_mdp_float('rvdw', 1.4)
    writer.set_nsteps(2000000)
```



## Make a gromacs top

Note: `TprParser` must be >= `0.1.51`. The top is not a full topology, such as `Restraint` is missing!

```python
from TprParser.TprMakeTop import make_top_from_tpr

make_top_from_tpr('md.tpr', 'out.top')
```



## Read gromacs edr file

Note: `TprParser` must be >= `0.1.57`. 

```python
from TprParser.EdrReader import EdrReader

# return a dictionary of all energies
# use gromacs unit, such as Energy -> KJ/mol, Length -> nm
energies = EdrReader('yourfile.edr').get_ene()

# available energies name
print(energies.keys())

# get vaules by available key
times = energies['Time']                # ps
temperature = energies['Temperature']   # K
```



## Other

Please see `TprReader` and `SimSettings` module annotation



# Cite

If `TprParser` is utilized in your work, please cite as follows in main text:

> Yujie Liu, TprParser, Version [xxx](), https://pypi.org/project/TprParser/



## TODO

* More parameters can be modified




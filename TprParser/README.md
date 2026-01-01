# Description

`TprParser` is a convenient Python module for reading and setting simulation parameters of gromacs `tpr` file. It does **NOT** rely on the GROMACS library and `only` requires C++ and pure Python environment.

This module mainly aimed to modify **atom property** of `tpr` and create a new tpr file (named `new.tpr`) after use any one `set_` method. 

The module supports many functions to get topology properties, such as atoms coordinates, velocity and force by `module.get_xvf(...)` function. 

Many properties can be set up by this module, such as total simulation time `nsteps`, simulation integrator interval `dt`, output control parameters (`nstxout, nstvout, etc.`) and temperature/pressure coupling parameters.



# Compatibility

GROMACS tpr version should between `4.0` to `2026`, too old tpr to be read by this module.



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
from TprParser.TprReader import TprReader	# import this module
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
# get paris (1-based index) and it's parameters
pairs = reader.get_nonbonded('pairs')

# get atomtypes lj parameters for each atoms [sigma, epsion], which crossbonding to reader.get_name('type')
pairs = reader.get_nonbonded('lj')[:, 1:]
```



## Get virtual sites parameters

NOTE: `TprParser` must be >= `0.1.60`

```python
# get virtual site type + paris (1-based index) + it's parameters
vsites = reader.get_vsites()
```



## Modify atom property

```python
newcoords = np.array([[1,2,3], [4,5,6], [...]], dtype=np.float32) # shape= N*3
# The step will create new.tpr that used newcoords
reader.set_xvf('x', newcoords)
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
    writer.set_nsteps(2000000)
```



## Make a gromacs top

Note: `TprParser` must be >= `0.1.51`. The top is not a full topology, such as `Virtual Site, Restraint` is missing!

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

* more parameters can be modified




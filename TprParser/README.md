# Description

`TprParser` is a convenient Python module for reading and setting simulation parameters in gromacs `tpr` file.

This module mainly aimed to modify **atom property** of `tpr` and create a new tpr file (named `new.tpr`) after use any one `set_` method. 

The module only supports get atoms coordinates, velocity and force if exist in `tpr` by `module.get_xvf(...)` function. 

However, many properties can be set up by module, such as total simulation time `nsteps`, simulation integrator interval `dt`, output control parameters (`nstxout, nstvout, etc.`) and temperature/pressure coupling parameters.

# Compatibility
GROMACS tpr version should between `4.0` to `2024`, too old tpr can not be read by this module.

# Install

* Requirements

  * `Python >= 3.8`

  * `Numpy`

  * `C++ compiler` (g++ for Linux, MSVC for Windows) supports `C++ 17` standard

* Install

  The module is installed by `pip` method:

  ```
  pip install TprParser -i https://pypi.org/simple
  ```

  Please **ALWAYS** install Latest version.

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
# ...
```



## Modify atom property

```python
newcoords = np.array([[1,2,3], [4,5,6], [...]], dtype=np.float32) # shape= N*3
# The step will create new.tpr that used newcoords
reader.set_xvf('x', newcoords)
```



# Modify system pressure

you can define a function do this work:

```Python
def Pressure(fname):
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



## Other

Please see `TprReader` module annotation


# Cite
If `TprParser` is utilized in your work, please cite as follows in main text:

> Yujie Liu, TprParser, Version [xxx](), https://pypi.org/project/TprParser/


## TODO

* More parameters can be modified
* Get More essential parameters 


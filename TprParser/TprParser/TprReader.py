from typing_extensions import TypeAlias
from typing import Literal
import numpy as np
import TprParser_
import shutil

class TprReader:
    """ @brief A wrapper of TprParser
        1. get atmic coordinates/velocity/force/box/mass/charge ... of tpr
        2. get full force field parameters for bonds/angles/dihedrals/impropers
        3. modify simulation nsteps/dt/integer/coordinates/velocity/force/box and save as new.tpr

        Parameters
        ---------
        fname: str
            The name of tpr file
        bGRO: bool, default ``False``
            If output a system gro file
        bMol2: bool, default ``False``
            If output a system mol2 file with correct bonds
        bCharge: bool, default ``False``
            If output a plain text file contained atom mass and charge
    """
    VecType: TypeAlias = Literal['x', 'X', 'v', 'V', 'f', 'F', 'box', 'BOX']
    VecType2: TypeAlias = Literal['m', 'M', 'q', 'Q']
    VecType3: TypeAlias = Literal['res', 'atom', 'type']
    VecType4: TypeAlias = Literal['resid', 'atnum']
    BondedType: TypeAlias = Literal['bonds', 'angles', 'dihedrals', 'impropers']
    NonBondedType: TypeAlias = Literal['pairs', 'lj', 'type']
    def __init__(self, fname, bGRO = False, bMol2 = False, bCharge = False) -> None:
        self.tprCapsule = TprParser_.load(fname, bGRO, bMol2, bCharge)
    
    def set_nsteps(self, nsteps):
        """ @brief set up nsteps of tpr, same as mdp

        Parameters
        ----------
        nsteps: the nsteps of simulation

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_nsteps(self.tprCapsule, nsteps)

    def set_dt(self, dt):
        """ @brief set up dt of tpr in ps, same as mdp

        Parameters
        ----------
        dt: the dt of simulation, ps

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_dt(self.tprCapsule, dt)

    def set_xvf(self, type:VecType, vec:np.array):
        """ @brief set up atomic coordinates/velocity/force/box of tpr

        Parameters
        ----------
        type: must be 'X', 'V', 'F', or 'BOX', represents atomic coordinates/velocity/force/box to set
        vec: a np.array(dtype=np.float32) of atom coordinates/velocity/force/box, 
        the dimension must be natoms * 3, except the box dimension is 9

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_xvf(self.tprCapsule, type, np.array(vec, dtype=np.float32).flatten())
    
    def set_pressure(self, epc, epct, tau_p, ref_p, compress, deform=np.zeros(9, dtype=np.float32)):
        """ @brief set up pressure coulping parts of tpr

        Parameters
        ----------
        epc: pressure coupling method, No, Berendsen, ParrinelloRahman, CRescale
        epct: pressure coupling type, Isotropic, SemiIsotropic, Anisotropic
        tau_p: the pressure coupling constant
        ref_p: a list of pressure in bar, the length must be 9
        compress: a list of compressibility in bar^-1, the length must be 9
        deform: optional, a list of deform value in nm/ps, the length must be 9, default all zero

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_pressure(self.tprCapsule, epc, epct, tau_p, ref_p, compress, deform)
    
    def set_temperature(self, etc, tau_t:list, ref_t:list):
        """ @brief set up temperature coulping parts of tpr

        Parameters
        ----------
        etc: temperature coupling method, No, Berendsen, NoseHoover, VRescale
        tau_t: the temperature coupling constant, the length must be same as old tpr
        ref_t: a list of temperature in K, the length must be same as old tpr

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_temperature(self.tprCapsule, etc, tau_t, ref_t)

    def set_mdp_integer(self, keyword:str, val:int):
        """ @brief set up integer keyword of tpr

        Parameters
        ----------
        keyword: the mdp keyword, nstlog, nstxout, nstvout, nstfout, nstenergy, nstxout_compressed,
        nsttcouple, nstpcouple, nstcalcenergy
        val: an int value for keyword

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_mdp_integer(self.tprCapsule, keyword, val)
    
    def get_prec(self):
        """ @brief get the precision of tpr

        Return
        ------
        return 4 is float, 8 is dobule
        """
        return TprParser_.get_prec(self.tprCapsule)
    
    def get_mdp_integer(self, keyword:str):
        """ @brief get integer keyword of tpr

        Parameters
        ----------
        keyword: the mdp keyword, nstlog, nstxout, nstvout, nstfout, nstenergy, nstxout_compressed,
        nsttcouple, nstpcouple, nstcalcenergy

        Returns
        -------
        return an int value for keyword
        """
        return TprParser_.get_mdp_integer(self.tprCapsule, keyword)
        
    def get_xvf(self, type:VecType) -> np.array:
        """ @brief get atomic coordinates/velocity/force/box from tpr if exist. 
        the unit is nm, nm/ps, kJ/mol/nm, nm

        Parameters
        ----------
        type: must be 'X', 'V', 'F', 'BOX', represents atomic coordinates/velocity/force/box to get

        Returns
        -------
        return a np.array(dtype=np.float32), the dimension is natoms * 3, except box is 3*3
        """
        vec = TprParser_.get_xvf(self.tprCapsule, type)
        return np.array(vec, np.float32).reshape(-1, 3)
    
    def get_mq(self, type:VecType2):
        """ @brief get atomic mass/charge from tpr

        Parameters
        ----------
        type: must be 'M', 'Q', represents atomic mass/charge to get

        Returns
        -------
        return a np.array(dtype=np.float32), the lengths is natoms
        """
        vec = TprParser_.get_xvf(self.tprCapsule, type)
        return np.array(vec, np.float32)
    
    def get_ivector(self, type:VecType4):
        """ @brief get resid/atomtypenumber from tpr

        Parameters
        ----------
        type: must be 'resid', 'atnum', represents resid/atomtypenumber to get

        Return
        ------
        return a np.array(dtype='<i'), the length is natoms for resid, the atomtypes for atnum
        """
        vec = TprParser_.get_ivector(self.tprCapsule, type)
        return np.array(vec, np.int32)

    def get_name(self, type:VecType3):
        """ @brief get resname/atomname/atomtype from tpr

        Parameters
        ----------
        type: must be 'res', 'atom', 'type', represents resname/atomname/atomtype to get

        Returns
        -------
        return a np.array(dtype='<U'), the length is natoms
        """
        vec = TprParser_.get_name(self.tprCapsule, type)
        return np.array(vec, dtype='<U')

    def get_bonded(self, type:BondedType):
        """ @brief get atom bonds/angles/dihedrals/impropers (1-based index) force field parameters from tpr if exist.

        Returns
        -------
        return a np.array(dtype=object), the length is the number of bonded. 
        For each bonded, composed of [atomid pairs + force field parameters], float precision error can be ignored

        Example:
        -------
        >>> bonds = reader.get_bonded('bonds')
        # print all information about the first bonds
        >>> print(bonds[0])     
        [1, 2, 1, 0.10100000351667404, 363171.1875, 0.10100000351667404, 363171.1875]
        # print atom index (1-based) of the first bond
        >>> print(bonds[0][:2]) 
        [1, 2]
        # print force field parameters of the first bond, includes functype+parameters
        >>> print(bonds[0][2:] )        
        [1, 0.10100000351667404, 363171.1875, 0.10100000351667404, 363171.1875]
        The first int value represents function type, corresponding to [ bonds ] function type in itp/top
        ---------------------------------------------
        >>> angles = reader.get_bonded('angles')
        # print all information about the first angle
        >>> print(angles[0]) 
        [1 5 6 1 109.5 418.3999938964844 109.5 418.3999938964844]
        # print atom index (1-based) of the first angle
        >>> print(angles[0][:3]) 
        [1 5 6]
        # print force field parameters
        >>> print(angles[0][3:]) 
        [1 109.5 418.3999938964844 109.5 418.3999938964844]
        The first int value corresponding to [ angles ] function type in itp/top
        ---------------------------------------------
        >>> dihedrals = reader.get_bonded('dihedrals')
        # print all information about the first dihedral
        >>> print(dihedrals[0]) 
        [1 5 7 8 9 0.0 0.6508399844169617 0.0 0.6508399844169617 3.0]
        # print atom index (1-based) of the first dihedral
        >>> print(dihedrals[0][:4]) 
        [1 5 6]
        # print force field parameters
        >>> print(dihedrals[0][4:]) 
        [9 0.0 0.6508399844169617 0.0 0.6508399844169617 3.0]
        The first int value corresponding to [ dihedrals ] function type in itp/top
        """
        bonded = TprParser_.get_bonded(self.tprCapsule, type)
        return np.array(bonded, dtype=object)
    
    def get_nonbonded(self, type:NonBondedType):
        
        """ @brief get pairs (1-based index)/LJ parameters of each atom/atomtype LJ information from tpr.
        
        type='pairs', the length is the number of nonbonded, composed of [atomid pairs + force field parameters] (ifunc=1)

        type='lj', the length is the number of atoms, composed of [force field parameters] (ifunc=3)

        type='type', the length is the number of [ atomtypes ], composed of [force field parameters] (ifunc=3)

        Returns
        -------
        return a np.array(dtype=object), the length is the number of nonbonded. 
        """
        nonbonded = TprParser_.get_nonbonded(self.tprCapsule, type)
        return np.array(nonbonded, dtype=object)


class SimSettings():
    """ @breif A wrapper of TprParser for setting multiple mdp parameters

    Parameters
    ---------
    See ``class TprParser`` all set_ methods

    Example
    -------
    >>> with SimSettings('input.tpr', 'output.tpr') as writer:
    >>>    writer.set_dt(0.001)
    >>>    writer.set_mdp_integer('nstxout', 100)
    
    # output.tpr 
    """  
    def __init__(self, fin, fout, bGRO=False, bMol2=False, bCharge=False) -> None:
        self.tempname = '_temp_.tpr'
        self.newname = 'new.tpr'
        self.fout = fout
        shutil.copy(fin, self.tempname) # copy src to temp.tpr

    def __movefile(self):
        """ @breif Move generated self.newname to self.tempname
        """
        shutil.move(self.newname, self.tempname)

    def set_dt(self, dt):
        reader = TprReader(self.tempname)
        reader.set_dt(dt)
        reader = None # only set None to relase handle, do not call __del__
        self.__movefile()
    
    def set_nsteps(self, nsteps):
        reader = TprReader(self.tempname)
        reader.set_nsteps(nsteps)
        reader = None
        self.__movefile()
    
    def set_mdp_integer(self, keyword: str, val: int):
        reader = TprReader(self.tempname)
        reader.set_mdp_integer(keyword, val)
        reader = None
        self.__movefile()

    def set_xvf(self, type, vec):
        reader = TprReader(self.tempname)
        reader.set_xvf(type, vec)
        reader = None
        self.__movefile()

    def set_pressure(self, epc, epct, tau_p, ref_p, compress, deform=np.zeros(9, dtype=np.float32)):
        reader = TprReader(self.tempname)
        reader.set_pressure(epc, epct, tau_p, ref_p, compress, deform)
        reader = None
        self.__movefile()

    def set_temperature(self, etc, tau_t: list, ref_t: list):
        reader = TprReader(self.tempname)
        reader.set_temperature(etc, tau_t, ref_t)
        reader = None
        self.__movefile()

    def __enter__(self):
        return self

    def __exit__(self, type, exec, tracback):
        try:
            shutil.move(self.tempname, self.fout)
        except:
            pass

from typing_extensions import TypeAlias
from typing import Literal
import numpy as np
import TprParser_

class TprReader:
    """ @brief A wrapper of TprParser_
        1. get atmic coordinates/velocity/force/mass/charge ... of tpr
        2. modify simulation nsteps/dt/integer/coordinates/velocity/force and save as new.tpr
    """
    VecType: TypeAlias = Literal['x', 'X', 'v', 'V', 'f', 'F', 'box', 'BOX']
    VecType2: TypeAlias = Literal['m', 'M', 'q', 'Q']
    VecType3: TypeAlias = Literal['res', 'atom']
    BondedType: TypeAlias = Literal['bonds', 'angles', 'dihedrals', 'impropers']
    def __init__(self, fname, bGRO = False, bMol2 = False, bCharge = False) -> None:
        # get internal object
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
        """ @brief set up atomic coordinates/velocity/force of tpr

        Parameters
        ----------
        type: must be 'X', 'V', or 'F', represents atomic coordinates/velocity/force to set
        vec: a np.array(dtype=np.float32) of atom coordinates/velocity/force, the dimension must be natoms * 3

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_xvf(self.tprCapsule, type, np.array(vec, dtype=np.float32).flatten())
    
    def set_pressure(self, epc, epct, tau_p, ref_p, compress):
        """ @brief set up pressure coulping parts of tpr

        Parameters
        ----------
        epc: pressure coupling method, No, Berendsen, ParrinelloRahman, CRescale
        epct: pressure coupling type, Isotropic, SemiIsotropic
        tau_p: the pressure coupling constant
        ref_p: a list of pressure in bar, the length must be 9
        compress: a list of compressibility in bar^-1, the length must be 9

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_pressure(self.tprCapsule, epc, epct, tau_p, ref_p, compress)
    
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
        keyword: the mdp keyword, nstlog, nstxout, nstvout, nstfout, nstenergy, nstxout_compressed, \
            nsttcouple, nstpcouple, nstcalcenergy
        val: an int value for keyword

        Returns
        -------
        return True if succeed
        """
        return TprParser_.set_mdp_integer(self.tprCapsule, keyword, val)
        
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
    
    def get_name(self, type:VecType3):
        """ @brief get resname/atomname from tpr

        Parameters
        ----------
        type: must be 'res', 'atom', represents resname/atomname to get

        Returns
        -------
        return a np.array(dtype='<U'), the lengths is natoms
        """
        vec = TprParser_.get_name(self.tprCapsule, type)
        return np.array(vec, dtype='<U')

    def get_bonded(self, type:BondedType):
        """ @brief get atom bonds/angles/dihedrals/impropers pairs (1-based index) from tpr if exist.

        Returns
        -------
        return a np.array(dtype=int), the length is the number of bonded
        """
        bonded = TprParser_.get_bonded(self.tprCapsule, type)
        dim = 2 if type=='bonds' else 3 if type=='angles' else 4
        return np.array(bonded, dtype=np.int32).reshape(-1, dim)

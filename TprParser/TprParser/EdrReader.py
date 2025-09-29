import EdrParser_

class EdrReader:
    """ @brief A wrapper of EdrParser
        Read all energies data from gromacs edr file

        Parameters
        ---------
        fname: str
            The name of edr file
    """
    def __init__(self, fname) -> None:
        self.edrCapsule = EdrParser_.load(fname)
    
    def get_ene(self):
        """ @brief get all energies data in gromacs unit, such as Time -> ps, Energy -> kJ/mol, Pressure -> bar, Length -> nm

        Returns
        -------
        return a dictionary of energies data
        """
        return EdrParser_.get_ene(self.edrCapsule)

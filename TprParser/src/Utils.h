#ifndef UTILS_H
#define UTILS_H

#include "define.h"

#include <vector>
#include <utility> // std::pair



//! \brief return bond function type id and force parameters. 
//! includes constraint derived from bonds
//! \return return std::pair(-1, {}) if failed
std::pair<int, std::vector<float>> get_bond_type(int ftype, const t_iparams* param);

//! \brief return angle function type id and force parameters,
//! \return return std::pair(-1, {}) if failed
std::pair<int, std::vector<float>> get_angle_type(int ftype, const t_iparams* param);

//! \brief return dihedral function type id and force parameters,
//! \return return std::pair(-1, {}) if failed
std::pair<int, std::vector<float>> get_dihedral_type(int ftype, const t_iparams* param);

//! \brief return impropers dihedral function type id and force parameters,
//! \return return std::pair(-1, {}) if failed
std::pair<int, std::vector<float>> get_improper_type(int ftype, const t_iparams* param);

#endif // !UTILS_H

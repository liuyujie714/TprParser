#ifndef UTILS_H
#define UTILS_H

#include "define.h"

#include <vector>
#include <utility> // std::pair



//! \brief return bond function type id and force parameters,
//! \return return type 0 if F_SETTLE, return -1 if failed
std::pair<int, std::vector<float>> get_bond_type(int Enum, const t_iparams* param);

//! \brief return angle function type id and force parameters,
//! \return return type 0 if F_SETTLE, return -1 if failed
int get_angle_type(int Enum);

//! \brief return dihedral function type id and force parameters,
//! \return return type 0 if F_SETTLE, return -1 if failed
int get_dihedral_type(int Enum);

#endif // !UTILS_H

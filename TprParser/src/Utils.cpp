#include "Utils.h"

std::pair<int, std::vector<float>> get_bond_type(int Enum, const t_iparams *param)
{
    std::vector<float> ffparam;
    switch (Enum)
    {
    case F_BONDS: 
        ffparam.push_back(param->harmonic.rA);
        ffparam.push_back(param->harmonic.krA);
        ffparam.push_back(param->harmonic.rB);
        ffparam.push_back(param->harmonic.krB);
        return std::make_pair(1, ffparam);
    case F_G96BONDS:
        ffparam.push_back(param->harmonic.rA);
        ffparam.push_back(param->harmonic.krA);
        ffparam.push_back(param->harmonic.rB);
        ffparam.push_back(param->harmonic.krB);
        return std::make_pair(2, ffparam);
    case F_MORSE:
        ffparam.push_back(param->morse.b0A);
        ffparam.push_back(param->morse.cbA);
        ffparam.push_back(param->morse.betaA);
        ffparam.push_back(param->morse.b0B);
        ffparam.push_back(param->morse.cbB);
        ffparam.push_back(param->morse.betaB);
        return std::make_pair(3, ffparam);
    case F_CUBICBONDS:
        ffparam.push_back(param->cubic.b0);
        ffparam.push_back(param->cubic.kb);
        ffparam.push_back(param->cubic.kcub);
        return std::make_pair(4, ffparam);
    case F_CONNBONDS:
        return std::make_pair(5, ffparam);
    case F_HARMONIC:
        ffparam.push_back(param->harmonic.rA);
        ffparam.push_back(param->harmonic.krA);
        ffparam.push_back(param->harmonic.rB);
        ffparam.push_back(param->harmonic.krB);
        return std::make_pair(6, ffparam);
    case F_FENEBONDS:
        ffparam.push_back(param->fene.bm);
        ffparam.push_back(param->fene.kb);
        return std::make_pair(7, ffparam);
    case F_TABBONDS:
        ffparam.push_back(param->tab.kA);
        ffparam.push_back(static_cast<float>(param->tab.table)); // int to float
        ffparam.push_back(param->tab.kB);
        return std::make_pair(8, ffparam);
    case F_TABBONDSNC:
        ffparam.push_back(param->tab.kA);
        ffparam.push_back(static_cast<float>(param->tab.table)); // int to float
        ffparam.push_back(param->tab.kB);
        return std::make_pair(9, ffparam);
    case F_RESTRBONDS:
        ffparam.push_back(param->restraint.lowA);
        ffparam.push_back(param->restraint.up1A);
        ffparam.push_back(param->restraint.up2A);
        ffparam.push_back(param->restraint.kA);
        ffparam.push_back(param->restraint.lowB);
        ffparam.push_back(param->restraint.up1B);
        ffparam.push_back(param->restraint.up2B);
        ffparam.push_back(param->restraint.kB);
        return std::make_pair(10, ffparam);
    // add settle for water
    case F_SETTLE:
        ffparam.push_back(param->settle.doh);
        ffparam.push_back(param->settle.dhh);
        return std::make_pair(0, ffparam);
    // 成键关系会转换成约束Constraint
    case F_CONSTR:
        ffparam.push_back(param->constr.dA); // 距离
        ffparam.push_back(param->constr.dB);
        return std::make_pair(1, ffparam);
    case F_CONSTRNC:
        ffparam.push_back(param->constr.dA); // 距离
        ffparam.push_back(param->constr.dB);
        return std::make_pair(2, ffparam);
    default:
        break;
    }
    return {};
}


// TODO
int get_angle_type(int Enum)
{
    switch (Enum)
    {
    case F_ANGLES:
        return 1;
    case F_G96ANGLES:
        return 2;
    case F_CROSS_BOND_BONDS:
        return 3;
    case F_CROSS_BOND_ANGLES:
        return 4;
    case F_UREY_BRADLEY:
        return 5;
    case F_QUARTIC_ANGLES:
        return 6;
    case F_TABANGLES:
        return 8;
    case F_LINEAR_ANGLES:
        return 9;
    case F_RESTRANGLES:
        return 10;
    default:
        break;
    }
    return -1;
}

// TODO
int get_dihedral_type(int Enum)
{
    switch (Enum)
    { 
    case F_PDIHS: // 周期性二面角多重
        //return 1;
        return 9;
    case F_IDIHS:
        return 2;
    case F_RBDIHS:
        return 3;
    case F_PIDIHS:
        return 4;
    case F_FOURDIHS:
        return 5;
    case F_TABDIHS:
        return 8;
    case F_RESTRDIHS:
        return 10;
    case F_CBTDIHS:
        return 11;
    default:
        break;
    }
    return -1;
}

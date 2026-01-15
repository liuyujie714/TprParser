#ifndef ENUM_H
#define ENUM_H

#if defined(_MSC_VER) || defined(_WIN32)
#    include <string.h>
#    define mystricmp _stricmp
#else
#    include <strings.h>
#    define mystricmp strcasecmp
#endif
#include <string>

#define STATIC_ASSERT_ENUM_STR(enum_type, arrstr)                                           \
    static_assert(static_cast<int>(enum_type::Count) == sizeof(arrstr) / sizeof(arrstr[0]), \
                  #enum_type " size is not matched")

enum class PbcType : int
{
    Xyz     = 0, //!< Periodic boundaries in all dimensions.
    No      = 1, //!< No periodic boundaries.
    XY      = 2, //!< Only two dimensions are periodic.
    Screw   = 3, //!< Screw.
    Unset   = 4, //!< The type of PBC is not set or invalid.
    Count   = 5,
    Default = Xyz
};

enum class SwapType : int
{
    No = 0,
    X,
    Y,
    Z
};

//! Cutoff scheme
enum class CutoffScheme : int
{
    Verlet,
    Group,
    Count,
    Default = Verlet
};

//< pressure coupling methods
enum class PressureCoupling : int
{
    No,
    Berendsen,
    ParrinelloRahman,
    Isotropic,
    Mttk,
    CRescale,
    Count
};
static const char* c_PressureCoupling[] =
    {"No", "Berendsen", "ParrinelloRahman", "Isotropic", "Mttk", "CRescale"};
STATIC_ASSERT_ENUM_STR(PressureCoupling, c_PressureCoupling);

//< pressure coupling type
enum class PressureCouplingType : int
{
    Isotropic,
    SemiIsotropic,
    Anisotropic,
    Count
};
static const char* c_PressureCouplingType[] = {"Isotropic", "SemiIsotropic", "Anisotropic"};
STATIC_ASSERT_ENUM_STR(PressureCouplingType, c_PressureCouplingType);

//< temperature coupling methods
enum class TemperatureCoupling : int
{
    No,
    Berendsen,
    NoseHoover,
    Yes,
    Andersen,
    AndersenMassive,
    VRescale,
    Count,
};
static const char* c_TemperatureCoupling[] =
    {"No", "Berendsen", "NoseHoover", "Yes", "Andersen", "AndersenMassive", "VRescale"};
STATIC_ASSERT_ENUM_STR(TemperatureCoupling, c_TemperatureCoupling);

// Integer mdp
enum class ParamsInteger : int
{
    nstlog,
    nstxout,
    nstvout,
    nstfout,
    nstenergy,
    nstxout_compressed,
    nsttcouple,
    nstpcouple,
    nstcalcenergy,
    nstlist,
    nstcomm,
    fourier_nx,
    fourier_ny,
    fourier_nz,
    userint1,
    userint2,
    userint3,
    userint4,
    cutoff_scheme,
    Count,
};
static const char* c_mdp_integer[] = {"nstlog",
                                      "nstxout",
                                      "nstvout",
                                      "nstfout",
                                      "nstenergy",
                                      "nstxout_compressed",
                                      "nsttcouple",
                                      "nstpcouple",
                                      "nstcalcenergy",
                                      "nstlist",
                                      "nstcomm",
                                      "fourier_nx",
                                      "fourier_ny",
                                      "fourier_nz",
                                      "userint1",
                                      "userint2",
                                      "userint3",
                                      "userint4",
                                      "cutoff_scheme"};
STATIC_ASSERT_ENUM_STR(ParamsInteger, c_mdp_integer);

// float mdp
enum class ParamsFloat : int
{
    dt,
    rlist,
    rvdw,
    rcoulomb,
    rvdw_switch,
    rcoulomb_switch,
    tau_p,
    verletbuf_tol,
    x_compression_precision,
    verletBufferPressureTolerance,
    epsilon_r,
    epsilon_rf,
    fourier_spacing,
    em_stepsize,
    em_tol,
    shake_tol,
    cos_accel,
    userreal1,
    userreal2,
    userreal3,
    userreal4,

    Count,
};
static const char* c_mdp_float[] = {"dt",
                                    "rlist",
                                    "rvdw",
                                    "rcoulomb",
                                    "rvdw_switch",
                                    "rcoulomb_switch",
                                    "tau_p",
                                    "verletbuf_tol",
                                    "x_compression_precision",
                                    "verletBufferPressureTolerance",
                                    "epsilon_r",
                                    "epsilon_rf",
                                    "fourier_spacing",
                                    "em_stepsize",
                                    "em_tol",
                                    "shake_tol",
                                    "cos_accel",
                                    "userreal1",
                                    "userreal2",
                                    "userreal3",
                                    "userreal4"};
STATIC_ASSERT_ENUM_STR(ParamsFloat, c_mdp_float);

// vector of tpr, X or V or F
enum class VecProps : int
{
    x,
    v,
    f,
    m,   // the mass of atoms
    q,   // the charge of atoms
    box, // the box vector
    ef,  // electric field
    Count
};
static const char* c_mdp_vector[] = {"x", "v", "f", "m", "q", "box", "ef"};
STATIC_ASSERT_ENUM_STR(VecProps, c_mdp_vector);

// int vector of tpr, such resid
enum class IVectorProps : int
{
    resid,
    atnum,     // atomtype number
    atomicnum, // atomic number
    Count
};
static const char* c_int_vector[] = {"resid", "atnum", "atomicnum"};
STATIC_ASSERT_ENUM_STR(IVectorProps, c_int_vector);

// vector of tpr, resname / atomname / atomtype name
enum class StringType : int
{
    res,
    atom,
    type, // atomtype name
    Count
};
static const char* c_name_vector[] = {"res", "atom", "type"};
STATIC_ASSERT_ENUM_STR(StringType, c_name_vector);


// type of bonded
enum class BondedType : int
{
    bonds,
    angles,
    dihedrals,
    impropers,
    cmaps,
    Count
};
static const char* c_bonded_type[] = {"bonds", "angles", "dihedrals", "impropers", "cmaps"};
STATIC_ASSERT_ENUM_STR(BondedType, c_bonded_type);

//< Non bonded type
enum class NonBondedType : int
{
    LJ,
    atomtype, // only LJ [ atomtypes ]
    LJ_14,    // that is [ pairs ]
    BH,       // only Buckingham [ atomtypes ]
    Count
};
static const char* c_nonbonded_type[] = {"lj", "type", "pairs", "bh"};
STATIC_ASSERT_ENUM_STR(NonBondedType, c_nonbonded_type);

//! convert array to string with ', should be ' prefix
template<int N>
static inline std::string arr_to_string(const char* (&arr)[N])
{
    std::string ret = ", should be ";
    for (int i = 0; i < N; i++)
    {
        ret += arr[i];
        if (i < N - 1) { ret += ", "; }
    }
    return ret;
}

//< check key words in a c_string array ignore case, return enum value if find, else return ENUM::Count
template<typename ENUM, const int count = static_cast<int>(ENUM::Count), int N>
static inline ENUM check_string(const char* str, const char* (&arr)[N])
{
    //! check length must be equal
    static_assert(N == count, "c_string length is not equal to enum length");
    for (int i = 0; i < count; i++)
    {
        if (!mystricmp(str, arr[i])) return static_cast<ENUM>(i);
    }
    return ENUM::Count;
}


//! Pulling algorithm.
enum class PullingAlgorithm : int
{
    Umbrella,
    Constraint,
    ConstantForce,
    FlatBottom,
    FlatBottomHigh,
    External,
    Count,
    Default = Umbrella
};

//! Control of pull groups
enum class PullGroupGeometry : int
{
    Distance,
    Direction,
    Cylinder,
    DirectionPBC,
    DirectionRelative,
    Angle,
    Dihedral,
    AngleAxis,
    Transformation,
    Count,
    Default = Distance
};


//! Enforced rotation group type.
enum class EnforcedRotationGroupType : int
{
    Iso,
    Isopf,
    Pm,
    Pmpf,
    Rm,
    Rmpf,
    Rm2,
    Rm2pf,
    Flex,
    Flext,
    Flex2,
    Flex2t,
    Count,
    Default = Iso
};

//! Rotation group fitting type
enum class RotationGroupFitting : int
{
    Rmsd,
    Norm,
    Pot,
    Count,
    Default = Rmsd
};

// for old tpr, such as gmx < 4.0, ref: https://github.com/gromacs/gromacs/commit/71e36ecc6806c373f452fa70ca24f43d15f15621
enum
{
    eelCUT,
    eelRF,
    eelGRF_NOTUSED,
    eelPME,
    eelEWALD,
    eelP3M_AD,
    eelPOISSON,
    eelSWITCH,
    eelSHIFT,
    eelUSER,
    eelGB_NOTUSED,
    eelRF_NEC_UNSUPPORTED,
    eelENCADSHIFT_NOTUSED,
    eelPMEUSER,
    eelPMESWITCH,
    eelPMEUSERSWITCH,
    eelRF_ZERO,
    eelNR
};
#define EEL_RF(e) \
    ((e) == eelRF || (e) == eelGRF_NOTUSED || (e) == eelRF_NEC_UNSUPPORTED || (e) == eelRF_ZERO)

enum
{
    eiMD,
    eiSteep,
    eiCG,
    eiBD,
    eiSD2_REMOVED,
    eiNM,
    eiLBFGS,
    eiTPI,
    eiTPIC,
    eiSD1,
    eiVV,
    eiVVAK,
    eiMimic,
    eiNR
};

#endif // !ENUM_H

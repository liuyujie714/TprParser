#ifndef UTILS_H
#define UTILS_H

#include <cinttypes>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <utility> // std::pair
#include <vector>

#include "define.h"
#include "FileSerializer.h"
#include "TprData.h"

/* \brief return bond function type id and force parameters.
 * includes constraint derived from bonds
 * \return return std::pair(-1, {}) if failed
 */
std::pair<int, std::vector<float>> get_bond_type(int ftype, const t_iparams* param);

/* \brief return angle function type id and force parameters,
 * \return return std::pair(-1, {}) if failed
 */
std::pair<int, std::vector<float>> get_angle_type(int ftype, const t_iparams* param);

/* \brief return dihedral function type id and force parameters,
 * \return return std::pair(-1, {}) if failed
 */
std::pair<int, std::vector<float>> get_dihedral_type(int ftype, const t_iparams* param);

/* \brief return impropers dihedral function type id and force parameters,
 * \return return std::pair(-1, {}) if failed
 */
std::pair<int, std::vector<float>> get_improper_type(int ftype, const t_iparams* param);

/* \brief return virtual sites function type id and force parameters,
 * \return return std::pair(-1, {}) if failed
 * \return return std::pair(vsiten.n, {}) only for F_VSITEN
 */
std::pair<int, std::vector<float>> get_vsite_type(int ftype, const t_iparams* param);

/* \brief return nonbonded (LJ/LJ-14) function type id and force parameters,
 * ifunc=1 is LJ-14, ifunc=3 is LJ, ifunc=2 is Buckingham
 * \return return std::pair(-1, {}) if failed
 */
std::pair<int, std::vector<float>> get_nonbonded_type(int ftype, const t_iparams* param);

/* \brief safely fopen */
FILE* efopen(const char* fname, const char* mod);

/* \brief A class to read applied forces from tpr file */
class AppliedForces
{
public:
    AppliedForces(const FileSerializer& tpr, std::unique_ptr<TprData>& data);

    //! 执行tpr解序列化操作
    bool deserialize();

public:
    //! 解序列化字典
    std::map<unsigned char, std::function<void(AppliedForces*)>> s_deserializers;
    //! 字符串 -> 浮点数组
    std::map<std::string, std::vector<float>> m_float;
    //! 字符串 -> Double数组 (强转为了float)
    std::map<std::string, std::vector<float>> m_double;
    //! 字符串 -> Int数组
    std::map<std::string, std::vector<int>> m_int;
    std::string                             m_name;    //! current string name
    unsigned char                           m_typeTag; //! current typeTag
    const FileSerializer&                   tpr_;      //! tpr reference
    std::unique_ptr<TprData>&               data_;     //! TprData reference
};

// clang-format off
// \brief Print vector, do nothing if vector is empty
template<typename T>
static inline void print_vec_release(const char* name, const std::vector<T>& vec, const int Ncol = 5, FILE* fp = stdout)
{
    if (vec.empty())
    {
        return;
    }

    fprintf(fp, "%-40s   ", name);
    if constexpr (std::is_same_v<T, int> || std::is_same_v<T, int64_t>)
    {
        fprintf(fp, "I");
    }
    else
    {
        fprintf(fp, "R");
    }
    fprintf(fp, "   N=%12zu\n", vec.size());

    int idx = 0;
    for (auto& i : vec)
    {
        if constexpr (std::is_same_v<T, int>)
        {
            fprintf(fp, "%12d", i);
        }
        else if constexpr (std::is_same_v<T, int64_t>)
        {
            fprintf(fp, "%12" PRId64, i);
        }
        else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
        {
            fprintf(fp, "%16.8E", i);
        }
        else
        {
            static_assert(dependent_false<T>, "Unsupported type in print_vec_release");
        }

        if (++idx % Ncol == 0)
        {
            fprintf(fp, "\n");
        }
    }

    if (idx % Ncol != 0)
    {
        fprintf(fp, "\n");
    }
}
// clang-format on


// clang-format off
// 特殊处理字符串
template<>
inline void print_vec_release(const char* name, const std::vector<unsigned char>& vec, const int Ncol, FILE* fp)
{
    if (vec.empty())
    {
        return;
    }

    fprintf(fp, "%-40s   ", name);
    const size_t width = 12;
    size_t N = vec.size() / width;
    if (vec.size() % width != 0)
    {
        N++;
    }
    fprintf(fp, "C   N=%12zu\n", N);

    int idx = 0;
    for (size_t i = 0; i < vec.size(); i += width)
    {
        size_t len = std::min(width, vec.size() - i);
        fprintf(fp, "%.*s", (int)len, vec.data() + i);
        // to match formchk output format
        for (size_t j = len; j < width; ++j)
        {
            fprintf(fp, " ");
        }

        if (++idx % Ncol == 0)
        {
            fprintf(fp, "\n");
        }
    }

    if (idx % Ncol != 0)
    {
        fprintf(fp, "\n");
    }
}
// clang-format on

//! Only used for debug
template<typename T>
static inline void print_vec_debug(const char* name, const std::vector<T>& vec, const int Ncol = 5)
{
#ifdef _DEBUG
    print_vec_release(name, vec, Ncol);
#endif // DEBUG
}

template<typename T>
static inline void print_vec_debug(const char* name, const T* vec, const int ndim = 9, const int Ncol = 3)
{
#ifdef _DEBUG
    std::vector<T> tempv;
    tempv.assign(vec, vec + ndim);
    print_vec_release(name, tempv, Ncol, stdout);
#endif // DEBUG
}

#endif // !UTILS_H

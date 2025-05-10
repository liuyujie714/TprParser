#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <utility> // std::pair
#include <vector>
#include <typeindex>
#include <map>
#include <string>
#include <functional>

union t_iparams;
class FileSerializer;

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

//! \brief return nonbonded (LJ/LJ-14) function type id and force parameters,
//! ifunc=1 is LJ-14, ifunc=3 is LJ
//! \return return std::pair(-1, {}) if failed
std::pair<int, std::vector<float>> get_nonbonded_type(int ftype, const t_iparams* param);

//! \brief safely fopen
FILE* efopen(const char* fname, const char* mod);


//! TODO:
class KeyValueTreeObj
{
public:

};

class KeyValueTreeArray
{
public:

};


template<typename T>
struct Deserializer;

template<>
struct Deserializer<KeyValueTreeObj>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<KeyValueTreeArray>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<std::string>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<bool>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<char>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<unsigned char>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<int>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<int64_t>
{
	static void deserialize()
	{

	}
};

template<>
struct Deserializer<float>
{
	static void deserialize()
	{

	}
};


template<>
struct Deserializer<double>
{
	static void deserialize()
	{

	}
};


struct Serializer
{
	unsigned char tag;
	std::function<void ()> deserialize;
};

#define SERIALIZER(tag, type) \
{ \
	std::type_index(typeid(type)), \
	{ \
		tag, &Deserializer<type>::deserialize \
	} \
}

static const std::map<std::type_index, Serializer> c_deserializers = {
	SERIALIZER('O', KeyValueTreeObj),
	SERIALIZER('A', KeyValueTreeArray),
	SERIALIZER('s', std::string),
	SERIALIZER('b', bool),
	SERIALIZER('c', char),
	SERIALIZER('u', unsigned char),
	SERIALIZER('i', int),
	SERIALIZER('l', int64_t),
	SERIALIZER('f', float),
	SERIALIZER('d', double),
};

/* \brief A class to read applied forces from tpr file
*/
class AppliedForces
{
public:
	AppliedForces(const FileSerializer& tpr);

	void deserialize();

private:
	const FileSerializer& tpr_;
};

#endif // !UTILS_H

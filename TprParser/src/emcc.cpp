#include <emscripten/bind.h>

#include "Reader.h"

class TprReaderWrapper : public TprReader
{
public:
    TprReaderWrapper(const std::string& fname, bool bGRO, bool bMol2, bool bCharge)
        : TprReader(fname.c_str(), bGRO, bMol2, bCharge)
    {
    }


    int                       get_precision() const { return TprReader::get_precision(); }
    const std::vector<float>& get_ef() const { return TprReader::get_ef(); }

    int get_mdp_integer(const std::string& prop) const
    {
        return TprReader::get_mdp_integer(prop.c_str());
    }

    const std::vector<float>& get_xvf(const std::string& type) const
    {
        return TprReader::get_xvf(type.c_str());
    }

    const std::vector<int>& get_ivector(const std::string& type) const
    {
        return TprReader::get_ivector(type.c_str());
    }

    const std::vector<std::string>& get_name(const std::string& type) const
    {
        return TprReader::get_name(type.c_str());
    }

    const std::vector<Bonded>& get_bonded(const std::string& type) const
    {
        return TprReader::get_bonded(type.c_str());
    }

    const std::vector<NonBonded>& get_nonbonded(const std::string& type) const
    {
        return TprReader::get_nonbonded(type.c_str());
    }
};

EMSCRIPTEN_BINDINGS(my_module)
{
    emscripten::class_<Bonded>("Bonded")
        .constructor<int, int, int, const std::vector<float>&>()
        .constructor<int, int, int, int, const std::vector<float>&>()
        .constructor<int, int, int, int, int, const std::vector<float>&>()
        .property("a", &Bonded::a)
        .property("b", &Bonded::b)
        .property("c", &Bonded::c)
        .property("d", &Bonded::d)
        .property("ff", &Bonded::ff)
        .property("ifunc", &Bonded::ifunc);

    emscripten::register_vector<int>("VecInt");
    emscripten::register_vector<float>("VecFloat");
    emscripten::register_vector<std::string>("VecString");
    emscripten::register_vector<Bonded>("VecBonded");

    emscripten::class_<TprReaderWrapper>("TprReader")
        .constructor<const std::string&, bool, bool, bool>()
        .function("get_precision", &TprReaderWrapper::get_precision)
        .function("get_ef", &TprReaderWrapper::get_ef)
        .function("get_mdp_integer", &TprReaderWrapper::get_mdp_integer)
        .function("get_xvf", &TprReaderWrapper::get_xvf)
        .function("get_ivector", &TprReaderWrapper::get_ivector)
        .function("get_name", &TprReaderWrapper::get_name)
        .function("get_bonded", &TprReaderWrapper::get_bonded)
        // .function("get_nonbonded", &TprReaderWrapper::get_nonbonded)
        ;
}

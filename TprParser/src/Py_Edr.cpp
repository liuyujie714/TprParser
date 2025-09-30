#ifndef _DEBUG

#    include "EdrReader.h"
#    include "Python.h"

//===========================EdrParser========================================

#    define MODULE_NAME "EdrParser_"

// free
static void destory_edr(PyObject* obj)
{
    delete (EdrReader*)PyCapsule_GetPointer(obj, MODULE_NAME);
}

static PyObject* new_edr(PyObject* self, PyObject* args)
{
    const char* fname = NULL;

    if (!PyArg_ParseTuple(args, "s", &fname)) { return NULL; }

    EdrReader* reader = NULL;
    try
    {
        reader = new EdrReader(fname);
    }
    catch (const std::exception& e)
    {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    return PyCapsule_New(reader, MODULE_NAME, destory_edr);
}

// get energy from edr - returns dictionary of energy data
static PyObject* get_ene(PyObject* self, PyObject* args)
{
    PyObject* capsule = NULL;

    if (!PyArg_ParseTuple(args, "O", &capsule)) { return NULL; }

    // get energy data from reader
    EdrReader* reader = static_cast<EdrReader*>(PyCapsule_GetPointer(capsule, MODULE_NAME));
    if (!reader)
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
        return NULL;
    }

    std::map<std::string, std::vector<double>> energy_data;
    try
    {
        energy_data = reader->get_ene();
    }
    catch (const std::exception& e)
    {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

    if (energy_data.empty())
    {
        PyErr_SetString(PyExc_RuntimeError, "No energy data found in edr");
        return NULL;
    }

    // 创建 Python 字典
    PyObject* dict = PyDict_New();
    if (!dict)
    {
        PyErr_SetString(PyExc_RuntimeError, "Cannot create dictionary");
        return NULL;
    }

    // 遍历能量数据并添加到字典中
    for (const auto& pair : energy_data)
    {
        const std::string&         key    = pair.first;
        const std::vector<double>& values = pair.second;

        // 创建 Python 列表来存储数值
        PyObject* value_list = PyList_New(values.size());
        if (!value_list)
        {
            PyErr_SetString(PyExc_RuntimeError, "Cannot create list for values");
            Py_DECREF(dict);
            return NULL;
        }

        // 将 double 值转换为 Python float 对象并添加到列表中
        for (Py_ssize_t i = 0; i < values.size(); i++)
        {
            PyObject* py_value = PyFloat_FromDouble(values[i]);
            if (!py_value)
            {
                PyErr_SetString(PyExc_RuntimeError, "Cannot convert double to Python float");
                Py_DECREF(value_list);
                Py_DECREF(dict);
                return NULL;
            }
            PyList_SET_ITEM(value_list, i, py_value);
        }

        // 创建 Python 字符串作为键
        PyObject* py_key = PyUnicode_FromString(key.c_str());
        if (!py_key)
        {
            PyErr_SetString(PyExc_RuntimeError, "Cannot create key string");
            Py_DECREF(value_list);
            Py_DECREF(dict);
            return NULL;
        }

        // 将键值对添加到字典中
        if (PyDict_SetItem(dict, py_key, value_list) < 0)
        {
            PyErr_SetString(PyExc_RuntimeError, "Cannot set dictionary item");
            Py_DECREF(py_key);
            Py_DECREF(value_list);
            Py_DECREF(dict);
            return NULL;
        }

        // 减少引用计数（PyDict_SetItem 会增加引用）
        Py_DECREF(py_key);
        Py_DECREF(value_list);
    }

    return dict;
}

static PyMethodDef edr_methods[] = {
    {"load", new_edr, METH_VARARGS, "Create a new EdrParser instance"},
    {"get_ene", get_ene, METH_VARARGS, "Get all energys from edr"},

    {NULL, NULL, 0, NULL}};

static struct PyModuleDef edr_module = {
    PyModuleDef_HEAD_INIT,
    MODULE_NAME, // m_name
    NULL,        // m_doc
    -1,          // m_size
    edr_methods  // m_methods
};

// EdrParser_  module name
PyMODINIT_FUNC PyInit_EdrParser_(void)
{
    return PyModule_Create(&edr_module);
}

#endif

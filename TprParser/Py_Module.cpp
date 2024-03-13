#include "Python.h"
#include "numpy/arrayobject.h"
#include "Reader.h"
#include <string.h>

// free 
static void destory_tpr(PyObject* obj)
{
	delete (TprReader*)PyCapsule_GetPointer(obj, "TprParser");
}

static PyObject* reader_new(PyObject* self, PyObject* args) 
{
	const char		* fname = nullptr;
	PyObject		* bGRO_obj = Py_False;
	PyObject		* bMol2_obj = Py_False;
	PyObject		* bCharge_obj = Py_False;

	if (!PyArg_ParseTuple(args, "s|OOO", &fname, &bGRO_obj, &bMol2_obj, &bCharge_obj))
	{
		return NULL;
	}

	bool bGRO = PyObject_IsTrue(bGRO_obj);
	bool bMol2 = PyObject_IsTrue(bMol2_obj);
	bool bCharge = PyObject_IsTrue(bCharge_obj);

	TprReader* reader = nullptr;
	try
	{
		reader = new TprReader(fname, bGRO, bMol2, bCharge);
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}
	return PyCapsule_New(reader, "TprParser", destory_tpr);
}

static PyObject* set_nsteps(PyObject* self, PyObject* args) 
{
	PyObject	* capsule = nullptr;
	int64_t		nsteps;
	if (!PyArg_ParseTuple(args, "OL", &capsule, &nsteps)) 
	{
		return nullptr;
	}

	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader) 
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	try
	{
		bool ret = reader->set_nsteps(nsteps);
		if (ret == TPR_FAILED)
		{
			PyErr_SetString(PyExc_RuntimeError, "set_nsteps faliled");
			return nullptr;
		}
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

static PyObject* set_dt(PyObject* self, PyObject* args) 
{
	PyObject	* capsule = nullptr;
	double      dt = 0.0;
	if (!PyArg_ParseTuple(args, "Od", &capsule, &dt)) 
	{
		return nullptr;
	}

	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader)
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	try
	{
		bool ret = reader->set_dt(dt);
		if (ret == TPR_FAILED)
		{
			PyErr_SetString(PyExc_RuntimeError, "set_dt faliled");
			return nullptr;
		}
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

//< get vector from given object, return nullptr if failed
static inline PyObject * get_vector(PyObject* vec_obj, std::vector<float> &vec)
{
	if (PyArray_Check(vec_obj))
	{
		// numpy to C array
		PyArrayObject* arr = (PyArrayObject*)(vec_obj);

		// 比较数据类型代码
		if (PyArray_TYPE(arr) != NPY_FLOAT32) {
			PyErr_SetString(PyExc_RuntimeError, "Python API Only support np.float32 array");
			return nullptr;
		}

		npy_intp* dims = PyArray_DIMS(arr);
		int ndim = PyArray_NDIM(arr);
		if (ndim != 1)
		{
			PyErr_SetString(PyExc_RuntimeError, "Input numpy dimension is not equal 1");
			return nullptr;
		}
		float* data = (float*)PyArray_DATA(vec_obj);
		vec.assign(data, data + dims[0]);
	}
	// if is list
	else if (PyList_Check(vec_obj))
	{
		Py_ssize_t size = PyList_Size(vec_obj);
		for (Py_ssize_t i = 0; i < size; i++)
		{
			PyObject* item = PyList_GetItem(vec_obj, i);
			//if (!PyFloat_Check(item))
			//{
			//	PyErr_SetString(PyExc_RuntimeError, "List must be float type");
			//	return nullptr;
			//}
			double value = PyFloat_AsDouble(item);
			if (PyErr_Occurred())
			{
				PyErr_SetString(PyExc_RuntimeError, "Some value of input vector can not be converted");
				return nullptr;
			}
			vec.push_back((float)(value)); // double to float
		}
	}
	else
	{
		PyErr_SetString(PyExc_RuntimeError, "Input vector must be numpy array or list");
		return nullptr;
	}
	// check array size
	if (vec.empty())
	{
		PyErr_SetString(PyExc_RuntimeError, "Empty input vector");
		return nullptr;
	}

	return vec_obj;
}

static PyObject* set_coordinates(PyObject* self, PyObject* args, PyObject* kwargs)
{
	PyObject* capsule = nullptr;
	PyObject* coords_obj = nullptr;

	static const char* keywords[] = { "capsule", "coords", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", (char**)keywords, &capsule, &coords_obj))
	{
		return nullptr;
	}

	// get coords
	std::vector<float> coords;
	if (!get_vector(coords_obj, coords))
	{
		return nullptr;
	}

	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader) 
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	try
	{
		bool ret = reader->set_coordinates(coords);
		if (ret == TPR_FAILED)
		{
			PyErr_SetString(PyExc_RuntimeError, "set_coordinates faliled");
			return nullptr;
		}
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

static PyObject* set_pressure(PyObject* self, PyObject* args, PyObject* kwargs)
{
	PyObject		* capsule = nullptr;
	const char		* epc = nullptr;
	const char		* epct = nullptr;
	float			tau_p;
	PyObject		* ref_p = nullptr;
	PyObject		* compress = nullptr;

	static const char* keywords[] = { "capsule", "epc", "epct", "tau_p", "ref_p", "compress", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OssfOO",
		(char**)keywords, &capsule, &epc, &epct, &tau_p, &ref_p, &compress))
	{
		return nullptr;
	}

	// get pressure
	std::vector<float> vec_press;
	if (!get_vector(ref_p, vec_press)) return nullptr;

	// get compress
	std::vector<float> vec_compress;
	if (!get_vector(compress, vec_compress)) return nullptr;

	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader)
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	try
	{
		bool ret = reader->set_pressure(epc, epct, tau_p, vec_press, vec_compress);
		if (ret == TPR_FAILED)
		{
			PyErr_SetString(PyExc_RuntimeError, "set_pressure faliled");
			return nullptr;
		}
	}
	catch (const std::exception& e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

static PyObject* set_temperature(PyObject* self, PyObject* args, PyObject* kwargs)
{
	PyObject* capsule = nullptr;
	const char* etc = nullptr;
	PyObject* ref_t = nullptr;
	PyObject* tau_t = nullptr;

	static const char* keywords[] = { "capsule", "etc", "tau_t", "ref_t", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OsOO",
		(char**)keywords, &capsule, &etc, &tau_t, &ref_t))
	{
		return nullptr;
	}

	// get tau_t
	std::vector<float> vec_tau;
	if (!get_vector(tau_t, vec_tau)) return nullptr;

	// get ref_t
	std::vector<float> vec_t;
	if (!get_vector(ref_t, vec_t)) return nullptr;

	// get handle
	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader)
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	try
	{
		bool ret = reader->set_temperature(etc, vec_tau, vec_t);
		if (ret == TPR_FAILED)
		{
			PyErr_SetString(PyExc_RuntimeError, "set_temperature faliled");
			return nullptr;
		}
	}
	catch (const std::exception& e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

static PyObject* get_coordinates(PyObject* self, PyObject* args)
{
	PyObject* capsule = nullptr;

	// get object handle
	if (!PyArg_ParseTuple(args, "O", &capsule))
	{
		return nullptr;
	}

	TprReader* reader = static_cast<TprReader*>(PyCapsule_GetPointer(capsule, "TprParser"));
	if (!reader)
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid capsule object");
		return nullptr;
	}

	std::vector<float> coords;
	try
	{
		coords = reader->get_coordinates();
	}
	catch (const std::exception& e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}
	if (coords.empty())
	{
		PyErr_SetString(PyExc_RuntimeError, "Can not find coords in tpr");
		return nullptr;
	}

	// coords to python list
	PyObject* list = PyList_New(coords.size());
	if (!list)
	{
		PyErr_SetString(PyExc_RuntimeError, "Can not new list for coords");
		return nullptr;
	}
	Py_ssize_t i = 0;
	for (const auto& coord : coords)
	{
		PyObject* value = PyFloat_FromDouble(static_cast<double>(coord));
		if (!value)
		{
			PyErr_SetString(PyExc_RuntimeError, "Can not convert coords to list");
			Py_DECREF(list); // free list
			return nullptr;
		}
		PyList_SET_ITEM(list, i++, value);
	}

	return list;
}

static PyMethodDef methods[] =
{
	{"load", reader_new, METH_VARARGS, "Create a new TprReader instance"},
	{"set_nsteps", set_nsteps, METH_VARARGS, "Set up nsteps"},
	{"set_dt", set_dt, METH_VARARGS, "Set up dt"},
	{"set_coordinates", (PyCFunction)set_coordinates, METH_VARARGS | METH_KEYWORDS, "Set up atomic coordinates"},
	{"set_pressure", (PyCFunction)set_pressure, METH_VARARGS | METH_KEYWORDS, "Set up pressure coupling parts"},
	{"set_temperature", (PyCFunction)set_temperature, METH_VARARGS | METH_KEYWORDS, "Set up temperature coupling parts"},

	{"get_coordinates", get_coordinates, METH_VARARGS, "Get coords from tpr"}, 
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef tpr_module = 
{
	PyModuleDef_HEAD_INIT,
	"TprParser",		//m_name
	NULL,				//m_doc
	-1,					//m_size
	methods				//m_methods
};

PyMODINIT_FUNC PyInit_TprParser(void) 
{
	// if is numpy
	import_array() // 使用numpy相关的函数时候必须先调用这个

	return PyModule_Create(&tpr_module);
}

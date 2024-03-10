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
		reader->set_nsteps(nsteps);
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
		reader->set_dt(dt);
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
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

	std::vector<float> coords;
	// if is numpy
	import_array() // 使用nymoy相关的函数时候必须先调用这个
	if (PyArray_Check(coords_obj))
	{
		// numpy to C array
		PyArrayObject* arr = (PyArrayObject*)(coords_obj);

		// 比较数据类型代码
		if (PyArray_TYPE(arr) != NPY_FLOAT32) {
			PyErr_SetString(PyExc_RuntimeError, "Only support np.float32 array");
			return nullptr;
		}

		npy_intp* dims = PyArray_DIMS(arr);
		int ndim = PyArray_NDIM(arr);
		if (ndim != 1)
		{
			PyErr_SetString(PyExc_RuntimeError, "Input numoy dimension is not equal 1");
			return nullptr;
		}
		float* data = (float*)PyArray_DATA(coords_obj);
		coords.assign(data, data + dims[0]);
	}
	// if is list
	else if (PyList_Check(coords_obj))
	{
		Py_ssize_t size = PyList_Size(coords_obj);
		for (Py_ssize_t i = 0; i < size; i++)
		{
			PyObject* item = PyList_GetItem(coords_obj, i);
			if (!PyFloat_Check(item))
			{
				PyErr_SetString(PyExc_RuntimeError, "List must be float type");
				return nullptr;
			}
			coords.push_back((float)(PyFloat_AsDouble(item)));
		}
	}
	else
	{
		PyErr_SetString(PyExc_RuntimeError, "Input coords must be numpy array or list");
		return nullptr;
	}
	// check array size
	if (coords.empty())
	{
		PyErr_SetString(PyExc_RuntimeError, "Empty input coords");
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
		reader->set_coordinates(coords);
	}
	catch (const std::exception&e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}

	Py_RETURN_TRUE;
}

static PyMethodDef methods[] =
{
	{"load", reader_new, METH_VARARGS, "Create a new TprReader instance"},
	{"set_nsteps", set_nsteps, METH_VARARGS, "Set up nsteps"},
	{"set_dt", set_dt, METH_VARARGS, "Set up dt"},
	{"set_coordinates", (PyCFunction)set_coordinates, METH_VARARGS | METH_KEYWORDS, "Set up atomic coordinates"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef tpr_module = 
{
	PyModuleDef_HEAD_INIT,
	"TprParser",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit_TprParser(void) 
{
	return PyModule_Create(&tpr_module);
}

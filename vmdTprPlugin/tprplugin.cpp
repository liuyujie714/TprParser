#define VMDPLUGIN_EXPORTS

#if (defined(_WIN32) || defined(__CYGWIN__)) && !defined(WIN32)
#define WIN32 // is needed for windows
#endif // !WIN32

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "molfile_plugin.h"
#include "../TprParser/src/Reader.h"

#ifndef M_PI
#define M_PI            3.14159265358979
#define M_PI_2          (M_PI/2.0)
#endif // !M_PI

//< safe use malloc, throw error if not allocated
template<typename T>
T* emalloc(size_t nelem)
{
	T* p = static_cast<T*>(malloc(sizeof(T) * nelem));
	if (!p)
	{
		throw std::runtime_error("Error! Can not allocated!");
	}
	return p;
}

struct TprPlugin
{
	TprPlugin() : tpr(NULL), from(NULL), to(NULL),
		angles(NULL), dihedrals(NULL), impropers(NULL), has_read(false)
	{
	}

	~TprPlugin()
	{
		if (tpr) {
			delete tpr; tpr = NULL;
		}
		if (from) {
			free(from); from = NULL;
		}
		if (to) {
			free(to); to = NULL;
		}
		if (angles) {
			free(angles); angles = NULL;
		}
		if (dihedrals) {
			free(dihedrals); dihedrals = NULL;
		}
		if (impropers) {
			free(impropers); impropers = NULL;
		}
	}

	TprReader	* tpr = NULL; // tpr handle
	//! some pointer to record address, for free memory
	int			* from, * to;
	int			*angles, *dihedrals, *impropers;
	bool		has_read = false; //< if read coords
};

#ifdef __cplusplus
extern "C" {
#endif // _cplusplus

// get tpr total atoms
static void* open_tpr_read(const char* fname, const char* ftype, int* natoms)
{
	TprPlugin * tprplugin = NULL;
	try
	{
		tprplugin = new TprPlugin;
		tprplugin->tpr = new TprReader(fname);
	}
	catch (const std::exception&e)
	{
		fprintf(stderr, "tprplugin) ERROR: %s\n", e.what());
		return NULL;
	}
	catch (...)
	{
		fprintf(stderr, "tprplugin) Unknown exception\n");
		return NULL;
	}

	// return natoms and tpr class
	*natoms = static_cast<int>(tprplugin->tpr->get_name("atom").size());
	// init value
	tprplugin->has_read = false;

	return tprplugin;
}

static void close_tpr_read(void* mydata)
{
	TprPlugin * tprplugin = static_cast<TprPlugin*>(mydata);
	delete tprplugin;
}

// get atom mass, charge, resname, atomname, atomtype
static int read_tpr_structure(void* mydata, int* optflags, molfile_atom_t* atoms)
{
	// 提供原子属性
	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_CHARGE;

	TprPlugin* tprplugin = static_cast<TprPlugin*> (mydata);
	auto tpr = tprplugin->tpr;
	std::vector<float> mass, charge;
	std::vector<int> resid;
	std::vector<std::string> resname, atomname, atype;

	try
	{
		mass		= tpr->get_xvf("m");
		charge		= tpr->get_xvf("q");
		resid		= tpr->get_ivector("resid");
		resname		= tpr->get_name("res");
		atomname	= tpr->get_name("atom");
		atype		= tpr->get_name("type");
	}
	catch (const std::exception&e)
	{
		fprintf(stderr, "tprplugin) ERROR: %s\n", e.what());
		return MOLFILE_ERROR;
	}

	for (size_t i = 0; i < resname.size(); i++)
	{
		auto atom = atoms + i; // to next pointer
		atom->mass		= mass[i];
		atom->charge	= charge[i];
		atom->resid		= resid[i] % 100000; // fix must be <= 99999 for gro format
		strcpy(atom->name, atomname[i].c_str());
		strcpy(atom->resname, resname[i].c_str());
		strcpy(atom->type, atype[i].c_str());
		atom->segid[0] = '\0';
		atom->chain[0] = '\0';
	}

	return MOLFILE_SUCCESS;
}

//< check if has velocity
static int read_tpr_timestep_metadata(void* mydata, molfile_timestep_metadata_t* meta) {

	TprPlugin* tprplugin = static_cast<TprPlugin*> (mydata);
	meta->count = -1;
	meta->has_velocities = 0;

	// determine if has velocity
	try
	{
		const auto &velocities = tprplugin->tpr->get_xvf("v");
		meta->has_velocities = !velocities.empty();
	}
	catch (...)
	{
		fprintf(stderr, "tprPlugin) No Velocities in tpr");
	}

	return MOLFILE_SUCCESS;
}

// read coords !
static int read_tpr_timestep(void* mydata, int natoms, molfile_timestep_t* ts)
{
	TprPlugin* tprplugin = static_cast<TprPlugin*>(mydata);

	if (tprplugin->has_read) return MOLFILE_ERROR;

	if (ts != NULL)
	{
		try
		{
			std::vector<float> coords;
			coords = tprplugin->tpr->get_xvf("x");
			for (size_t i = 0; i < 3L * natoms; i++)
			{
				ts->coords[i] = coords[i] * 10; // to A
			}
		}
		catch (const std::exception&e)
		{
			fprintf(stderr, "tprplugin) ERROR: %s\n", e.what());
			return MOLFILE_ERROR;
		}

		//velocity if exist
		try
		{
			std::vector<float> velocity;
			velocity = tprplugin->tpr->get_xvf("v");
			for (size_t i = 0; i < 3L * natoms; i++)
			{
				ts->velocities[i] = velocity[i] * 10; // to A/ps
			}
		}
		catch (const std::exception& e)
		{
			fprintf(stderr, "tprplugin) Warning: %s\n", e.what());
		}


		//box if exist
		try
		{
			std::vector<float> box = tprplugin->tpr->get_xvf("box");

			ts->A = sqrt(box[0] * box[0] + box[1] * box[1] + box[2] * box[2]) * 10;
			ts->B = sqrt(box[3] * box[3] + box[4] * box[4] + box[5] * box[5]) * 10;
			ts->C = sqrt(box[6] * box[6] + box[7] * box[7] + box[8] * box[8]) * 10;

			if (ts->A <= 0 || ts->B <= 0 || ts->C <= 0)
			{
				ts->A = ts->B = ts->C = 0;
				ts->alpha = ts->beta = ts->gamma = 0;
			}
			else
			{
				ts->gamma = acos((box[0] * box[3] + box[1] * box[4] + box[2] * box[5]) * 100 / (ts->A * ts->B)) * 90.0 / M_PI_2;
				ts->beta  = acos((box[0] * box[6] + box[1] * box[7] + box[2] * box[8]) * 100 / (ts->A * ts->C)) * 90.0 / M_PI_2;
				ts->alpha = acos((box[3] * box[6] + box[4] * box[7] + box[5] * box[8]) * 100 / (ts->B * ts->C)) * 90.0 / M_PI_2;
			}
		}
		catch (const std::exception& e)
		{
			fprintf(stderr, "tprplugin) Warning: %s\n", e.what());
		}
	}

	tprplugin->has_read = true;
	return MOLFILE_SUCCESS;
}

static int read_tpr_bonds(void* mydata, int* nbonds, int** from, int** to, float** bondorder,
	int** bondtype, int* nbondtypes, char*** bondtypename)
{
	*nbonds = 0;
	*from = NULL;
	*to = NULL;
	*bondorder = NULL;
	*bondtype = NULL;
	*nbondtypes = 0;
	*bondtypename = NULL;

	TprPlugin* tprplugin = static_cast<TprPlugin*>(mydata);
	try
	{
		const std::vector<Bonded> &bonds = tprplugin->tpr->get_bonded("bonds");
		int nb = static_cast<int>(bonds.size());
		*nbonds = nb;
		*from = tprplugin->from = emalloc<int>(nb);
		*to = tprplugin->to = emalloc<int>(nb);
		for (int i = 0; i < nb; i++)
		{
			(*from)[i]     = bonds[i].a;
			(*to)[i]       = bonds[i].b;
		}
	}
	catch (...)
	{
		fprintf(stderr, "tprplugin) No bonds found!\n");
	}

	return MOLFILE_SUCCESS;
}

// all angles
static int read_tpr_angles(void* mydata, int* numangles, int** angles, int** angletypes,
	int* numangletypes, char*** angletypenames, int* numdihedrals,
	int** dihedrals, int** dihedraltypes, int* numdihedraltypes,
	char*** dihedraltypenames, int* numimpropers, int** impropers,
	int** impropertypes, int* numimpropertypes, char*** impropertypenames,
	int* numcterms, int** cterms, int* ctermcols, int* ctermrows)
{
	/* initialize data to zero */
	*numangles = 0;
	*angles = NULL;
	*angletypes = NULL;
	*numangletypes = 0;
	*angletypenames = NULL;
	*numdihedrals = 0;
	*dihedrals = NULL;
	*dihedraltypes = NULL;
	*numdihedraltypes = 0;
	*dihedraltypenames = NULL;
	*numimpropers = 0;
	*impropers = NULL;
	*impropertypes = NULL;
	*numimpropertypes = 0;
	*impropertypenames = NULL;
	*numcterms = 0;
	*cterms = NULL;
	*ctermrows = 0;
	*ctermcols = 0;

	// angles
	TprPlugin* tprplugin = static_cast<TprPlugin*>(mydata);
	try
	{
		const std::vector<Bonded>& ang = tprplugin->tpr->get_bonded("angles");
		int nA = static_cast<int>(ang.size());
		*numangles = nA;
		*angles = tprplugin->angles = emalloc<int>(nA * 3L);

		for (int i = 0; i < nA; i++)
		{
			(*angles)[3L * i]     = ang[i].a;
			(*angles)[3L * i + 1] = ang[i].b;
			(*angles)[3L * i + 2] = ang[i].c;
		}
	}
	catch (...)
	{
		fprintf(stderr, "tprplugin) No angles found!\n");
	}

	// dihedrals
	try
	{
		const std::vector<Bonded>& dih = tprplugin->tpr->get_bonded("dihedrals");
		int nD = static_cast<int>(dih.size());
		*numdihedrals = nD;
		*dihedrals = tprplugin->dihedrals = emalloc<int>(nD * 4L);
		for (int i = 0; i < nD; i++)
		{
			(*dihedrals)[4L * i]     = dih[i].a;
			(*dihedrals)[4L * i + 1] = dih[i].b;
			(*dihedrals)[4L * i + 2] = dih[i].c;
			(*dihedrals)[4L * i + 3] = dih[i].d;
		}
	}
	catch (...)
	{
		fprintf(stderr, "tprplugin) No dihedrals found!\n");
	}


	// impropers
	try
	{
		const std::vector<Bonded>& improper = tprplugin->tpr->get_bonded("impropers");
		int nImp = static_cast<int>(improper.size());
		*numimpropers = nImp;
		*impropers = tprplugin->impropers = emalloc<int>(nImp * 4L);
		for (int i = 0; i < nImp; i++)
		{
			(*impropers)[4L * i]     = improper[i].a;
			(*impropers)[4L * i + 1] = improper[i].b;
			(*impropers)[4L * i + 2] = improper[i].c;
			(*impropers)[4L * i + 3] = improper[i].d;
		}
	}
	catch (...)
	{
		fprintf(stderr, "tprplugin) No impropers found!\n");
	}

	return MOLFILE_SUCCESS;
}


static molfile_plugin_t tpr_plugin;
VMDPLUGIN_API int VMDPLUGIN_init()
{
	//< gromacs tpr plugin
	memset(&tpr_plugin, 0, sizeof(molfile_plugin_t));
	tpr_plugin.abiversion = vmdplugin_ABIVERSION;
	tpr_plugin.type = MOLFILE_PLUGIN_TYPE;
	tpr_plugin.name = "tpr";
	tpr_plugin.prettyname = "Gromacs TPR format"; // open name
	tpr_plugin.author = "Yujie Liu";
	tpr_plugin.majorv = 1;
	tpr_plugin.minorv = 1;
	tpr_plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
	tpr_plugin.filename_extension = "tpr";
	tpr_plugin.open_file_read = open_tpr_read;
	tpr_plugin.read_structure = read_tpr_structure;
	tpr_plugin.read_timestep_metadata = read_tpr_timestep_metadata;
	tpr_plugin.read_next_timestep = read_tpr_timestep;
	tpr_plugin.read_bonds = read_tpr_bonds;
	tpr_plugin.read_angles = read_tpr_angles;
	tpr_plugin.close_file_read = close_tpr_read;

	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void* v, vmdplugin_register_cb cb)
{
	(*cb)(v, (vmdplugin_t*)&tpr_plugin);
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLIGIN_fini()
{
	return VMDPLUGIN_SUCCESS;
}

#ifdef __cplusplus
}
#endif // _cplusplus

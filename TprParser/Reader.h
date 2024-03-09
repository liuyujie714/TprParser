#ifndef READER_H
#define	 READER_H

#include "define.h"
#include "Bytes.h"
#include <stdio.h>
#include <vector>
#include <array>

using vecI2D = std::vector<std::vector<int>>;
using vecF2D = std::vector<std::vector<float>>;
using vecU2D = std::vector<std::vector<unsigned short>>;

enum class PbcType : int
{
	Xyz = 0, //!< Periodic boundaries in all dimensions.
	No = 1, //!< No periodic boundaries.
	XY = 2, //!< Only two dimensions are periodic.
	Screw = 3, //!< Screw.
	Unset = 4, //!< The type of PBC is not set or invalid.
	Count = 5,
	Default = Xyz
};

struct TprData
{
	// clear zero
	TprData() : symtab(nullptr), lambda(0.0f)
	{
		prec = filever = vergen = natoms = ngtc = fep_state = 0;
		symtablen = nmoltypes = nmolblock = 0;
		bIr = bTop = bX = bV = bF = bBox = bInter = false;
	}

	int					prec; //< the precision of tpr, 4 or 8
	int					filever; //< the version of file format, fver
	int					vergen; //< the verions of generation code, fgen
	int					natoms; //< the total natoms
	int					ngtc; //< The number of temperature coupling groups.
	int					fep_state; //< fep state
	float				lambda; //< lambda
	bool				bIr; //< if has ir 
	bool				bTop; //< if has top 
	bool				bX; //< if has coordinates 
	bool				bV; //< if has velocity
	bool				bF; //< if has force
	bool				bBox; //< if has box 
	bool				bInter; //< if has inter-molecular bonds
	float				box[DIM * DIM] = {0}; //< box size
	char				* symtab;//< symb name, truncate to 8 characters
	int					symtablen, nmoltypes, nmolblock;

	std::vector<int>	atomsinmol;
	std::vector<int>	resinmol;
	std::vector<int>	molnames;
	std::vector<int>	molbtype;
	std::vector<int>	molbnmol;
	std::vector<int>	molbnatoms; // 每个单分子有多少个原子构成
	vecF2D				charges;
	vecF2D				masses;
	vecI2D				resids;
	std::vector<int>	trueresids; // actually residues number in tpr
	vecI2D				ptypes;
	vecU2D				types;
	vecI2D				atomnameids;
	vecI2D				atomtypeids;
	vecI2D				resnames;
	vecI2D				atomicnumbers;

	// mdp parameters
	struct IR
	{
		PbcType				pbc = PbcType::Unset; //< which pbc type
		bool				pbcmol; //< periodic-molecules
		int64_t				nsteps; // the number of simulation steps 
		int64_t				init_step; // simulation init steps
		int					simulation_part;
		int					nstcalcenergy;
		int					cutoff_scheme; // int to enum
		int					nstlist;
		int					nstcomm;
		int					comm_mode; // int to enum, 0=Linear, 1=Angular
		int					nstcgsteep; // Number of steps after which a steepest descents step is done while doing cg
		int					nbfgscorr; // Number of corrections to the Hessian to keep
		int					nstlog; // number of log steps
		int					nstxout; // number of trr coordinates steps
		int					nstvout; // number of velocity steps
		int					nstfout; // number of force steps
		int					nstenergy; // number of energy output steps
		int					nstxout_compressed; // number of xtc coordinates steps
		double				init_t = 0.0; // init time, ps
		double				dt = 0.0; // time steps, ps

		float				x_compression_precision; /// precision of xtc coordinates
		float				verletbuf_tol; // tolerance of verlet buffer
		float				verletBufferPressureTolerance;
		float				rlist;
		int					coulombtype; // int to enum, 0=Cut, 1=RF, 3=Pme
		int					coulomb_modifier; // int to enum, 0=PotShiftVerletUnsupported, 1=PotShift, 2=None
		float				rcoulomb_switch;
		float				rcoulomb;
		int					vdwtype;// int to enum, 0=Cut, 1=Switch,2=Shift, ...
		int					vdw_modifier; // int to enum


		float				rvdw_switch;
		float				rvdw;
		int					eDispCorr;
		float				epsilon_r;
		float				epsilon_rf;
		float				tabext;

		bool				implicit_solvent = false; // if has implicit solvent

		float				fourier_spacing;
		int					nkx;
		int					nky;
		int					nkz;
		int					pme_order;
		float				ewald_rtol;
		float				ewald_rtol_lj;
		int					ewald_geometry; // int to enum, 0=3D, 1=3DC
		float				epsilon_surface; 
		int					ljpme_combination_rule; // int to enum, 0=Geom, 1=LB
		bool				bContinuation; 
		// int to enum, 0=No,1=Berendsen,2=NoseHoover,3=Yes,4=Andersen,5=AndersenMassive
		// 6=VRescale
		int					etc;

		int					nsttcouple;
		int					nstpcouple;
		int					epc; // int to enum
		int					epct; // int to enum
		float				tau_p;
		float				ref_p[DIM * DIM] = { 0 };
		float				compress[DIM * DIM] = { 0 };
		float				posres_com[DIM] = { 0 };
		float				posres_comB[DIM] = { 0 };
		int					refcoord_scaling; // int to enum,0=No,1=All,2=Com

		float				shake_tol; // tolerance of shake
		int					efep; // int to enum
		int					n_lambda = 0; //The number of foreign lambda points

		bool				bSimTemp;// if has simulation temperature
		int					eSimTempScale; // enum to int
		float				simtemp_high;
		float				simtemp_low;

		bool				bExpanded = false; // Whether expanded ensembles are used

		// em
		float				em_stepsize;
		float				em_tol;
		int64_t				ld_seed = 0;

		// deform
		float				deform[DIM * DIM] = { 0 };
		float				cos_accel = 0;
		int					userint1, userint2, userint3, userint4;
		float				userreal1, userreal2, userreal3, userreal4;
	} ir;


	struct
	{
		vecI2D			    interactionlist[F_NRE];
		std::vector<int>	nr[F_NRE];
	} 
	ilist,					// 分子相互作用列表
	inter_molecular_ilist;	// 全局指定的分子间相互作用

	// 原子属性
	struct
	{
		std::vector<float>			x; //< coordinates
		std::vector<float>			v; //< velocity
		std::vector<float>			f; //< force
		std::vector<std::string>	atomname;
		std::vector<std::string>	resname;
		std::vector<int>			resid;
		std::vector<float>			mass;
		std::vector<float>			charge;
	} atoms;

	// bonds
	std::vector<std::array<int, 2>> bonds;
	// angles
	std::vector<std::array<int, 3>>	angles;

	// mdp属性位置, 所有变量都必须初始化为0
	struct
	{
		long		nsteps = 0; //< the started nsteps position in tpr
		long		dt = 0; //< the started dt position in tpr
		long		x = 0; //< the started atom coordinates position in tpr
	} property;
};

class TprReader
{
public:
	const char* fout_ = nullptr;

public:
	// data_ 不能用memset清零含有模板类的结构体
	TprReader(
		const char *fname, 
		bool bGRO = false, 
		bool bMol2 = false, 
		bool bCharge = false
	) 
		: tpr_(fname, "rb"), data_(new TprData), fout_("new.tpr"), 
		bGRO_(bGRO), bMol2_(bMol2), bCharge_(bCharge)
	{
		if (tpr_header() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_header()");
		}
		if (tpr_body() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_body()");
		}
		if (tpr_mtop() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_mtop()");
		}
		if (tpr_xvf() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_xvf()");
		}
		if (tpr_chargemass() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_chargemass()");
		}
		if (tpr_bonds() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_bonds()");
		}
		if (tpr_angles() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for tpr_angles()");
		}
		if (do_ir() != TPR_SUCCESS)
		{
			throw std::runtime_error("error for do_ir()");
		}
	}

	~TprReader()
	{
		if (data_->symtab) delete [] data_->symtab;
		if (data_) delete data_;

		msg("End of TprReader\n");
	}

	// read header
	bool tpr_header();

	// read body of tpr
	bool tpr_body();

	// read mtop
	bool tpr_mtop();

	//< read coodinates, velocity and force of atoms
	bool tpr_xvf();

	//< dump charges and mass
	bool tpr_chargemass();

	//< dump bonds of tpr
	bool tpr_bonds();

	//< dump angles of tpr
	bool tpr_angles();

	//< do_ir
	bool do_ir();

public:
	//< change tpr file nsteps
	void set_nsteps(int64_t nsteps);

	//< change tpr file dt (ps)
	void set_dt(double dt);

	//< change tpr atomic coordinates
	template<typename T>
	void set_coordinates(std::vector<T> &coords)
	{
		// check if has coordinates of tpr
		if (!data_->bX)
		{
			throw std::runtime_error("Input tpr has not coordinates information");
		}

		// check vector size 
		if (coords.size() != data_->natoms * DIM)
		{
			throw std::runtime_error("Input vector size is not equal to natoms * 3");
		}

		// check data type float or double, must be same as data_->prec
		if (sizeof(T) != data_->prec)
		{
			throw std::runtime_error("Input data type is not equal to data_->precision: " + std::to_string(data_->prec));
		}

		FileSerializer  newtpr(fout_, "wb");
		long            fsize = 0;
		const char* buffer = tpr_.get_file_buffer(&fsize);
		// 原始位置不为0
		if (fsize && data_->property.x)
		{
			// write nsteps before 
			if (newtpr.fwrite_(buffer, data_->property.x * sizeof(char), 1) != 1)
			{
				throw std::runtime_error("fwrite_ error in set_coordinates before");
			}

			// write new coordinates
			newtpr.do_vector(coords.data(), (int)coords.size(), data_->prec);

			// write coordinates after
			size_t size = coords.size() * sizeof(T);
			long len = fsize - data_->property.x - (long)size;
			if (newtpr.fwrite_(&buffer[data_->property.x + size], len * sizeof(char), 1) != 1)
			{
				throw std::runtime_error("fwrite_ error in set_coordinates after");
			}
		}
	}

private:
	//< read forcefield parameters
	bool tpr_readff();

	//< read parameters
	bool do_iparams(int ftype, t_iparams * iparams, int filever, int prec);

	//< moltype dump
	bool do_atoms();

	//< atomtype
	bool do_atomtypes();

	//< read cmap
	bool do_cmap();

	//< read groups
	bool do_groups();

	//< do_ilists
	bool do_ilists(int ntype, std::vector<int>(&nr)[F_NRE], vecI2D(&interactionlist)[F_NRE]);

	//< do_fepvals
	bool do_fepvals();

private:
	FileSerializer			tpr_;
	TprData					*data_;
	std::vector<t_iparams>  iparams_; // 力场参数
	bool					bGRO_ = false; //< if write a gro
	bool					bMol2_ = false; //< if write a mol2 whith bonds
	bool					bCharge_ = false; //< if write atomic charge and mass to file
};

#endif // !READER_H

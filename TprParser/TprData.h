#ifndef TPRDATA_H
#define TPRDATA_H

#include "define.h"
#include "Enum.h"
#include <vector>
#include <string>
#include <array>

using vecI2D = std::vector<std::vector<int>>;
using vecF2D = std::vector<std::vector<float>>;
using vecU2D = std::vector<std::vector<unsigned short>>;

//< set up box rel
static inline void do_box_rel(int ndim, const float deform[DIM][DIM], float box_rel[DIM][DIM], float b[DIM][DIM], bool bInit)
{
	for (int d = YY; d <= ZZ; ++d)
	{
		for (int d2 = XX; d2 < ndim; ++d2)
		{
			/* We need to check if this box component is deformed
			 * or if deformation of another component might cause
			 * changes in this component due to box corrections.
			 */
			if (deform[d][d2] == 0
				&& !(d == ZZ && d2 == XX && deform[d][YY] != 0 && (b[YY][d2] != 0 || deform[YY][d2] != 0)))
			{
				if (bInit)
				{
					box_rel[d][d2] = b[d][d2] / b[XX][XX];
				}
				else
				{
					b[d][d2] = b[XX][XX] * box_rel[d][d2];
				}
			}
		}
	}
}


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
	float				box[DIM * DIM] = { 0 }; //< box size
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
	struct
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
		int					epc; // Pressure coupling algorithm，int to enum, 0=No, 1=Berendsen, 2=ParrinelloRahman,5=CRescale
		int					epct; //Pressure coupling type, int to enum, 0=Isotropic, 1=SemiIsotropic
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


		int					ngacc = 0; // 加速组个数
		int					ngfrz = 0; // 冻结组个数
		int					ngener = 0; // 能量组个数

		// 控温部分，维度ngtc
		int					ngtc = 0; // 控温组数目
		int					nhchainlength = 1;
		std::vector<float>  nrdf; // 每个组自由度
		std::vector<float>  ref_t; // 每个组参考温度
		std::vector<int>    annealing; // 每个组模拟退火类型, enum to int
		std::vector<int>    anneal_npoints; // 每个组模拟退火点数
		vecF2D				anneal_time; // 每个组模拟退火时间点
		vecF2D				anneal_temp; // 每个组模拟退火温度点
		std::vector<float>	tau_t;  // 每个组模拟退火时间常数

		std::vector<std::array<int, DIM>> nFreeze; // 每个组在三个方向时候被冻结，维度ngfrz
		std::vector<std::array<float, DIM>> acceleration; // 每个组在三个方向时候被冻结，维度ngacc
		std::vector<int>	egp_flags; // 能量组每对之间的Exclusions/tables，维度ngener*ngener

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
		// 压力设置参数位置
		struct
		{
			long	box_rel = 0; //< the started vector position for preserve box shape, is DIM*DIM vector
			long	epc = 0; //< the started pressure coupling method position
			long	epct = 0; //< the started pressure coupling type position
			long	tau_p = 0; //< the started tau_p position
			long	ref_p = 0; //< the started ref pressure value position, is DIM*DIM vector
			long	compress = 0; //< the started compressibility value position, is DIM*DIM vector

			//< return False if get all parameters
			bool empty() const
			{
				return !(box_rel && epc && epct && tau_p && ref_p && compress);
			}
		} press;
		
		// 温度设置参数位置
		struct 
		{
			long	g_ngtc = 0; //< the started number of temperature coupling group position in tpr header, g_ngtc==ir->ngtc
			long	etc = 0; //< the started temperature coupling type position, enum to int, 0=No, 1=Berendsen,2=NoseHoover,6=VRescale
			long	ngtc = 0; //< the started number of temperature coupling group position
			long    nhchainlength = 0; //< the Nose-Hoover chain length if use Nose-Hoover temperature coupling
			long	ref_t = 0; //< the started ref temperature position, is ir->ngtc vector
			long	tau_t = 0; //< the started temperature coupling constant position, is ir->ngtc vector

			//< return True if can not read temperature position due to pull code, AWH have not yet finish in tpr reader (TODO)
			bool empty() const
			{
				return !(ref_t && tau_t && etc && ngtc && g_ngtc);
			}
		} temperature;
		
		// 单个整数属性mdp设置位置
		struct
		{
			long			nstlog = 0; // started 
			long			nstxout = 0;
			long			nstvout = 0;
			long			nstfout = 0;
			long			nstenergy = 0;
			long			nstxout_compressed = 0;

			// must be data_->filever >= 71
			long			nsttcouple = 0;
			long			nstpcouple = 0;

			//must be data_->filever >= 67
			long			nstcalcenergy = 0;

			//< return True if can not read any one position
			bool empty() const
			{
				return !(nstlog && nstxout && nstvout && nstfout && nstenergy && nstxout_compressed && nsttcouple && nstpcouple && nstcalcenergy);
			}
		} integer;
	} property;
};


#endif // !TPRDATA_H

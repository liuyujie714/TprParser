#ifndef READER_H
#define	 READER_H

#include "define.h"
#include "Bytes.h"
#include <stdio.h>
#include <vector>

using vecI2D = std::vector<std::vector<int>>;
using vecF2D = std::vector<std::vector<float>>;
using vecU2D = std::vector<std::vector<unsigned short>>;

// angles
struct t_angle {
	t_angle(int a_, int b_, int c_) : a(a_), b(b_), c(c_) {}
	int a, b, c;
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

	// 分子相互作用列表
	struct
	{
		vecI2D			    interactionlist[F_NRE];
		std::vector<int>	nr[F_NRE];
	} ilist;
	// 全局指定的分子间相互作用
	struct
	{
		vecI2D				interactionlist[F_NRE];
		std::vector<int>	nr[F_NRE];
	} inter_molecular_ilist;

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
	std::vector<std::pair<int, int>> bonds;
	// angles
	std::vector<t_angle>			angles;
};

class TprReader
{
public:
	TprReader(const char *fname) : tpr_(fname, "rb")
	{
		if (tpr_header() != TPR_SUCCESS)
		{
			msg("error for tpr_header()\n");
		}
		if (tpr_body() != TPR_SUCCESS)
		{
			msg("error for tpr_body()\n");
		}
		if (tpr_mtop() != TPR_SUCCESS)
		{
			msg("error for tpr_mtop()\n");
		}
		if (tpr_xvf() != TPR_SUCCESS)
		{
			msg("error for tpr_xvf()\n");
		}
		if (tpr_chargemass() != TPR_SUCCESS)
		{
			msg("error for tpr_chargemass()\n");
		}
		if (tpr_bonds() != TPR_SUCCESS)
		{
			msg("error for tpr_bonds()\n");
		}
		if (tpr_angles() != TPR_SUCCESS)
		{
			msg("error for tpr_angles()\n");
		}
	}

	~TprReader()
	{
		if (data_->symtab) delete [] data_->symtab;
		if (data_) delete data_;

		msg("End of reading\n");
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

private:
	TpxSerializer			tpr_;
	TprData					*data_;
	std::vector<t_iparams>  iparams_; // 力场参数
};


#endif // !READER_H

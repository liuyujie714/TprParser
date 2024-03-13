#ifndef READER_H
#define	 READER_H

#include "define.h"
#include "Bytes.h"
#include "TprData.h"
#include <string.h>

#define INSERT_POS(prop) data_->property.prop = tpr_.ftell_()

class TprReader
{
public:
	const char* fout_ = nullptr;

public:
	// data_ 不能用memset清零含有模板类的结构体
	TprReader(
		const char* fname,
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
		if (data_->symtab) delete[] data_->symtab;
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
	bool set_nsteps(int64_t nsteps);

	//< change tpr file dt (ps)
	bool set_dt(double dt);

	//< change tpr atomic coordinates
	template<typename T>
	bool set_coordinates(std::vector<T>& coords)
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

		long            fsize = 0;
		const char* buffer = tpr_.get_file_buffer(&fsize);
		// 原始位置不为0
		if (fsize && data_->property.x)
		{
			FileSerializer  newtpr(fout_, "wb");

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

			return TPR_SUCCESS;
		}

		return TPR_FAILED;
	}

	//< set pressure coupling parts
	bool set_pressure(
		const char* method,
		const char* type,
		float tau_p,
		std::vector<float>& ref_p,
		std::vector<float>& compress);

	//< set temperature coupling parts. Have not yet set groups name
	bool set_temperature(
		const char* method,
		std::vector<float>& tau_t,
		std::vector<float>& ref_t);

	//< get coords
	const std::vector<float> &get_coordinates() const
	{
		// check if has coordinates of tpr
		if (!data_->bX)
		{
			throw std::runtime_error("Input tpr has not coordinates information");
		}

		return data_->atoms.x;
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

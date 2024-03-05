#include "Reader.h"
#include <set>
#include <algorithm>

bool TprReader::tpr_header()
{
	// read the first int at the first of tpr
	int tempint;
	if (!tpr_.do_int(&tempint)) return TPR_FAILED;

	// read string contains gmx version
	char filever[MAX_LEN];
	if (!tpr_.xdr_string(filever, MAX_LEN)) return TPR_FAILED;
	msg("gmx version: %s\n", filever);

	// read precision int
	if (!tpr_.do_int(&data_->prec)) return TPR_FAILED;
	msg("gmx precision: %s\n", data_->prec == sizeof(float) ? "float" : "double");
	if (data_->prec != sizeof(float) && data_->prec != sizeof(double))
	{
        throw std::runtime_error("TpxSerializer unsupports precision: " + std::to_string(data_->prec));
	}
	
	return TPR_SUCCESS;
}

static void print_box(const char *name, float *arr)
{
#ifdef DEBUG
    msg(name);
    for (int i = 0; i < 9; i++) printf("%f ", arr[i]);
    puts("");
#endif // DEBUG
}

bool TprReader::tpr_body()
{
	// read file foramt version of tpr
	if (!tpr_.do_int(&data_->filever)) return TPR_FAILED;
	msg("File Format Version: %d\n", data_->filever);
	/* This is for backward compatibility with development versions 77-79
	 * where the tag was, mistakenly, placed before the generation,
	 * which would cause a segv instead of a proper error message
	 * when reading the topology only from tpx with <77 code.
	 */
	if (data_->filever >= 77 && data_->filever <= 79)
	{
		char release[MAX_LEN];
		if (!tpr_.xdr_string(release, MAX_LEN)) return TPR_FAILED;
		msg("%s\n", release);
	}
	if (!tpr_.do_int(&data_->vergen)) return TPR_FAILED;
	msg("file generator: %d\n", data_->vergen);

	// string ?
	if (data_->filever >= 81)
	{
		char buf[MAX_LEN];
		int tempint;
		// 前4个字节未使用
		if (!tpr_.do_int(&tempint)) return TPR_FAILED;

		if (!tpr_.xdr_string(buf, MAX_LEN)) return TPR_FAILED;
		msg("%s\n", buf);
	}

	// natoms and ngtc
	if (!tpr_.do_int(&data_->natoms)) return TPR_FAILED;
	if (!tpr_.do_int(&data_->ngtc)) return TPR_FAILED;
	msg("natoms= %d, ngtc= %d\n", data_->natoms, data_->ngtc);

	// fep state and lambda
	if (data_->filever < 62)
	{
		int tempint;
		float tempreal;
		if (!tpr_.do_int(&tempint)) return TPR_FAILED;
		if (!tpr_.do_real(&tempreal, data_->prec)) return TPR_FAILED;
	}
	if (data_->filever >= 79)
	{
		// fep state
		if (!tpr_.do_int(&data_->fep_state)) return TPR_FAILED;
		msg("fep_state= %d\n", data_->fep_state);
	}
	// lambda 
	if (!tpr_.do_real(&data_->lambda, data_->prec)) return TPR_FAILED;
	msg("lambda= %f\n", data_->lambda);
	// bool type
	if (!tpr_.do_bool(&data_->bIr)) return TPR_FAILED;
	if (!tpr_.do_bool(&data_->bTop)) return TPR_FAILED;
	if (!tpr_.do_bool(&data_->bX)) return TPR_FAILED;
	if (!tpr_.do_bool(&data_->bV)) return TPR_FAILED;
	if (!tpr_.do_bool(&data_->bF)) return TPR_FAILED;
	if (!tpr_.do_bool(&data_->bBox)) return TPR_FAILED;
	msg("bIr= %d, bTop= %d, bX= %d, bV= %d, bF= %d, bBox= %d\n",
		data_->bIr ? 1 : 0,
		data_->bTop ? 1 : 0,
		data_->bX ? 1 : 0,
		data_->bV ? 1 : 0,
		data_->bF ? 1 : 0,
		data_->bBox ? 1 : 0
	);

	if (data_->filever >= tpxv_AddSizeField && data_->vergen >= 27)
	{
		int64_t fsize;
		if (!tpr_.do_int64(&fsize)) return TPR_FAILED;
		msg("Size of tpr body= %lld bytes\n", fsize);
	}
	if (data_->vergen > 28) data_->bIr = false; // This can only happen if TopOnlyOK=TRUE

	// read box size
	if (data_->bBox)
	{
		tpr_.do_vector(data_->box, 9, data_->prec);
        print_box("box= ", data_->box);

		// Relative box vectors characteristic of the box shape, used to to preserve that box shape
		if (data_->filever >= 51)
		{
			float box_rel[9] = { 0 };
			tpr_.do_vector(box_rel, 9, data_->prec);
            print_box("box_rel= ", box_rel);
		}

		// Box velocities for Parrinello-Rahman P-coupling
		float boxv[9] = { 0 };
		tpr_.do_vector(boxv, 9, data_->prec);
        print_box("boxv= ", boxv);

		if (data_->filever < 56)
		{
			float dump[9] = { 0 };
			tpr_.do_vector(dump, 9, data_->prec);
		}
	}

	// 温度耦合组
	if (data_->ngtc > 0)
	{
		float* temparr = new float[data_->ngtc];
		if (data_->filever < 69)
		{
			if (!tpr_.do_vector(temparr, data_->ngtc, data_->prec)) return TPR_FAILED;
		}
		//These used to be the Berendsen tcoupl_lambda's
		if(!tpr_.do_vector(temparr, data_->ngtc, data_->prec)) return TPR_FAILED;
		delete temparr;
	}

	return TPR_SUCCESS;
}

bool TprReader::tpr_mtop()
{
	//do_mtop starts here, which starts by reading the symtab (do_symtab)
	if (!tpr_.do_int(&data_->symtablen)) return TPR_FAILED;
    msg("symtablen= %d\n", data_->symtablen);
	data_->symtab = new char[SAVELEN * data_->symtablen];
	// clear data
	memset(data_->symtab, 0, sizeof(char) * SAVELEN * data_->symtablen);

	// 原子类型名称和组名
	for (int i = 0; i < data_->symtablen; i++)
	{
		tpr_.tpr_save_string(&data_->symtab[SAVELEN * i], data_->vergen);

		// print symb
		if constexpr (0)
		{
			const char* start = &data_->symtab[SAVELEN * i];
			while (*start)
			{
				printf("%c", *start++);
			}
			printf("\n");
		}
	}

	// temp int
	int tempint;
	if (!tpr_.do_int(&tempint)) return TPR_FAILED;
	msg("tempint= %d\n", tempint);

	// read forcefiled parameters
	if (!tpr_readff()) return TPR_FAILED;

    // read type of molecules
    if (!tpr_.do_int(&data_->nmoltypes)) return TPR_FAILED;
    msg("nmoltypes= %d\n", data_->nmoltypes);
    if (!do_atoms()) return TPR_FAILED;

    
    // 保存分子和原子信息到atoms结构体中
    data_->atoms.atomname.resize(data_->natoms);
    data_->atoms.resname.resize(data_->natoms);
    data_->atoms.resid.resize(data_->natoms);
    data_->atoms.mass.resize(data_->natoms);
    data_->atoms.charge.resize(data_->natoms);
    unsigned int idx = 0;
    int startedresindex = 1;
    for (int i = 0; i < data_->nmolblock; i++)
    {
        int m = data_->molbtype[i];
        for (int j = 0; j < data_->molbnmol[i]; j++)
        {
            std::set<int> residx;
            for (int k = 0; k < data_->molbnatoms[i]; k++)
            {
                int resind = data_->resids[m][k];
                data_->atoms.atomname[idx]  = &data_->symtab[SAVELEN * data_->atomnameids[m][k]];
                data_->atoms.resname[idx]   = &data_->symtab[SAVELEN * data_->resnames[m][resind]];

                // 此处的残基编号有问题，当tpr中不连续时候处理不了
                residx.insert(resind);
                data_->atoms.resid[idx]     = resind + startedresindex;

                data_->atoms.mass[idx]      = data_->masses[m][k];
                data_->atoms.charge[idx]    = data_->charges[m][k];
                idx++;
            }
            startedresindex += static_cast<int>(residx.size());
        }
    }

	return TPR_SUCCESS;
}

bool TprReader::tpr_xvf()
{
    if (data_->bX)
    {
        data_->atoms.x.resize(data_->natoms * DIM);// 3N 
        if (!tpr_.do_vector(data_->atoms.x.data(), data_->natoms * DIM, data_->prec)) return TPR_FAILED;
    }
    if (data_->bV)
    {
        data_->atoms.v.resize(data_->natoms * DIM);// 3N 
        if (!tpr_.do_vector(data_->atoms.v.data(), data_->natoms * DIM, data_->prec)) return TPR_FAILED;
    }
    if (data_->bF)
    {
        data_->atoms.f.resize(data_->natoms * DIM);// 3N 
        if (!tpr_.do_vector(data_->atoms.f.data(), data_->natoms * DIM, data_->prec)) return TPR_FAILED;
    }

    // write a gro
    if (data_->bX)
    {
        FILE* fp = fopen("dump.gro", "w");

        fprintf(fp, "MOL\n%d\n", data_->natoms);
        for (int i = 0; i < data_->natoms; i++)
        {
            fprintf(fp, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f",
                data_->atoms.resid[i] % 100000, data_->atoms.resname[i].c_str(), data_->atoms.atomname[i].c_str(), (i + 1) % 100000,
                data_->atoms.x[3L * i + 0], data_->atoms.x[3L * i + 1], data_->atoms.x[3L * i + 2]);

            if (data_->bV)
            {
                fprintf(fp, "%8.4f%8.4f%8.4f", data_->atoms.v[3L * i + 0], data_->atoms.v[3L * i + 1], data_->atoms.v[3L * i + 2]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
            data_->box[0], data_->box[4], data_->box[8],
            data_->box[1], data_->box[2], data_->box[3],
            data_->box[5], data_->box[6], data_->box[7]);
        fclose(fp);
    }

    return TPR_SUCCESS;
}

bool TprReader::tpr_chargemass()
{
    FILE* fp = fopen("chgmass.dat", "w");
    for (int i = 0; i < data_->atoms.atomname.size(); i++)
    {
        fprintf(fp, "%5s %10.6f %10.6f\n", data_->atoms.atomname[i].c_str(),
            data_->atoms.charge[i], data_->atoms.mass[i]);
    }
    fclose(fp);

    return TPR_SUCCESS;
}

bool TprReader::tpr_bonds()
{
    // bonds type
    const int interactions[] = {
        F_BONDS, F_G96BONDS, F_MORSE, F_CUBICBONDS, F_CONNBONDS, F_HARMONIC, F_FENEBONDS,
        F_CONSTR, F_CONSTRNC, F_TABBONDS, F_TABBONDSNC, F_SETTLE
    };
    constexpr int nBonds = asize(interactions);

    int aoffset = 0;
    for (int i = 0; i < data_->nmolblock; i++)
    {
        int mtype = data_->molbtype[i];
        for (int j = 0; j < data_->molbnmol[i]; j++)
        {
            for (int k = 0; k < nBonds; k++)
            {
                int type = interactions[k];
                // settle 
                if (type == F_SETTLE)
                {
                    // settle algorithm for water molecules
                    for (int m = 0; m < data_->ilist.nr[type][mtype] / 2; m++)
                    {
                        // OW-HW1
                        data_->bonds.push_back(std::make_pair<int, int>(1 + aoffset, 2 + aoffset));
                        // OW-HW2
                        data_->bonds.push_back(std::make_pair<int, int>(1 + aoffset, 3 + aoffset));
                    }
                }
                else
                {
                    for (int m = 0; m < data_->ilist.nr[type][mtype] / 3; m++)
                    {
                        int a = 3 * m + 1;
                        int b = 3 * m + 2;
                        data_->bonds.push_back(std::make_pair<int, int>(
                            1 + data_->ilist.interactionlist[type][mtype][a] + aoffset,
                            1 + data_->ilist.interactionlist[type][mtype][b] + aoffset)
                        );
                    }
                }
            }
            aoffset += data_->atomsinmol[mtype];
        }
    }

    // inter-molecular bonds
    if (data_->bInter)
    {
        // use global atom index
        for (int m = 0; m < data_->inter_molecular_ilist.nr[F_HARMONIC][0] / 3; m++)
        {
            int a = 3 * m + 1;
            int b = 3 * m + 2;
            data_->bonds.push_back(std::make_pair<int, int>(
                1 + data_->inter_molecular_ilist.interactionlist[F_HARMONIC][0][a],
                1 + data_->inter_molecular_ilist.interactionlist[F_HARMONIC][0][b])
            );
        }
    }


    // write a mol2 format
    FILE* fp = fopen("dump.mol2", "w");
    fprintf(fp, "@<TRIPOS>MOLECULE\nMOL\n%d %d 1 0 0\nSMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", data_->natoms, (int)(data_->bonds.size()));
    for (int i = 0; i < data_->natoms; i++)
    {
        fprintf(fp, "%3d %5s %8.4f %8.4f %8.4f %c %5d %5s %8.4f\n",
            i + 1, data_->atoms.atomname[i].c_str(),
            data_->atoms.x[3L * i + 0] * 10.0, data_->atoms.x[3L * i + 1] * 10.0, data_->atoms.x[3L * i + 2] * 10.0,
            data_->atoms.atomname[i].c_str()[0], 1, data_->atoms.resname[i].c_str(), data_->atoms.charge[i]);
    }
    fprintf(fp, "@<TRIPOS>BOND\n");
    int i = 1;
    for (auto& bond : data_->bonds)
    {
        fprintf(fp, "%5d %5d %5d %5d\n", i++, bond.first, bond.second, 1);
    }
    fclose(fp);

    return TPR_SUCCESS;
}

bool TprReader::tpr_angles()
{
    // angles type
    const int interactions[] = {
        F_ANGLES, F_G96ANGLES, F_CROSS_BOND_BONDS, F_CROSS_BOND_ANGLES, F_UREY_BRADLEY,
        F_QUARTIC_ANGLES, F_RESTRANGLES, F_TABANGLES, F_SETTLE
    };
    constexpr int nAngles = asize(interactions);

    int aoffset = 0;
    for (int i = 0; i < data_->nmolblock; i++)
    {
        int mtype = data_->molbtype[i];
        for (int j = 0; j < data_->molbnmol[i]; j++)
        {
            for (int k = 0; k < nAngles; k++)
            {
                int type = interactions[k];
                // settle 
                if (type == F_SETTLE)
                {
                    // settle algorithm for water molecules
                    for (int m = 0; m < data_->ilist.nr[type][mtype] / 4; m++)
                    {
                        // HW1 - OW - HW2
                        data_->angles.push_back(t_angle(2 + aoffset, 1 + aoffset, 3 + aoffset));
                    }
                }
                else
                {
                    for (int m = 0; m < data_->ilist.nr[type][mtype] / 4; m++)
                    {
                        int a = 4 * m + 1;
                        int b = 4 * m + 2;
                        int c = 4 * m + 3;
                        data_->angles.push_back(
                            t_angle(
                                1 + data_->ilist.interactionlist[type][mtype][a] + aoffset,
                                1 + data_->ilist.interactionlist[type][mtype][b] + aoffset,
                                1 + data_->ilist.interactionlist[type][mtype][c] + aoffset
                            )
                        );
                    }
                }
            }
            aoffset += data_->atomsinmol[mtype];
        }
    }

    // TODO
    // 1. dihedrals, impdihedral...
    // 2. inter-molecular angles, dihedrals, imp...

    return TPR_SUCCESS;
}

bool TprReader::tpr_readff()
{
    int		            atnr, ntypes;
	double              reppow = 12.0;
	float	            fudge = 0.5;
    std::vector<int>    functype;

	if (!tpr_.do_int(&atnr)) return TPR_FAILED;
	if (!tpr_.do_int(&ntypes)) return TPR_FAILED;
	msg("ntypes= %d\n", ntypes);

	// 函数类型
    functype.resize(ntypes);
	for (int i = 0; i < ntypes; i++)
	{
		if (!tpr_.do_int(&functype[i])) return TPR_FAILED;
	}
	if (data_->filever >= 66) if (!tpr_.do_double(&reppow)) return TPR_FAILED;
	if (!tpr_.do_real(&fudge, data_->prec)) return TPR_FAILED;
	msg("fudge= %f\n", fudge);
	
	// 调整所有函数类型
    iparams_.resize(ntypes);
	for (int i = 0; i < ntypes; i++)
	{
		for (int j = 0; j < NFTUPD; j++)
		{
			if (data_->filever < ftupd[j].fvnr && functype[i] >= ftupd[j].ftype)
			{
				functype[i] += 1;
			}
		}
		// 读力场参数
		if (!do_iparams(functype[i], &iparams_[i], data_->filever, data_->prec)) return TPR_FAILED;
	}

	return TPR_SUCCESS;
}

bool TprReader::do_iparams(int ftype, t_iparams* iparams, int filever, int prec)
{
    int         idum;
    float       rdum;

    switch (ftype)
    {
    case F_ANGLES:
    case F_G96ANGLES:
    case F_BONDS:
    case F_G96BONDS:
    case F_HARMONIC:
    case F_IDIHS:
        tpr_.do_real(&iparams->harmonic.rA, prec);
        tpr_.do_real(&iparams->harmonic.krA, prec);
        tpr_.do_real(&iparams->harmonic.rB, prec);
        tpr_.do_real(&iparams->harmonic.krB, prec);
        if ((ftype == F_ANGRES || ftype == F_ANGRESZ))
        {
            /* Correct incorrect storage of parameters */
            iparams->pdihs.phiB = iparams->pdihs.phiA;
            iparams->pdihs.cpB = iparams->pdihs.cpA;
        }
        break;
    case F_RESTRANGLES:
        tpr_.do_real(&iparams->harmonic.rA, prec);
        tpr_.do_real(&iparams->harmonic.krA, prec);
        break;
    case F_LINEAR_ANGLES:
        tpr_.do_real(&iparams->linangle.klinA, prec);
        tpr_.do_real(&iparams->linangle.aA, prec);
        tpr_.do_real(&iparams->linangle.klinB, prec);
        tpr_.do_real(&iparams->linangle.aB, prec);
        break;
    case F_FENEBONDS:
        tpr_.do_real(&iparams->fene.bm, prec);
        tpr_.do_real(&iparams->fene.kb, prec);
        break;

    case F_RESTRBONDS:
        tpr_.do_real(&iparams->restraint.lowA, prec);
        tpr_.do_real(&iparams->restraint.up1A, prec);
        tpr_.do_real(&iparams->restraint.up2A, prec);
        tpr_.do_real(&iparams->restraint.kA, prec);
        tpr_.do_real(&iparams->restraint.lowB, prec);
        tpr_.do_real(&iparams->restraint.up1B, prec);
        tpr_.do_real(&iparams->restraint.up2B, prec);
        tpr_.do_real(&iparams->restraint.kB, prec);
        break;
    case F_TABBONDS:
    case F_TABBONDSNC:
    case F_TABANGLES:
    case F_TABDIHS:
        tpr_.do_real(&iparams->tab.kA, prec);
        tpr_.do_int(&iparams->tab.table);
        tpr_.do_real(&iparams->tab.kB, prec);
        break;
    case F_CROSS_BOND_BONDS:
        tpr_.do_real(&iparams->cross_bb.r1e, prec);
        tpr_.do_real(&iparams->cross_bb.r2e, prec);
        tpr_.do_real(&iparams->cross_bb.krr, prec);
        break;
    case F_CROSS_BOND_ANGLES:
        tpr_.do_real(&iparams->cross_ba.r1e, prec);
        tpr_.do_real(&iparams->cross_ba.r2e, prec);
        tpr_.do_real(&iparams->cross_ba.r3e, prec);
        tpr_.do_real(&iparams->cross_ba.krt, prec);
        break;
    case F_UREY_BRADLEY:
        tpr_.do_real(&iparams->u_b.thetaA, prec);
        tpr_.do_real(&iparams->u_b.kthetaA, prec);
        tpr_.do_real(&iparams->u_b.r13A, prec);
        tpr_.do_real(&iparams->u_b.kUBA, prec);
        if (filever >= 79)
        {
            tpr_.do_real(&iparams->u_b.thetaB, prec);
            tpr_.do_real(&iparams->u_b.kthetaB, prec);
            tpr_.do_real(&iparams->u_b.r13B, prec);
            tpr_.do_real(&iparams->u_b.kUBB, prec);
        }
        else
        {
            iparams->u_b.thetaB = iparams->u_b.thetaA;
            iparams->u_b.kthetaB = iparams->u_b.kthetaA;
            iparams->u_b.r13B = iparams->u_b.r13A;
            iparams->u_b.kUBB = iparams->u_b.kUBA;
        }
        break;
    case F_QUARTIC_ANGLES:
        tpr_.do_real(&iparams->qangle.theta, prec);
        tpr_.do_vector(iparams->qangle.c, 5, data_->prec);
        break;
    case F_BHAM:
        tpr_.do_real(&iparams->bham.a, prec);
        tpr_.do_real(&iparams->bham.b, prec);
        tpr_.do_real(&iparams->bham.c, prec);
        break;
    case F_MORSE:
        tpr_.do_real(&iparams->morse.b0A, prec);
        tpr_.do_real(&iparams->morse.cbA, prec);
        tpr_.do_real(&iparams->morse.betaA, prec);
        if (filever >= 79)
        {
            tpr_.do_real(&iparams->morse.b0B, prec);
            tpr_.do_real(&iparams->morse.cbB, prec);
            tpr_.do_real(&iparams->morse.betaB, prec);
        }
        else
        {
            iparams->morse.b0B = iparams->morse.b0A;
            iparams->morse.cbB = iparams->morse.cbA;
            iparams->morse.betaB = iparams->morse.betaA;
        }
        break;
    case F_CUBICBONDS:
        tpr_.do_real(&iparams->cubic.b0, prec);
        tpr_.do_real(&iparams->cubic.kb, prec);
        tpr_.do_real(&iparams->cubic.kcub, prec);
        break;
    case F_CONNBONDS: break;
    case F_POLARIZATION: 
        tpr_.do_real(&iparams->polarize.alpha, prec); 
        break;
    case F_ANHARM_POL:
        tpr_.do_real(&iparams->anharm_polarize.alpha, prec);
        tpr_.do_real(&iparams->anharm_polarize.drcut, prec);
        tpr_.do_real(&iparams->anharm_polarize.khyp, prec);
        break;
    case F_WATER_POL:
        tpr_.do_real(&iparams->wpol.al_x, prec);
        tpr_.do_real(&iparams->wpol.al_y, prec);
        tpr_.do_real(&iparams->wpol.al_z, prec);
        tpr_.do_real(&iparams->wpol.rOH, prec);
        tpr_.do_real(&iparams->wpol.rHH, prec);
        tpr_.do_real(&iparams->wpol.rOD, prec);
        break;
    case F_THOLE_POL:
        tpr_.do_real(&iparams->thole.a, prec);
        tpr_.do_real(&iparams->thole.alpha1, prec);
        tpr_.do_real(&iparams->thole.alpha2, prec);
        if (filever < tpxv_RemoveTholeRfac)
        {
            float noRfac = 0;
            tpr_.do_real(&noRfac, prec);
        }

        break;
    case F_LJ:
        tpr_.do_real(&iparams->lj.c6, prec);
        tpr_.do_real(&iparams->lj.c12, prec);
        break;
    case F_LJ14:
        tpr_.do_real(&iparams->lj14.c6A, prec);
        tpr_.do_real(&iparams->lj14.c12A, prec);
        tpr_.do_real(&iparams->lj14.c6B, prec);
        tpr_.do_real(&iparams->lj14.c12B, prec);
        break;
    case F_LJC14_Q:
        tpr_.do_real(&iparams->ljc14.fqq, prec);
        tpr_.do_real(&iparams->ljc14.qi, prec);
        tpr_.do_real(&iparams->ljc14.qj, prec);
        tpr_.do_real(&iparams->ljc14.c6, prec);
        tpr_.do_real(&iparams->ljc14.c12, prec);
        break;
    case F_LJC_PAIRS_NB:
        tpr_.do_real(&iparams->ljcnb.qi, prec);
        tpr_.do_real(&iparams->ljcnb.qj, prec);
        tpr_.do_real(&iparams->ljcnb.c6, prec);
        tpr_.do_real(&iparams->ljcnb.c12, prec);
        break;
    case F_PDIHS:
    case F_PIDIHS:
    case F_ANGRES:
    case F_ANGRESZ:
        tpr_.do_real(&iparams->pdihs.phiA, prec);
        tpr_.do_real(&iparams->pdihs.cpA, prec);
        tpr_.do_real(&iparams->pdihs.phiB, prec);
        tpr_.do_real(&iparams->pdihs.cpB, prec);
        tpr_.do_int(&iparams->pdihs.mult);
        break;
    case F_RESTRDIHS:
        tpr_.do_real(&iparams->pdihs.phiA, prec);
        tpr_.do_real(&iparams->pdihs.cpA, prec);
        break;
    case F_DISRES:
        tpr_.do_int(&iparams->disres.label);
        tpr_.do_int(&iparams->disres.type);
        tpr_.do_real(&iparams->disres.low, prec);
        tpr_.do_real(&iparams->disres.up1, prec);
        tpr_.do_real(&iparams->disres.up2, prec);
        tpr_.do_real(&iparams->disres.kfac, prec);
        break;
    case F_ORIRES:
        tpr_.do_int(&iparams->orires.ex);
        tpr_.do_int(&iparams->orires.label);
        tpr_.do_int(&iparams->orires.power);
        tpr_.do_real(&iparams->orires.c, prec);
        tpr_.do_real(&iparams->orires.obs, prec);
        tpr_.do_real(&iparams->orires.kfac, prec);
        break;
    case F_DIHRES:
        if (filever < 82)
        {
            tpr_.do_int(&idum);
            tpr_.do_int(&idum);
        }
        tpr_.do_real(&iparams->dihres.phiA, prec);
        tpr_.do_real(&iparams->dihres.dphiA, prec);
        tpr_.do_real(&iparams->dihres.kfacA, prec);
        if (filever >= 82)
        {
            tpr_.do_real(&iparams->dihres.phiB, prec);
            tpr_.do_real(&iparams->dihres.dphiB, prec);
            tpr_.do_real(&iparams->dihres.kfacB, prec);
        }
        else
        {
            iparams->dihres.phiB = iparams->dihres.phiA;
            iparams->dihres.dphiB = iparams->dihres.dphiA;
            iparams->dihres.kfacB = iparams->dihres.kfacA;
        }
        break;
    case F_POSRES:
        tpr_.do_vector(iparams->posres.pos0A, DIM, data_->prec);
        tpr_.do_vector(iparams->posres.fcA, DIM, data_->prec);
        tpr_.do_vector(iparams->posres.pos0B, DIM, data_->prec);
        tpr_.do_vector(iparams->posres.fcB, DIM, data_->prec);
        break;
    case F_FBPOSRES:
        tpr_.do_int(&iparams->fbposres.geom);
        tpr_.do_vector(iparams->fbposres.pos0, DIM, data_->prec);
        tpr_.do_real(&iparams->fbposres.r, prec);
        tpr_.do_real(&iparams->fbposres.k, prec);
        break;
    case F_CBTDIHS: 
        tpr_.do_vector(iparams->cbtdihs.cbtcA, NR_CBTDIHS, data_->prec);
        break;
    case F_RBDIHS:
        // Fall-through intended
    case F_FOURDIHS:
        /* Fourier dihedrals are internally represented
         * as Ryckaert-Bellemans since those are faster to compute.
         */
        tpr_.do_vector(iparams->rbdihs.rbcA, NR_RBDIHS, data_->prec);
        tpr_.do_vector(iparams->rbdihs.rbcB, NR_RBDIHS, data_->prec);
        break;
    case F_CONSTR:
    case F_CONSTRNC:
        tpr_.do_real(&iparams->constr.dA, prec);
        tpr_.do_real(&iparams->constr.dB, prec);
        break;
    case F_SETTLE:
        tpr_.do_real(&iparams->settle.doh, prec);
        tpr_.do_real(&iparams->settle.dhh, prec);
        break;
    case F_VSITE1: 
        break; // VSite1 has 0 parameters
    case F_VSITE2:
    case F_VSITE2FD:
        tpr_.do_real(&iparams->vsite.a, prec); 
        break;
    case F_VSITE3:
    case F_VSITE3FD:
    case F_VSITE3FAD:
        tpr_.do_real(&iparams->vsite.a, prec);
        tpr_.do_real(&iparams->vsite.b, prec);
        break;
    case F_VSITE3OUT:
    case F_VSITE4FD:
    case F_VSITE4FDN:
        tpr_.do_real(&iparams->vsite.a, prec);
        tpr_.do_real(&iparams->vsite.b, prec);
        tpr_.do_real(&iparams->vsite.c, prec);
        break;
    case F_VSITEN:
        tpr_.do_int(&iparams->vsiten.n);
        tpr_.do_real(&iparams->vsiten.a, prec);
        break;
    case F_GB12_NOLONGERUSED:
    case F_GB13_NOLONGERUSED:
    case F_GB14_NOLONGERUSED:
        // Implicit solvent parameters can still be read, but never used
        //if (serializer->reading())
        {
            if (filever < 68)
            {
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
            }
            if (filever < tpxv_RemoveImplicitSolvation)
            {
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
                tpr_.do_real(&rdum, prec);
            }
        }
        break;
    case F_CMAP:
        tpr_.do_int(&iparams->cmap.cmapA);
        tpr_.do_int(&iparams->cmap.cmapB);
        break;
    default:
        msg("Unknown function type %d", ftype);
    }

	return TPR_SUCCESS;
}

bool TprReader::do_atoms()
{
    const int       n = data_->nmoltypes;
    float           rdum;
    int             idum, idum2;
    unsigned char   ucdum;
    unsigned short  usdum;

    data_->atomsinmol.resize(n);
    data_->molnames.resize(n);
    data_->resinmol.resize(n);
    data_->charges.resize(n);
    data_->masses.resize(n);
    data_->ptypes.resize(n);
    data_->types.resize(n);
    data_->resids.resize(n);
    data_->atomnameids.resize(n);
    data_->atomtypeids.resize(n);
    data_->atomicnumbers.resize(n);
    data_->resnames.resize(n);
    for (int i = 0; i < F_NRE; i++)
    {
        data_->ilist.interactionlist[i].resize(n);
        data_->ilist.nr[i].resize(n);
    }

    // read each mol
    for (int i = 0; i < n; i++)
    {
        // 分子名称长度
        if (!tpr_.do_int(&data_->molnames[i])) return TPR_FAILED;
        msg("data_->molnames[i]= %d\n", data_->molnames[i]);

        // 每个moltype的原子数目
        if (!tpr_.do_int(&data_->atomsinmol[i])) return TPR_FAILED;
        // 每个moltype的残基数目
        if (!tpr_.do_int(&data_->resinmol[i])) return TPR_FAILED;

        // allocate for 2D vector
        data_->charges[i].resize(data_->atomsinmol[i]);
        data_->masses[i].resize(data_->atomsinmol[i]);
        data_->types[i].resize(data_->atomsinmol[i]);
        data_->ptypes[i].resize(data_->atomsinmol[i]);
        data_->resids[i].resize(data_->atomsinmol[i]);
        data_->atomicnumbers[i].resize(data_->atomsinmol[i]);
        for (int j = 0; j < data_->atomsinmol[i]; j++)
        {
            // 读原子质量，电荷
            if (!tpr_.do_real(&data_->masses[i][j], data_->prec)) return TPR_FAILED;
            if (!tpr_.do_real(&data_->charges[i][j], data_->prec)) return TPR_FAILED;
            if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED; //mB
            if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED; //qB

            // for gmx>=2020，short只读2字节
            if (!tpr_.do_ushort(&data_->types[i][j], data_->vergen)) return TPR_FAILED;
            if (!tpr_.do_ushort(&usdum, data_->vergen)) return TPR_FAILED;

            if (!tpr_.do_int(&data_->ptypes[i][j])) return TPR_FAILED;
            if (!tpr_.do_int(&data_->resids[i][j])) return TPR_FAILED;
            if (data_->filever >= 52)
            {
                if (!tpr_.do_int(&data_->atomicnumbers[i][j])) return TPR_FAILED;
            }
            //msg("data_->types[i][j]= %d\n", data_->types[i][j]);
            //msg("data_->ptypes[i][j]= %d\n", data_->ptypes[i][j]);
            //msg("data_->resids[i][j]= %d\n", data_->resids[i][j]);
            //msg("data_->atomicnumbers[i][j]= %d\n", data_->atomicnumbers[i][j]);
        }

        // allocated for 2D vector
        data_->atomnameids[i].resize(data_->atomsinmol[i]);
        data_->atomtypeids[i].resize(data_->atomsinmol[i]);
        if (!tpr_.do_vector(data_->atomnameids[i].data(), data_->atomsinmol[i], data_->prec)) return TPR_FAILED;
        if (!tpr_.do_vector(data_->atomtypeids[i].data(), data_->atomsinmol[i], data_->prec)) return TPR_FAILED;
        // typeB
        for (int j = 0; j < data_->atomsinmol[i]; j++)
        {
            if (!tpr_.do_int(&idum)) return TPR_FAILED;
        }

        //read residues
        data_->resnames[i].resize(data_->resinmol[i]);
        for (int j = 0; j < data_->resinmol[i]; j++)
        {
            if (!tpr_.do_int(&data_->resnames[i][j])) return TPR_FAILED;
            //msg("data_->resnames[i][j]= %d\n", data_->resnames[i][j]);

            if (data_->filever >= 63)
            {
                // true Residue number 
                if (!tpr_.do_int(&idum)) return TPR_FAILED;
                data_->trueresids.push_back(idum);

                // gmx >= 2020, 只读一字节
                if (!tpr_.do_uchar(&ucdum, data_->vergen)) return TPR_FAILED;
            }
            else
            {
                data_->resnames[i][j] += 1;
            }
        }

        // do_ilists
        msg("do_ilists\n");
        if (!do_ilists(i, data_->ilist.nr, data_->ilist.interactionlist)) return TPR_FAILED;

        // charge groups parts
        if (!tpr_.do_int(&idum)) return TPR_FAILED;
        std::vector<int> temp(idum + 1); // need +1
        if (!tpr_.do_vector(temp.data(), idum + 1, data_->prec)) return TPR_FAILED;
        // doListOfLists
        if (!tpr_.do_int(&idum)) return TPR_FAILED;
        if (!tpr_.do_int(&idum2)) return TPR_FAILED;
        temp.resize(idum + 1); // need +1
        if (!tpr_.do_vector(temp.data(), idum + 1, data_->prec)) return TPR_FAILED;
        temp.resize(idum2); // not need +1
        if (!tpr_.do_vector(temp.data(), idum2, data_->prec)) return TPR_FAILED;
    } 


    // do molblock
    if (!tpr_.do_int(&data_->nmolblock)) return TPR_FAILED;
    msg("nmolblock= %d\n", data_->nmolblock);
    data_->molbtype.resize(data_->nmolblock);
    data_->molbnmol.resize(data_->nmolblock);
    data_->molbnatoms.resize(data_->nmolblock);
    for (int i = 0; i < data_->nmolblock; i++)
    {
        if (!tpr_.do_int(&data_->molbtype[i])) return TPR_FAILED;
        if (!tpr_.do_int(&data_->molbnmol[i])) return TPR_FAILED;
        if (!tpr_.do_int(&data_->molbnatoms[i])) return TPR_FAILED;
        msg("data_->molbtype[i]= %d\n", data_->molbtype[i]);
        msg("data_->molbnmol[i]= %d\n", data_->molbnmol[i]);
        msg("data_->molbnatoms[i]= %d\n", data_->molbnatoms[i]);

        // posres
        if (!tpr_.do_int(&idum)) return TPR_FAILED;  //posres_xA
        msg("posres_xA= %d\n", idum);
        if (idum > 0)
        {
            std::vector<float> temp(idum * DIM);
            if (!tpr_.do_vector(temp.data(), idum * DIM, data_->prec)) return TPR_FAILED;
        }

        if (!tpr_.do_int(&idum)) return TPR_FAILED;  //posres_xB
        msg("posres_xB= %d\n", idum);
        if (idum > 0)
        {
            std::vector<float> temp(idum * DIM);
            if (!tpr_.do_vector(temp.data(), idum * DIM, data_->prec)) return TPR_FAILED;
        }
    }
    // 体系全局原子数
    if (!tpr_.do_int(&idum)) return TPR_FAILED;
    msg("The number of atoms= %d\n", idum);
    
    // inter-molecularbonds
    if (data_->filever >= tpxv_IntermolecularBondeds)
    {
        // for gmx>=2020, read 1 byte
        if (!tpr_.do_bool(&data_->bInter, data_->vergen))  return TPR_FAILED;
        msg("bInter= %d\n", data_->bInter ? 1 : 0);
        if (data_->bInter)
        {
            // allocated
            for (int i = 0; i < F_NRE; i++)
            {
                data_->inter_molecular_ilist.interactionlist[i].resize(1);
                data_->inter_molecular_ilist.nr[i].resize(1);
            }
            do_ilists(0, data_->inter_molecular_ilist.nr, data_->inter_molecular_ilist.interactionlist);
        }
    }

    if (data_->filever < tpxv_RemoveAtomtypes)
    {
        if (!do_atomtypes()) return TPR_FAILED;
    }

    if (data_->filever >= 65)
    {
        if (!do_cmap()) return TPR_FAILED;
    }

    if (!do_groups()) return TPR_FAILED;

    if (data_->filever >= tpxv_StoreNonBondedInteractionExclusionGroup)
    {
        int64_t intermolecularExclusionGroupSize;
        if (!tpr_.do_int64(&intermolecularExclusionGroupSize)) return TPR_FAILED;
        std::vector<int> temp(intermolecularExclusionGroupSize);
        if (!tpr_.do_vector(temp.data(), static_cast<int>(intermolecularExclusionGroupSize))) return TPR_FAILED;
    }

    return TPR_SUCCESS;
}

bool TprReader::do_atomtypes()
{
    int     nr;

    if (!tpr_.do_int(&nr)) return TPR_FAILED;
    if (data_->filever < tpxv_RemoveImplicitSolvation)
    {
        std::vector<float> temp(nr);
        if (!tpr_.do_vector(temp.data(), nr, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_vector(temp.data(), nr, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_vector(temp.data(), nr, data_->prec)) return TPR_FAILED;
    }
    std::vector<int> atomnumbers(nr);
    if (!tpr_.do_vector(atomnumbers.data(), nr, data_->prec)) return TPR_FAILED;

    if (data_->filever >= 60 && data_->filever < tpxv_RemoveImplicitSolvation)
    {
        std::vector<float> temp(nr);
        if (!tpr_.do_vector(temp.data(), nr, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_vector(temp.data(), nr, data_->prec)) return TPR_FAILED;
    }

    return TPR_SUCCESS;
}

bool TprReader::do_cmap()
{
    int         ngrid, gridspace;
    float       rdum;

    if (!tpr_.do_int(&ngrid)) return TPR_FAILED;
    if (!tpr_.do_int(&gridspace)) return TPR_FAILED;
    msg("ngrid= %d, gridspace= %d\n", ngrid, gridspace);

    for (int i = 0; i < ngrid * gridspace * gridspace; i++)
    {
        if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED;
        if (!tpr_.do_real(&rdum, data_->prec)) return TPR_FAILED;
    }

    return TPR_SUCCESS;
}

bool TprReader::do_groups()
{
    //do_grps
    int     idum;
    for (int i = 0; i < egcNR; i++)
    {
        if (!tpr_.do_int(&idum)) return TPR_FAILED;

        std::vector<int> temp(idum);
        if (!tpr_.do_vector(temp.data(), idum, data_->prec)) return TPR_FAILED;
    }

    if (!tpr_.do_int(&idum)) return TPR_FAILED;
    msg("number of group names= %d\n", idum);
    for (int i = 0; i < idum; i++)
    {
        int j;
        if (!tpr_.do_int(&j)) return TPR_FAILED;
    }

    for (int i = 0; i < egcNR; i++)
    {
        if (!tpr_.do_int(&idum)) return TPR_FAILED;
        if (idum != 0)
        {
            // for gmx >= 2020，uchar只读1字节
            std::vector<unsigned char> temp(idum);
            tpr_.do_vector(temp.data(), idum, data_->prec, data_->vergen);
        }
    }

    return TPR_SUCCESS;
}

bool TprReader::do_ilists(int ntype, std::vector<int> (&nr)[F_NRE], vecI2D (&interactionlist)[F_NRE])
{
    for (int i = 0; i < F_NRE; i++)
    {
        bool bClear = false;
        for (int k = 0; k < NFTUPD; k++)
        {
            if ((data_->filever < ftupd[k].fvnr) && (i == ftupd[k].ftype))
            {
                bClear = true;
            }
        }

        if (bClear)
        {
            nr[i][ntype] = 0;
            interactionlist[i][ntype].clear();
        }
        else
        {
            // get the number of interactions
            if (!tpr_.do_int(&nr[i][ntype])) return TPR_FAILED;
            if (nr[i][ntype] == 0) continue; // empty

            // allocated memory
            interactionlist[i][ntype].resize(nr[i][ntype]);
            if (!tpr_.do_vector(interactionlist[i][ntype].data(), nr[i][ntype], data_->prec)) return TPR_FAILED;

            if constexpr (0)
            {
                msg("%d ", nr[i][ntype]);
                for (int j = 0; j < nr[i][ntype]; j++)
                {
                    printf("%d ", interactionlist[i][ntype][j]);
                }
                printf("\n");
            }
        }
    }

    return TPR_SUCCESS;
}

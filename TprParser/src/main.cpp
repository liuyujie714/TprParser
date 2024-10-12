/*
* update: 
* 2024.06.24 - support set up deform 
* 2024.07.07 - support read electric field
* 2024.09.30 - fix empty lj parameters
*/

#include <iostream>
#include "Reader.h"

int main(int argc, char *argv[])
{
	try
	{
		//TprReader reader("test/em.tpr");
		//TprReader reader("test/large_2021_aa_posres.tpr");
		//TprReader reader("test/double_2023.tpr");
		//TprReader reader("test/benchMEM.tpr");
		//TprReader reader("test/md_cg.tpr");
		//TprReader reader("test/2lyz_gmx_4.0.tpr");
		//TprReader reader("test/inter-md.tpr");
		//TprReader reader("test/md2024.tpr");
		//TprReader reader("test/2023demo.tpr");
		//TprReader reader("test/2022.tpr");
		//TprReader reader("test/2lyz_gmx_2021.tpr");
		//TprReader reader("test/2020.4_gra.tpr");
		//TprReader reader("test/annealing.tpr");
		TprReader reader("test/npt2025-dev.tpr");

		// 电场测试
		{
			////TprReader reader("test/inter-md.tpr");
			////TprReader reader("test/elec.tpr");
			//TprReader reader("test/elecxyz.tpr");
			////TprReader reader("test/elec5.1.2.tpr");
			//auto &ef = reader.get_ef();
			//std::cout << "Electric field:\n";
			//for (auto& e : ef) {
			//	std::cout << e << "\t";
			//}

			//// write new tpr
			//std::vector<float> efnew{
			//	2, 8, 6, 4,
			//	0, 0.225, 0, 0,
			//	5, 0, 0, 0.11
			//};
			//reader.set_xvf("ef", efnew);
		}

		// 修改MD总步数(模拟时长)
		{
			//TprReader reader("test/md.tpr");
			//reader.set_nsteps(10000000);
		}

		// 修改MD步长（dt, ps）
		{
			//TprReader reader("test/md.tpr");
			//reader.set_dt(0.001); // 0.001 ps
		}

		// 修改原子坐标
		{
			//TprReader reader("test/oneWater.tpr");
			//std::vector<float> newCoords = {
			//	1.901f,   2.939f,   2.362f,
			//	2.000f,   2.928f,   2.353f,
			//	1.861f,   2.853f,   2.394f
			//};
			//reader.set_coordinates(newCoords);
		}

		// 分开控温
		{
			//TprReader reader("test/nvt_WAT.tpr");
			//std::vector<float> ref_t = { 200.0f };
			//std::vector<float> tau_t = { 2.0f };
			//reader.set_temperature("NoseHoover", tau_t, ref_t);
		}

		{
			//TprReader reader("test/Inter-2019.6.tpr");
		}

		{
			//TprReader reader("test/CO2_LineAngle.tpr", false, false);
			//auto ret = reader.get_nonbonded("pairs");
		}

		{
			//TprReader reader("test/2lyz_gmx_2021_double.tpr");
		}

		// 无LJ参数测试
		{
			//TprReader reader("test/extra-interactions-2018.tpr");
			//auto &ret = reader.get_nonbonded("lj");
			//std::cout << ret[0].ifunc << std::endl;
		}
	}
	catch (const std::exception&e)
	{
		fprintf(stderr, "%s\n", e.what());
	}

	return 0;
}


#include "Reader.h"

int main(int argc, char *argv[])
{
	try
	{
		TprReader reader1("test/em.tpr");
		//TprReader reader1("test/large_2021_aa_posres.tpr");
		//TprReader reader2("test/double_2023.tpr");
		//TprReader reader3("test/benchMEM.tpr");
		//TprReader reader4("test/md_cg.tpr");
		//TprReader reader5("test/2lyz_gmx_4.0.tpr");
		//TprReader reader6("test/inter-md.tpr");
		//TprReader reader7("test/md2024.tpr");
		//TprReader reader8("test/2023demo.tpr");
		//TprReader reader9("test/2022.tpr");
		//TprReader reader10("test/2lyz_gmx_2021.tpr");
		//TprReader reader11("test/2020.4_gra.tpr");
		//TprReader reader12("test/2lyz_gmx_2021_double.tpr");
	}
	catch (const std::exception&e)
	{
		msg("%s\n", e.what());
	}

	return 0;
}


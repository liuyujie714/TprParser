#include "Reader.h"

int main(int argc, char *argv[])
{
	try
	{
		//TprReader reader("inter-md.tpr");
		//TprReader reader(argv[1]);
		TprReader reader("md2024.tpr");
		//TprReader reader("2023demo.tpr");
		//TprReader reader("2022.tpr");
		//TprReader reader("2lyz_gmx_2021.tpr");
		//TprReader reader("2020.4_gra.tpr");
		//TprReader reader2("2020.4_gra.tpr");
	}
	catch (const std::exception&e)
	{
		printf("%s\n", e.what());
	}

	return 0;
}


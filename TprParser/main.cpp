#include "Reader.h"

int main(int argc, char *argv[])
{
	try
	{
		//TprReader reader("test/inter-md.tpr");
		TprReader reader1("test/md2024.tpr");
		//TprReader reader2("test/2023demo.tpr");
		//TprReader reader3("test/2022.tpr");
		//TprReader reader4("test/2lyz_gmx_2021.tpr");
		//TprReader reader5("test/2020.4_gra.tpr");
		//TprReader reader6("test/2020.4_gra.tpr");
	}
	catch (const std::exception&e)
	{
		msg("%s\n", e.what());
	}

	return 0;
}


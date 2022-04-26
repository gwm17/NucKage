#include "Simulator.h"
#include "Utils/Timer.h"
#include "Tests/tests.h"

int main(int argc, char** argv)
{
	if(false) //turn on/off testing
	{
		NucKage::EnergyLossTest();
		return 0;
	}
	
	int nthreads;
	std::string role;
	if(argc < 2)
	{
		return 1;
	}
	else if (argc == 2)
	{
		nthreads = 1;
		role = argv[1];
	}
	else if (argc == 3)
	{
		nthreads = std::stoi(argv[1]);
		role = argv[2];
	}
	else
	{
		return 1;
	}

	
	
	NucKage::Simulator* sim = NucKage::CreateSimulator(nthreads);
	sim->LoadConfig(role);
	NucKage::Timer stopwatch("WholeProgram");
	sim->Run();
	std::cout<<"Program duration: "<<stopwatch.ElapsedMilliseconds()<<" ms"<<std::endl;

	return 0;
}
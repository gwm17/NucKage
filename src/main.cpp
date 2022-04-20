#include "Simulator.h"
#include "Utils/Timer.h"
#include "Tests/tests.h"

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		return 1;
	}

	if(false) //turn on/off testing
	{
		NucKage::EnergyLossTest();
		return 0;
	}
	
	NucKage::Simulator* sim = NucKage::CreateSimulator(3);
	sim->LoadConfig(argv[1]);
	NucKage::Timer stopwatch("WholeProgram");
	sim->Run();
	std::cout<<"Program duration: "<<stopwatch.ElapsedMilliseconds()<<" ms"<<std::endl;

	return 0;
}
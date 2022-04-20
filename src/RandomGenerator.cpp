#include "RandomGenerator.h"

namespace NucKage {

	RandomGenerator::RandomGenerator() 
	{
		std::random_device rd;
		rng.seed(rd());
	}

	RandomGenerator::~RandomGenerator() {}
}
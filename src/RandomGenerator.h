#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>
#include <thread>
#include <mutex>

namespace NucKage {

	class RandomGenerator {
	public:
		RandomGenerator();
		~RandomGenerator();
		
		inline std::mt19937& GetGenerator() { return rng; }
		inline static RandomGenerator& GetInstance()
		{
			thread_local RandomGenerator s_generator; //implicitly static
			return s_generator;
		}

	private:
		std::mt19937 rng;
	};

}

#endif
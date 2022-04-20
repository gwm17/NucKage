#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <string>
#include "ReactorChain.h"
#include "Detectors/DetectorArray.h"
#include "ThreadPool.h"
#include "RootPlotter.h"

namespace NucKage {

	class Simulator
	{
	public:
		Simulator();
		Simulator(int nthreads);
		~Simulator();

		void LoadConfig(const std::string& filename);
		void Run();
		void GeneratePlots();

		inline static Simulator& GetInstance() { return *s_instance; }

	private:
		void RunChain(int index);
		static Simulator* s_instance;

		std::string m_outputFile;
		std::atomic<uint64_t> m_samples;
		bool m_initFlag;

		std::vector<ReactorChain> m_chains;
		std::vector<ChainResult> m_results;
		DetectorArray m_array;
		RootPlotter m_plotter;

		ThreadPool m_pool;
	};

	Simulator* CreateSimulator(int nthreads);
}

#endif
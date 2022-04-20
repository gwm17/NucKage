#ifndef REACTOR_CHAIN_H
#define REACTOR_CHAIN_H

#include "Reactor.h"
#include "RandomGenerator.h"
#include "EnergyLoss/Target.h"

namespace NucKage {

	struct SamplingParameters
	{
		double meanBeamKE=0.0;
		double sigmaBeamKE=0.0;
		double meanEx=0.0;
		double sigmaEx=0.0;
	};

	struct ChainResult
	{
		int chainID=-1;
		std::vector<ReactorProducts> products;
	};

	class ReactorChain
	{
	public:
		ReactorChain();
		~ReactorChain();
		void AddReactor(const std::vector<int>& Z, const std::vector<int>& A, const SamplingParameters& params);
		void SetTarget(const std::vector<int>& ZT, const std::vector<int>& stoich, double thickness);
		void BindTarget();
		bool VerifyChain();
		inline const int GetChainID() const { return m_result.chainID; }
		ChainResult& GenerateProducts();

	private:
		static int s_globalChainID;
		std::vector<Reactor> m_reactors;
		ChainResult m_result;
		Target m_target;

		std::uniform_real_distribution<double> m_thetaDist;
		std::uniform_real_distribution<double> m_phiDist;
		std::uniform_real_distribution<double> m_targetDist;

		std::vector<std::normal_distribution<double>> m_beamDists;
		std::vector<std::normal_distribution<double>> m_exDists;
	};
}

#endif
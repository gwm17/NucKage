#include "ReactorChain.h"
#include <iostream>

namespace NucKage {

	int ReactorChain::s_globalChainID = 0;
	ReactorChain::ReactorChain() :
		m_thetaDist(-1.0, 1.0), m_phiDist(0.0, 2.0*M_PI), m_targetDist(0.0, 1.0)
	{
		++s_globalChainID;
		m_result.chainID = s_globalChainID;
	}

	ReactorChain::~ReactorChain() {}

	void ReactorChain::AddReactor(const std::vector<int>& Z, const std::vector<int>& A, const SamplingParameters& params)
	{
		m_reactors.emplace_back(Z, A);
		m_result.products.push_back(ReactorProducts());
		m_beamDists.emplace_back(params.meanBeamKE, params.sigmaBeamKE);
		m_exDists.emplace_back(params.meanEx, params.sigmaEx);
	}

	void ReactorChain::SetTarget(const std::vector<int>& ZT, const std::vector<int>& stoich, double thickness)
	{
		m_target.SetParameters(ZT, stoich, thickness);
		BindTarget();
	}

	void ReactorChain::BindTarget()
	{
		for(auto& reactor : m_reactors)
			reactor.BindTarget(&m_target);
	}

	bool ReactorChain::VerifyChain()
	{
		if(m_reactors.size() == 0)
			return false;
		else if(m_reactors.size() == 1)
			return true;

		bool result = true;
		std::cout<<m_reactors[0].GetEquation();
		for(size_t i=1; i<m_reactors.size(); i++)
		{
			std::cout<<"->"<<m_reactors[i].GetEquation();
			const Nucleus& prevResidual = m_reactors[i-1].GetResidual();
			const Nucleus& currentTarget = m_reactors[i].GetTarget();
			if(currentTarget.Z != prevResidual.Z || currentTarget.A != prevResidual.A)
			{
				result = false;
				break;
			}
		}
		std::cout<<std::endl;
		return result;
	}

	//assumes all steps short enough that occur at roughly same reaction location
	ChainResult& ReactorChain::GenerateProducts()
	{
		ReactionParameters params;
		RandomGenerator& generator = RandomGenerator::GetInstance();
		params.targetFraction = m_targetDist(generator.GetGenerator()); //determine location of rxn
		for(size_t i=0; i<m_reactors.size(); i++)
		{
			i != 0 ? m_reactors[i].SetTarget4Vector(m_reactors[i-1].GetResidual().pvector) : void();
			switch(m_reactors[i].GetType())
			{
				case Reactor::Type::Reaction:
				{
					params.beamEnergy = m_beamDists[i](generator.GetGenerator());
					params.excitationEnergy = m_exDists[i](generator.GetGenerator());
					params.thetaCM = std::acos(m_thetaDist(generator.GetGenerator()));
					params.phiCM = m_phiDist(generator.GetGenerator());
					break;
				}
				case Reactor::Type::Decay:
				{
					params.excitationEnergy = m_exDists[i](generator.GetGenerator());
					params.thetaCM = std::acos(m_thetaDist(generator.GetGenerator()));
					params.phiCM = m_phiDist(generator.GetGenerator());
					break;
				}
				case Reactor::Type::None:
				{
					std::cerr<<"ERR -- Reactor of None type encountered in ReactorChain::GenerateProducts()"<<std::endl;
					return m_result;
				}
			}
			m_result.products[i] = m_reactors[i].GenerateProducts(params);
		}

		return m_result;
	}
}
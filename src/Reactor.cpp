#include "Reactor.h"
#include <iostream>

namespace NucKage {

	Reactor::Reactor(const std::vector<int>& Z, const std::vector<int>& A) :
		m_target(nullptr)
	{
		if(Z.size() == 2 && A.size() == 2)
			m_type = Type::Decay;
		else if(Z.size() == 3 && A.size() == 3)
			m_type = Type::Reaction;
		else
			m_type = Type::None;

		if(m_type == Type::Decay)
		{
			for(size_t i=0; i<Z.size(); i++)
			{
				m_reactants.emplace_back(Z[i], A[i]);
			}
			int AR = A[0] - A[1];
			int ZR = Z[0] - Z[1];
			m_reactants.emplace_back(ZR, AR);
			if(ZR < 0 || AR < 1 || ZR > AR)
			{
				m_reactants.clear();
				m_type = Type::None;
			}			
		}
		else if (m_type == Type::Reaction)
		{
			for(size_t i=0; i<Z.size(); i++)
			{
				m_reactants.emplace_back(Z[i], A[i]);
			}
			int AR = A[0] + A[1] - A[2];
			int ZR = Z[0] + Z[1] - Z[2];
			m_reactants.emplace_back(ZR, AR);
			if(ZR < 0 || AR < 1 || ZR > AR)
			{
				m_reactants.clear();
				m_type = Type::None;
			}
		}
	}

	Reactor::~Reactor() {}

	const std::string Reactor::GetEquation() const
	{
		switch(m_type)
		{
			case Type::Decay:
				return m_reactants[0].symbol + "->" + m_reactants[1].symbol + "+" + m_reactants[2].symbol;
			case Type::Reaction:
				return m_reactants[0].symbol + "(" + m_reactants[1].symbol + "," + m_reactants[2].symbol + ")" + m_reactants[3].symbol;
			case Type::None:
				return "None";
		}
		return "None";
	}

	ReactorProducts Reactor::GenerateProducts(const ReactionParameters& params)
	{
		switch(m_type)
		{
			case Type::Decay: return CalculateDecay(params);
			case Type::Reaction: return CalculateReaction(params);
			case Type::None: return ReactorProducts();
		}

		return ReactorProducts();
	}

	ReactorProducts Reactor::CalculateReaction(const ReactionParameters& params)
	{
		ReactorProducts prods;
		prods.reactorName = GetEquation();
		double beamE, beamPz;
		if(m_target) //if target do energy loss
		{
			double beamKE = params.beamEnergy;
			beamKE -= m_target->GetEnergyLossFractionalDepth(m_reactants[1].Z, m_reactants[1].A, beamKE, 0.0, params.targetFraction);
			beamPz = std::sqrt(beamKE*(beamKE + 2.0*m_reactants[1].mass));
			beamE = beamKE + m_reactants[1].mass;
		}
		else
		{
			beamPz = std::sqrt(params.beamEnergy*(params.beamEnergy + 2.0 * m_reactants[1].mass));
			beamE = params.beamEnergy + m_reactants[1].mass;
		}
		m_reactants[1].pvector.SetPxPyPzE(0.,0.,beamPz,beamE);

	
	
		double Q = m_reactants[0].pvector.M() + m_reactants[1].mass - (m_reactants[2].mass + m_reactants[3].mass + params.excitationEnergy);
	
		double Ethresh = -Q*(m_reactants[2].mass+m_reactants[3].mass)/(m_reactants[2].mass + m_reactants[3].mass - m_reactants[1].mass);
		if(params.beamEnergy < Ethresh)
		{
			return prods;
		}
		
		auto parent = m_reactants[0].pvector + m_reactants[1].pvector;
		auto boost = parent.BoostVector();
		parent.Boost(-1.0*boost);
		double ejectE_cm = (std::pow(m_reactants[2].mass, 2.0) - std::pow(m_reactants[3].mass + params.excitationEnergy, 2.0) +
							std::pow(parent.E(),2.0))/(2.0*parent.E());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - std::pow(m_reactants[2].mass, 2.0));
		m_reactants[2].pvector.SetPxPyPzE(ejectP_cm*std::sin(params.thetaCM)*std::cos(params.phiCM),
										  ejectP_cm*std::sin(params.thetaCM)*std::sin(params.phiCM),
										  ejectP_cm*std::cos(params.thetaCM),
										  ejectE_cm);
		m_reactants[2].pvector.Boost(boost);
		m_reactants[3].pvector = m_reactants[0].pvector + m_reactants[1].pvector - m_reactants[2].pvector;

		if(m_target) // ejectile energy loss
		{
			double ejectKE = m_reactants[2].pvector.E() - m_reactants[2].pvector.M();
			double percent_depth = params.targetFraction;
			double ejectTheta = m_reactants[2].pvector.Theta();
			double ejectPhi = m_reactants[2].pvector.Phi();
			if(ejectTheta < M_PI/2.0) //forwards through target (other fraction)
				percent_depth = 1.0 - percent_depth;
			ejectKE -= m_target->GetEnergyLossFractionalDepth(m_reactants[2].Z, m_reactants[2].A, ejectKE, ejectTheta, params.targetFraction);
			double ejectP = std::sqrt(ejectKE*(ejectKE + 2.0*m_reactants[2].pvector.M()));
			double ejectE = ejectKE + m_reactants[2].pvector.M();
			m_reactants[2].pvector.SetPxPyPzE(ejectP*std::sin(ejectTheta)*std::cos(ejectPhi),
											  ejectP*std::sin(ejectTheta)*std::sin(ejectPhi),
											  ejectP*std::cos(ejectTheta),
											  ejectE);
		}

		prods.target = m_reactants[0];
		prods.projectile = m_reactants[1];
		prods.ejectile = m_reactants[2];
		prods.residual = m_reactants[3];

		return prods;
	}

	ReactorProducts Reactor::CalculateDecay(const ReactionParameters& params)
	{
		ReactorProducts prods;
		prods.reactorName = GetEquation();
		double Q = m_reactants[0].pvector.M() - m_reactants[1].mass - m_reactants[2].mass;
		if(Q < 0)
		{
			return prods;
		}
	
		auto boost = m_reactants[0].pvector.BoostVector();
		m_reactants[0].pvector.Boost(-1.0*boost);
		double ejectE_cm = (m_reactants[1].mass*m_reactants[1].mass - m_reactants[2].mass*m_reactants[2].mass + 
							m_reactants[0].pvector.E()*m_reactants[0].pvector.E())/(2.0*m_reactants[0].pvector.E());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - m_reactants[1].mass*m_reactants[1].mass);
	
		m_reactants[1].pvector.SetPxPyPzE(ejectP_cm*std::sin(params.thetaCM)*std::cos(params.phiCM),
										  ejectP_cm*std::sin(params.thetaCM)*std::sin(params.phiCM),
										  ejectP_cm*std::cos(params.thetaCM),
										  ejectE_cm);
	
		m_reactants[0].pvector.Boost(boost);
		m_reactants[1].pvector.Boost(boost);
	
		m_reactants[2].pvector = m_reactants[0].pvector - m_reactants[1].pvector;

		
		if(m_target) //ejectile energy loss
		{
			double ejectKE = m_reactants[1].pvector.E() - m_reactants[1].pvector.M();
			double percent_depth = params.targetFraction;
			double ejectTheta = m_reactants[1].pvector.Theta();
			double ejectPhi = m_reactants[1].pvector.Phi();
			if(ejectTheta < M_PI/2.0) //forwards through target (other fraction)
				percent_depth = 1.0 - percent_depth;
			ejectKE -= m_target->GetEnergyLossFractionalDepth(m_reactants[1].Z, m_reactants[1].A, ejectKE, ejectTheta, params.targetFraction);
			double ejectP = std::sqrt(ejectKE*(ejectKE + 2.0*m_reactants[1].pvector.M()));
			double ejectE = ejectKE + m_reactants[1].pvector.M();
			m_reactants[1].pvector.SetPxPyPzE(ejectP*std::sin(ejectTheta)*std::cos(ejectPhi),
											  ejectP*std::sin(ejectTheta)*std::sin(ejectPhi),
											  ejectP*std::cos(ejectTheta),
											  ejectE);
		}

		prods.target = m_reactants[0];
		prods.ejectile = m_reactants[1];
		prods.residual = m_reactants[2];

		return prods;
	}
}
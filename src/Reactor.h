#ifndef REACTOR_H
#define REACTOR_H

#include <vector>
#include <iostream>
#include "Nucleus.h"
#include "EnergyLoss/Target.h"

namespace NucKage {

	struct ReactionParameters
	{
		double beamEnergy; //MeV
		double thetaCM; //radians
		double phiCM; //radians
		double excitationEnergy;
		double targetFraction; //percentage from upstream to downstream (percent of target travelled by beam to get to rxn)
	};

	struct ReactorProducts
	{
		std::string reactorName;
		Nucleus target;
		Nucleus projectile;
		Nucleus ejectile;
		Nucleus residual;
	};

	class Reactor
	{
	public:
		enum class Type
		{
			Reaction,
			Decay,
			None
		};

		Reactor(const std::vector<int>& Z, const std::vector<int>& A);
		~Reactor();

		inline const Type GetType() const { return m_type; }
		inline void SetTarget4Vector(const TLorentzVector& vec) { if(m_reactants.size() > 0) m_reactants[0].pvector = vec; }
		inline void ResetTarget4Vector() { if(m_reactants.size() > 0) m_reactants[0].pvector.SetPxPyPzE(0.,0.,0.,m_reactants[0].mass); }
		inline const Nucleus& GetTarget() const { if(m_reactants.size() > 0) return m_reactants[0]; else return m_blank;}
		inline const Nucleus& GetResidual() const 
		{ 
			if(m_reactants.size() > 0)
			{
				if(m_type == Type::Reaction)
				{
					return m_reactants[3];
				}
				else if(m_type == Type::Decay)
					return m_reactants[2];
			} 
			return m_blank;
		}
		const std::string GetEquation() const;

		ReactorProducts GenerateProducts(const ReactionParameters& params);

		inline void BindTarget(Target* target) { m_target = target; }

	private:
		ReactorProducts CalculateReaction(const ReactionParameters& params);
		ReactorProducts CalculateDecay(const ReactionParameters& params);

		Type m_type;
		std::vector<Nucleus> m_reactants;
		Target* m_target; //Not owned by reactor! do not delete

		Nucleus m_blank;
	};
}

#endif
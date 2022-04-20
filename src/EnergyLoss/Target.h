/*

Target.h
A basic target unit for use in the Mask environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include "EnergyLoss.h"

namespace NucKage {

	class Target
	{

	public:
		Target();
	 	Target(const std::vector<int>& z, const std::vector<int>& stoich, double thick);
	 	~Target();

	 	void SetParameters(const std::vector<int>& z, const std::vector<int>& stoich, double thick);
	 	double GetEnergyLossTotal(int zp, int ap, double startEnergy, double angle);
	 	double GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double angle);
	 	double GetEnergyLossFractionalDepth(int zp, int ap, double startEnergy, double angle, double percent_depth);
	 	double GetReverseEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double angle, double percent_depth);

	 	inline const EnergyLoss::Parameters& GetParameters() const { return m_params; }
	 	inline const double GetTotalThickness() const { return m_totalThickness; }
	 	inline const bool IsValid() { return m_isValid; }
	
	private:
		EnergyLoss::Parameters m_params;
		double m_totalThickness;
		bool m_isValid;
	};

}

#endif 

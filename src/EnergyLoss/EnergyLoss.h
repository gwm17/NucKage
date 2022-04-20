/*

EnergyLoss.h
Code for calculating the energy loss of a charged, massive particle through an arbitrary medium.
Based on code written by D.W. Visser at Yale for the original SPANC. Uses energy loss calulations
described by Ziegler in various SRIM textbooks.

Written by G.W. McCann Aug. 2020

*/

#ifndef ENERGYLOSS_H
#define ENERGYLOSS_H

#include <vector>

namespace NucKage {

	namespace EnergyLoss {

		struct Parameters
		{
			int ZP;
			double massP;
			std::vector<int> ZT;
			std::vector<double> composition; //percent composition
			double energy;
			double thickness;
		};
	
		//Main integration functions
		double GetEnergyLoss(const Parameters& params);
		double GetReverseEnergyLoss(const Parameters& params);
		
		//Helpers
		double GetTotalStoppingPower(const Parameters& params, double current_energy);
		double GetElectronicStoppingPower(const Parameters& params, double current_energy);
		double GetNuclearStoppingPower(const Parameters& params, double current_energy);
		double Hydrogen_dEdx_Low(double ePerU, int z);
		double Hydrogen_dEdx_Med(double ePerU, int z);
		double Hydrogen_dEdx_High(double ePerU, double massP, double energy, int z);
		double CalculateEffectiveChargeRatio(double ePerU, int zp, int z);

	}
}

#endif

#include "EnergyLoss.h"
#include "EnergyLossConstants.h"
#include <cstdlib>
#include <cmath>

namespace NucKage {

	namespace EnergyLoss {

		double GetEnergyLoss(const Parameters& params)
		{
	
			if(params.thickness == 0.0 || params.energy == 0.0 || params.ZP == 0)
				return 0.0;


			double energyFinal = params.energy;
			double xTraversed = 0;
			double xStep = 0.25*params.thickness; //initial step in x
			double energyStep = GetTotalStoppingPower(params, energyFinal)*xStep/1000.0; //initial step in e
			double energyThreshold = 0.05*params.energy;
		
			int depth=0;
		
			while(true)
			{
				//If intial guess of step size is too large, shrink until in range
				if(energyStep/energyFinal > maxFracStep && depth < maxDepth)
				{
					depth++;
					xStep *= 0.5;
					energyStep = GetTotalStoppingPower(params, energyFinal)*xStep/1000.0;
				}
				else if((xStep + xTraversed) >= params.thickness)
				{ //last chunk
				  	xStep = params.thickness - xTraversed; //get valid portion of last chunk
				  	energyFinal -= GetTotalStoppingPower(params, energyFinal)*xStep/1000.0;
					if(energyFinal <= energyThreshold)
						return params.energy;
					else
						break;
				} 
				else if(depth == maxDepth)
				{
					return params.energy;
				}
				else
				{
					xTraversed += xStep;
					energyStep = GetTotalStoppingPower(params, energyFinal)*xStep/1000.0;
					energyFinal -= energyStep;
					if(energyFinal <= energyThreshold)
						return params.energy;
				}
			}
			return params.energy - energyFinal;
		}
	
		double GetReverseEnergyLoss(const Parameters& params)
		{
			double energyInitial = params.energy;
			double xTraversed = 0.0;
			double xStep = 0.25*params.thickness; //initial step in x
			double energyStep = GetTotalStoppingPower(params, energyInitial)*xStep/1000.0; //initial step in E
		
			while(true)
			{
				if(energyStep/energyInitial > maxFracStep)
				{
					xStep *= 0.5;
					energyStep = GetTotalStoppingPower(params, energyInitial)*xStep/1000.0;
				}
				else if (xTraversed+xStep > params.thickness)
				{
					xStep = params.thickness - xTraversed;
					energyInitial += GetTotalStoppingPower(params, energyInitial)*xStep/1000.0;
					break;
				}
				else
				{
					xTraversed += xStep;
					energyStep = GetTotalStoppingPower(params, energyInitial)*xStep/1000.0;
					energyInitial += energyStep;
				}
			}
		
			return energyInitial-params.energy;
		}
	
		/*Wrapper function for aquiring total stopping (elec + nuc)*/
		double GetTotalStoppingPower(const Parameters& params, double current_energy)
		{
			if(params.ZP == 0)
				return GetNuclearStoppingPower(params, current_energy);
		
			return GetElectronicStoppingPower(params, current_energy)+GetNuclearStoppingPower(params, current_energy);
		}
	
		/*Charge rel to H*/
		double GetElectronicStoppingPower(const Parameters& params, double current_energy)
		{
			//Wants in units of keV
			current_energy *= 1000.0;
			double ePerU = current_energy/params.massP;
			std::vector<double> values;
			if(ePerU > maxHEperU)
			{
				return 0.0;
			}
			else if (ePerU > 1000.0)
			{
				for(auto& z: params.ZT)
					values.push_back(Hydrogen_dEdx_High(ePerU, params.massP, current_energy, z));
			}
			else if (ePerU > 10.0)
			{
				for(auto& z: params.ZT)
					values.push_back(Hydrogen_dEdx_Med(ePerU, z));
			}
			else if (ePerU > 0.0)
			{
				for(auto& z: params.ZT)
					values.push_back(Hydrogen_dEdx_Low(ePerU, z));
			}
			else
			{
				return 0.0;
			}
		
			if(values.size() == 0)
				return 0.0; //bad
		
			if(params.ZP > 1)
			{ //not hydrogen, need to account for effective charge
				for(unsigned int i=0; i<values.size(); i++)
					values[i] *= CalculateEffectiveChargeRatio(ePerU, params.ZP, params.ZT[i]);
			}
			
			double stopping_total = 0;
			double conversion_factor = 0;
			for(size_t i=0; i< params.ZT.size(); i++)
			{
				conversion_factor += params.composition[i]*naturalMassList[params.ZT[i]];
				stopping_total += values[i]*params.composition[i];
			}
			stopping_total *= avogadro/conversion_factor;
		
			return stopping_total;
		}
	
		//Returns units of keV/(ug/cm^2)
		double GetNuclearStoppingPower(const Parameters& params, double energy)
		{
			energy *= 1000.0;
			double stopping_total = 0.0;
			double sn, x, epsilon, conversion_factor, massT;
			for(size_t i=0; i<params.ZT.size(); i++)
			{
				massT = naturalMassList[params.ZT[i]];
				x = (params.massP + massT) * std::sqrt(std::pow(params.ZP, 2.0/3.0) + std::pow(params.ZT[i], 2.0/3.0));
				epsilon = 32.53*massT*energy/(params.ZP*params.ZT[i]*x);
				sn = 8.462*(0.5*std::log(1.0+epsilon)/(epsilon+0.10718*std::pow(epsilon, 0.37544)))*params.ZP*params.ZT[i]*params.massP/x;
				conversion_factor = avogadro/massT;
				stopping_total += sn*conversion_factor*params.composition[i];
			}
		
			return stopping_total;
		}
	
		double Hydrogen_dEdx_Low(double ePerU, int z)
		{
			return std::sqrt(ePerU)*hydrogenCoefficients[z][0];
		}
	
		double Hydrogen_dEdx_Med(double ePerU, int z)
		{
			double x = hydrogenCoefficients[z][1]*std::pow(ePerU, 0.45);
			double y = hydrogenCoefficients[z][2]/ePerU * std::log(1.0+hydrogenCoefficients[z][3]/ePerU+hydrogenCoefficients[z][4]*ePerU);
			return x*y/(x+y);
		}
	
		double Hydrogen_dEdx_High(double ePerU, double massP, double energy, int z)
		{
			energy /= 1000.0; //back to MeV for ease of beta calc
			double beta_sq = energy * (energy+2.0*massP/mev2u)/std::pow(energy+massP/mev2u, 2.0);
			double alpha = hydrogenCoefficients[z][5]/beta_sq;
			double epsilon = hydrogenCoefficients[z][6]*beta_sq/(1.0-beta_sq) - beta_sq - hydrogenCoefficients[z][7];
			for(int i=1; i<5; i++)
				epsilon += hydrogenCoefficients[z][7+i]*std::pow(std::log(ePerU), i);
		
			return alpha * std::log(epsilon);
		}
	
		/*Charge rel to H*/
		double CalculateEffectiveChargeRatio(double ePerU, int zp, int z)
		{
			double z_ratio;
			if(zp == 2)
			{
				double ln_epu = std::log(ePerU);
				double gamma = 1.0+(0.007+0.00005*z)*std::exp(-std::pow(7.6-ln_epu,2.0));
				double alpha = 0.7446 + 0.1429*ln_epu + 0.01562*std::pow(ln_epu, 2.0) - 0.00267*std::pow(ln_epu,3.0)
								+ 1.338E-6*std::pow(ln_epu,8.0);
				z_ratio = gamma*(1.0-std::exp(-alpha))*2.0;
			}
			else if (zp == 3)
			{
				double ln_epu = std::log(ePerU);
				double gamma = 1.0+(0.007+0.00005*z)*std::exp(-std::pow(7.6-ln_epu,2.0));
				double alpha = 0.7138+0.002797*ePerU+1.348E-6*std::pow(ePerU, 2.0);
				z_ratio = gamma*(1-std::exp(-alpha))*3.0;
			}
			else
			{
				double B = 0.886*std::pow(ePerU/25.0, 0.5)/std::pow(zp, 2.0/3.0);
				double A = B + 0.0378*std::sin(M_PI/2.0*B);
				z_ratio = (1.0 - std::exp(-A)*(1.034-0.1777*std::exp(-0.08114*zp)))*zp;
			}
			return z_ratio*z_ratio; //for stopping power uses ratio sq. 
		}
	}
}
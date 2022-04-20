#include <iostream>
#include "EnergyLoss/Target.h"
#include <cmath>
#include "Utils/Timer.h"

namespace NucKage {
	void EnergyLossTest()
	{
		Target target;
	
		std::vector<int> Z = {6};
		std::vector<int> A = {12};
		std::vector<int> S = {1};
	
		double thickness = 100.0; //ug/cm^2
	
		target.SetParameters(Z, S, thickness);
	
		int ZP = 2;
		int AP = 4;
		double proj_angle = 0.0; //radians
		double energy = 10.0; //MeV
		double eloss_result;
		float duration;
	
		std::cout<<"------------EnergyLoss Unit Tests---------------"<<std::endl;
		std::cout<<"Target given components: "<<std::endl;
		for(size_t i=0; i<Z.size(); i++)
			std::cout<<"Z: "<<Z[i]<<" S: "<<S[i]<<std::endl;
		std::cout<<"Target given thickness: "<<thickness<<" ug/cm^2"<<std::endl;
		auto& params = target.GetParameters();
		std::cout<<"Target recieved parameters: "<<std::endl;
		for(size_t i=0; i<Z.size(); i++)
			std::cout<<"Z: "<<params.ZT[i]<<" pecent comp: "<<params.composition[i]<<std::endl;
		std::cout<<"Target recieved thickness: "<<target.GetTotalThickness()<<std::endl;
		std::cout<<"Projectile Z: "<<ZP<<" A: "<<AP<<" energy: "<<energy<<" MeV"<<" angle: "<<proj_angle*180.0/M_PI<<" deg"<<std::endl;
		std::cout<<"Testing whole target energy loss..."<<std::endl;
		Timer stopwatch("ELossTimer");
		eloss_result = target.GetEnergyLossTotal(ZP, AP, energy, proj_angle);
		duration  = stopwatch.ElapsedMilliseconds();
		std::cout<<"Calculated energy loss of "<<eloss_result<<" MeV in "<<duration<<" ms"<<std::endl;
		std::cout<<"Testing half target energy loss..."<<std::endl;
		stopwatch.Restart();
		eloss_result = target.GetEnergyLossFractionalDepth(ZP, AP, energy, proj_angle, 0.5);
		duration  = stopwatch.ElapsedMilliseconds();
		std::cout<<"Calculated energy loss of "<<eloss_result<<" MeV in "<<duration<<" ms"<<std::endl;
		std::cout<<"Testing whole target reverse energy loss..."<<std::endl;
		stopwatch.Restart();
		eloss_result = target.GetReverseEnergyLossTotal(ZP, AP, energy, proj_angle);
		duration  = stopwatch.ElapsedMilliseconds();
		std::cout<<"Calculated energy loss of "<<eloss_result<<" MeV in "<<duration<<" ms"<<std::endl;
		std::cout<<"Testing reverse half target energy loss..."<<std::endl;
		stopwatch.Restart();
		eloss_result = target.GetReverseEnergyLossFractionalDepth(ZP, AP, energy, proj_angle, 0.5);
		duration  = stopwatch.ElapsedMilliseconds();
		std::cout<<"Calculated energy loss of "<<eloss_result<<" MeV in "<<duration<<" ms"<<std::endl;
		std::cout<<"------------------------------------------------"<<std::endl;
	}

}
#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <TLorentzVector.h>
#include "MassLookup.h"

namespace NucKage {

	struct Nucleus
	{
		Nucleus() {};

		Nucleus(int z, int a)
		{
			Z = z;
			A = a;
			MassLookup& table = MassLookup::GetInstance();
			mass = table.FindMass(Z, A);
			symbol = table.FindSymbol(Z, A);
			pvector.SetPxPyPzE(0.,0.,0.,mass);
		}

		void SetIsotope(int z, int a)
		{
			Z = z;
			A = a;
			MassLookup& table = MassLookup::GetInstance();
			mass = table.FindMass(Z, A);
			symbol = table.FindSymbol(Z, A);
			pvector.SetPxPyPzE(0.,0.,0.,mass);
		}

		int Z=0;
		int A=0;
		double mass=0.0;
		std::string symbol="";
		TLorentzVector pvector;

		//Detector crap
		bool detected=false;
		std::string detectorName="";
		int detectorID=-1;
		int detectorFrontChannel=-1;
		int detectorBackChannel=-1;
		double rho=0.0;
	};
}

#endif
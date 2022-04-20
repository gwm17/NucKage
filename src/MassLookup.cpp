/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#include "MassLookup.h"
#include <iostream>

namespace NucKage {

	MassLookup* MassLookup::s_instance = new MassLookup();
	
	MassLookup::MassLookup()
	{
		std::ifstream massfile("etc/mass.txt");
		if(massfile.is_open())
		{
			std::string junk, A, element;
			int Z;
			double atomicMassBig, atomicMassSmall, isotopicMass;
			getline(massfile,junk);
			getline(massfile,junk);
			while(massfile>>junk)
			{
				massfile>>Z>>A>>element>>atomicMassBig>>atomicMassSmall;
				isotopicMass = (atomicMassBig + atomicMassSmall*1e-6 - Z*electron_mass)*u_to_mev;
				std::string key = "("+std::to_string(Z)+","+A+")";
				massTable[key] = isotopicMass;
				elementTable[Z] = element;
			}
		}
		else
		{
			std::cerr<<"ERR -- Unable to open mass file, check etc directory"<<std::endl;
		}
	}
	
	MassLookup::~MassLookup() {}
	
	//Returns nuclear mass in MeV
	double MassLookup::FindMass(int Z, int A)
	{
		std::lock_guard<std::mutex> guard(m_massMutex);
		std::string key = "("+std::to_string(Z)+","+std::to_string(A)+")";
		auto data = massTable.find(key);
		if(data == massTable.end())
		{
			std::cerr<<"WARN -- Unable to find mass of (Z,A)=("<<Z<<","<<A<<")."<<std::endl;
			return 0.0;
		}
	
		return data->second;
	}
	
	//returns element symbol
	std::string MassLookup::FindSymbol(int Z, int A)
	{
		std::lock_guard<std::mutex> guard(m_massMutex);
		auto data = elementTable.find(Z);
		if(data == elementTable.end())
		{
			std::cerr<<"WARN -- Unable to find symbol of (Z,A)=("<<Z<<","<<A<<")."<<std::endl;
			return "";
		}
	
		std::string fullsymbol = std::to_string(A) + data->second;
		return fullsymbol;
	}

}
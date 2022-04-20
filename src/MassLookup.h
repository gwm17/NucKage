/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

Converted to true singleton to simplify usage -- Aug. 2021 GWM
*/
#ifndef MASS_LOOKUP_H
#define MASS_LOOKUP_H

#include <fstream>
#include <string>
#include <unordered_map>
#include <mutex>

namespace NucKage {

	class MassLookup
	{
	public:
		MassLookup();
		~MassLookup();
		double FindMass(int Z, int A);
		std::string FindSymbol(int Z, int A);
	
		inline static MassLookup& GetInstance() { return *s_instance; }
	
	private:
		static MassLookup* s_instance;
		std::unordered_map<std::string, double> massTable;
		std::unordered_map<int, std::string> elementTable;

		std::mutex m_massMutex;
	
		//constants
		static constexpr double u_to_mev = 931.4940954;
		static constexpr double electron_mass = 0.000548579909;
	    
	};

}

#endif

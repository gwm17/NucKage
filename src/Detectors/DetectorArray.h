#ifndef DETECTOR_ARRAY_H
#define DETECTOR_ARRAY_H

#include "SabreDetector.h"
#include "ReactorChain.h"
#include "FocalPlaneDetector.h"

namespace NucKage {

	class DetectorArray
	{
	public:
		DetectorArray();
		~DetectorArray();

		inline void SetFocalPlane(const FocalPlaneDetector::Parameters& params) { m_focalPlane.SetParameters(params); }
		void SetSabre(bool draw); //SABRE is always SABRE, no user input. Simply turn it on or off.

		void ProcessData(ChainResult& data);
		void MakeSabreFile(const std::string& name);
		void TestSabre();
		inline double GetFPSolidAngle() { return m_focalPlane.CalculateSolidAngle(); }

	private:
		static constexpr double s_deg2rad = M_PI/180.0;
		std::vector<SabreDetector> m_sabre;
		FocalPlaneDetector m_focalPlane;
	};

}

#endif
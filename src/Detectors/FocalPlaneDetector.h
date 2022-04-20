#ifndef FOCAL_PLANE_DETECTOR_H
#define FOCAL_PLANE_DETECTOR_H

#include "Nucleus.h"
#include "TRotation.h"

namespace NucKage {

	class FocalPlaneDetector
	{
	public:
		struct Parameters
		{
			Parameters(double theta, double bf) :
				angle(theta*M_PI/180.0), bfield(bf)
			{
			}
			double angle; //rad
			double bfield;  //kG
		};
		FocalPlaneDetector();
		FocalPlaneDetector(const Parameters& params);
		~FocalPlaneDetector();

		void SetParameters(const Parameters& params);
		double CalculateSolidAngle();

		void CheckNucleus(Nucleus& nucleus);

	private:
		bool PassesAperture(double theta, double phi);

		//inline double FullPhi(double phi) { return phi >= 0.0 ? phi : 2.0*M_PI+phi; }

		//FP Always has DetectorID of 99
		static constexpr int s_detectorID = 99;

		/*Geometry*/
		//Aperture coords when at spsTheta = 0.0 (centered on +z)
		static constexpr double s_apertureLeftX = -0.0023;
		static constexpr double s_apertureRightX = 0.0024;
		static constexpr double s_apertureLeftZ = 0.1731;
		static constexpr double s_apertureRightZ = 0.1842;
		static constexpr double s_apertureTopY = -0.0071;
		static constexpr double s_apertureBottomY = 0.0071;
		static constexpr double s_apertureVertZ = 0.1546;

		static constexpr double s_rhoMin = 69.5; //cm
		static constexpr double s_rhoMax = 83.5; //cm

		//Physics
		static constexpr double s_qbrho2p = 0.299792458; //C*1e-9 (kG*cm->MeV/c)

		//data
		double m_bfield;
		double m_angle;
		TRotation m_spsRotation;
	};
}

#endif
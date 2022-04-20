#include "FocalPlaneDetector.h"
#include "RandomGenerator.h"

namespace NucKage {

	FocalPlaneDetector::FocalPlaneDetector() :
		m_bfield(0.0), m_angle(0.0)
	{
		m_spsRotation.RotateY(0.0);
	}

	FocalPlaneDetector::FocalPlaneDetector(const Parameters& params) :
		m_bfield(params.bfield), m_angle(params.angle)
	{
		m_spsRotation.RotateY(-params.angle); //rotate back to SPS 0 deg (90 deg in x-z plane)
	}

	FocalPlaneDetector::~FocalPlaneDetector() {}

	void FocalPlaneDetector::SetParameters(const Parameters& params)
	{
		m_bfield = params.bfield;
		m_angle = params.angle;
		m_spsRotation.RotateY(-params.angle); //rotate back to SPS 0 deg
	}

	bool FocalPlaneDetector::PassesAperture(double theta, double phi)
	{
		if(theta == 0.0 || theta > M_PI/2.0) //Calculation doesn't explicitly check forward vs. backward
				return false;

		TVector3 unitVec(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
		auto rotated_unitVec = m_spsRotation*unitVec;

		double R_at_verticalSlit = s_apertureVertZ/std::cos(rotated_unitVec.Theta());
		double Y_at_verticalSlit = R_at_verticalSlit*std::sin(rotated_unitVec.Theta())*std::sin(rotated_unitVec.Phi());
		if(Y_at_verticalSlit > s_apertureBottomY || Y_at_verticalSlit < s_apertureTopY)
			return false;

		double R_at_horizLeft = s_apertureLeftZ/std::cos(rotated_unitVec.Theta());
		double X_at_horizLeft = R_at_horizLeft*std::sin(rotated_unitVec.Theta())*std::cos(rotated_unitVec.Phi());
		double R_at_horizRight = s_apertureRightZ/std::cos(rotated_unitVec.Theta());
		double X_at_horizRight = R_at_horizRight*std::sin(rotated_unitVec.Theta())*std::cos(rotated_unitVec.Phi());

		if(X_at_horizRight > s_apertureRightX || X_at_horizLeft < s_apertureLeftX)
			return false;

		return true;
	}

	double FocalPlaneDetector::CalculateSolidAngle()
	{
		uint64_t samples = 10e6;
		double theta, phi;
		std::uniform_real_distribution<double> costheta_dist(-1.0, 1.0);
		std::uniform_real_distribution<double> phi_dist(0.0, 2.0*M_PI);
		RandomGenerator& gen = RandomGenerator::GetInstance();
		uint64_t passed=0;
		for(uint64_t i=0; i<samples; i++)
		{
			theta = std::acos(costheta_dist(gen.GetGenerator()));
			phi = phi_dist(gen.GetGenerator());
			if(PassesAperture(theta, phi))
				passed++;
		}
		return double(passed)/double(samples)*4.0*M_PI;
	}

	void FocalPlaneDetector::CheckNucleus(Nucleus& nucleus)
	{
		if(!PassesAperture(nucleus.pvector.Theta(), nucleus.pvector.Phi()))
		{
			nucleus.detected = false;
			nucleus.detectorName = "";
			return;
		}

		double rho = nucleus.pvector.P()/(s_qbrho2p*m_bfield*nucleus.Z);
		if(rho < s_rhoMin || rho > s_rhoMax)
		{
			nucleus.detected = false;
			nucleus.detectorName = "";
			return;
		}

		nucleus.detected = true;
		nucleus.detectorName = "FocalPlane";
		nucleus.detectorID = s_detectorID;
		nucleus.detectorFrontChannel = -1;
		nucleus.detectorBackChannel = -1;
		nucleus.rho = rho;
	}
}
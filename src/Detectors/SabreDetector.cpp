/*

  Class which represents a single MMM detector in the SABRE array at FSU. Origial code by KGH, re-written by
  GWM.

  Distances in meters, angles in radians.

  The channel arrays have four points, one for each corner. The corners are
  as follows, as if looking BACK along beam (i.e. from the target's pov):

  0---------------------1
  |                     |
  |                     |      x
  |                     |      <-----
  |                     |      		|
  |                     |      		|
  3---------------------2      		y
                               (z is hence positive along beam direction) 

  The channel numbers, also as looking back from target pov, are:

  >> rings are 0 -- 15 from inner to outer:

    15 -------------------
    14 -------------------
    13 -------------------
       .
       .
       .
     2 -------------------
     1 -------------------
     0 -------------------

  >> wedges are 0 -- 7 moving counterclockwise:

      7 6 ... 1 0
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |


  >> Note that the detector starts centered on the x-axis (central phi = 0) untilted, and then is rotated to wherever the frick
  	 it is supposed to go; phi = 90 is centered on y axis, pointing down towards the bottom of the scattering chamber

  -- gwm, Dec 2020; based on the og code from kgh

*/


#include "SabreDetector.h"
#include <iostream>

namespace NucKage {

	void PrintMatrix(TRotation& mat)
	{
		std::cout<<"XX: "<<mat.XX()<<" XY: "<<mat.XY()<<" XZ: "<<mat.XZ()<<std::endl;
		std::cout<<"YX: "<<mat.YX()<<" YY: "<<mat.YY()<<" YZ: "<<mat.YZ()<<std::endl;
		std::cout<<"ZX: "<<mat.ZX()<<" ZY: "<<mat.ZY()<<" ZZ: "<<mat.ZZ()<<std::endl;
	}

	SabreDetector::SabreDetector() :
		m_phiCentral(0.0), m_tilt(0.0), m_translation(0.,0.,0.), m_norm_flat(0,0,1.0), m_drawingFlag(true), m_channelSmear(0.0, 1.0),
		m_detectorID(-1)
	{
		m_YRot.RotateY(-1.0*m_tilt);
		m_ZRot.RotateZ(m_phiCentral);
	
		//Initialize the coordinate arrays
		m_ringCoords_flat.resize(s_nRings);
		m_ringCoords_tilt.resize(s_nRings);
		m_wedgeCoords_flat.resize(s_nWedges);
		m_wedgeCoords_tilt.resize(s_nWedges);
		for(int i=0; i<s_nRings; i++)
		{
			m_ringCoords_flat[i].resize(4);
			m_ringCoords_tilt[i].resize(4);
		}
		for(int i=0; i<s_nWedges; i++)
		{
			m_wedgeCoords_flat[i].resize(4);
			m_wedgeCoords_tilt[i].resize(4);
		}
	
		m_deltaR_flat = s_Router - s_Rinner;
		m_deltaR_flat_ring = m_deltaR_flat/s_nRings;
	
		CalculateCorners();
	}
	
	SabreDetector::SabreDetector(const Parameters& params) :
		m_phiCentral(params.phiCenter), m_tilt(params.tilt), m_translation(0., 0., params.zOffset), m_norm_flat(0,0,1.0), 
		m_drawingFlag(params.drawing), m_channelSmear(0.0, 1.0), m_detectorID(params.detID)
	{
		m_YRot.RotateY(-1.0*m_tilt); //clockwise rotation
		m_ZRot.RotateZ(m_phiCentral);
	
		//Initialize coordinate arrays
		if(m_drawingFlag)
		{
			m_ringCoords_flat.resize(s_nRings);
			m_ringCoords_tilt.resize(s_nRings);
			m_wedgeCoords_flat.resize(s_nWedges);
			m_wedgeCoords_tilt.resize(s_nWedges);
			for(int i=0; i<s_nRings; i++)
			{
				m_ringCoords_flat[i].resize(4);
				m_ringCoords_tilt[i].resize(4);
			}
			for(int i=0; i<s_nWedges; i++)
			{
				m_wedgeCoords_flat[i].resize(4);
				m_wedgeCoords_tilt[i].resize(4);
			}
		
			m_deltaR_flat = s_Router - s_Rinner;
			m_deltaR_flat_ring = m_deltaR_flat/s_nRings;
			m_deltaPhi_flat_wedge = s_deltaPhi_flat/s_nWedges;
		
			CalculateCorners();
		}
	}
	
	SabreDetector::~SabreDetector() {}
	
	void SabreDetector::CalculateCorners() {
	
		double x0, x1, x2, x3;
		double y0, y1, y2, y3;
		double z0, z1, z2, z3;
	
		//Generate flat ring corner coordinates
		for(int i=0; i<s_nRings; i++)
		{
			x0 = (s_Rinner + m_deltaR_flat_ring*(i+1))*std::cos(-s_deltaPhi_flat/2.0);
			y0 = (s_Rinner + m_deltaR_flat_ring*(i+1))*std::sin(-s_deltaPhi_flat/2.0);
			z0 = 0.0;
			m_ringCoords_flat[i][0].SetXYZ(x0, y0, z0);
	
			x1 = (s_Rinner + m_deltaR_flat_ring*(i))*std::cos(-s_deltaPhi_flat/2.0);
			y1 = (s_Rinner + m_deltaR_flat_ring*(i))*std::sin(-s_deltaPhi_flat/2.0);
			z1 = 0.0;
			m_ringCoords_flat[i][1].SetXYZ(x1, y1, z1);
	
			x2 = (s_Rinner + m_deltaR_flat_ring*(i))*std::cos(s_deltaPhi_flat/2.0);
			y2 = (s_Rinner + m_deltaR_flat_ring*(i))*std::sin(s_deltaPhi_flat/2.0);
			z2 = 0.0;
			m_ringCoords_flat[i][2].SetXYZ(x2, y2, z2);
	
			x3 = (s_Rinner + m_deltaR_flat_ring*(i+1))*std::cos(s_deltaPhi_flat/2.0);
			y3 = (s_Rinner + m_deltaR_flat_ring*(i+1))*std::sin(s_deltaPhi_flat/2.0);
			z3 = 0.0;
			m_ringCoords_flat[i][3].SetXYZ(x3, y3, z3);
		}
	
		//Generate flat wedge corner coordinates
		for(int i=0; i<s_nWedges; i++)
		{
			x0 = s_Router * std::cos(-s_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
			y0 = s_Router * std::sin(-s_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
			z0 = 0.0;
			m_wedgeCoords_flat[i][0].SetXYZ(x0, y0, z0);
	
			x1 = s_Rinner * std::cos(-s_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
			y1 = s_Rinner * std::sin(-s_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
			z1 = 0.0;
			m_wedgeCoords_flat[i][1].SetXYZ(x1, y1, z1);
	
			x2 = s_Rinner * std::cos(-s_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
			y2 = s_Rinner * std::sin(-s_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
			z2 = 0.0;
			m_wedgeCoords_flat[i][2].SetXYZ(x2, y2, z2);
	
			x3 = s_Router * std::cos(-s_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
			y3 = s_Router * std::sin(-s_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
			z3 = 0.0;
			m_wedgeCoords_flat[i][3].SetXYZ(x3, y3, z3);
		}
	
		//Generate tilted rings
		for(int i=0; i<s_nRings; i++)
		{
			for(int j=0; j<4; j++)
				m_ringCoords_tilt[i][j] = TransformToTiltedFrame(m_ringCoords_flat[i][j]);
		}
	
		//Generate tilted wedges
		for(int i=0; i<s_nWedges; i++)
		{
			for(int j=0; j<4; j++)
				m_wedgeCoords_tilt[i][j] = TransformToTiltedFrame(m_wedgeCoords_flat[i][j]);
		}
	}

	void SabreDetector::CheckNucleus(Nucleus& nucleus)
	{
		auto result = GetTrajectoryRingWedge(nucleus.pvector.Theta(), nucleus.pvector.Phi());
		if(result.first == -1 || result.second == -1)
		{
			//std::cout<<"Here: first: "<<result.first<<" second: "<<result.second<<std::endl;
			nucleus.detected = false;
			nucleus.detectorName = "";
		}
		else
		{
			
			nucleus.detected = true;
			nucleus.detectorName = "SABRE";
			nucleus.detectorID = m_detectorID;
			nucleus.detectorFrontChannel = result.first;
			nucleus.detectorBackChannel = result.second;
		}
	}
	
	/*
		Given a unit vector (R=1, theta, phi) which corresponds to some particle's trajectory,
		determine whether that particle will intersect with this SABRE detector. If it does calculate
		the coordinates of the hit. The equation is as follows:
		
		Rz(eta)*Ry(psi)*Flat_vector(R', theta'=PI/2, phi') + translation = Tilted_vector(R, theta, phi)
	
		Where Rz is the ZRotation, Ry is the YRotation, F_vector is the vector of the coordinates in the flat detector frame,
		and Tilted_vector is the vector of the hit coordinates in the tilted frame. The theta and phi of the the Tilted_vector correspond
		to the input arguments of the function.
	
		It checks to deterime whether or not the particle hits within the borders (read: edges) of the SABRE detector, and does not account for
		the spacing between rings and wedges.
	
		!NOTE: This currently only applies to a configuration where there is no translation in x & y. The math becomes significantly messier in these cases.
		Also, don't use tan(). It's behavior near PI/2 makes it basically useless for these.
	*/
	TVector3 SabreDetector::GetTrajectoryCoordinates(double theta, double phi)
	{
		if(m_translation.X() != 0.0 || m_translation.Y() != 0.0)
			return TVector3();
	
		//Calculate the *potential* phi in the flat detector
		double phi_numerator = std::cos(m_tilt)*(std::sin(phi)*std::cos(m_phiCentral) - std::sin(m_phiCentral)*std::cos(phi));
		double phi_denominator = std::cos(m_phiCentral)*std::cos(phi) + std::sin(m_phiCentral)*std::sin(phi);
		double phi_flat = std::atan2(phi_numerator, phi_denominator);
		if(phi_flat < 0) phi_flat += M_PI*2.0;
	
		//Calculate the *potential* R in the flat detector
		double r_numerator = m_translation.Z()*std::cos(phi)*std::sin(theta);
		double r_denominator = std::cos(phi_flat)*std::cos(m_phiCentral)*std::cos(m_tilt)*std::cos(theta) - std::sin(phi_flat)*std::sin(m_phiCentral)*std::cos(theta) - std::cos(phi_flat)*std::sin(m_tilt)*std::cos(phi)*std::sin(theta);
		double r_flat = r_numerator/r_denominator;
	
		//Calculate the distance from the origin to the hit on the detector
		double R_to_detector = (r_flat*std::cos(phi_flat)*std::sin(m_tilt) + m_translation.Z())/std::cos(theta);
		double xhit = R_to_detector*std::sin(theta)*std::cos(phi);
		double yhit = R_to_detector*std::sin(theta)*std::sin(phi);
		double zhit = R_to_detector*std::cos(theta);
	
	
		//Check to see if our flat coords fall inside the flat detector
		if(IsInside(r_flat, phi_flat))
			return TVector3(xhit, yhit, zhit);
		else
			return TVector3();
	}
	
	/*
		Given a unit vector (R=1, theta, phi) which corresponds to some particle's trajectory,
		determine whether that particle will intersect with this SABRE detector. If it does determine 
		which ring and wedge the hit occurs in. The equation is as follows:
		
		Rz(eta)*Ry(psi)*Flat_vector(R', theta'=PI/2, phi') + translation = Tilted_vector(R, theta, phi)
	
		Where Rz is the ZRotation, Ry is the YRotation, F_vector is the vector of the coordinates in the flat detector frame,
		and Tilted_vector is the vector of the hit coordinates in the tilted frame. The theta and phi of the the Tilted_vector correspond
		to the input arguments of the function.
	
		Then using the flat coordinate R' and phi' determine which ring/wedge channels are hit. For precision purposes, the channel is not calculated, but
		rather found using comparisions. This method accounts for the spacing between rings and wedges.
	
		!NOTE: This currently only applies to a configuration where there is no translation in x & y. The math becomes significantly messier in these cases.
		Also, don't use tan(). It's behavior near PI/2 makes it basically useless for these.
	*/
	std::pair<int, int> SabreDetector::GetTrajectoryRingWedge(double theta, double phi)
	{
		phi = phi < 0 ? 2.0*M_PI + phi : phi;
		if(m_translation.X() != 0.0 || m_translation.Y() != 0.0)
			return std::make_pair(-1, -1);
	
		//Calculate the *potential* phi in the flat detector
		double phi_numerator = std::cos(m_tilt)*(std::sin(phi)*std::cos(m_phiCentral) - std::sin(m_phiCentral)*std::cos(phi));
		double phi_denominator = std::cos(m_phiCentral)*std::cos(phi) + std::sin(m_phiCentral)*std::sin(phi);
		double phi_flat = std::atan2(phi_numerator, phi_denominator);
		if(phi_flat < 0) phi_flat += M_PI*2.0;
	
		//Calculate the *potential* R in the flat detector
		double r_numerator = m_translation.Z()*std::cos(phi)*std::sin(theta);
		double r_denominator = std::cos(phi_flat)*std::cos(m_phiCentral)*std::cos(m_tilt)*std::cos(theta) - std::sin(phi_flat)*std::sin(m_phiCentral)*std::cos(theta) - std::cos(phi_flat)*std::sin(m_tilt)*std::cos(phi)*std::sin(theta);
		double r_flat = r_numerator/r_denominator;
	
		//Calculate the distance from the origin to the hit on the detector
		//double R_to_detector = (r_flat*std::cos(phi_flat)*std::sin(m_tilt) + m_translation.GetZ())/std::cos(theta);
	
		//Check to see if our flat coords fall inside the flat detector
		if(IsInside(r_flat, phi_flat))
		{
			int ringch, wedgech;
			if(phi_flat > M_PI) phi_flat -= 2.0*M_PI; //Need phi in terms of [-deltaPhi_flat/2, deltaPhi_flat/2]
			for(int i=0; i<s_nRings; i++) {
				if(IsRingTopEdge(r_flat, i) || IsRingBottomEdge(r_flat, i))
				{ //If it falls in the interstrip spacing, kill it
					ringch = -1;
					break;
				}
				else if(IsRing(r_flat, i))
				{
					ringch = i;
					break;
				}
			}
			for(int i=0; i<s_nWedges; i++)
			{
				if(IsWedgeTopEdge(phi_flat, i) || IsWedgeBottomEdge(phi_flat, i))
				{ //If it falls in the interstrip spacing, kill it
					wedgech = -1;
					break;
				}
				else if(IsWedge(phi_flat, i))
				{
					wedgech = i;
					break;
				}
			}
			return std::make_pair(ringch, wedgech);
		}
		else
		{
			return std::make_pair(-1,-1);
		}
	}
	
	/*
		Given a ring/wedge of this SABRE detector, calculate the coordinates of a hit.
		Currently gives a point in the *center* of the pixel. Better method would be to
		randomly wiggle the point within the pixel. Method intended for use with data, or
		to smear out simulated data to mimic real data.
	*/
	TVector3 SabreDetector::GetHitCoordinates(int ringch, int wedgech)
	{
		if(!CheckRingChannel(ringch) || !CheckWedgeChannel(wedgech))
			return TVector3();

		RandomGenerator& gen = RandomGenerator::GetInstance();
		double ringSmear = m_channelSmear(gen.GetGenerator());
		double wedgeSmear = m_channelSmear(gen.GetGenerator());
		double r_center  = s_Rinner + (ringch + ringSmear)*m_deltaR_flat_ring;
		double phi_center = -s_deltaPhi_flat/2.0 + (wedgech + wedgeSmear)*m_deltaPhi_flat_wedge;
		double x = r_center*std::cos(phi_center);
		double y = r_center*std::sin(phi_center);
		double z = 0;
	
		TVector3 hit(x, y, z);

		return TransformToTiltedFrame(hit);
	}

}
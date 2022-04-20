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

  -- GWM, Dec 2020; based on the og code from kgh

*/

#ifndef SABREDETECTOR_H
#define SABREDETECTOR_H

#include <vector>
#include <cmath>

#include "TVector3.h"
#include "TRotation.h"

#include "RandomGenerator.h"
#include "Nucleus.h"

namespace NucKage {

	class SabreDetector
	{
	public:
		struct Parameters
		{
			Parameters(double phi, double tilt, double z, bool draw, int id) :
				phiCenter(phi*M_PI/180.0), tilt(tilt*M_PI/180.0), zOffset(z), drawing(draw), detID(id)
			{
			}
			double phiCenter=0.0;
			double tilt=0.0;
			double zOffset=0.0;
			bool drawing=false;
			int detID=-1;
		};
	
		SabreDetector();
		SabreDetector(const Parameters& params);
		~SabreDetector();

		void CheckNucleus(Nucleus& nucleus);
	
		/*Return coordinates of the corners of each ring/wedge in SABRE*/
		inline TVector3 GetRingFlatCoords(int ch, int corner) { return m_drawingFlag && CheckRingLocation(ch, corner) ? m_ringCoords_flat[ch][corner] : TVector3(); }
		inline TVector3 GetWedgeFlatCoords(int ch, int corner) { return m_drawingFlag && CheckWedgeLocation(ch, corner) ? m_wedgeCoords_flat[ch][corner] : TVector3(); }
		inline TVector3 GetRingTiltCoords(int ch, int corner) { return m_drawingFlag && CheckRingLocation(ch, corner) ? m_ringCoords_tilt[ch][corner] : TVector3(); }
		inline TVector3 GetWedgeTiltCoords(int ch, int corner) { return m_drawingFlag && CheckWedgeLocation(ch, corner) ? m_wedgeCoords_tilt[ch][corner] : TVector3(); }
	
		TVector3 GetTrajectoryCoordinates(double theta, double phi);
		std::pair<int, int> GetTrajectoryRingWedge(double theta, double phi);
		TVector3 GetHitCoordinates(int ringch, int wedgech);

		inline const int GetDetectorID() const { return m_detectorID; }
	
		/*Basic getters*/
		inline TVector3 GetNormTilted() { return TransformToTiltedFrame(m_norm_flat); }
	
	
	private:
	
		/*Class constants*/
		static constexpr int s_nRings = 16;
		static constexpr int s_nWedges = 8;
		static constexpr double s_deg2rad = M_PI/180.0;
		static constexpr double s_Rinner = 0.0326;
		static constexpr double s_Router = 0.1351;
		static constexpr double s_deltaPhi_flat = 54.4 * M_PI/180.0;
		/*These are implicitly the width of the spacing between detector active strips*/
		static constexpr double position_tol = 0.0001; //0.1 mm position tolerance
		static constexpr double angular_tol = 0.1*M_PI/180.0; // 0.1 degree angular tolerance
	
		void CalculateCorners();
	
		/*Performs the transformation to the tilted,rotated,translated frame of the SABRE detector*/
		inline TVector3 TransformToTiltedFrame(TVector3& vector) { return (vector.Transform(m_YRot)).Transform(m_ZRot) + m_translation; }
	
		/*Determine if a given channel/corner combo is valid*/
		inline bool CheckRingChannel(int ch) { return (ch<s_nRings && ch>=0) ? true : false; }
		inline bool CheckWedgeChannel(int ch) { return (ch<s_nWedges && ch >=0) ? true : false; }
		inline bool CheckCorner(int corner) { return (corner < 4 && corner >=0) ? true : false; }
		inline bool CheckRingLocation(int ch, int corner) { return CheckRingChannel(ch) && CheckCorner(corner); }
		inline bool CheckWedgeLocation(int ch, int corner) { return CheckWedgeChannel(ch) && CheckCorner(corner); }
	
		/*
			For all of the calculations, need a limit precision to determine if values are actually equal or not
			Here the approx. size of the strip spacing is used as the precision.
		*/
		inline bool CheckPositionEqual(double val1,double val2) { return fabs(val1-val2) > position_tol ? false : true; }
		inline bool CheckAngleEqual(double val1,double val2) { return fabs(val1-val2) > angular_tol ? false : true; }
	
		/*Determine if a hit is within the bulk detector*/
		inline bool IsInside(double r, double phi)
		{ 
			double phi_1 = s_deltaPhi_flat/2.0;
			double phi_2 = M_PI*2.0 - s_deltaPhi_flat/2.0;
			return (((r > s_Rinner && r < s_Router) || CheckPositionEqual(r, s_Rinner) || CheckPositionEqual(r, s_Router))
					  && (phi > phi_2 || phi < phi_1 || CheckAngleEqual(phi, phi_1) || CheckAngleEqual(phi, phi_2)));
		}
	
		/*
			For a given radius/phi are you inside of a given ring/wedge channel,
			or are you on the spacing between these channels
		*/
		inline bool IsRing(double r, int ringch)
		{
			double ringtop = s_Rinner + m_deltaR_flat_ring*(ringch + 1);
			double ringbottom = s_Rinner + m_deltaR_flat_ring*(ringch);
			return (r>ringbottom && r<ringtop); 
		}
	
		inline bool IsRingTopEdge(double r, int ringch)
		{
			double ringtop = s_Rinner + m_deltaR_flat_ring*(ringch + 1);
			return CheckPositionEqual(r, ringtop); 
		}
	
		inline bool IsRingBottomEdge(double r, int ringch)
		{
			double ringbottom = s_Rinner + m_deltaR_flat_ring*(ringch);
			return CheckPositionEqual(r, ringbottom); 
		}
	
		inline bool IsWedge(double phi, int wedgech)
		{
			double wedgetop = -s_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech+1);
			double wedgebottom = -s_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech);
			return ((phi>wedgebottom && phi<wedgetop));
		}
	
		inline bool IsWedgeTopEdge(double phi, int wedgech)
		{
			double wedgetop = -s_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech+1);
			return CheckAngleEqual(phi, wedgetop);
		}
	
		inline bool IsWedgeBottomEdge(double phi, int wedgech)
		{
			double wedgebottom = -s_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech);
			return CheckAngleEqual(phi, wedgebottom);
		}


		/*Class data*/
		double m_phiCentral, m_tilt;
		TVector3 m_translation;
		TRotation m_YRot;
		TRotation m_ZRot;
		double m_deltaR_flat, m_deltaR_flat_ring, m_deltaPhi_flat_wedge;
		TVector3 m_norm_flat;
		bool m_drawingFlag;
		int m_detectorID;

		std::uniform_real_distribution<double> m_channelSmear;
	
		std::vector<std::vector<TVector3>> m_ringCoords_flat, m_wedgeCoords_flat;
		std::vector<std::vector<TVector3>> m_ringCoords_tilt, m_wedgeCoords_tilt;
	
	};
}

#endif

#include "DetectorArray.h"
#include <fstream>

namespace NucKage {

	DetectorArray::DetectorArray() {}

	DetectorArray::~DetectorArray() {}

	void DetectorArray::SetSabre(bool draw)
	{
		m_sabre.emplace_back(SabreDetector::Parameters(306.0, 40.0, -0.1245, draw, 0));
		m_sabre.emplace_back(SabreDetector::Parameters(18.0, 40.0, -0.1245, draw, 1));
		m_sabre.emplace_back(SabreDetector::Parameters(234.0, 40.0, -0.1245, draw, 2));
		m_sabre.emplace_back(SabreDetector::Parameters(162.0, 40.0, -0.1245, draw, 3));
		m_sabre.emplace_back(SabreDetector::Parameters(90.0, 40.0, -0.1245, draw, 4));
	}

	void DetectorArray::ProcessData(ChainResult& data)
	{
		for(int i=0; i<(data.products.size() -1); i++)
		{
			Nucleus& particle = data.products[i].ejectile;
			m_focalPlane.CheckNucleus(particle);
			if(particle.detected)
				continue;
			
			for(auto& sabdet : m_sabre)
			{
				sabdet.CheckNucleus(particle);
				if(particle.detected)
					break;
			}
		}

		//Last reaction in the chain look at both residual and ejectile

		Nucleus& particle = data.products[data.products.size()-1].ejectile;
		m_focalPlane.CheckNucleus(particle);
		if(!particle.detected)
		{
			for(auto& sabdet : m_sabre)
			{
				sabdet.CheckNucleus(particle);
				if(particle.detected)
					break;
			}
		}
		particle = data.products[data.products.size()-1].residual;
		m_focalPlane.CheckNucleus(particle);
		if(!particle.detected)
		{
			for(auto& sabdet : m_sabre)
			{
				sabdet.CheckNucleus(particle);
				if(particle.detected)
					break;
			}
		}
	}

	void DetectorArray::MakeSabreFile(const std::string& name)
	{
		std::ofstream output(name);
		TVector3 vec;
		for(auto& det : m_sabre)
		{
			for(int i=0; i<16; i++)
			{
				for(int j=0; j<4; j++)
				{
					vec = det.GetRingTiltCoords(i, j);
					output<<vec.X()<<" "<<vec.Y()<<" "<<vec.Z();
				}
			}
			for(int i=0; i<8; i++)
			{
				for(int j=0; j<4; j++)
				{
					vec = det.GetWedgeTiltCoords(i, j);
					output<<vec.X()<<" "<<vec.Y()<<" "<<vec.Z();
				}
			}
		}
		output.close();
	}

	void DetectorArray::TestSabre()
	{
		TVector3 vec;
		std::pair<int, int> chans;
		for(auto& det : m_sabre)
		{
			for(int i=0; i<16; i++)
			{
				for(int j=0; j<8; j++)
				{
					vec = det.GetHitCoordinates(i, j);
					chans = det.GetTrajectoryRingWedge(vec.Theta(), vec.Phi());
					if(chans.first != i || chans.second != j)
					{
						std::cout<<"ERR -- Misidentified ring/wedge in SABRE "<<det.GetDetectorID()<<std::endl;
						std::cout<<"Given ring:wedge "<<i<<":"<<j<<" found "<<chans.first<<":"<<chans.second<<std::endl;
					}
				}
			}
		}
	}
}
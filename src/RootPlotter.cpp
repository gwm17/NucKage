#include "RootPlotter.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>

namespace NucKage {

	RootPlotter::RootPlotter() :
		m_queueSize(0)
	{
		TH1::AddDirectory(kFALSE);
	}

	RootPlotter::RootPlotter(const std::string& name) :
		m_queueSize(0)
	{
		TH1::AddDirectory(kFALSE);
		Open(name);
	}

	RootPlotter::~RootPlotter()
	{
		if(m_file && m_openFlag)
		{
			Close();
		}
	}

	void RootPlotter::Open(const std::string& name)
	{
		m_file = TFile::Open(name.c_str(), "RECREATE");
		if(m_file && m_file->IsOpen())
			m_openFlag = true;
		else
			m_openFlag = false;
	}

	void RootPlotter::Close()
	{
		std::lock_guard<std::mutex> guard(m_rootMutex);
		for(auto& iter : m_map)
			iter.second->Write();
		m_file->Close();
		m_openFlag = false;
	}

	void RootPlotter::MyFill2D(const Histo2DParams& params, double valueX, double valueY)
	{
		auto h = std::static_pointer_cast<TH2>(m_map[params.name]);
		if(h) {
			h->Fill(valueX, valueY);
		} else {
			h = std::make_shared<TH2F>(params.name.c_str(), params.name.c_str(), params.binsX, params.minX, params.maxX, params.binsY, params.minY, params.maxY);
			h->Fill(valueX, valueY);
			m_map[params.name] = h;
		}
	}

	void RootPlotter::MyFill1D(const Histo1DParams& params, double valueX)
	{
		auto h = std::static_pointer_cast<TH1>(m_map[params.name]);
		if(h) {
			h->Fill(valueX);
		} else {
			h = std::make_shared<TH1F>(params.name.c_str(), params.name.c_str(), params.binsX, params.minX, params.maxX);
			h->Fill(valueX);
			m_map[params.name] = h;
		}
	}

	void RootPlotter::MyFillGraph(const std::string& name, double valueX, double valueY, int color)
	{
		auto g = std::static_pointer_cast<TGraph>(m_map[name]);
		if(g)
		{
			g->SetPoint(g->GetN(), valueX, valueY);
		}
		else
		{
			g = std::make_shared<TGraph>(1, &valueX, &valueY);
			g->SetName(name.c_str());
			g->SetTitle(name.c_str());
			g->SetMarkerColor(color);
			m_map[name] = g;
		}
		
	}


	void RootPlotter::PlotData()
	{
		if(!IsOpen() || m_queueSize == 0)
			return;

		ChainResult data = PopData();
		Histo1DParams h1pars;
		std::string graph_name;
		for(auto& result : data.products)
		{
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.target.pvector.Theta()*s_rad2deg, result.target.pvector.E()-result.target.pvector.M(), 2);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.target.pvector.Phi())*s_rad2deg, result.target.pvector.E()-result.target.pvector.M(), 2);
			h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_Ex;E_{x} (MeV);counts";
			h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
			MyFill1D(h1pars, result.target.pvector.M()-result.target.mass);

			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.ejectile.pvector.Theta()*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.ejectile.pvector.Phi())*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
			
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.residual.pvector.Theta()*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 5);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.residual.pvector.Phi())*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 5);
			h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_Ex;E_{x} (MeV);counts";
			h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
			MyFill1D(h1pars, result.residual.pvector.M()-result.residual.mass);

			if(!result.projectile.symbol.empty())
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.projectile.pvector.Theta()*s_rad2deg, result.projectile.pvector.E()-result.projectile.pvector.M(), 3);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.projectile.pvector.Phi())*s_rad2deg, result.projectile.pvector.E()-result.projectile.pvector.M(), 3);
				h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KE;KE (MeV);counts";
				h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
				MyFill1D(h1pars, result.projectile.pvector.E()-result.projectile.pvector.M());
			}

			if(result.ejectile.detected)
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvTheta_detect;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.ejectile.pvector.Theta()*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvPhi_detect;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.ejectile.pvector.Phi())*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
				if(result.ejectile.detectorName == "FocalPlane")
				{
					h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_rho;#rho (cm);counts";
					h1pars.binsX = 1400, h1pars.minX = 69.5, h1pars.maxX = 83.5;
					MyFill1D(h1pars, result.ejectile.rho);
				}
			}

			if(result.residual.detected)
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvTheta_detect;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.residual.pvector.Theta()*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 4);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvPhi_detect;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.residual.pvector.Phi())*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 4);
				if(result.residual.detectorName == "FocalPlane")
				{
					h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_rho;#rho (cm);counts";
					h1pars.binsX = 1400, h1pars.minX = 69.5, h1pars.maxX = 83.5;
					MyFill1D(h1pars, result.residual.rho);
				}
			}
			
		}
		m_queueSize--;
	}

	//Testing single thread only! Not thread safe
	void RootPlotter::PlotData(const ChainResult& data)
	{
		if(!IsOpen())
			return;

		Histo1DParams h1pars;
		std::string graph_name;
		for(auto& result : data.products)
		{
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.target.pvector.Theta()*s_rad2deg, result.target.pvector.E()-result.target.pvector.M(), 2);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.target.pvector.Phi())*s_rad2deg, result.target.pvector.E()-result.target.pvector.M(), 2);
			h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.target.symbol+"_Ex;E_{x} (MeV);counts";
			h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
			MyFill1D(h1pars, result.target.pvector.M()-result.target.mass);

			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.ejectile.pvector.Theta()*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.ejectile.pvector.Phi())*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
			
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, result.residual.pvector.Theta()*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 5);
			graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
			MyFillGraph(graph_name, FullPhi(result.residual.pvector.Phi())*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 5);
			h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_Ex;E_{x} (MeV);counts";
			h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
			MyFill1D(h1pars, result.residual.pvector.M()-result.residual.mass);

			if(!result.projectile.symbol.empty())
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KEvTheta;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.projectile.pvector.Theta()*s_rad2deg, result.projectile.pvector.E()-result.projectile.pvector.M(), 3);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KEvPhi;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.projectile.pvector.Phi())*s_rad2deg, result.projectile.pvector.E()-result.projectile.pvector.M(), 3);
				h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.projectile.symbol+"_KE;KE (MeV);counts";
				h1pars.binsX = 300, h1pars.minX = 0.0, h1pars.maxX = 30.0;
				MyFill1D(h1pars, result.projectile.pvector.E()-result.projectile.pvector.M());
			}

			if(result.ejectile.detected)
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvTheta_detect;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.ejectile.pvector.Theta()*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_KEvPhi_detect;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.ejectile.pvector.Phi())*s_rad2deg, result.ejectile.pvector.E()-result.ejectile.pvector.M(), 4);
				if(result.ejectile.detectorName == "FocalPlane")
				{
					h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.ejectile.symbol+"_rho;#rho (cm);counts";
					h1pars.binsX = 1400, h1pars.minX = 69.5, h1pars.maxX = 83.5;
					MyFill1D(h1pars, result.ejectile.rho);
				}
			}

			if(result.residual.detected)
			{
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvTheta_detect;#theta_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, result.residual.pvector.Theta()*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 4);
				graph_name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_KEvPhi_detect;#phi_{Lab}(deg);KE (MeV)";
				MyFillGraph(graph_name, FullPhi(result.residual.pvector.Phi())*s_rad2deg, result.residual.pvector.E()-result.residual.pvector.M(), 4);
				if(result.residual.detectorName == "FocalPlane")
				{
					h1pars.name = "Chain_"+std::to_string(data.chainID)+"_Rxn_"+result.reactorName+"_Nuc_"+result.residual.symbol+"_rho;#rho (cm);counts";
					h1pars.binsX = 1400, h1pars.minX = 69.5, h1pars.maxX = 83.5;
					MyFill1D(h1pars, result.residual.rho);
				}
			}
			
		}
	}
}
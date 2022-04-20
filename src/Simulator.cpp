#include "Simulator.h"
#include <iostream>
#include <fstream>
#include <future>

namespace NucKage {

	Simulator* Simulator::s_instance = nullptr;

	Simulator* CreateSimulator(int nthreads)
	{
		return new Simulator(nthreads);
	}

	Simulator::Simulator() :
		m_outputFile(""), m_initFlag(false), m_samples(0), m_pool(1)
	{
		if(s_instance)
		{
			std::cerr<<"WARN -- Attempting to create a second Simulator instance!"<<std::endl;
			return;
		}

		s_instance = this;
	}

	Simulator::Simulator(int nthreads) :
		m_outputFile(""), m_initFlag(false), m_samples(0), m_pool(nthreads)
	{
		if(s_instance)
		{
			std::cerr<<"WARN -- Attempting to create a second Simulator instance!"<<std::endl;
			return;
		}

		s_instance = this;
	}

	Simulator::~Simulator()
	{
	}

	void Simulator::LoadConfig(const std::string& filename)
	{
		m_initFlag = false;
		std::ifstream input(filename);
		if(!input.is_open())
		{
			std::cerr<<"WARN -- Unable to open input configuration file: "<<filename<<std::endl;
			return;
		}

		std::string junk;
		input>>junk;
		uint64_t temp;
		if(junk == "begin_simulator")
		{
			input>>m_outputFile;
			input>>temp;
			m_samples = temp;
		}
		else
		{
			std::cerr<<"Bad input file, incorrect config in file: "<<filename<<std::endl;
			return;
		}

		std::cout<<"Output: "<<m_outputFile<<" samples: "<<m_samples<<std::endl;

		SamplingParameters temp_params;
		std::vector<int> Z, A, stoich;
		int a;
		double thickness;
		while(input>>junk)
		{
			if(junk == "begin_reactorchain")
			{
				ReactorChain temp_chain;
				while(input>>junk)
				{
					if(junk == "begin_reactor")
					{
						Z.clear();
						A.clear();
						input>>temp_params.meanEx;
						input>>temp_params.sigmaEx;
						input>>temp_params.meanBeamKE;
						input>>temp_params.sigmaBeamKE;
						while(input>>junk)
						{
							if(junk == "begin_nuclei")
								continue;
							else if(junk == "end_nuclei")
								break;
							else
							{
								Z.push_back(std::stoi(junk));
								input>>a;
								A.push_back(a);
							}
						}
						temp_chain.AddReactor(Z, A, temp_params);
					}
					else if(junk == "begin_target")
					{
						Z.clear();
						stoich.clear();
						input>>thickness;
						while(input>>junk)
						{
							if(junk == "begin_elements")
								continue;
							else if(junk == "end_elements")
								break;
							else
							{
								Z.push_back(std::stoi(junk));
								input>>a;
								stoich.push_back(a);
							}
						}
						temp_chain.SetTarget(Z, stoich, thickness);
					}
					else if(junk == "end_target")
						continue;
					else if(junk == "end_reactor")
						continue;
					else if(junk == "end_reactorchain")
						break;
					else
					{
						std::cerr<<"Bad input file, incorrect config in file, unexpected option "<<junk<<" in reactor chain "<<filename<<std::endl;
						return;
					}
				}
				m_chains.push_back(temp_chain);
			}
			else if(junk == "begin_detectorarray")
			{
				while(input >> junk)
				{
					if(junk == "focalplane")
					{
						double theta, b;
						input>>theta>>b;
						m_array.SetFocalPlane({theta, b});
					}
					else if(junk == "sabre")
					{
						m_array.SetSabre(true);
					}
					else if(junk == "end_detectorarray")
						break;
					else
					{
						std::cerr<<"Bad input file, incorrect config in file, unexpected option in detector array "<<filename<<std::endl;
						return;
					}
				}
			}
			else if(junk == "end_simulator")
				break;
			else
			{
				std::cerr<<"Bad input file, incorrect config in file, unexpected option in simulator "<<filename<<std::endl;
				return;
			}
		}

		m_initFlag = true;

		input.close();
	}

	void Simulator::Run()
	{
		if(!m_initFlag)
		{
			std::cerr<<"ERR -- Simulator not properly initialized!"<<std::endl;
			return;
		}
		
		for(auto& chain : m_chains)
		{
			if(!chain.VerifyChain())
			{
				std::cerr<<"ERR -- Invalid chain with id "<<chain.GetChainID()<<std::endl;
				return;
			}
			chain.BindTarget();
		}

		m_plotter.Open(m_outputFile);
		if(!m_plotter.IsOpen())
		{
			std::cerr<<"ERR -- Unable to open file "<<m_outputFile<<std::endl;
			return;
		}

		for(int i=0; i<m_chains.size(); i++)
			m_pool.PushJob({std::bind(&Simulator::RunChain, std::ref(*this), std::placeholders::_1), i});
		while(true)
		{
			if(m_pool.IsFinished() && m_plotter.GetQueueSize() == 0)
			{
				break;
			}
			else
				m_plotter.PlotData();
		}
		std::cout<<std::endl;
		m_pool.Shutdown();
		std::cout<<"Thread pool shutdown"<<std::endl;
		m_plotter.Close();
		std::cout<<"Data written to file"<<std::endl;
	}

	void Simulator::RunChain(int index)
	{
		uint64_t flush_val = m_samples * 0.1;
		int flush_count = 0;
		int count=0;
		std::vector<ChainResult> results;
		ChainResult this_result;
		results.reserve(m_samples);
		for(uint64_t i=0; i<m_samples; i++)
		{
			count++;
			if(count == flush_val)
			{
				flush_count++;
				count = 0;
				std::cout<<"\rPercent simulated: "<<flush_count*10<<"%"<<std::flush;
			}
			this_result = m_chains[index].GenerateProducts();
			m_array.ProcessData(this_result);
			m_plotter.PushData(this_result);
		}
	}
}
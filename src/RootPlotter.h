#ifndef ROOT_PLOTTER_H
#define ROOT_PLOTTER_H

#include <string>
#include <vector>
#include "ReactorChain.h"
#include <TFile.h>
#include <unordered_map>
#include <mutex>
#include <atomic>
#include <queue>
#include <memory>

namespace NucKage {

	class RootPlotter
	{
	public:
		RootPlotter();
		RootPlotter(const std::string& name);
		~RootPlotter();
		inline bool IsOpen() { return m_openFlag; }
		inline void PushData(const ChainResult& data)
		{
			std::lock_guard<std::mutex> guard(m_rootMutex);
			m_queue.push(data);
			m_queueSize++;
		}

		inline uint64_t GetQueueSize()
		{
			return m_queueSize;
		}
		void PlotData();//const ChainResult& data);
		void PlotData(const ChainResult& data); //For single thread testing
		void Close();
		void Open(const std::string& name);

	private:
		static constexpr double s_rad2deg = 180.0/M_PI;

		struct Histo2DParams
		{
			std::string name;
			int binsX=0;
			int binsY=0;
			double minY=0.0;
			double minX=0.0;
			double maxX=0.0;
			double maxY=0.0;
		};
		struct Histo1DParams
		{
			std::string name;
			int binsX=0;
			double minX=0.0;
			double maxX=0.0;
		};

		inline double FullPhi(double phi) { return phi >= 0.0 ? phi : 2.0*M_PI+phi; }
		inline ChainResult PopData() //Do not decrement queue size here, will cause early exit of program
		{
			std::lock_guard<std::mutex> guard(m_rootMutex);
			auto result = m_queue.front();
			m_queue.pop();
			return result;
		}

		void MyFill2D(const Histo2DParams& params, double valueX, double valueY);
		void MyFill1D(const Histo1DParams& params, double valueX);
		void MyFillGraph(const std::string& name, double valueX, double valueY, int color);

		std::queue<ChainResult, std::deque<ChainResult>> m_queue;

		std::atomic<bool> m_openFlag;
		std::atomic<uint64_t> m_queueSize; //Important! Do not use actual queue size, we want to hold off until the data is processed
		TFile* m_file;
		std::unordered_map<std::string, std::shared_ptr<TObject>> m_map;

		std::mutex m_rootMutex;
	};

}

#endif
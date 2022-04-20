#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <functional>
#include <queue>
#include <condition_variable>


namespace NucKage {

	using JobFunction = std::function<void(int)>;
	struct Job
	{
		JobFunction func;
		int argument;
	};

	class ThreadPool
	{
	public:
		ThreadPool(int nthreads) :
			m_isStopped(false), m_numberRunning(0), m_queueSize(0), m_initShutdown(false)
		{
			for(int i=0; i<nthreads; i++)
			{
				m_pool.push_back(std::thread(std::bind(&ThreadPool::ExecuteJob, std::ref(*this))));
			}
		}
		~ThreadPool()
		{
			if(!m_isStopped)
				Shutdown();
		}

		void PushJob(Job job)
		{
			{
				std::unique_lock<std::mutex> guard(m_poolMutex);
				m_queue.push(job);
				m_queueSize++;
			}
			m_wakeCondition.notify_one();
		}

		bool IsFinished()
		{
			return m_numberRunning == 0 && m_queueSize == 0;
		}

		void Shutdown()
		{
			m_initShutdown = true;
			m_wakeCondition.notify_all();
			for(auto& thread : m_pool)
			{
				thread.join();
			}

			m_pool.clear();
			m_isStopped = true;
		}

	private:
		void ExecuteJob()
		{
			while(true)
			{
				Job job;
				{
					std::unique_lock<std::mutex> guard(m_poolMutex);
					m_wakeCondition.wait(guard, [this](){ 
						return (!m_queue.empty() || m_initShutdown); 
					});
					if(m_initShutdown)
						return;
	
					job = m_queue.front();
					m_queue.pop();
				}

				//Change number running first so that no crash 
				m_numberRunning++;
				m_queueSize--;

				job.func(job.argument);
				m_numberRunning--;
			}
			
		}

		std::queue<Job, std::deque<Job>> m_queue;
		std::vector<std::thread> m_pool;
		std::mutex m_poolMutex;
		std::condition_variable m_wakeCondition;

		bool m_isStopped;
		std::atomic<int> m_numberRunning;
		std::atomic<int> m_queueSize;
		std::atomic<bool> m_initShutdown;
	};
}

#endif
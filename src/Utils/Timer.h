#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace NucKage {

	class Timer
	{
	public:
		Timer(const char* name) :
			m_name(name)
		{
			Restart();
		}

		~Timer()
		{
		}

		void Restart()
		{
			m_startTime = Clock::now();
		}


		float ElapsedMilliseconds()
		{
			float duration = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - m_startTime).count()*1.0e-6;
			return duration;
		}

	private:
		using Time = std::chrono::steady_clock::time_point;
		using Clock = std::chrono::steady_clock;

		const char* m_name;
		Time m_startTime;
	};
}

#endif
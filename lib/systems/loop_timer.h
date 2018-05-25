#ifndef LOOPTIMER_H
#define LOOPTIMER_H
#include <chrono>

class LoopTimer {
public:
  void loop_start();
  void loop_pause();
  void loop_end();
  double average_loop_period();
  double loop_period();

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> t0;
  std::chrono::duration<double> dt_loop;
  double dt_average = 0;
  int N = 0;
  bool loop_start_ = false;
  bool loop_pause_ = false;
};

#endif // LOOPTIMER_H

#include "loop_timer.h"
#include <stdexcept>

void LoopTimer::loop_start() {
  if (loop_start_ && !loop_pause_) {
    throw std::runtime_error("Cannot start loop without ending the loop");
  }
  if (!loop_pause_) {
    dt_loop = std::chrono::duration<double>::zero();
  }
  loop_start_ = true;
  loop_pause_ = false;
  t0 = std::chrono::high_resolution_clock::now();
}

void LoopTimer::loop_pause() {
  if (loop_pause_) {
    throw std::runtime_error("Loop already paused");
  }
  if (!loop_start_) {
    throw std::runtime_error("Loop did not start");
  }
  dt_loop = std::chrono::duration<double>(
      std::chrono::high_resolution_clock::now() - t0);
  loop_pause_ = true;
}

void LoopTimer::loop_end() {
  if (!loop_start_) {
    throw std::runtime_error("Loop did not start");
  }
  dt_loop = dt_loop + std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - t0);
  if (N == 0) {
    dt_average = dt_loop.count();
    N = 1;
  } else {
    dt_average = (N * dt_average + dt_loop.count()) / (++N);
  }
  loop_pause_ = false;
  loop_start_ = false;
}

double LoopTimer::average_loop_period() { return dt_average; }

double LoopTimer::loop_period() { return dt_loop.count(); }

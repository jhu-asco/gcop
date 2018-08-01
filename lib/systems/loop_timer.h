#ifndef LOOPTIMER_H
#define LOOPTIMER_H
#include <chrono>

/**
 * @brief The LoopTimer class
 *
 * Helper class to find the time taken in a loop
 */
class LoopTimer {
public:
  /**
   * @brief loop_start
   *
   * Starts a timer. Resets internal t0
   * if the timer is called after loop_end
   * If timer is called after loop_pause, does not reset t0
   * If neither loop_pause or loop_end are called, the
   * function throw error. To begin with the state is
   * assumed to be in loop_end state.
   */
  void loop_start();
  /**
   * @brief loop_pause
   *
   * Temporarily pause the loop timer to stop time counting
   */
  void loop_pause();
  /**
   * @brief loop_end
   *
   * End the loop timer and update the internal average loop time
   */
  void loop_end();
  /**
   * @brief average_loop_period
   * @return  The average time period of a loop in seconds
   */
  double average_loop_period();
  /**
   * @brief loop_period
   * @return  The current loop period in seconds
   */
  double loop_period();

private:
  /**
   * @brief t0 The start time for a loop
   */
  std::chrono::time_point<std::chrono::high_resolution_clock> t0;
  /**
   * @brief dt_loop The time elapsed in current loop
   */
  std::chrono::duration<double> dt_loop;
  /**
   * @brief dt_average The average time taken in a loop
   */
  double dt_average = 0;
  /**
   * @brief N Number of loops
   */
  int N = 0;
  /**
   * @brief loop_start_
   *
   * Flag to indicate whether loop started
   */
  bool loop_start_ = false;
  /**
   * @brief loop_pause_
   *
   * Flag to indicate whether loop paused
   */
  bool loop_pause_ = false;
};

#endif // LOOPTIMER_H

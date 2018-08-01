#include "loop_timer.h"
#include "gtest/gtest.h"


TEST(LoopTimer, simple_loop) {
  LoopTimer loop_timer;
  for(int i = 0; i < 20; ++i) {
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_end();
  }
  ASSERT_NEAR(loop_timer.average_loop_period(), 0.005, 1e-3);
}

TEST(LoopTimer, single_pause) {
  LoopTimer loop_timer;
  for(int i = 0; i < 20; ++i) {
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_pause();
    usleep(5000);
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_end();
  }
  ASSERT_NEAR(loop_timer.average_loop_period(), 0.01, 1e-3);
}

TEST(LoopTimer, multi_pause) {
  LoopTimer loop_timer;
  for(int i = 0; i < 20; ++i) {
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_pause();
    usleep(5000); // Not counted
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_pause();
    usleep(5000);//Not counted
    loop_timer.loop_start();
    usleep(5000);
    loop_timer.loop_end();
  }
  ASSERT_NEAR(loop_timer.average_loop_period(), 0.015, 1e-3);
}

TEST(LoopTimer, start_twice){
  LoopTimer loop_timer;
  loop_timer.loop_start();
  ASSERT_THROW(loop_timer.loop_start(), std::runtime_error);
}

TEST(LoopTimer, pause_twice){
  LoopTimer loop_timer;
  loop_timer.loop_start();
  loop_timer.loop_pause();
  ASSERT_THROW(loop_timer.loop_pause(), std::runtime_error);
}

TEST(LoopTimer, end_before_start){
  LoopTimer loop_timer;
  ASSERT_THROW(loop_timer.loop_end(), std::runtime_error);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

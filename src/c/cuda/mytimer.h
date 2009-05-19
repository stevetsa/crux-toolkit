#include <time.h>
struct my_timer {
  long start_counts;
  long stop_counts;
  struct timespec start;
  struct timespec stop;
  struct timespec total;


};
#ifdef CUDA_NVCC
//extern "C" {
#endif
void my_timer_reset(struct my_timer* timer);
void my_timer_start(struct my_timer* timer);
void my_timer_stop(struct my_timer* timer);
void my_timer_print(struct my_timer* timer);
struct my_timer* new_my_timer();

void diff(struct timespec start, struct timespec end, struct timespec* ans);
void add(struct timespec* ans, struct timespec add);
#ifdef CUDA_NVCC
//}
#endif

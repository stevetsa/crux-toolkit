#include "mytimer.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NANOSECOND 1000000000

void my_timer_reset(struct my_timer* timer) {
  timer -> start_counts = 0;
  timer -> stop_counts = 0;
  memset(&(timer -> start),0,sizeof(struct timespec));
  memset(&(timer -> stop),0,sizeof(struct timespec));
  memset(&(timer -> total),0,sizeof(struct timespec));


}
void my_timer_start(struct my_timer* timer) {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&(timer -> start));
  timer -> start_counts++;
}
void my_timer_stop(struct my_timer* timer) {
  struct timespec temp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&(timer -> stop));
  timer -> stop_counts++;
  diff(timer -> start, timer -> stop, &temp);
  //printf("temp.sec:%ld\n",temp.tv_sec);
  //printf("temp.nsec:%ld\n",temp.tv_nsec);


  add(&(timer -> total), temp);

  //printf("total.sec:%ld\n",timer -> total.tv_sec);
  //printf("total.nsec:%ld\n", timer -> total.tv_nsec);


}
void my_timer_print(struct my_timer* timer) {
  printf("Number of start:%ld\n",timer -> start_counts);
  printf("Number of stop:%ld\n", timer -> stop_counts);
  printf("Total Sec:%ld\n",timer -> total.tv_sec);
  printf("Total NSec:%ld\n", timer -> total.tv_nsec);
  printf("Average Sec:%lf\n", (double)timer -> total.tv_sec / (double)timer -> start_counts);
  printf("Average NSec:%lf\n",(double)timer -> total.tv_nsec / (double)timer -> start_counts);
}

struct my_timer* new_my_timer() {
  struct my_timer* ans = malloc(sizeof(struct my_timer));
  my_timer_reset(ans);
  return ans;
}

void diff(struct timespec start, struct timespec end, struct timespec* ans) {
  if ((end.tv_nsec - start.tv_nsec) < 0) {
    ans -> tv_sec = end.tv_sec - start.tv_sec - 1;
    ans -> tv_nsec = NANOSECOND + end.tv_nsec - start.tv_nsec;
  }
  else
    {
      ans -> tv_sec = end.tv_sec - start.tv_sec;
      ans -> tv_nsec = end.tv_nsec - start.tv_nsec;
    }
}

void add(struct timespec* ans, struct timespec add) {
  ans -> tv_nsec += add.tv_nsec;
  if (ans -> tv_nsec > NANOSECOND) {
    ans -> tv_nsec -= NANOSECOND;
    ans -> tv_sec ++;
  }
  ans -> tv_sec += add.tv_sec;

}

/*
int main(void) {
  int temp;

  int i;
  struct my_timer timer;
  struct my_timer timer2;
  my_timer_reset(&timer);
  my_timer_reset(&timer2);
  my_timer_start(&timer);
  for (i=0;i<1000;i++) {
    my_timer_start(&timer2);
    temp += temp;
    my_timer_stop(&timer2);
  }
  my_timer_stop(&timer);


  printf("Timer1:\n");
  my_timer_print(&timer);
  printf("Timer2:\n");
  my_timer_print(&timer2);
}
*/

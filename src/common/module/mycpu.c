//#define _GNU_SOURCE

#include <sched.h>

/*
 * Find the core the thread belongs to
 */

int mycpu_ ()
{
  /* int sched_getcpu(void); */
  int cpu;
  cpu = sched_getcpu();
  return cpu;
}
int mycpu() { return mycpu_(); }

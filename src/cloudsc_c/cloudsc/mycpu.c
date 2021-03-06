/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

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

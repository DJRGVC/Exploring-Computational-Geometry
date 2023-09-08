/*
Laura Toma
*/



#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "rtimer.h"

char *
rt_sprint(char *buf, Rtimer rt) {
  if(rt_w_useconds(rt) == 0) {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    rt_u_useconds(rt)/1000000,
	    100.0*rt_u_useconds(rt)/rt_w_useconds(rt),
	    rt_s_useconds(rt)/1000000,
	    100.0*rt_s_useconds(rt)/rt_w_useconds(rt),
	    rt_w_useconds(rt)/1000000,
	    100.0*(rt_u_useconds(rt)+rt_s_useconds(rt)) / rt_w_useconds(rt));
  }
  return buf;
}



/* prints the timer divided by k --- used when k experiments kept in
   the same timer */
char *
rt_sprint_average(char *buf, Rtimer rt, int k) {

  assert(k>0);
  if(rt_w_useconds(rt) == 0) {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    rt_u_useconds(rt)/1000000/k,
	    100.0*rt_u_useconds(rt)/rt_w_useconds(rt),
	    rt_s_useconds(rt)/1000000/k,
	    100.0*rt_s_useconds(rt)/rt_w_useconds(rt),
	    rt_w_useconds(rt)/1000000/k,
	    100.0*(rt_u_useconds(rt)+rt_s_useconds(rt)) / rt_w_useconds(rt));
  }
  return buf;
}


/* to be called after rt_stop_and_accumulate to print the total time */
/* assumes rt.tu_usec, rt.ts_usec, rt.tw_usec have been computed
 */
char *
rt_sprint_total(char *buf, Rtimer rt) {
  if(rt.tw_usec == 0) {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    sprintf(buf, "[%4.2fu (%.0f%%) %4.2fs (%.0f%%) %4.2f %.1f%%]",
	    rt.tu_usec/1000000,
	    100.0*rt.tu_usec/rt.tw_usec,
	    rt.ts_usec/1000000,
	    100.0*rt.ts_usec/rt.tw_usec,
	    rt.tw_usec/1000000,
	    100.0*(rt.tu_usec+rt.ts_usec) / rt.tw_usec);
  }
  return buf;
}

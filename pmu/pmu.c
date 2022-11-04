/*
 * Performance Monitor Unit (PMU)
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * libpmu is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * libpmu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#include "pmu.h"
#include <stdarg.h>

#ifdef __cplusplus
extern "C"{
#endif

/********************************
 *    Performance Evaluation    *
 ********************************/

static uint64_t cpucycles_overhead(void) {
  uint64_t t0, t1, overhead = -1LL;
  unsigned int i;
  for(i=0;i<100000;i++) {
    t0 = rdtsc();
    __asm__ volatile ("");
    t1 = rdtsc();
    if(t1 - t0 < overhead)
      overhead = t1 - t0;
  }
  return overhead;
}

static int cmp_uint64(const void *a, const void *b) {
  if(*(uint64_t *)a < *(uint64_t *)b) return -1;
  if(*(uint64_t *)a > *(uint64_t *)b) return  1;
  return 0;
}

static uint64_t median(uint64_t *t, size_t ntests) {
  qsort(t,ntests,sizeof(uint64_t),cmp_uint64);
  if(ntests&1)
    return t[ntests/2];
  else
    return (t[ntests/2-1]+t[ntests/2])/2;
}

static uint64_t average(uint64_t *t, size_t ntests) {
  size_t i;
  uint64_t acc=0;
  for(i=0;i<ntests;i++)
    acc += t[i];
  return acc/ntests;
}

static void strf_timediff(char *__restrict s, unsigned long diff)
{
  const long int K=1000, M=1000000, G=1000000000;
  time_t  sec = diff/G;
  time_t nsec = diff-sec*G;
  if (sec==0)
    goto zerosec;
  else
    goto nonzerosec;
zerosec:
  if (nsec < K)
    sprintf(s, "%lins", nsec);
  else if (nsec < M) /* 1e3 <= diff < 1e6 */
    sprintf(s, "%lius %03lins", nsec/K, nsec%K);
  else               /* 1e6 <= diff < 1e9 */
    sprintf(s, "%lims %03lius %03lins", nsec/M, (nsec%M)/K, nsec%K);
  return;
nonzerosec:
  if (sec < 60)
    sprintf(s, "%lis %03lims %03lius %03lins", sec, nsec/M, (nsec%M)/K, nsec%K);
  else if (sec < 3600)
    sprintf(s, "%lim %02lis %03lims %03lius %03lins",
      sec/60, sec%60, nsec/M, (nsec%M)/K, nsec%K);
  else
    sprintf(s, "%lih %02lim %02lis %03lims %03lius %03lins",
      sec/3600, (sec%3600)/60, sec%60, nsec/M, (nsec%M)/K, nsec%K);
  return;
}

/**
 * print_cpuinfo - Print CPU info
 * 
 * The CPU info refers to the CPU cycles per tick.
 * 
 * @s: title
 * @cpu: cpu cycles recoder
 * @ntests: number of tests
 */
void print_cpuinfo(const char *__restrict s, uint64_t *cpu, size_t ntests)
{
  if(ntests < 2) {
    fprintf(stderr, "ERROR: Need a least two cycle counts!\n");
    return;
  }

  static uint64_t overhead = -1;
  if(overhead  == (uint64_t)-1)
    overhead = cpucycles_overhead();

  ntests--;
  for(size_t i=0; i<ntests; ++i)
    cpu[i] = cpu[i+1] - cpu[i] - overhead;

  printf("\033[4mcpu info\033[0m %s:\n", s);
  printf("\\_ freq:    %f GHz\n",rdcpufreq());
  printf("\\_ median:  %lu cycles/ticks\n", median(cpu, ntests));
  printf("\\_ average: %lu cycles/ticks\n", average(cpu, ntests));
}

/**
 * print_clkinfo - Print clock info (in nanoseconds)
 * 
 * The time info should be obtained by 
 * - timespec_get(&clk, TIME_UTC), or
 * - clock_gettime(CLOCK_MONOTONIC, &clk).
 * This uses <time.h>.
 * 
 * @s: title
 * @clk: clock recoder
 * @ntests: number of tests
 */
void print_clkinfo(const char *__restrict s, struct timespec *clk, size_t ntests)
{
  if(ntests < 2) {
    fprintf(stderr, "ERROR: Need a least two cycle counts!\n");
    return;
  }

  ntests--;
  unsigned long diff[ntests];
  for (size_t i=0; i<ntests; i++){
    time_t  sec = clk[i+1].tv_sec  - clk[i].tv_sec;
    time_t nsec = clk[i+1].tv_nsec - clk[i].tv_nsec;
    const unsigned long G=1000000000;
    if (nsec < 0) {
      sec--;
      nsec += G;
    }
    diff[i] = sec*G + nsec;
  }

  printf("\033[4mclock info\033[0m %s:\n", s);

  char str[80];
  strf_timediff(str, median(diff, ntests));
  printf("\\_ median:  %s\n", str);
  strf_timediff(str, average(diff, ntests));
  printf("\\_ average: %s\n", str);
}

/**
 * print_meminfo - Print memory info
 * 
 * The memory info refers to the maximum resident set size (rss) in kilobyte.
 * This uses <sys/resource.h>.
 * 
 * @s: title
 * @mem: memory recoder
 * @ntests: number of tests
 */
void print_meminfo(const char *__restrict s, uint64_t *mem, size_t ntests)
{
  if(ntests < 2) {
    fprintf(stderr, "ERROR: Need a least two cycle counts!\n");
    return;
  }

  ntests--;
  for(size_t i=0; i<ntests; ++i)
    mem[i] = mem[i+1] - mem[i];

  const uint64_t octkilo = 8*1024;

  printf("\033[4mmemory info\033[0m %s:\n", s);
  uint64_t m_mid = median(mem, ntests); /* memory median */
  if (m_mid<octkilo)
    printf("\\_ median:  %lu Kb\n", m_mid);
  else
    printf("\\_ median:  %lu Kb = %f Mb\n", m_mid, (double)m_mid/octkilo);
  uint64_t m_aver = average(mem, ntests); /* memory average */
  if (m_aver<octkilo)
    printf("\\_ average: %lu Kb\n", m_aver);
  else
    printf("\\_ average: %lu Kb = %f Mb\n", m_aver, (double)m_aver/octkilo);
}


/*******************
 *    Unit Test    *
 *******************/

struct procinfo {
  char func_name[80];
  bool fail_flag;
  size_t task_count;
  size_t fail_count;
  struct timespec time_test[2];
  struct timespec time_case[2];
  struct rusage rusage_test[2];
  struct rusage rusage_case[2];
};

/** process status */
static struct procinfo ps;

/** 
 * strfps - Format timespec and rusage into string.
 * 
 * @str: the returned string that store the info.
 * @timespec: {tv_sec. tv_nsec}.
 * @rusage: need ru_maxrss (maximum resident set size) (in kilobytes).
 * 
 * @return: Formatted string of the form "01m 04s 763.911ms 13.024Mb used"
 */
static void strfps(char *str, struct timespec *ts, struct rusage *ru)
{
  /* time diff */
  char str_time[80];
  time_t  s = ts[1].tv_sec  - ts[0].tv_sec;
  time_t ns = ts[1].tv_nsec - ts[0].tv_nsec;
  const unsigned long G=1000000000;
  if (ns < 0) {
    s--;
    ns += G;
  }
  strf_timediff(str_time, s*G+ns);
  /* maxrss diff */
  char str_maxrss[80];
  double diff_maxrss = ru[1].ru_maxrss - ru[0].ru_maxrss;
  const int octkilo = 8*1024;
  if (diff_maxrss<octkilo)
    sprintf(str_maxrss, "%.0fKb", diff_maxrss);
  else
    sprintf(str_maxrss, "%.3fMb", diff_maxrss/octkilo);
  sprintf(str, "%s %s used", str_time, str_maxrss);
}

/** 
 * pmu_test_begin - Start to track the ps of a test task.
 * 
 * @func_name: The name the current test function: __func__
 * @usage: current resource usage (will use user time + max rss)
 */
void pmu_test_begin(const char *__restrict func_name)
{
  timespec_get(&ps.time_test[0], TIME_UTC);
  getrusage(RUSAGE_SELF, &ps.rusage_test[0]);
  char hms[9]; /* 8 chars for "hh:mm:ss", and 1 for '\0' */
  strftime(hms, 9, "%T", localtime(&ps.time_test[0].tv_sec));
  memset(ps.func_name, 0x00, 80);
  strcpy(ps.func_name, func_name);
  ps.fail_flag = 0;
  ps.fail_count = 0;
  ps.task_count = 0;
  printf("\033[32m[===============]\033[0m\n");
  printf("\033[32m[%s BEGIN ]\033[0m \033[1m%s\033[0m\n", hms, ps.func_name);
  printf("[---------------]\n");
}

/** 
 * pmu_test_do - Start to track the ps of an individual test case.
 * 
 * @case_name: user specified test case name
 * @usage: current resource usage (will use user time + max rss)
 */
void pmu_test_do(const char *__restrict case_name)
{
  timespec_get(&ps.time_case[0], TIME_UTC);
  getrusage(RUSAGE_SELF, &ps.rusage_case[0]);
  char hms[9]; /* 8 chars for "hh:mm:ss", and 1 for '\0' */
  strftime(hms, 9, "%T", localtime(&ps.time_case[0].tv_sec));
  ps.fail_flag = 0;
  ps.task_count++;
  printf("\033[32m[%s CASE  ]\033[0m %s\n", hms, case_name);
}

void pmu_test_placeholder()
{
  printf("\033[32m[               ]\033[0m ");
}

/** 
 * pmu_test_done - Use when an individual test case is done.
 * 
 * @usage: current resource usage (will use user time + max rss)
 */
void pmu_test_done()
{
  timespec_get(&ps.time_case[1], TIME_UTC);
  getrusage(RUSAGE_SELF, &ps.rusage_case[1]);
  char hms[9]; /* 8 chars for "hh:mm:ss", and 1 for '\0' */
  strftime(hms, 9, "%T", localtime(&ps.time_case[1].tv_sec));
  char str[BUFSIZ]={0};
  strfps(str, ps.time_case, ps.rusage_case);
  if (ps.fail_flag) {
    printf("\033[31m[%s FAILED]\033[0m %s.\n", hms, str);
    ps.fail_count++;
  } else {
    printf("\033[32m[%s     OK]\033[0m %s.\n", hms, str);
  }
  printf("[---------------]\n");
}

/** 
 * pmu_test_end - Use when a test task is finished.
 * 
 * @usage: current resource usage (will use user time + max rss)
 */
void pmu_test_end()
{
  timespec_get(&ps.time_test[1], TIME_UTC);
  getrusage(RUSAGE_SELF, &ps.rusage_test[1]);
  char hms[9]; /* 8 chars for "hh:mm:ss", and 1 for '\0' */
  strftime(hms, 9, "%T", localtime(&ps.time_test[1].tv_sec));
  char str[BUFSIZ]={0};
  strfps(str, ps.time_test, ps.rusage_test);
  if (ps.fail_count) {
    printf("\033[31m[%s FAILED]\033[0m \033[1m%s\033[0m: "
      "%lu/%lu subtask(s) failed (%s).\n", hms, ps.func_name,
      ps.fail_count, ps.task_count, str);
    printf("\033[31m[===============]\033[0m\n\n");
  } else {
    printf("\033[32m[%s PASSED]\033[0m \033[1m%s\033[0m: "
      "%lu/%lu subtask(s) passed (%s).\n", hms, ps.func_name,
      ps.task_count, ps.task_count, str);
    printf("\033[32m[===============]\033[0m\n\n");
  }
}

int pmu_test_reterr()
{
  if (ps.fail_count)
    return 1;
  else
    return 0;
}

void pmu_test_pile(const char *__restrict file,
  const char *__restrict func, const int line)
{
#if 0
  char *note;
  va_list ap;
  va_start (ap, s) ;
  *s++;
  note = va_arg(ap, char *);
  //sprintf(note, s, ap);
  //vfprintf (stdout, s, ap);
  printf("%s\n",note);
  va_end (ap);
#endif
  printf("\033[33mPile\033[0m in %s:\033[1m%s:%i\033[0m.\n", file, func, line);
}


/*********************
 *    Assertation    *
 *********************/

void pmu_assert_fail(const char *__restrict expr,
  const char *__restrict file, const char *__restrict func, const int line)
{
  ps.fail_flag=1;
  printf("\033[31mFail\033[0m in %s:\033[1m%s:%i\033[0m: Assertion `%s`.\n",
    file, func, line, expr);
}

#ifdef __cplusplus
}
#endif

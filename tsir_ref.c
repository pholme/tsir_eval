// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018/2020)

// this code uses a standard approach (going throught the contacts in time
// order)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdint.h>

#define NAVG 100000 // number of runs for averages

#define S(x) (g.rtime[(x)] == UINT_MAX)
#define I(x) ((g.rtime[(x)] > now) && (g.rtime[(x)] < UINT_MAX))
#define R(x) (g.rtime[(x)] <= now)

// auxiliary macro
#define SQ(x) ((x) * (x))

typedef struct GLOBALS {
	// INPUT PARAMETERS
	double a; // -1/nu where nu = recovery rate (input nu in units of duration, but internally in units of time steps)
	unsigned int ubeta;
	// RECOVERY TIMES
	unsigned int *rtime;
	// NETWORK SPECS
	unsigned int n, nc, dur;
	// OUTBREAK STATS
	unsigned int ni;
	// FOR RNG
} GLOBALS;

typedef struct CONTACT {
	unsigned int left, right, t;
} CONTACT;

// misc.c

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();

GLOBALS g;
CONTACT *c;
uint64_t pcg_state;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// giving exponential random numbers with a mean reciprocal of the recovery rate

unsigned int exptime () {
	uint32_t r = pcg_32();

	if (r == 4294967295u) return 0;

	return g.a * log((r + 1) / 4294967296.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// using the bisection method to find the first contact

unsigned int start_index (unsigned int now) {
	unsigned int lo = 0, mid, hi = g.nc - 1;

	while (lo < hi) {
		mid = (lo + hi) >> 1;
		if (c[mid].t > now) hi = mid;
		else lo = mid + 1;
	}

	return lo;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random starting node

void sir () {
	unsigned int i, me, you, now;
	unsigned int t0, tx; // source time, extinction time
	
	// initialize
	for (i = 0; i < g.n; i++) g.rtime[i] = UINT_MAX;

	// get & infect the source
	t0 = pcg_32_bounded(g.dur);
	g.rtime[pcg_32_bounded(g.n)] = tx = t0 + exptime();
	g.ni = 1;

	// run the outbreak
	for (i = start_index(t0); i < g.nc; i++) {
		me = c[i].left;
		you = c[i].right;
		now = c[i].t;
		if (now > tx) break;
		if (I(me) && S(you)) {
			if (pcg_32() < g.ubeta) { // infecting
				g.ni++;
				g.rtime[you] = now + exptime();
				if (tx < g.rtime[you]) tx = g.rtime[you];
			}
		} else if (S(me) && I(you)) {
			if (pcg_32() < g.ubeta) { // infecting
				g.ni++;
				g.rtime[me] = now + exptime();
				if (tx < g.rtime[me]) tx = g.rtime[me];
			}
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling i/o

int main (int argc, char *argv[]) {
	unsigned int i, a[NAVG];
	FILE *fp;
	struct timespec t0, t1;
	
	// initialize parameters
	g.ubeta = UINT_MAX * atof(argv[1]);
	
	// read the network
	if (3 != scanf("%u %u %u\n", &g.n, &g.nc, &g.dur)) {
		fprintf(stderr, "input error 1\n");
		exit(1);
	}

	c = malloc(g.nc * sizeof(CONTACT));
	
	for (i = 0; i < g.nc; i++) {
		if (3 != scanf("%u %u %u\n", &c[i].left, &c[i].right, &c[i].t)) {
			fprintf(stderr, "input error 2\n");
			exit(2);
		}
	}

	pcg_state = strtoull(argv[3], NULL, 10); // argv[3] is the RNG state

	g.a = -(g.dur / atof(argv[2])); // argv[4] is nu in units of the duration of the data

	g.rtime = malloc(g.n * sizeof(unsigned int));

	// run the simulations and summing for averages
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
	for (i = 0; i < NAVG; i++) {
		sir();
		a[i] = g.ni;
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);

	// print result
	printf("%g\n", ((t1.tv_sec - t0.tv_sec) + 1.0e-9 * (t1.tv_nsec - t0.tv_nsec)) / NAVG);
	for (i = 0; i < NAVG; i++) printf("%u\n", a[i]);
	
	// cleaning up
	free(g.rtime); free(c);
	 
	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

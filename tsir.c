// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on temporal networks by Petter Holme (2018/2020)

#include "tsir.h"

GLOBALS g;
NODE *n;
uint64_t pcg_state;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine first localizes the first contact later than 'now' in t
// then picks the contact that can infect (a chain of bernouilli trials)
// among the rest of the contacts. It returns the time of the infecting contact

unsigned int contagious_contact (unsigned int *t, unsigned int nt, unsigned int now) {
	unsigned int lo = 0, mid, hi = nt - 1;

	if (t[hi] <= now) return END; // no need to search further bcoz t is sorted. Note that the bisection search depends on this line.

	// the actual bisection search
	while (lo < hi) {
		mid = (lo + hi) >> 1; //  (lo & hi) + ((lo ^ hi) >> 1) would make it possible to go to max 2^32 lo and hi
		if (t[mid] > now) hi = mid;
		else lo = mid + 1;
	}

	// get a random contact
	hi += g.rnd2inx[pcg_16()];

	if (hi >= nt) return NONE; // if the contact is too late, skip it

	// return the time of the contact
	return t[hi];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine does the book keeping for an infection event

void infect () {
	unsigned int i, you, t, me = g.heap[1];
	unsigned int now = n[me].time, duration = exptime();

	del_root(); // take the newly infected off the heap

	if (duration > 0) { // if the duration is zero, no one else can be infected
		n[me].time += duration;

		// go through the neighbors of the infected node . .
		for (i = 0; i < n[me].deg; i++) {
			you = n[me].nb[i];
			if (S(you)) { // if you is S, you can be infected
				// find the infection time of you
				t = contagious_contact(n[me].t[i], n[me].nc[i], now);
				if (t == END) break; // bcoz the sorting of nbs, we can break

				// if the infection time is before when me gets infected,
				// and (if it was already listed for infection) before the
				// previously listed infection event, then list it
				if ((t <= n[me].time) && (t < n[you].time)) {
					n[you].time = t; // set you's infection time
					if (n[you].heap == NONE) { // if not listed before, then extend the heap
						g.heap[++g.nheap] = you;
						n[you].heap = g.nheap;
					}
					up_heap(n[you].heap); // this works bcoz the only heap relationship that can be violated is the one between you and its parent
				}
			}
		}
	}

	g.i[g.ni++] = me; // to get the outbreak size
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random starting node

void sir () {
	unsigned int i, source;
	
	g.ni = 0;
	
	// get & infect the source

	source = pcg_32_bounded(g.n);
	n[source].time = pcg_32_bounded(g.dur);
	n[source].heap = 1;
	g.heap[g.nheap = 1] = source;

	// run the outbreak
	while (g.nheap) infect();

	// clean
	for (i = 0; i < g.ni; i++) n[g.i[i]].heap = n[g.i[i]].time = NONE;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling i/o

int main (int argc, char *argv[]) {
	unsigned int i, j, a[NAVG];
	double d, x;
	struct timespec t0, t1;
	
	pcg_state = strtoull(argv[3], NULL, 10); // argv[3] is the RNG state

	// read network
	read_data();

	// initialize parameters
	d = atof(argv[1]);
	if (d < 1.0) {
		d = 1.0 / log(1.0 - d); // argv[1] is beta
		for (i = 0; i < 0x10000; i++) {
			x = d * log((i + 1) / 65536.0);
			g.rnd2inx[i] = (x > USHRT_MAX) ? USHRT_MAX : x;
		}
	} else for (i = 0; i < 0x10000; i++) g.rnd2inx[i] = 0;

	g.a = -(g.dur / atof(argv[2])); // argv[2] is nu in units of the duration of the data

	// allocating the heap (N + 1) because it's indices are 1,...,N
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));
	g.i = calloc(g.n, sizeof(unsigned int));

	// initialize
	for (i = 0; i < g.n; i++) n[i].heap = n[i].time = NONE;
	
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
	for (i = 0; i < g.n; i++) {
		for (j = 0; j < n[i].deg; j++) free(n[i].t[j]);
		free(n[i].nb);
		free(n[i].nc);
		free(n[i].t);
	}
	free(n); free(g.heap); free(g.i);
	 
	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

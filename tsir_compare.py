# Code for comparing the event-driven and straightforward simulations of SIR on temporal networks

from sys import argv
from random import getrandbits
from subprocess import Popen, PIPE, STDOUT
from scipy.stats import mannwhitneyu
import networkx as nx
import numpy as np

nrun = 10

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
# can s be converted to a floating point number?

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
# creating the network representation for the straightforward method

def txt_ref (a,G):

	nwk = str(G.number_of_nodes()) + ' ' + str(a.shape[0]) + ' ' + str(a[-1,2]) + '\n'
	for i in range(a.shape[0]):
		nwk += str(a[i,0]) + ' ' + str(a[i,1]) + ' ' + str(a[i,2]) + '\n'

	return nwk

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
# creating the network representation for the event-driven method

def txt (a,G):

	tdic = {}
	
	for i in range(a.shape[0]):
		u = a[i,0]
		v = a[i,1]
		tme = a[i,2]
		if u < v:
			e = (u,v)
		else:
			e = (v,u)

		if e in tdic:
			tdic[e].append(tme)
		else:
			tdic[e] = [tme]

	ntxt = [''] * G.number_of_nodes()
	tmax = a[-1,2]

	deg = [0] * G.number_of_nodes()

	for (e,t) in sorted(tdic.items(), key=lambda x: (-x[1][-1],len(x[1]))):
		v = e[0]
		deg[v] += 1
		ntxt[v] += str(e[1]) + ' ' + str(len(t)) + '\n'
		for tme in t:
			ntxt[v] += str(tme) + '\n'

		v = e[1]
		deg[v] += 1
		ntxt[v] += str(e[0]) + ' ' + str(len(t)) + '\n'
		for tme in t:
			ntxt[v] += str(tme) + '\n'

	nwk = str(G.number_of_nodes()) + ' ' + str(tmax) + '\n'
	for i in range(G.number_of_nodes()):
		nwk += str(deg[i]) + '\n' + ntxt[i]

	return nwk

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
# generates a temporal network

def gen_tn (n, z, c):

	G = nx.fast_gnp_random_graph(n,z/n)

	links = [(u,v) for (u,v) in G.edges()]

	nl = G.number_of_edges()

	contacts = []

	for i in range(nl):
		tnow = np.random.exponential()
		while tnow < c:
			contacts.append((tnow,links[i][0],links[i][1]))
			tnow += np.random.exponential()

	contacts.sort()

	scalefac = (2**31 - 1.0) / contacts[-1][0]

	A = np.empty((len(contacts),3),dtype=np.int)

	for i in range(A.shape[0]):
		A[i,0] = contacts[i][1]
		A[i,1] = contacts[i][2]
		A[i,2] = int(contacts[i][0] * scalefac)

	return A,G
			
##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
# running the C code

def run_sim (cmnd, nwk, beta, nu):

	p = Popen([cmnd, str(beta), str(nu), str(getrandbits(64))], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
	o = p.communicate(input=bytes(nwk,encoding='utf-8'))[0].decode()

	# interpret and report the output
	a = o.split('\n')
	if not is_number(a[0]):
		print(o)
		exit(1)
		
	b = np.array(a[1:-1])

	return float(a[0]), b.astype(np.int)

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ## }

if __name__ == '__main__':

	if len(argv) != 6:
		print('usage: python3 tsir_compare.py [n] [z] [c] [beta] [nu]')
		exit(1)

	n = int(argv[1])
	z = float(argv[2])
	c = float(argv[3])
	beta = float(argv[4])
	nu = float(argv[5])

	st0 = 0.0
	st1 = 0.0
	v0 = np.array([])
	v1 = np.array([])

	for i in range(nrun):
		A,G = gen_tn(n, z, c)

		t1,s1 = run_sim('./tsir', txt(A,G), beta, nu)
		st1 += t1
		v1 = np.concatenate((v1,s1))

		t0,s0 = run_sim('./tsir_ref', txt_ref(A,G), beta, nu)
		st0 += t0
		v0 = np.concatenate((v0,s0))

	
	st, p = mannwhitneyu(s0,s1)

	print('Relative speed-up:', st0/st1)
	print('P-value of the Mann-Whitney U test:', p)
	print('Mean outbreak size:', np.mean(np.concatenate((v0,v1))))

##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##

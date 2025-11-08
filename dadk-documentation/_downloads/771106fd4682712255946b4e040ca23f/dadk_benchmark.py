import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from memory_profiler import memory_usage
import networkx as nx
from collections import defaultdict
import itertools

from pyqubo import Array, Constraint, Placeholder
from dadk.BinPol import *
import dimod

###############################################################################################################
# DWAVE SDK: based on https://github.com/dwavesystems/dwave-networkx/blob/main/dwave_networkx/algorithms/tsp.py
###############################################################################################################

def tsp_dwave(n_city, lagrange=None, weight='weight', missing_edge_weight=None):
    """Return the QUBO with ground states corresponding to a minimum TSP route.

    If :math:`|G|` is the number of nodes in the graph, the resulting qubo will have:

    * :math:`|G|^2` variables/nodes
    * :math:`2 |G|^2 (|G| - 1)` interactions/edges

    Parameters
    ----------
    G : NetworkX graph
        A complete graph in which each edge has a attribute giving its weight.

    lagrange : number, optional (default None)
        Lagrange parameter to weight constraints (no edges within set)
        versus objective (largest set possible).

    weight : optional (default 'weight')
        The name of the edge attribute containing the weight.
    
    missing_edge_weight : number, optional (default None)
        For bi-directional graphs, the weight given to missing edges.
        If None is given (the default), missing edges will be set to
        the sum of all weights.

    Returns
    -------
    QUBO : dict
       The QUBO with ground states corresponding to a minimum travelling
       salesperson route. The QUBO variables are labelled `(c, t)` where `c`
       is a node in `G` and `t` is the time index. For instance, if `('a', 0)`
       is 1 in the ground state, that means the node 'a' is visted first.

    """
    t0 = time.time()
    G = nx.Graph()
    G.add_weighted_edges_from((u, v, 10) for u, v in nx.complete_graph(n_city).edges())

    N = G.number_of_nodes()

    if lagrange is None:
        # If no lagrange parameter provided, set to 'average' tour length.
        # Usually a good estimate for a lagrange parameter is between 75-150%
        # of the objective function value, so we come up with an estimate for 
        # tour length and use that.
        if G.number_of_edges() > 0:
            lagrange = G.size(weight=weight) * G.number_of_nodes() / G.number_of_edges()
        else:
            lagrange = 2

    # calculate default missing_edge_weight if required
    if missing_edge_weight is None:
        # networkx method to calculate sum of all weights
        missing_edge_weight = G.size(weight=weight)

    # some input checking
    if N in (1, 2):
        msg = "graph must have at least 3 nodes or be empty"
        raise ValueError(msg)

    # Creating the QUBO
    Q = defaultdict(float)

    # Constraint that each row has exactly one 1
    for node in G:
        for pos_1 in range(N):
            Q[((node, pos_1), (node, pos_1))] -= lagrange
            for pos_2 in range(pos_1 + 1, N):
                Q[((node, pos_1), (node, pos_2))] += 2.0 * lagrange

    # Constraint that each col has exactly one 1
    for pos in range(N):
        for node_1 in G:
            Q[((node_1, pos), (node_1, pos))] -= lagrange
            for node_2 in set(G) - {node_1}:
                # QUBO coefficient is 2*lagrange, but we are placing this value 
                # above *and* below the diagonal, so we put half in each position.
                Q[((node_1, pos), (node_2, pos))] += lagrange

    # Objective that minimizes distance
    for u, v in itertools.combinations(G.nodes, 2):
        for pos in range(N):
            nextpos = (pos + 1) % N

            # going from u -> v
            try:
                value = G[u][v][weight]
            except KeyError:
                value = missing_edge_weight

            Q[((u, pos), (v, nextpos))] += value

            # going from v -> u
            try:
                value = G[v][u][weight]
            except KeyError:
                value = missing_edge_weight

            Q[((v, pos), (u, nextpos))] += value

    bqm = dimod.BinaryQuadraticModel(Q, dimod.BINARY)
    t1 = time.time()
    return t1 - t0

def measure_dwave(n_cities_list):
    times = {}
    for n in n_cities_list:
        mem_start = memory_usage(max_usage=True)
        max_memory = memory_usage((tsp_dwave, (n, 1, 'weight', None)), max_usage=True)
        times[n] = [tsp_dwave(n, 1, 'weight', None), max_memory - mem_start]
    return times

#####################################################################################################################
# PYQUBO SDK: based on https://github.com/recruit-communications/pyqubo/blob/master/benchmark/benchmark.py
#####################################################################################################################

def tsp_pyqubo(n_city):
    t0 = time.time()
    x = Array.create('c', (n_city, n_city), 'BINARY')
    use_for_loop = True

    # Constraint not to visit more than two cities at the same time.
    time_const = 0.0
    for i in range(n_city):
        # If you wrap the hamiltonian by Const(...), this part is recognized as constraint
        time_const += Constraint((sum(x[i, j] for j in range(n_city)) - 1) ** 2, label="time{}".format(i))

    # Constraint not to visit the same city more than twice.
    city_const = 0.0
    for j in range(n_city):
        city_const += Constraint((sum(x[i, j] for i in range(n_city)) - 1) ** 2, label="city{}".format(j))

    # distance of route
    feed_dict = {}

    if use_for_loop:
        distance = 0.0
        for i in range(n_city):
            for j in range(n_city):
                for k in range(n_city):
                    # we set the constant distance
                    d_ij = 10
                    distance += d_ij * x[k, i] * x[(k + 1) % n_city, j]

    else:
        distance = []
        for i in range(n_city):
            for j in range(n_city):
                for k in range(n_city):
                    # we set the constant distance
                    d_ij = 10
                    distance.append(d_ij * x[k, i] * x[(k + 1) % n_city, j])
        distance = sum(distance)

    # Construct hamiltonian
    A = Placeholder("A")
    Q = distance
    H = Q + city_const + time_const

    # Compile model
    t1 = time.time()
    model = H.compile()
    t2 = time.time()
    qubo, offset = model.to_qubo(index_label=False)
    t3 = time.time()

    # total, express_time, compile_time, to_qubo_time
    return [t3 - t0, t1 - t0, t2 - t1, t3 - t2]

def measure_pyqubo(n_cities_list):
    times = {}
    for n in n_cities_list:
        mem_start = memory_usage(max_usage=True)
        max_memory = memory_usage((tsp_pyqubo, (n,)), max_usage=True)
        times[n] = tsp_pyqubo(n) + [max_memory - mem_start]
    return times

#####################################################################################################################
# DADK SDK: based on tutorial
#####################################################################################################################
# 
def tsp_dadk(n_city):
    t0 = time.time()

    # square array x with constant bits is created as VarShapeSet            
    var_shape_set = VarShapeSet(BitArrayShape('x', (n_city, n_city), one_hot=OneHot.two_way))

    # -P_time-------------------------------------------------------------------------------
    # -P_city-------------------------------------------------------------------------------
    constraints = BinPol(var_shape_set)
    for t in range(n_city):
        constraints.add_exactly_1_bit_on(bits=[('x', t, c) for c in range(n_city)])
    for c in range(n_city):
        constraints.add_exactly_1_bit_on(bits=[('x', t, c) for t in range(n_city)])

    # -P_distance---------------------------------------------------------------------------
    objective = BinPol(var_shape_set)
    for t in range(n_city):
        for c_0 in range(n_city):
            for c_1 in range(n_city):
                if c_0 != c_1:
                    objective.add_term(10, ('x', t, c_0), ('x', (t + 1) % n_city, c_1))

    # -HQ-----------------------------------------------------------------------------------
    qubo = BinPol(var_shape_set)
    qubo.add(objective, 1)
    qubo.add(constraints, 1)

    t1 = time.time()
    return t1 - t0

def measure_dadk(n_cities_list):
    times = {}
    for n in n_cities_list:
        mem_start = memory_usage(max_usage=True)
        max_memory = memory_usage((tsp_dadk, (n,)), max_usage=True)
        times[n] = [tsp_dadk(n), max_memory - mem_start]
    return times

#####################################################################################################################
# measurements
#####################################################################################################################

n_cities_list = [5, 10, 25, 50, 75, 100, 150, 200, 250, 300]
times_pyqubo = measure_pyqubo(n_cities_list)
times_dwave = measure_dwave(n_cities_list)
times_dadk = measure_dadk(n_cities_list)

#####################################################################################################################
# plot
#####################################################################################################################

fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1, 1], 'hspace': 0})

ax[0].plot(times_pyqubo.keys(), [v[0] for v in times_pyqubo.values()], color='blue', marker='x', linestyle='-', label='pyQubo')
ax[0].plot(times_dwave.keys(), [v[0] for v in times_dwave.values()], color='orange', marker='o', linestyle='-', label='dimod')
ax[0].plot(times_dadk.keys(), [v[0] for v in times_dadk.values()], color='red', marker='+', linestyle='-', label='dadk')

ax[1].plot(times_pyqubo.keys(), [v[-1] for v in times_pyqubo.values()], color='blue', marker='x', linestyle='-')
ax[1].plot(times_dwave.keys(), [v[-1] for v in times_dwave.values()], color='orange', marker='o', linestyle='-')
ax[1].plot(times_dadk.keys(), [v[-1] for v in times_dadk.values()], color='red', marker='+', linestyle='-')

def forward(x):
    return x ** 2

def inverse(x):
    return np.sqrt(x)

secax = ax[0].secondary_xaxis('top', functions=(forward, inverse))

secax.set_xlabel('Number of bits')

ax[1].set_xlabel('Number of cities')
ax[0].set_ylabel('Time to QUBO [s]')
ax[1].set_ylabel('Memory [MiB]')

ax[0].set_title(f'Comparison of QUBO SDKs for TSP QUBO')
ax[0].legend(loc='best')

fig.savefig(f"dadk_comp.png", dpi=300)

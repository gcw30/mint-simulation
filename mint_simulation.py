#! /usr/bin/env python3

print('Starting program')

import networkx as nx, numpy as np, random, openpyxl as xl
from networkx.algorithms import bipartite

random.seed(0)
np.random.seed(0)

### VARIABLES ###
k = 1                    # Shape of gamma distribution (k=1 suggested by Esty 2011)
theta = 450              # Die life in minutes (25hrs suggested by Carter 1983)
a = 1.5                  # How much longer obverses are to last than reverses
t = 600                  # Work period length in minutes
D = 20                   # How many obverse dies to generate in graph
W = 1                    # Number of workstations
odb = 1                  # Number of obverses in die box
rdb = 1                  # Number of reverses in die box
o_loss = 0               # Percentage of obverse dies to delete from full chart
r_loss = 0               # Percentage of reverse dies to delete from full chart
e_loss = 0               # Percentage of edges to delete from full chart
i = 10000                # Iterations

assert odb >= W, 'Must have more obverse dies in die box than available workstations'
assert rdb >= W, 'Must have more reverse dies in die box than available workstations'

def get_new_obv(obv):
    """Create a new reverse die."""
    obv += 1
    master_obv = obv
    obv_life = np.random.gamma(k, (a*theta))
    die_lifetimes['Obverses'].setdefault('O' + str(obv), obv_life)
    return master_obv, 'O' + str(obv)

def get_new_rev(rev):
    """Create a new reverse die."""
    rev += 1
    master_rev = rev
    rev_life = np.random.gamma(k, theta)
    die_lifetimes['Reverses'].setdefault(master_rev, rev_life)
    return master_rev, rev

def pop_odb(master_obv, obvs):
    """Create obverse dies for available spaces in die box"""
    to_add = odb - len(obvs)
    for y in range(to_add):
        obvs.append(None)
        master_obv, obvs[-1] = get_new_obv(master_obv)
    return master_obv, obvs

def pop_rdb(master_rev, revs):
    """Create reverse dies for available spaces in die box"""
    to_add = rdb - len(revs)
    for z in range(to_add):
        revs.append(None)
        master_rev, revs[-1] = get_new_rev(master_rev)
    return master_rev, revs

def create_graph():
    """Create a die chart based on input variables"""
    t_values = {}
    obvs = []
    master_obv = 0
    revs = []
    master_rev = 0
    G = nx.Graph()
    master_obv, obvs = pop_odb(master_obv, obvs)
    for o in obvs:
        G.add_node(o, label=o)
        t_values.setdefault(o, 0)
    master_rev, revs = pop_rdb(master_rev, revs)
    for r in revs:
        G.add_node(r, label=str(r))
        t_values.setdefault(r, 0)
    while master_obv < D:
        T = 0
        stations = []
        obv_in_use = []
        rev_in_use = []
        e = np.random.choice(len(obvs), W, replace=False)
        f = np.random.choice(len(revs), W, replace=False)
        for x in range(W):
            stations.append([obvs[e[x]],revs[f[x]]])
            obv_in_use.append(obvs[e[x]])
            rev_in_use.append(revs[f[x]])

        while T != t:
            for workstation, pair in enumerate(stations):
                obv_idx = pair[0]
                rev_idx = pair[1]
                if pair in G.edges:
                    G.edges[pair]['t_count'] += 1
                else:
                    G.add_edge(obv_idx, rev_idx, t_count=1)
                t_values[obv_idx] += 1
                t_values[rev_idx] += 1
                if t_values[rev_idx] >= die_lifetimes['Reverses'][rev_idx]:
                    revs.remove(rev_idx)
                    rev_in_use.remove(rev_idx)
                    av_revs = list(set(revs) - set(rev_in_use))
                    if av_revs == []:
                        master_rev, revs = pop_rdb(master_rev, revs)
                        for r in revs:
                            G.add_node(r, label=str(r))
                            t_values.setdefault(r, 0)
                    av_revs = list(set(revs) - set(rev_in_use))
                    f = np.random.choice(len(av_revs))
                    stations[workstation][1] = av_revs[f]
                    rev_in_use.append(av_revs[f])
                if t_values[obv_idx] >= die_lifetimes['Obverses'][obv_idx]:
                    obvs.remove(obv_idx)
                    obv_in_use.remove(obv_idx)
                    av_obvs = list(set(obvs) - set(obv_in_use))
                    if av_obvs == []:
                        master_obv, obvs = pop_odb(master_obv, obvs)
                        for o in obvs:
                            G.add_node(o, label=o)
                            t_values.setdefault(o, 0)
                    av_obvs = list(set(obvs) - set(obv_in_use))
                    e = np.random.choice(len(av_obvs))
                    stations[workstation][0] = av_obvs[e]
                    obv_in_use.append(av_obvs[e])
            T += 1
    return G

def disruption(G, o_loss, r_loss, e_loss):
    """Delete nodes and edges from graph according to set parameters"""
    obv_nodes = {n for n, d in G.nodes(data=True) if d['label'].startswith('O')}
    rev_nodes = set(G) - obv_nodes
    count_obvs_to_delete = int(len(obv_nodes) * (o_loss/100))
    count_revs_to_delete = int(len(rev_nodes) * (r_loss/100))

    # Create probabilities for obv deletion that are inverse of obverse die lifetimes
    o_ratios = []
    o_probs = []
    for v in list(die_lifetimes['Obverses'].values()):
        o_ratios.append(1/v)
    for r in o_ratios:
        o_probs.append(r/sum(o_ratios))
    o_probs = tuple(o_probs)

    # Create probabilities for rev deletion that are inverse of reverse die lifetimes
    r_ratios = []
    r_probs = []
    for v in list(die_lifetimes['Reverses'].values()):
        r_ratios.append(1/v)
    for r in r_ratios:
        r_probs.append(r/sum(r_ratios))
    r_probs = tuple(r_probs)

    # Choose nodes to delete
    n_to_delete = list(np.random.choice(list(die_lifetimes['Obverses'].keys()), size=count_obvs_to_delete, replace=False, p=o_probs))
    for x in list(np.random.choice(list(die_lifetimes['Reverses'].keys()), size=count_revs_to_delete, replace=False, p=r_probs)):
        n_to_delete.append(x)

    G.remove_nodes_from(n_to_delete)

    # Create probabilities for edge deletion that are inverse of t_counts
    e_ratios = []
    e_probs = []
    for u, v, t_count in G.edges.data('t_count'):
        e_ratios.append(1/t_count)
    for r in e_ratios:
        e_probs.append(r/sum(e_ratios))
    e_probs = tuple(e_probs)

    # Choose edges to delete
    count_edges_to_delete = int(nx.number_of_edges(G) * (e_loss/100))
    e_to_delete = []
    for e in list(np.random.choice(len(G.edges()), size=count_edges_to_delete, replace=False, p=e_probs)):
        e_to_delete.append(list(G.edges())[e])

    G.remove_edges_from(e_to_delete)

    return G

def analysis(G):
    """Calculate properties of a given cie chart"""
    obv_nodes = {n for n, d in G.nodes(data=True) if d['label'].startswith('O')}
    rev_nodes = set(G) - obv_nodes

    obv_count = len(obv_nodes)
    rev_count = len(rev_nodes)
    edge_count = nx.number_of_edges(G)

    r_deg, o_deg = bipartite.degrees(G, obv_nodes)

    o_deg_max = max(x[1] for x in o_deg)
    r_deg_max = max(x[1] for x in r_deg)

    o_deg_min = min(x[1] for x in o_deg)
    r_deg_min = min(x[1] for x in r_deg)

    o_deg_mean = float((sum(x[1] for x in o_deg))/len(o_deg))
    r_deg_mean = float((sum(x[1] for x in r_deg))/len(r_deg))

    o_av_cluster = bipartite.average_clustering(G, obv_nodes)
    r_av_cluster = bipartite.average_clustering(G, rev_nodes)

    connect = nx.is_connected(G)

    c = sorted(nx.connected_components(G), key = len, reverse=True)
    min_comp = nx.number_of_edges(G.subgraph((c[-1])))
    max_comp = nx.number_of_edges(G.subgraph((c[0])))

    """
    If the graph is connected, calculate distance metrics for the whole graph.
    If the graph is not connected, calculate distance metrics for the largest
    connected component.
    """
    if connect:
        spl = nx.average_shortest_path_length(G)
        diam = nx.diameter(G)
    else:
        spl = nx.average_shortest_path_length(G.subgraph(c[0]))             # finds average shortest path length of largest connected component
        diam = nx.diameter(G.subgraph(c[0]))

    density = nx.density(G)

    planar = nx.check_planarity(G)[0]

    return obv_count, rev_count, edge_count, o_deg_min, o_deg_mean, o_deg_max, r_deg_min, r_deg_mean, r_deg_max, o_av_cluster, r_av_cluster, connect, min_comp, max_comp, spl, diam, density, planar

# Set up output Excel file
wb = xl.Workbook()
wb.active.title = 'Variables'
wb.create_sheet(index=1, title='Output')
sheet1 = wb['Variables']
sheet2 = wb['Output']

sheet1['A2'] = 'Gamma distribution shape'
sheet1['A3'] = 'Gamma distribution scale'
sheet1['A4'] = 'a (scale factor)'
sheet1['A5'] = 't (length of working period)'
sheet1['A6'] = 'Obverse dies generated'
sheet1['A7'] = 'Number of workstations'
sheet1['A8'] = 'Obverses in die box'
sheet1['A9'] = 'Reverses in die box'
sheet1['A10'] = 'Obverse loss (%)'
sheet1['A11'] = 'Reverse loss (%)'
sheet1['A12'] = 'Edges loss (%)'
sheet1['A13'] = 'Iterations'

sheet1['B2'] = k
sheet1['B3'] = theta
sheet1['B4'] = a
sheet1['B5'] = t
sheet1['B6'] = D
sheet1['B7'] = W
sheet1['B8'] = odb
sheet1['B9'] = rdb
sheet1['B10'] = o_loss
sheet1['B11'] = r_loss
sheet1['B12'] = e_loss
sheet1['B13'] = i

sheet1.column_dimensions['A'].width = 30

sheet2['A1'] = 'Iteration'
sheet2['B1'] = 'Revs/Obv'
sheet2['C1'] = 'Obv. min. degree'
sheet2.column_dimensions['C'].width = 17
sheet2['D1'] = 'Obv. mean degree'
sheet2.column_dimensions['D'].width = 17
sheet2['E1'] = 'Obv. max. degree'
sheet2.column_dimensions['E'].width = 17
sheet2['F1'] = 'Rev. min. degree'
sheet2.column_dimensions['F'].width = 17
sheet2['G1'] = 'Rev. mean degree'
sheet2.column_dimensions['G'].width = 17
sheet2['H1'] = 'Rev. max. degree'
sheet2.column_dimensions['H'].width = 17
sheet2['I1'] = 'Obv. av. clustering'
sheet2.column_dimensions['I'].width = 20
sheet2['J1'] = 'Rev. av. clustering'
sheet2.column_dimensions['J'].width = 20
sheet2['K1'] = 'Connected?'
sheet2['L1'] = 'Edges in smallest component'
sheet2.column_dimensions['L'].width = 24
sheet2['M1'] = 'Edges in largest component'
sheet2.column_dimensions['M'].width = 24
sheet2['N1'] = 'Av. Shortest Path Length'
sheet2.column_dimensions['N'].width = 20
sheet2['O1'] = 'Diameter'
sheet2['P1'] = 'Density'
sheet2['Q1'] = 'Planar?'
sheet2['R1'] = 'Die combo ratio'
sheet2.column_dimensions['R'].width = 15

# Create die charts and record properties to output file
for j in range(0, i):
    die_lifetimes = {'Obverses':{}, 'Reverses':{}}
    G = create_graph()
    if o_loss != 0 or r_loss != 0 or e_loss != 0:
        G = disruption(G, o_loss, r_loss, e_loss)
    G.remove_nodes_from(list(nx.isolates(G)))
    obv_count, \
    rev_count, \
    edge_count, \
    o_deg_min, \
    o_deg_mean, \
    o_deg_max, \
    r_deg_min, \
    r_deg_mean, \
    r_deg_max, \
    o_av_cluster, \
    r_av_cluster, \
    connect, \
    min_comp, \
    max_comp, \
    spl, \
    diam, \
    density, \
    planar = analysis(G)
    sheet2['A' + str(j + 2)] = j + 1
    sheet2['B' + str(j + 2)] = rev_count/obv_count
    sheet2['C' + str(j + 2)] = o_deg_min
    sheet2['D' + str(j + 2)] = o_deg_mean
    sheet2['E' + str(j + 2)] = o_deg_max
    sheet2['F' + str(j + 2)] = r_deg_min
    sheet2['G' + str(j + 2)] = r_deg_mean
    sheet2['H' + str(j + 2)] = r_deg_max
    sheet2['I' + str(j + 2)] = o_av_cluster
    sheet2['J' + str(j + 2)] = r_av_cluster
    sheet2['K' + str(j + 2)] = connect
    sheet2['L' + str(j + 2)] = min_comp
    sheet2['M' + str(j + 2)] = max_comp
    sheet2['N' + str(j + 2)] = spl
    sheet2['O' + str(j + 2)] = diam
    sheet2['P' + str(j + 2)] = density
    sheet2['Q' + str(j + 2)] = planar
    sheet2['R' + str(j + 2)] = edge_count/(obv_count + rev_count)

wb.save('output.xlsx')
print('End of Program')

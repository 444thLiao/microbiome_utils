import numpy as np
import pandas as pd
import networkx as nx
from sklearn.neighbors import KNeighborsClassifier as knc
from networkx.readwrite import json_graph
from nets2 import visualize
from itertools import count
distance_matrix_file = '/home/liaoth/data2/16s/171027_16s/16s_pipeliens/XK/final_result/analysis/beta/unifrac.txt'
dm = pd.read_csv(distance_matrix_file, sep='\t', header=0, index_col=0)
assert list(dm.columns) == list(dm.index)

# MST part
nodes_n = dm.shape[0]
G = nx.from_numpy_matrix(dm.values)
nx.set_node_attributes(G, 'samples', dict(zip(range(nodes_n), [[it] for it in dm.index])))
names = dict(zip(range(nodes_n), dm.index))
G = nx.relabel_nodes(G, names)
mst = nx.minimum_spanning_tree(G)
mst_json = json_graph.node_link_data(mst)
for each in mst_json['links']:
    each['type'] = 'mst'

#mst_json['links'] = []
# if you uncomment upper line, this will also display MST lines or it will only display KNN line.


# KNN part
# KNN config
nearest_num = 1
###
M = knc(weights='distance', metric='precomputed')
M.fit(dm.values.tolist(), list(dm.index))

query_dict = {}
for _idx, name in enumerate(mst_json['nodes']):
    name = name['id']
    query_dict[name] = _idx

for _idx, name in enumerate(list(dm.index)):
    temp = M.kneighbors(np.array(dm.values.tolist()[_idx]).reshape(1, -1), n_neighbors=nearest_num + 1)
    for num in range(nearest_num):
        links = {'source': query_dict[name],
                 'target': query_dict[list(dm.index)[temp[1][0][num + 1]]],
                 'weight': temp[0][0][num + 1],
                 'type': 'knn',
                 'rank': num}
        mst_json['links'].append(links)

visualize(mst_json,
          path_html='mst_vis.html',
          graph_charge=-200,
          graph_gravity=0.05,
          graph_link_distance=10)


import numpy as np
from sklearn.manifold import MDS
from sklearn.neighbors import KNeighborsClassifier as knc
import plotly,os
from plotly import graph_objs as go
import networkx as nx
from utils import read_data
import seaborn as sns

#Default Config (DO NOT CHANGE.)
otherinfo = None
dot_color = None
new_labels = None
MST_color = ('rgb(205, 12, 24)')
MST_width = 8
KNN_line_mode = 'dash'
KNN_width = 4
# Config
metadata = './example/MDS/sample_info.csv'
distance_file = './example/MDS/unifrac_binary.txt'
col_need = 'class'  # which columns in metadata you need to parse name in distance to.
otherinfo_col = 'person' # which columns you need to grouped, in order to use color to
output_dir = './example/MDS/drawn/'

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

distance = read_data(distance_file)
names = list(distance.index)

if metadata:
    sample_info = read_data(metadata)
    new_labels = list(
        sample_info.loc[names, col_need])  # change id to new name. New name is assign at one col in metadata
    otherinfo = list(sample_info.loc[names, otherinfo_col])  # extract person name or id.

# dot color assign
if otherinfo:
    person_ids = set(otherinfo)
    colors = sns.color_palette("Set1", n_colors=len(person_ids), desat=.5).as_hex()
    cc = dict(zip(person_ids, colors))
    dot_color = [cc[label] for label in otherinfo]

if not dot_color:
    dot_color = '#000000'


# build MDS matrix, it will be the basic.
mds = MDS(n_components=2, dissimilarity='precomputed')
mds.fit(distance.values.tolist())
vals = mds.embedding_
# vals is a big dataframe which is storged coordination of MDS's matrix. Each sample is a row of x and y.


# collect all traces.
draw_data = []

# MST
G = nx.from_numpy_matrix(np.array(distance))
# Using distance to construct a graph in order to cal MST.
G = nx.relabel_nodes(G, {_i: list(distance.columns)[_i] for _i in range(len(G.nodes()))})
# Above graph it isn't have any labels. we need to assign its.
mst_G = nx.mst.minimum_spanning_tree(G)
# generate a new graph which is MST version of ori G.
for edge in mst_G.edges():
    draw_data.append(go.Scatter(x=[vals[names.index(edge[0]), 0], vals[names.index(edge[1]), 0]],
                                y=[vals[names.index(edge[0]), 1], vals[names.index(edge[1]), 1]],
                                mode='lines',
                                showlegend=False,
                                line=dict(color=MST_color,
                                          width=MST_width)
                                )
                     )


# KNN
M = knc(n_neighbors=3, weights='distance', metric='precomputed')
M.fit(distance.values.tolist(), list(distance.index))

# same as before, Model doen't have any label, so below if using whole row to represent sample.
for _idx in range(len(names)):
    temp = M.kneighbors(distance.values.tolist()[_idx], n_neighbors=2)[1]
    # 2 is mean choose closest one sample, total is 2, but need to consider itself is the nearest.

    # temp is a complex structure object.
    # last [1:] mean except itself, left -> right is distance nearest -> far.
    for _x in temp[0][1:]:
        if type(dot_color) == list:
            current_color = dot_color[_idx]
        else:
            current_color = ['#000000']
        draw_data.append(go.Scatter(x=[vals[_idx, 0], vals[_x, 0]],
                                    y=[vals[_idx, 1], vals[_x, 1]],
                                    showlegend=False,
                                    mode='lines',
                                    line=dict(width=KNN_width,
                                              color=current_color,
                                              dash=KNN_line_mode)))

### just for add legend,so it doen't need to loop.
draw_data.append(go.Scatter(x=[0], y=[0], showlegend=True, mode='lines', name='KNN(K Nearest Neighbor)',
                            line=dict(width=KNN_width, color='#000000', dash=KNN_line_mode)))
draw_data.append(go.Scatter(x=[0],
                            y=[0],
                            mode='lines',
                            showlegend=True,
                            name='MST(minimum spanning tree)',
                            line=dict(color=MST_color,
                                      width=MST_width)))

# MDS
# Draw text in each dot.
if new_labels:
    name_text = new_labels
else:
    name_text = names

if type(dot_color) != list:
    current_color = [dot_color] * len(vals[:, 0])
else:
    current_color = dot_color
draw_data.append(go.Scatter(x=vals[:, 0], y=vals[:, 1], text=name_text, marker=dict(size=28, color=current_color),
                            textfont=dict(size=20), textposition='top left', mode='markers+text',
                            name='Sample projection'))

layout = go.Layout(title='MDS projection with MST and KNN(Unifrac binary)',
                   font=dict(size=20))
plotly.offline.plot(dict(data=draw_data, layout=layout), filename=output_dir + '/mds-data.html')

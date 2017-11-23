import plotly
from plotly import figure_factory as ff
from plotly import graph_objs as go
from self_16s.utils import read_data
from sklearn.manifold import MDS
import itertools
from scipy.stats.stats import ranksums
def product_dot(samples1,samples2,name):
    repeated = []
    return_dot_y = []
    for s1, s2 in itertools.product(samples1, samples2):
        if s1 != s2 and sorted((s1, s2)) not in repeated:
            return_dot_y.append(distance.loc[s1, s2])
            repeated.append(sorted((s1, s2)))
            all_x.append(name)
    return return_dot_y

metric = 'unifrac'
distance = read_data('/home/liaoth/data2/16s/171027_16s/16s_pipeliens/NY/analysis/routine_analysis_TAN/beta/%s.txt' % metric)

all_sample = list(distance.index)



def group(sample,rename=False):
    if not rename:
        if '-N' in sample or 'YN' in sample:
            return True
        else:
            return False
    else:
        if '-N' in sample:
            new_s = sample.replace('-N', '-T')
        elif 'YN' in sample:
            new_s = sample.replace('YN','YT')
        else:
            new_s = None
        return new_s

all_N = [s for s in all_sample if group(s)]
all_T = [group(s,rename=True) for s in all_N]

all_N_idx = [all_sample.index(s) for s in all_N]
all_T_idx = [all_sample.index(s) for s in all_T]

all_x = []

normal_dot = product_dot(all_N,all_N,'normal')


paired_dot = []
for _N,_T in zip(all_N,all_T):
    all_x.append('Paired')
    paired_dot.append(distance.loc[_N,_T])

tumor_dot = product_dot(all_T,all_T,'tumor')

draw_data = []
#violin_df = df(index=range(len(all_x)))
#fig = ff.create_violin(violin_df, data_header='distance', group_header='group')
draw_data.append(go.Box(x=all_x,y = normal_dot+tumor_dot+paired_dot,boxpoints='all',showlegend=False))
layout = go.Layout(title='Pairwise distance boxplot<Br>%s' % metric ,
                   font=dict(size=20),
                   xaxis=dict(title='factor(group)'),
                   yaxis=dict(title='distance'))

print ranksums(normal_dot,tumor_dot).pvalue,'N-T'
print ranksums(normal_dot,paired_dot).pvalue,'N-P'
print ranksums(tumor_dot,paired_dot).pvalue,'T-P'
plotly.offline.plot(dict(data=draw_data,layout=layout),
                    filename='/home/liaoth/data2/16s/171027_16s/16s_pipeliens/NY/analysis/drawn/pairwise_distance/%s.html' % metric)
                  #  filename='/home/liaoth/data2/16s/171027_16s/16s_pipeliens/NY/analysis/drawn/pairwise_marker/pairwise %s' % fea)


# unifrac
# 1.33274017035e-05 N-T
# 0.490240956372 N-P
# 0.0507292976232 T-P

# unifrac_binary
# 0.254693432503 N-T
# 0.181198765269 N-P
# 0.258897696146 T-P
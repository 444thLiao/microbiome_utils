from __future__ import print_function
import pandas as pd
import numpy as np
import csv
import re,os
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from boruta import BorutaPy
import plotly.graph_objs as go
import plotly.offline as ply
from ..utils import read_data

def parse_tax(tax):
    pat = '^([d,p,c,o,f,g]):(.+)\(([\.\d]+)\)$'
    m = re.match(pat, tax)
    level, name, confidence = m.groups()
    name = name.strip('"')
    return level, name, float(confidence)


def get_otu_tax(infile, otus,c=0.8,filter=None):
    """
    to-do: allow for specification of rank levels to be featured
    :param infile:
    :param c:
    :param filter: 'g' for genus only
    :return:
    """
    otu_tax = defaultdict(list)
    with open(infile) as infh:
        reader = csv.reader(infh, delimiter='\t')
        for row in reader:
            otu_id, tax_predicted, strand, tax_assigned = row
            if otu_id in otus:
                if strand != '+':
                    print('Warning: OTU in minus strand: %s' % (otu_id,))
                tax_ranks = [parse_tax(tax) for tax in tax_predicted.split(',')]
                tax_ranks_prune = [tax_rank for tax_rank in tax_ranks if tax_rank[-1] >= c]
                for i in range(len(tax_ranks_prune)):
                    tax_key = ','.join(['%s:%s' % (it[0], it[1]) for it in tax_ranks_prune[:i+1]])
                    otu_tax[tax_key].append(otu_id)
                # append OTUs (leaves of the taxonomy tree)
                tax_key_otu = ','.join([tax_key, otu_id])
                otu_tax[tax_key_otu].append(otu_id)
    if filter:
        new_otu_tax = defaultdict(list)
        for key in otu_tax.keys():
            if key.split(',')[-1].startswith('%s:' % filter):
                new_otu_tax[key.split(',')[-1]] += otu_tax[key]
        return new_otu_tax
    return otu_tax


def get_tax_profile(otu_tab, otu_tax):
    """
    otu_tab should take a form of rows(samples) by cols (OTUs)
    :param otu_tab:
    :param otu_tax:
    :return:
    """
    tax_profile = pd.DataFrame()
    for tax_name, otus in otu_tax.items():
        tax_profile[tax_name] = otu_tab.loc[:,otus].sum(axis=1)
    return tax_profile

def run_boruta(tax_file,otu_tab,group_method,filter=None,fn=None,output_result_dir=None,subset_otu=None,max_depth=5,max_iter=1000,is_normalized=False):
    """

    :param tax_file:
    :param otu_tab:
    :param group_method: a function to group the sample, receive each sample,return a group label.
    :param filter: receive a tax level abbre. Such as, if you want genus level only without upper or lower,you can type 'g'.
    :param fn: output html absoulute filename.
    :param subset_otu: samples name you want to use only
    :return:
    """

    otu_tab = read_data(otu_tab)
    # check
    if not otu_tab.index.values[0].startswith('OTU'):
        otu_tab = otu_tab.T

    if subset_otu:
        otu_tab = otu_tab.loc[:,subset_otu]
        otu_tab = otu_tab[otu_tab.sum(1)!=0]

    # transpose into samples (rows) by OTUs (cols)
    otu_tab = otu_tab.T
    if not is_normalized:
        # normalization into relative abundance
        otu_tab = otu_tab.div(otu_tab.sum(axis=1), axis=0)

    otu_tab = otu_tab.loc[:, otu_tab.sum(0) != 0]
    otus = otu_tab.columns.values.tolist()

    otu_tax = get_otu_tax(tax_file, otus, filter=filter)
    # propagate with tax profiles
    try:
        tax_tab = get_tax_profile(otu_tab, otu_tax)
    except:
        import pdb;pdb.set_trace()
    # tax_tab
    # get sample metadata
    samples = list(tax_tab.index)
    groups = [_ for _ in map(group_method, samples)]


    ###########################################################################################
    ### boruta ################################################################################
    ###########################################################################################
    X = tax_tab.values
    print(X.shape)
    y = groups
    features = tax_tab.columns
    max_depth = max_depth
    max_iter = max_iter


    # define random forest classifier, with utilising all cores and
    # sampling in proportion to y labels
    # [tree pruning!!!] highly recommend using pruned trees with a depth between 3-7
    rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=max_depth, random_state=123)
    # define Boruta feature selection method
    feat_selector = BorutaPy(rf, n_estimators='auto', verbose=1, random_state=123, max_iter=max_iter)
    # find all relevant features
    try:
        feat_selector.fit(X, y)
    except:
        import pdb;pdb.set_trace()
    # check selected features
    #print feat_selector.support_
    # check ranking of features
    #print feat_selector.ranking_
    # call transform() on X to filter it down to selected features
    # X_filtered = feat_selector.transform(X)
    print("confirmed features:")
    print(features[feat_selector.support_])
    print("tentative features:")
    print(features[feat_selector.support_weak_])
    # find out the trend of importances.
    fea_status = []
    if len(set(y)) == 2:
        a,b = tuple(set(y))
        for fea in list(features):
            med_a = np.median(tax_tab.loc[ [i for i,v in zip(tax_tab.index.values.tolist(),y) if v == a] ,fea])
            med_b = np.median(tax_tab.loc[ [i for i,v in zip(tax_tab.index.values.tolist(),y) if v == b] ,fea])
            if med_a >= med_b:
                fea_status.append('(%s > %s)' % (a,b))
            else:
                fea_status.append('(%s < %s)' % (a, b))
        features = np.array([i + v for i,v in zip(features.values.tolist(),fea_status)])
    # get feature importances
    # forest = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=max_depth,
    #                                 random_state=123, n_estimators=max_iter)
    # forest.fit(X, y)
    importances = feat_selector._get_imp(X, y)
    std = np.std([tree.feature_importances_ for tree in feat_selector.estimator.estimators_], axis=0)
    # reverse sorted index
    colors = np.empty(len(features), dtype=object)
    colors.fill('blue')
    colors[feat_selector.support_] = 'red'
    colors[feat_selector.support_weak_] = 'yellow'
    indices = np.argsort(importances)
    selected_indices = indices[-40:]
    print(selected_indices)
    trace = go.Bar(y=features[selected_indices],
                   x=importances[selected_indices],
                   marker=dict(color=colors[selected_indices]),
                   error_x=dict(visible=True, arrayminus=std[selected_indices]),
                   orientation='h',
                   )
    layout = go.Layout(title="Feature importances",
                       margin=go.Margin(l=800))
    fig = go.Figure(data=[trace], layout=layout)
    if output_result_dir:
        if not os.path.isdir(output_result_dir):
            os.makedirs(output_result_dir)
        tax_tab.to_csv(output_result_dir+'/input_data.tab',sep='\t')
        tax_tab.loc[:,feat_selector.support_].to_csv(output_result_dir+'/confirmed_fetures_data.tab',sep='\t')
        tax_tab.loc[:, feat_selector.support_weak_].to_csv(output_result_dir + '/tentative_fetures_data.tab', sep='\t')
        tmp = pd.DataFrame(index=features[indices],columns=['importances'])
        tmp.loc[:,'importances'] = importances[indices]
        tmp.to_csv(output_result_dir + '/fetures_importances.tab', sep='\t')
    if fn:
        ply.plot(fig,filename=fn)
    else:
        ply.plot(fig)





if __name__ == '__main__':
    tax_file = '/home/liaoth/data2/16s/171027_16s/16s_pipeliens/NY/analysis/regular_analysis/sintax.txt'
    otu_tax = get_otu_tax(tax_file)
    otu_tab = pd.read_csv('/home/liaoth/temp.csv', sep='\t', index_col=0, header=0)
    # transpose into samples (rows) by OTUs (cols)
    otu_tab = otu_tab.T

    # otu_tab = otu_tab.loc[set(otu_tab.index).difference(set(['TAN-3-N-2', 'TAN-1-N-2'])),:]
    # otu_tab = otu_tab.loc[[i for i in list(otu_tab.index) if 'NY' not in i],:]
    otu_tab = otu_tab.loc[:, otu_tab.sum(0) != 0]

    # normalization into relative abundance
    otu_tab = otu_tab.div(otu_tab.sum(axis=1), axis=0)

    # propagate with tax profiles
    tax_tab = get_tax_profile(otu_tab, otu_tax)

    # get sample metadata
    def get_meta_data(sample):
        if sample.endswith('-T'):
            return 'T'
        elif sample.startswith('NYT'):
            return 'T'
        else:
            return 'N'

    samples = list(tax_tab.index)
    groups = map(get_meta_data, samples)

    ###########################################################################################
    ### boruta ################################################################################
    ###########################################################################################
    X = tax_tab.values
    print(X.shape)
    y = groups
    features = tax_tab.columns
    max_depth = 5
    max_iter = 500



    # define random forest classifier, with utilising all cores and
    # sampling in proportion to y labels
    # [tree pruning!!!] highly recommend using pruned trees with a depth between 3-7
    rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=max_depth, random_state=123)
    # define Boruta feature selection method
    feat_selector = BorutaPy(rf, n_estimators='balanced', verbose=1, random_state=123, max_iter=max_iter)
    # find all relevant features
    feat_selector.fit(X, y)
    # check selected features
    #print feat_selector.support_
    # check ranking of features
    #print feat_selector.ranking_
    # call transform() on X to filter it down to selected features
    # X_filtered = feat_selector.transform(X)
    print("confirmed features:")
    print(features[feat_selector.support_])
    print("tentative features:")
    print(features[feat_selector.support_weak_])

    # get feature importances
    forest = RandomForestClassifier(n_jobs=-1, class_weight='auto', max_depth=max_depth,
                                    random_state=123, n_estimators=max_iter)
    forest.fit(X, y)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    # reverse sorted index
    colors = np.empty(len(features), dtype=object)
    colors.fill('blue')
    colors[feat_selector.support_] = 'red'
    colors[feat_selector.support_weak_] = 'yellow'
    indices = np.argsort(importances)
    selected_indices = indices[-40:]
    print(selected_indices)
    trace = go.Bar(y=features[selected_indices],
                   x=importances[selected_indices],
                   marker=dict(color=colors[selected_indices]),
                   error_x=dict(visible=True, arrayminus=std[selected_indices]),
                   orientation='h',
                   )
    layout = go.Layout(title="Feature importances",
                       margin=go.Margin(l=800))
    fig = go.Figure(data=[trace], layout=layout)
    ply.plot(fig)
             #filename='data/feature.importances.html')



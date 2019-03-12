from __future__ import print_function
import argparse
import os
import pandas
import pickle
import random
from collections import defaultdict, Counter
from multiprocessing import Pool,Manager
from utils import read_data
import plotly
import plotly.graph_objs as go
import tqdm
from pandas import DataFrame as df
from skbio import TreeNode
from skbio.diversity import alpha_diversity, get_alpha_diversity_metrics


# tree_file = '/home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otus.tree'
# Otu_table = '/home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otu_raw_filterd_e5.txt'
# metadata = '/home/liaoth/data2/hospital_seven/PG17H11010002Amp01/sample_info.csv'
# step = 500
# repeat = 100
# output_dir = '/home/liaoth/data2/hospital_seven/PG17H11010002Amp01/drawn/'
# col_need = 'class'  # which columns in 'metadata' you need to parse name in distance to.

# python /home/liaoth/project/16s_pipelines/self_16s/draw_PD.py -t /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otus.tree -i /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otu_raw_filterd_e5.txt -m /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/sample_info.csv -o /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/drawn/


# Subsampling
def subsampling(ori_df, num, prebuild=None):
    """
    :param ori_df: row is OTU,cols is sample.
    :param num:
    :return:
    """
    if not prebuild:
        # offer a prebuilded dict.(For speed)
        sample_bucket = defaultdict(list)
        # if not, init a dict for storge all sample's OTU.
        for s in list(ori_df.columns):
            for r in list(ori_df.index):
                sample_bucket[s] += [r] * ori_df.loc[r, s]
                # each sample: a list include 'OTU' and weighted by it num in matrix.
    else:
        sample_bucket = prebuild

    sub_df = df(index=ori_df.columns, columns=ori_df.index,
                data=0)  # init a subsampled matrix, row is sample, col is OTU.

    bucket = {sample: Counter([random.choice(sample_bucket[sample]) for _ in range(num)])
              for sample in sample_bucket.keys()}
    # Use twice comprehension to speed subsampling process.  sample as key, otu count as col and val.
    sub_df.update(df(bucket).fillna(0).T)
    # first it will full of NAN, so fill nan with 0. And it need to transpose.
    return sub_df.astype(int)


def subsampling2(ori_df, num, prebuild=None):
    """
    :param ori_df: row is OTU,cols is sample.
    :param num:
    :return:
    """
    if not prebuild:
        sample_bucket = defaultdict(list)
        for s in list(ori_df.columns):
            for r in list(ori_df.index):
                sample_bucket[s] += df([r] * ori_df.loc[r, s])
    else:
        sample_bucket = prebuild
        for key in sample_bucket.keys():
            if type(sample_bucket[key]) != df:
                sample_bucket[key] = df(sample_bucket[key])

    sub_df = df(index=ori_df.columns, columns=ori_df.index, data=0)
    for sample in sample_bucket.keys():
        bucket = list(sample_bucket[sample].sample(num).iloc[:, 0])
        sub_df.update(df(Counter(bucket), index=[sample]))
    return sub_df.astype(int)


def generate_step(step,maxiumn):
    "For generate unequal step for draw more smooth rarefaction curve"
    generate_range = []
    small_bin = range(0,step+1,100)[1:]
    for idx,small_max in enumerate(small_bin[1:]):
        for s_step in range(small_bin[idx],small_max+1,int(small_bin[idx]/10)):
            generate_range.append(s_step)
    return generate_range+list(range(step,maxiumn+1,step))

def multiprocess_subsample(bucket,ori_otu, num, prebuild_dict):
    bucket.extend([subsampling(ori_otu, num, prebuild=prebuild_dict)])

def diversity_ana(metric,subsample,ids,**kwargs):
    if metric == 'faith_pd':
        each = alpha_diversity('faith_pd', subsample, ids=ids,
                               otu_ids=kwargs['otu_ids'], tree=kwargs['tree'])
    elif metric == 'shannon':
        each = alpha_diversity('shannon', subsample, ids=ids)
    elif metric == 'observed_otus':
        each = alpha_diversity('observed_otus', subsample, ids=ids)
    else:
        try:
            each = alpha_diversity(metric, subsample, ids=ids)
        except:
            print('Metric you can use is listed below: \n' + '\n'.join(get_alpha_diversity_metrics()))
            exit()
    return each

# Start cal alpha diversity
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    If you want to do example. Like """)
    parser.add_argument('-t', "--tree", help="Otus tree in newick format.")
    parser.add_argument('-i', "--input", help="Otus table")
    parser.add_argument('-m', "--metadata", help="Metadata to parse otus sample id into normal or other group name.")
    parser.add_argument('-o', "--output", help="Output dir, File name dosen't need to set.")

    parser.add_argument('-M', "--metric", help="Metric you want to use.[default: %(default)s]", default='faith_pd')
    parser.add_argument('-s', "--step", help="Subsampling step [default: %(default)s]", default=500)
    parser.add_argument('-r', "--repeat", help="Subsampling iteration num [default: %(default)s]", default=100)
    parser.add_argument('-thread', help="subsampleing process. [default: %(default)s]", default=25)
    parser.add_argument('-cols',
                        help="which columns in 'metadata' you need to parse name in distance to. [default: %(default)s]",
                        default='class')
    parser.add_argument('--list_metric', action='store_true',
                        help="List all metric name you could use.")
    parser.add_argument('--output_fig',
                        help="If you need to remodify figure, you can assign file path, we will pickle it into file.")

    args = parser.parse_args()
    threads = args.thread
    tree_file = args.tree
    Otu_table = args.input
    metadata = args.metadata
    output_dir = args.output
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if args.list_metric:
        print('Metric you can use is listed below: \n' + '\n'.join(get_alpha_diversity_metrics()))
        exit()
    metric = args.metric
    if ',' in metric:
        metric = metric.split(',')
    else:
        metric = [metric]
    step = args.step
    repeat = args.repeat
    col_need = args.cols

    ori_otu = read_data(Otu_table)
    if metadata:
        sample_info = read_data(metadata)
    if tree_file:
        tree = TreeNode.read(tree_file, format='newick')
    else:
        tree = None

    # Above all is parse the argument.
    print('Start building the query dict')
    prebuild_dict = defaultdict(list)
    for s in list(ori_otu.columns):
        prebuild_dict[s] = sum(map(lambda y: [y[0]] * y[1], zip(ori_otu.index, list(ori_otu.loc[:, s]))), [])
        # for r in list(ori_otu.index):
        #     prebuild_dict[s] += [r] * ori_otu.loc[r, s]

    print('Down the query dict. Start the progressBar.')
    # Build the query dict.
    # each sample: a list include 'OTU' and weighted by it num in matrix.

    pbar_total = sum([repeat * 5 * len(metric) if i < step else repeat * len(metric) for i in
                      generate_step(step, min(ori_otu.sum(0)))])
    # with tqdm.tqdm(pbar_total) as pbar:
    #pbar = progressbar.ProgressBar().start()
    df_b = defaultdict(list)
    for _idx, num in enumerate(tqdm.tqdm(generate_step(step,min(ori_otu.sum(0))))):
    # for _idx, num in enumerate(generate_step(step,min(ori_otu.sum(0)))):
        # subsampling step.
        current = defaultdict(list)

        manager = Manager()
        subsampled = manager.list()
        p = Pool(int(threads))

        if num < step:
            repeat_times = repeat *5
        else:
            repeat_times = repeat

        for _ in range(repeat_times):
            p.apply_async(multiprocess_subsample,
                               args=(subsampled,ori_otu, num, prebuild_dict))
        p.close()
        p.join()

        for _m in metric:
            current[_m] = [diversity_ana(_m, temp.values.tolist(), list(ori_otu.columns), otu_ids=list(ori_otu.index),
                                         tree=tree) for temp in subsampled]

            current[_m] = sum(current[_m])
            if subsampled:
                current[_m] = current[_m] / len(subsampled) # one 'step' is run 'repeat' time for smooth the curve.
            else:
                import pdb;pdb.set_trace()
                current[_m] = 0

            df_b[_m].append(df(current[_m]).T) # collect all the DataFrame which have the values in each step.

    result = {}
    for _m in df_b.keys():
        result[_m] = pandas.concat(df_b[_m]) # sum them up.
        result[_m].index = generate_step(step,min(ori_otu.sum(0))) # assign the index.

    #Drawing part
    for _m in result.keys():
        datas = []
        if metadata:
            for col in sorted(list(result[_m].columns)):
                # Draw each sample curve. col is sample. row is the step(reasd num)
                datas.append(go.Scatter(x=[0] + list(result[_m].loc[:, col].index),
                                        y=[0] + list(result[_m].loc[:, col]),
                                        mode='lines+markers',
                                        name=sample_info.loc[col, col_need],
                                        hoverinfo='text+name',
                                        marker=dict(size=3),
                                        line=dict(shape='spline')
                                        ))
        else:
            for col in sorted(list(result[_m].columns)):
                datas.append(go.Scatter(x=[0] + list(result[_m].loc[:, col].index),
                                        y=[0] + list(result[_m].loc[:, col]),
                                        mode='lines+markers',
                                        name=col,
                                        hoverinfo='text+name',
                                        marker=dict(size=3),
                                        line=dict(shape='spline')
                                        ))

        layout = go.Layout(
            title='Rarefaction curve(%s)' % _m.title(),
            xaxis = dict(title = 'Reads count'),
            yaxis = dict(title = '%s' % _m.title())
        )
        fig = go.Figure(data=datas, layout=layout)
        # put data and layout config into one.

        fn = output_dir + '/rarefaction_curve(%s)' % _m.title()
        # construct a dir to storge result[_m].
        if args.output_fig:
        # In case to reuse the data, we can use package pickle to storge the fig into file.
        # You can use pickle.load to reuse the fig data instead of plot it.
            pickle.dump(fig, args.output_fig)

        if not os.path.isfile(fn):
            plotly.offline.plot(fig, filename=fn)
        else:
            button = input(
                "There are existed html in your dir. If you want to overlap this. Enter 'Y/y',or enter a new filepath(include dir and file name).")
            if button.upper() == 'Y':
                plotly.offline.plot(fig, filename=fn)
            else:
                if os.path.isdir(os.path.dirname(button)):
                    plotly.offline.plot(fig, filename=button)

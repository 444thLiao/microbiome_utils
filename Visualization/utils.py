import pandas as pd
from pandas import DataFrame as df
from collections import defaultdict
import numpy as np

def read_data(fp):
    try:
        first = pd.read_csv(fp)
        if len(first.columns)==1 or len(first.columns) == 0 :
            first = pd.read_table(fp,index_col=0)
        if len(first.columns)==0:
            print 'Unresolved delimiter.'
    except:
        first = pd.read_excel(fp)

    return first



def get_otu_table(infile):
    """
    Receive a OTU table and return a
    :param infile: absolute file path
    :return: a tuple with two dict. first is {sample_name:{OTU:reads_count in this sample}},
                                    second is {OTU:{sample_name:otu reads count in this sample}}
    """
    samples = defaultdict(dict)
    otus = {}

    with open(infile) as infh:
        header = infh.next()
        assert header.startswith('#')
        header = header.strip()
        header = header.lstrip('#')
        sample_ids = header.split('\t')[1:]

        for line in infh:
            line = line.strip()
            if line:
                its = line.split('\t')
                otu_id = its[0]
                assert otu_id not in otus

                counts = its[1:]
                counts = [int(count) for count in counts]
                assert len(sample_ids) == len(counts)
                obs = zip(sample_ids,counts)
                obs = [it for it in obs if it[1]>0]

                otus[otu_id] = dict(obs)
                for it in obs:
                    samples[it[0]][otu_id] = it[1]
    return samples, otus

def get_tax(infile):
    """
    Receive a sintax.txt file return a dict.
    :param infile:
    :return: A dict {OTU:a string of taxonomy} string is like "d:Bacteria(1.0000),p:"Proteobacteria"(1.0000),c:Alphaproteobacteria(1.0000),o:Caulobacterales(1.0000),f:Caulobacteraceae(1.0000),g:Brevundimonas(0.9900)"
    """
    otu_tax = {}
    with open(infile) as infh:
        for line in infh:
            line = line.strip()
            if line:
                its = line.split('\t')
                if len(its) != 4:
                    continue
                otu_id = its[0].split(';size')[0]
                tax = its[1]
                strand = its[2]
                if strand == '+':
                    tax_pass = its[3]
                    assert otu_id not in otu_tax
                    otu_tax[otu_id] = tax
                else:
                    pass
    return otu_tax


def get_rank(tax):
    """
    Receive a string like "d:Bacteria(1.0000),p:"Proteobacteria"(1.0000),c:Alphaproteobacteria(1.0000),o:Caulobacterales(1.0000),f:Caulobacteraceae(1.0000),g:Brevundimonas(0.9900)"
    :param tax:
    :return: a dict {tax_rank:(tax_name, confidence)} like.{'d':('Bacteria',1.0000),'p':('Proteobacteria',1.0000)......}
    """
    ranks = tax.split(',')
    info = {}
    for it in ranks:
        rank, name = it.split(':')
        name, confidence = name.split('(')
        confidence = confidence.rstrip(')')
        assert rank not in info
        info[rank] = (name, float(confidence))
    return info


def tax_anno(otu_tax_file,otu_tab_file,output_anno,tax = 'f'):
    otu_tax = get_tax(otu_tax_file)
    otu_tab = defaultdict(list)
    with open(otu_tab_file) as infh:
        header = infh.next()
        header = header.rstrip()
        for line in infh:
            line = line.strip()
            if line:
                its = line.split('\t')
                otu_id = its[0]
                otu_abd = its[1:]
                otu_abd = np.asarray([int(it) for it in otu_abd])
                try:
                    rank_info = get_rank(otu_tax[otu_id])
                    rank = rank_info.get(tax)
                    if rank and rank[1] >= 0.8:
                        otu_tab[rank[0]].append(otu_abd)
                    else:
                        otu_tab['unclassfied(%s)' % tax].append(otu_abd)
                        #print rank_info
                except:
                    #import pdb;pdb.set_trace()
                    print 'missing %s, maybe self deleted it' % otu_id

    outfile = output_anno
    with open(outfile, 'w') as outfh:
        print >>outfh, header
        for tax, abds in otu_tab.items():
            abds = np.asarray(abds)
            abd = np.sum(abds, axis=0)
            outs = [tax] + list(abd)
            outs = [str(it) for it in outs]
            print >>outfh, '\t'.join(outs)

def norm_otu(OTU_table_with_tax,outfile):
    samples, otus = get_otu_table(OTU_table_with_tax)
    #outfile = '/home/liaoth/data2/16s/170821_170519_combined_16s/XK_analysis/otutab_raw_filtered_5k_e5_genus.anno.txt'
    with open(outfile, 'w') as outfh:
        otu_ids = otus.keys()
        otu_ids.sort()
        header = ['sample_id']+otu_ids
        print >>outfh, '\t'.join(header)
        sample_ids = samples.keys()
        sample_ids.sort()

        for sample_id in sample_ids:
            outs = [sample_id]
            obs = [samples[sample_id].get(otu_id,0) for otu_id in otu_ids]
            total = sum(obs)
            obs_relative = [it*1.0/total for it in obs]
            outs = outs + obs_relative
            outs = [str(it) for it in outs]
            print >>outfh, '\t'.join(outs)
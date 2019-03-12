from __future__ import print_function
import plotly,random,os
from plotly import graph_objs as go
from utils import read_data
import seaborn as sns



def parse_data(tax,tax_otus_index,tax_otu_sub,colors):
    bucket = []
    if tax != 'family':
        judge = False
    else:
        judge = True

    for idx,_i in enumerate(list(tax_otus_index)):
        # loop it and raw each OTU(taxed)
        bucket.append(go.Bar(
            x = tax_otu_sub.columns,
            y = tax_otu_sub.loc[_i,:]/tax_otu_sub.sum(0), # Normalization with each sample.
            visible = judge,
            name=_i,
            # name will display at legend.
            marker=dict(color=colors[idx],
                        line=dict(width=1,color='#FFFFFF')
                        # use white line to separate the bar, it can use it to make it more clear.
                        )
        ))
    return bucket

# def sort_way(taxotu_columns):
#     return sorted(list(taxotu_columns))

def sort_way(taxotu_columns):
    sort_for = [(k.rpartition('_')[-1],k) if '.' not in k else (k.rpartition('_')[-1].split('.')[-1],k) for k in list(taxotu_columns)]
    result = [k for i,k in sorted(sort_for)]
    return result

def generate_html():
    for tax_otutab in tax_otutabs:
        tax_names = os.path.basename(tax_otutab).split('_')[-1].split('.')[0]

        tax_otu = read_data(tax_otutab)
        #tax_otu = tax_otu.loc[tax_otu.sum(0)!=0,:]
        # tax_otu_sub = tax_otu_sub.loc[:, sorted(list(tax_otu_sub.columns))]
        # new_cols = [i for i in sorted(list(tax_otu.columns)) if '-N' in i] + [i for i in sorted(list(tax_otu.columns)) if '-T' in i]
        #tax_otu = tax_otu.loc[:,new_cols]
        tax_otu.loc[:,'sum_all'] = tax_otu.sum(1)
        # sum all samples reads in each OTU(taxed).
        tax_otu.sort_values('sum_all',inplace=True)
        # Use total reads to sort in order to make big OTU been drawn above.
        tax_otu_sub = tax_otu.iloc[:,:-1]
        # deleted the sum_all columns.
        tax_otu_sub = tax_otu_sub.loc[:, sort_way(tax_otu.columns)]
        # sort the columns according to the header.
        tax_otus_index = tax_otu_sub.index

        #processing lot of colors
        data = []
        colors = sns.husl_palette(len(tax_otus_index))
        random.shuffle(colors)
        # shuffle the colors to make it divergence
        colors = colors.as_hex()

        if tax_names == 'family':
            family_dis = parse_data(tax_names,tax_otus_index,tax_otu_sub,colors)
        elif tax_names == 'genus':
            genus_dis = parse_data(tax_names, tax_otus_index, tax_otu_sub, colors)
        elif tax_names == 'phylum':
            phylum_dis = parse_data(tax_names, tax_otus_index, tax_otu_sub, colors)
        elif tax_names == 'OTU':
            OTU_dis = parse_data(tax_names, tax_otus_index, tax_otu_sub, colors)
        else:
            print(tax_names)
            #exit()
    family_shows = [True] * len(family_dis) + [False] * len(genus_dis) + [False] * len(phylum_dis)
    genus_shows =  [False] * len(family_dis) + [True] * len(genus_dis) + [False] * len(phylum_dis)
    phylum_shows = [False] * len(family_dis) + [False] * len(genus_dis) + [True] * len(phylum_dis)
    OTU_shows = [False] * len(family_dis) + [False] * len(genus_dis) + [False] * len(phylum_dis)
    updatemenus = list([
        dict(type="buttons",
             active=1,
             buttons=list([
                 dict(label = 'phylum',
                     method = 'restyle',
                     args = ['visible', phylum_shows]),
                dict(label = 'family',
                     method = 'restyle',
                     args = ['visible', family_shows]),
                dict(label = 'genus',
                     method = 'restyle',
                     args = ['visible', genus_shows]),
                    # dict(label = 'OTU',
                    #  method = 'restyle',
                    #  args = ['visible', OTU_shows])
            ]),
        )
    ])


    layout = go.Layout(
        title='Taxonomic distribution Plus',
        barmode='stack',
        #height=2000,
        updatemenus=updatemenus

    )

    fig = go.Figure(data=family_dis+genus_dis+phylum_dis, layout=layout)
    plotly.offline.plot(fig, filename=input_dir+'/sepcies_distribution')

if __name__ == '__main__':
    import sys
    # Config
    metadata = None
    input_dir = '/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/analysis_storge/routine_analysis/taxonomy_report'
    base_name = 'otu_norm_filtered_25k_s2'
    tax_otutabs = ['%s/%s_family.txt' % (input_dir, base_name),
                   '%s/%s_genus.txt' % (input_dir, base_name),
                   '%s/%s_phylum.txt' % (input_dir, base_name),
                   # '%s/%s_OTU.txt' % (input_dir, base_name)
                   ]

    tax_levels = [os.path.basename(i).split('_')[-1].split('.')[0] for i in tax_otutabs]

    generate_html()
    # if len(sys.argv) >= 2:
    #     input_dir = os.path.abspath(sys.argv[-1])

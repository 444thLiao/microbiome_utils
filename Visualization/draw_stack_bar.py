import plotly,random,os
from plotly import graph_objs as go
from utils import read_data
import seaborn as sns

#Config
metadata = './example/stack_bar/sample_info.csv'
tax_otutab = './example/stack_bar/otu_raw_filterd_e5_genus.txt'
output_dir = './example/stack_bar/drawn/'
col_need = 'class'  # which columns in metadata you need to parse name in distance to.


tax_otu = read_data(tax_otutab)
tax_otu.loc[:,'sum_all'] = tax_otu.sum(1)
# sum all samples reads in each OTU(taxed).
tax_otu.sort_values('sum_all',inplace=True)
# Use total reads to sort in order to make big OTU been drawn above.
tax_otu_sub = tax_otu.iloc[:,:-1]
# deleted the sum_all columns.
tax_otu_sub = tax_otu_sub.loc[:,sorted(list(tax_otu_sub.columns))]
# sort the columns according to the header.
tax_otus_index = tax_otu_sub.index

if metadata and col_need:
    sample_info = read_data(metadata)
    new_names = list(sample_info.loc[list(tax_otu_sub.columns),col_need])
    tax_otu_sub.columns = new_names
    tax_otu_sub = tax_otu_sub.loc[:,sorted(new_names)]
# if metadata and col_need, we rename the header and sort it.


#processing lot of colors
data = []
colors = sns.husl_palette(len(tax_otus_index))
random.shuffle(colors)
# shuffle the colors to make it divergence

colors = colors.as_hex()
for idx,_i in enumerate(list(tax_otus_index)):
    # loop it and raw each OTU(taxed)
    data.append(go.Bar(
        x = tax_otu_sub.columns,
        y = tax_otu_sub.loc[_i,:]/tax_otu_sub.sum(0), # Normalization with each sample.
        name=_i,
        # name will display at legend.
        marker=dict(color=colors[idx],
                    line=dict(width=1,color='#FFFFFF')
                    # use white line to separate the bar, it can use it to make it more clear.
                    )
    ))

layout = go.Layout(
    title='Species distribution (%s)' % os.path.basename(tax_otutab),
    barmode='stack',
    height=2000
)
fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig, filename=output_dir+'/sepcies_distribution')


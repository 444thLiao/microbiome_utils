import os
from utils import tax_anno,norm_otu
import argparse
USEARCH = '/usr/bin/usearch'
RDP_DB = '~/data2/rdp_16s_v16.fa'
draw_PD = "../Visualization/draw_PD.py"
def regular_analysis(OTU_table,rep_fa,outputdir,draw_pd = False):
    """
    OTU_TABLE is normalize and filtered.
    :param OTU_table:
    :param rep_fa:
    :param outputdir:
    :return:
    """
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    os.makedirs(outputdir+'/beta')
    os.makedirs(outputdir+'/alpha')
    os.makedirs(outputdir + '/rarefaction')
    os.makedirs(outputdir + '/taxonomy_report')
    print '%s -cluster_agg %s -treeout %s' % (USEARCH,
                                                  rep_fa,
                                                  os.path.join(outputdir,'otus.tree'))
    os.system('%s -cluster_agg %s -treeout %s' % (USEARCH,
                                                  rep_fa,
                                                  os.path.join(outputdir,'otus.tree')))
    os.system('%s -alpha_div %s -output %s' % (USEARCH,
                                               OTU_table,
                                               os.path.join(outputdir+'/alpha', 'alpha.txt')))
    os.system('%s -beta_div %s -tree %s -filename_prefix %s' % (USEARCH,
                                                                OTU_table,
                                                                os.path.join(outputdir, 'otus.tree'),
                                                                outputdir + '/beta/'))
    os.system("%s -sintax %s -db %s -strand both -tabbedout %s -sintax_cutoff 0.8" % (USEARCH,
                                                                                      rep_fa,
                                                                                      RDP_DB,
                                                                                      os.path.join(outputdir,'sintax.txt')))
    os.system("%s -alpha_div_rare %s -output %s/rare.txt" % (USEARCH,
                                                             OTU_table,
                                                             outputdir + '/rarefaction'))
    for taxs in [('g','genus'),('p','phylum'),('f','family')]:
        tax_summary_out_file = os.path.join(outputdir + '/taxonomy_report', '%s_%s.txt' % (os.path.basename(OTU_table).split('.')[0],taxs[1]))

        tax_anno(os.path.join(outputdir,'sintax.txt'),
                 OTU_table,
                 tax_summary_out_file,
                 tax=taxs[0])

        norm_otu(tax_summary_out_file,
                 tax_summary_out_file.replace('.txt','.norm.txt'))
    if draw_pd:
        os.system("python %s -t %s -i %s -o %s -M shannon,observed_otus,faith_pd" % (draw_PD,
                                                   os.path.join(outputdir, 'otus.tree'),
                                                   OTU_table,outputdir + '/rarefaction'))

    os.system('cp %s %s' % (OTU_table,outputdir+'/'))
    os.system('cp %s %s' % (rep_fa, outputdir + '/'))

# regular_analysis('/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/otu_norm_filtered_40k_s2.tab',
#                  '/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/otus_rep.fa',
#                  '/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/routine_analysis/')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='otutab', type=str, required=True,
                        help="OTU table path. Make sure you have filtered.")
    parser.add_argument('-fa',  dest='fasta', type=str, required=True,
                        help="OTU represent sequence fasta file path.")
    parser.add_argument('-o', '--output', dest='odir', type=str, required=True,
                        help="output Dir path")
    parser.add_argument('-rc', dest='rarefaction', action='store_true',
                        help="if you want to draw a rarefaction curve(which will take long time.But it will have a progress bar to display its.)")
    args = parser.parse_args()

    otu_tab = args.otutab
    fasta = args.fasta
    odir = args.odir
    draw_rc = args.rarefaction

    regular_analysis(otu_tab,fasta,odir,draw_pd=draw_rc)

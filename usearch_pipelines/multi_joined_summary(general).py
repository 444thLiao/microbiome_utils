"""
170905
written by: TianHua Liao

This script is tend to used after multiplt_joined.
There have four purpose in this script.
    1.remove unjoined fastq file incase mixed in after pipelines.
    2.rename the complex sample name into simpify one.
    3.remove the tags from fastq_screen tools in case affect after pipelines.
    4.calculate the reads maintain_ratio after joined process.

"""
import glob,os,re,argparse
from pandas import DataFrame as df

#samplename_pattern = '(XK-[0-9]+F-?[0-9]{0,3}).*_R[12]'
#samplename_pattern = re.findall('((NYN|NYT|TAN)-[0-9]+-?[NT]{0,1}-?[0-9]{0,2}).*_R[12]',path)[0][0]

def each_join_summary(dir_path):
    joined_file = dir_path+'/fastqjoin.join.fastq'
    unjoined_file1 = dir_path+'/fastqjoin.un1.fastq'
    unjoined_file2 = dir_path + '/fastqjoin.un2.fastq'

    j_reads_num = int(os.popen("grep -c '^+$' %s" % joined_file).read())
    unj_reads_num1 = int(os.popen("grep -c '^+$' %s" % unjoined_file1).read())
    unj_reads_num2 = int(os.popen("grep -c '^+$' %s" % unjoined_file2).read())
    if unj_reads_num1 != unj_reads_num2:
        print 'WARNING: R1 and R2 has different reads num.'
    return float(j_reads_num)/(j_reads_num+(unj_reads_num1+unj_reads_num2)/2),j_reads_num,j_reads_num+(unj_reads_num1+unj_reads_num2)/2

def joined_summary(joined_result_dir,output_csv,samplename_pattern):
    """
    :param joined_result_dir: the dir which is multiple_join_paired_ends.py outputed without last '/'.
    :param output_csv: the summary file path you want to storge.
    :return: None
    """
    sample_names = []
    sample_joined_summary = []
    joined_reads = []
    total_reads = []
    for each in glob.glob(joined_result_dir+'/*'):
        if os.path.isdir(each):
            sample_name = os.path.basename(each)
            simply_sn = re.findall(samplename_pattern, sample_name)[0]
            sample_names.append(simply_sn)
            sample_joined_summary.append(each_join_summary(each)[0]*100)
            joined_reads.append(each_join_summary(each)[1])
            total_reads.append(each_join_summary(each)[2])
    result = df(data={'sample name':sample_names,'joined per':sample_joined_summary,
                      'joined read num':joined_reads,
                      'ori total read num':total_reads})
    with open(output_csv,'w') as f1:
        result.loc[:,['sample name','joined read num','ori total read num','joined per']].to_csv(f1,index=False)

def rename_delete(joined_result_dir,samplename_pattern):
    os.system('rm %s/*/*.un[12].fastq' % joined_result_dir)
    for each in glob.glob(joined_result_dir + '/*'):
        if os.path.isdir(each):
            sample_name = os.path.basename(each)
            simplify_sn = re.findall(samplename_pattern, sample_name)[0]
            try:
                os.rename(joined_result_dir+'/'+sample_name,joined_result_dir+'/'+simplify_sn)
            except:
                print each
            with open(joined_result_dir + '/' + simplify_sn+'/fastqjoin.join.fastq') as data:
                buffer_str = ''
                for i in data:
                    if i.startswith('@') and '#FQST' in i:
                        buffer_str+=i.rpartition('#FQST')[0]+'\n'
                    else:
                        buffer_str += i
            with open(joined_result_dir + '/' + simplify_sn + '/%s.fastq' % simplify_sn,'w') as f1:
                f1.write(buffer_str)
    os.system("rm %s" % joined_result_dir + '/*/'+'/fastqjoin.join.fastq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='jrd', type=str, required=True,
                        help="joined process output path")
    parser.add_argument('-o', '--output', dest='oc', type=str, required=True,
                        help="csv output file path")
    parser.add_argument('-pm', dest='pm', type=str, required=False,
                        help="project name")
    parser.add_argument('-pattern', dest='pattern', type=str, required=False,default='(.*)',
                        help="name pattern will storged in summary and renamed with. [Default] is the ori name.")
    parser.add_argument('-updated', dest='updated_or_not', action='store_true',
                        help="if you want to updated the filename and delete unjoined file.")
    args = parser.parse_args()
    joined_result_dir = os.path.abspath(args.jrd)
    pattern = args.pattern
    output_csv = os.path.abspath(args.oc)


    if not os.path.isfile(output_csv):
        # import pdb;pdb.set_trace()
        test = raw_input("Rename from %s to %s. \nIf you sure about this,please enter Y/y:\n If you want to change this re pattern please use -pm." % (glob.glob(joined_result_dir + '/*')[0],re.findall(pattern, glob.glob(joined_result_dir + '/*')[0])[0]))
        if test.upper() == 'Y':
            joined_summary(joined_result_dir,output_csv,pattern)
        else:
            print 'Exit now.'
            exit(0)
    else:
        print 'Output file already exist, please make sure and delete it.'
        pass
    if args.updated_or_not:
        rename_delete(joined_result_dir,pattern)
    else:
        print 'Processing completed, exit now...'
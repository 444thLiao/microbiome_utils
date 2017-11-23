import os,glob,argparse,re
from pandas import DataFrame as df


#For joined screened accessment

def each_join_summary(ori_dir,aft_dir,suffix):
    ori = os.path.join(ori_dir,suffix)
    after = os.path.join(aft_dir,suffix)

    ori_reads_num = int(os.popen("zgrep -c '^+$' %s/*.fastq" % ori).read())
    aft_reads_num = int(os.popen("zgrep -c '^+$' %s.tagged_filter.fastq" % after).read())
    # if ori_reads_num == 0:
    #     print ori,after
    return float(aft_reads_num)/(ori_reads_num),aft_reads_num,ori_reads_num

def reads_num_summary(raw_dir,screened_dir,output_csv):
    """
    :param joined_result_dir: the dir which is multiple_join_paired_ends.py outputed without last '/'.
    :param output_csv: the summary file path you want to storge.
    :return: None
    """
    sample_names = []
    sample_processed_summary = []
    processed_reads = []
    ori_reads = []
    for each in glob.glob(raw_dir+'/*'):
        if os.path.isdir(each):
            sample_name = os.path.basename(each)
            sample_names.append(sample_name)
            r_nums = each_join_summary(raw_dir,screened_dir,sample_name)
            sample_processed_summary.append(r_nums[0]*100)
            processed_reads.append(r_nums[1])
            ori_reads.append(r_nums[2])
    result = df(data={'sample name':sample_names,'screened per':sample_processed_summary,
                      'screened read num':processed_reads,
                      'ori read num':ori_reads})
    #return result
    with open(output_csv,'w') as f1:
        result.loc[:,['sample name','screened read num','ori read num','screened per']].to_csv(f1,index=False)


def rename_delete(screened_dir):
    for each in glob.glob(screened_dir + '/*.tagged_filter.fastq'):
        sample_name = os.path.basename(each)
        simplify_sn = sample_name.replace('.tagged_filter','')
        os.rename(each,screened_dir+'/'+simplify_sn)
        with open(screened_dir + '/' + simplify_sn) as data:
            buffer_str = ''
            for i in data:
                if i.startswith('@') and '#FQST' in i:
                    buffer_str+=i.rpartition('#FQST')[0]+'\n'
                else:
                    buffer_str += i
        with open(screened_dir + '/' + simplify_sn,'w') as f1:
            f1.write(buffer_str)
    os.system("rm %s" % screened_dir+'/*.tagged.fastq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--joined', dest='jo', type=str, required=True,
                        help="multiple joined process output path")
    parser.add_argument('-s', '--screened', dest='scr', type=str, required=True,
                        help="screened process output path")
    parser.add_argument('-o', '--output', dest='oc', type=str, required=True,
                        help="csv output file path")
    parser.add_argument('-updated', dest='updated_or_not', action='store_true',
                        help="if you want to updated the filename and delete tagged file.")
    args = parser.parse_args()

    screened_result_dir = args.scr
    joined_result_dir = args.jo
    output_csv = args.oc

    if not os.path.isfile(output_csv):
        reads_num_summary(joined_result_dir,screened_result_dir,output_csv)

    if args.updated_or_not:
        print 'Starting delete and rename......'
        rename_delete(screened_result_dir)
        print 'Completeing.'
    else:
        print 'Processing completed, exit now...'
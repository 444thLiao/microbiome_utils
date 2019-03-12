

import glob,os,tqdm

raw_data_dir = '/home/liaoth/data_bank/t2d_from_HMP/16s_raw_data/'
all_file = glob.glob(raw_data_dir+'*.fastq.gz')
all_id = list(set([_i.split('/')[-1].split('_')[0] for _i in all_file]))


for sample_id in tqdm.tqdm(all_id):
    paths = sorted(glob.glob('%s/%s_*' % (raw_data_dir,sample_id)))
    des_path = os.path.dirname(paths[0])

    file_names = [os.path.basename(path).replace('.fastq.gz','') for path in paths]

    odir = raw_data_dir+'trimed_20/'
    if not os.path.isdir(odir):
        os.makedirs(odir)
    os.system(
'''java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100 '''.format(base_in=des_path,input1=file_names[0],input2 = file_names[1],base_out=odir,output=odir+'%s.log' % sample_id))




from multiprocessing import Process

def run_trim(*args):
    des_path, file_names, odir, sample_id = args
    # print '''java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100 '''.format(
    #         base_in=des_path, input1=file_names[0], input2=file_names[1], base_out=odir,
    #         output=odir + '%s.log' % sample_id)
    os.system(
        '''java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100 '''.format(
            base_in=des_path, input1=file_names[0], input2=file_names[1], base_out=odir,
            output=odir + '%s.log' % sample_id))

# raw_data_dir = '/home/liaoth/data2/16s/171027_16s/raw_data/second_splited/'
raw_data_dir = '/home/liaoth/temp_/ordered/'
all_file = glob.glob(raw_data_dir+'*.fastq.gz')
all_id = list(set([_i.split('/')[-1].split('_')[0] for _i in all_file]))

for sample_id in all_id:
    paths = sorted(glob.glob('%s/%s_*' % (raw_data_dir,sample_id)))
    des_path = os.path.dirname(paths[0])

    file_names = [os.path.basename(path).replace('.fastq.gz','') for path in paths]

    odir = raw_data_dir+'trimed_20/'
    if not os.path.isdir(odir):
        os.makedirs(odir)
    p = Process(target=run_trim, args=(des_path,file_names,odir,sample_id))
    p.start()
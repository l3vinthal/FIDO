import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import subprocess
import time

from env import *

def read_fasta(filename):
  output = {}
  with open(filename, 'r') as file:
    for line in file:
      if ">" in line:
        key = line.replace('\n','')
        output[key] = ''
      else:
        output[key] += line.replace('\n','')
  return output

def insert_newlines(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

def fasta_dict_to_file(fasta_dic, output_filename):
    with open(output_filename,'w') as output:
        for i in fasta_dic:
            try:
                #title = i.replace(">",'').replace("sp","").split("|")[1]   
                #output.write(">"+title+"\n"+fasta_dic[i]+"\n")
                
                #Use these next lines for NCBI sequence data.
                title = i.replace(">",'').replace("sp","").split(" ")[0]   
                output.write(">"+title+"\n"+fasta_dic[i]+"\n")

            except:
                title = i.replace(">",'').replace("sp","").split("|")[0]  
                output.write(">"+title+"\n"+fasta_dic[i]+"\n")                
            
def filter_seqs(fasta_filename, output_dir, shortest = MIN_SEQ_LEN, longest = MAX_SEQ_LEN):
    data = read_fasta(fasta_filename)
    output = {}
    removed = 0
    kept = 0
    for i in data:
        if len(data[i]) > shortest and len(data[i]) < longest and not('X' in data[i]):
            output[i] = insert_newlines(data[i].replace('\n',''), every=80)  
            kept += 1
        else:
            removed += 1
    print ('Sequences kept:', kept, 'Sequences removed:', removed)
    print ('Average sequence length:', np.average([len(i) for i in list(output.values())]).astype('int'))
    fasta_dict_to_file(output, 
                       output_filename=output_dir+"step_1_filtered.fasta",
                      )
    return output_dir+"step_1_filtered.fasta"
    
def blast_filter(fasta_filename, ref_fasta, output_dir,
                pident_min = MIN_IDEN,
                pident_max = MAX_IDEN,
                qcovhsp_min = QUERY_COVERAGE,
                blast_bin = ''):
    
    # Create the output directory for BLAST data    
    _ = subprocess.run(['mkdir ' + output_dir + 'blastp_data'], stderr=subprocess.STDOUT,shell=True)
    
    # Create a BLAST database from the input FASTA file
    _ = subprocess.run([blast_bin+'makeblastdb -in '+fasta_filename+' -dbtype prot -out ' + output_dir+"blastp_data/blast_fasta_db"], stderr=subprocess.STDOUT,shell=True)
    
    # Run BLASTP with the specified parameters
    _ = subprocess.run([blast_bin+'blastp -db ' + output_dir + 'blastp_data/blast_fasta_db -out ' + output_dir +'blastp_data/blast_alignments.log -query ' + ref_fasta + ' -outfmt "10 sseqid sacc evalue pident nident qcovhsp qcovs" -num_alignments 1000000'], stderr=subprocess.STDOUT,shell=True)
    _ = subprocess.run([blast_bin+'blastp -db ' + output_dir + 'blastp_data/blast_fasta_db -out ' + output_dir +'blastp_data/blast_alignments_verbose.log -query ' + ref_fasta], stderr=subprocess.STDOUT,shell=True)   
    
    # Parse BLASTP output and filter sequences    
    seqs = read_fasta(fasta_filename)  
    filtered_seqs = {}
    with open(output_dir+'blastp_data/blast_alignments.log') as alignment_log:
        for line in alignment_log:
            line = line.replace('\n','').split(',')
            if (float(line[3]) > pident_min) and (float(line[3]) < pident_max) and (float(line[5]) > qcovhsp_min):
                filtered_seqs[line[0]] = seqs[">" + line[0].replace('>','')]   
    #Extract output.  

    print ('Alignment complete. Writing to ' + output_dir+"blastp_filtered.fasta")
    
    fasta_dict_to_file(filtered_seqs, 
                   output_filename=output_dir+"step_2_blastp_filtered.fasta",
                  )
    
    return output_dir+"step_2_blastp_filtered.fasta"

def run_mmseq(fasta_filename, output_dir,
             min_seq_id = 0.7):

    # Create the output directory for BLAST data    
    _ = subprocess.run(['mkdir ' + output_dir + 'mmseq_hmmer_data'], stderr=subprocess.STDOUT,shell=True)
    _ = subprocess.run(['mmseqs easy-cluster ' + fasta_filename + ' ' + output_dir+'mmseq_hmmer_data/step_3_temp_mmseq tmp --min-seq-id ' + str(min_seq_id)], stderr=subprocess.STDOUT,shell=True)
    print ('Mseq2 clustering complete. Writing data to ' + output_dir)
    return output_dir + 'mmseq_hmmer_data/step_3_temp_mmseq_rep_seq.fasta', output_dir + 'mmseq_hmmer_data/step_3_temp_mmseq_cluster.tsv'

def clustalo(fasta_filename, output_dir, clustalo):
    # Extract the directory path from the FASTA filename

    _ = subprocess.run([clustalo + ' -i ' + fasta_filename + ' --dealign -o ' + output_dir + 'step_4_rep_seq_aligned.fasta --outfmt=fasta -force'], stderr=subprocess.STDOUT,shell=True)
    return output_dir + 'step_4_rep_seq_aligned.fasta'
    
def hmmer_build_and_align(rep_fasta_filename, output_dir, fasta_full_db, hmmer_bin):
    # Extract the directory path from the FASTA filename

    _ = subprocess.run([hmmer_bin + '/hmmbuild ' + output_dir+'rep_seq_hmm ' +  rep_fasta_filename], stderr=subprocess.STDOUT,shell=True)    
    _ = subprocess.run([hmmer_bin + '/hmmalign --outformat afa -o ' + output_dir + 'full_db_alignment.fasta ' + output_dir + 'rep_seq_hmm ' + fasta_full_db], stderr=subprocess.STDOUT,shell=True)     
    print ('HMM-based alignment of fasta database is complete.')
    return output_dir + 'full_db_alignment.fasta'

def pad_seqs(fasta_dic, uppercase=False):

    output = pd.DataFrame(fasta_dic.items())
    if uppercase:
        output[1] = output[1]#.str.upper()
    output = output.rename(columns={0: "accession", 1: "sequence"})
    max_len = np.max([len(i) for i in output.sequence.values])
    padded_seqs = [''.join(i+'-'*(max_len-len(i))) for i in output.sequence.values]
    output['sequence']=padded_seqs
    return output, max_len
    
def build_dataset(fasta_filename, clust_file, output_dir):
    data = read_fasta(fasta_filename)
    output, max_len = pad_seqs(data)

    assert max_len == np.max([len(i) for i in output.sequence.values]), 'Incorrect length detected in preprocessing.build_dataset. Check padding'

    output['cluster_ID'] = ''

    clusters = pd.read_table(clust_file, header=None)
    clusters = clusters.rename(columns={0: "cluster", 1: "member"})
    define_cluster_parent = {}
    for i in clusters.iterrows():
        define_cluster_parent[i[1]['member']] = i[1]['cluster']
    
    for i,j in enumerate(output.iterrows()):
        output.loc[i,'cluster_ID'] = define_cluster_parent[output.loc[i,'accession'].replace('>','')]
    output.to_csv(output_dir+'final_alignment_with_clust_ids.csv')

   
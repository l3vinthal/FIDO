import src.preprocessing as pp
import argparse
import subprocess 
from env import *

def main():
    parser = argparse.ArgumentParser(prog='preprocessor',
                                    description='Preprocess fasta file for deep learning application')
    parser.add_argument('-in', '--input', type=str, help='Input fasta file', required=True)
    parser.add_argument('-ref', '--reference', type=str, help='Input reference fasta to align to', required=True)
    #parser.add_argument('-co', '--clustalomega_dir', type=str, help='Directory of Clustal Omega executable', required=True)
    #parser.add_argument('-hmmer', '--hmmer_dir', type=str, help='Bin of hmmer executables', required=True)
    parser.add_argument('-out', '--output', type=str, help='Output directory', required=True)    
    parser.add_argument('-min_seq_id', '--min_seq_id', type=str, help='Minimum sequence identity for clustering using MMseq2', required=False)
    parser.add_argument('-proj_name', '--proj_name', type=str, help='Project name for files (Optional)', required=False)
    parser.add_argument('-add_ref', dest='ref_file', help='Add reference sequence from fasta file to representative sequences after MMseq clustering', required=False)
    args = parser.parse_args()

    try:
        _ = subprocess.run(['mkdir ' + args.output], stderr=subprocess.STDOUT,shell=True)

        #if args.min_seq_id
        fasta_filename = pp.blast_filter(pp.filter_seqs(args.input, args.output), args.reference, args.output)

        pp.add_seq_to_fasta(args.ref_file, fasta_filename)

        rep_seq_filename, cluster_filename = pp.run_mmseq(fasta_filename, args.output, min_seq_id = 0.9)

        pp.add_seq_to_fasta(args.ref_file, rep_seq_filename)
        
        rep_seq_aligned_filename = pp.clustalo(rep_seq_filename, args.output, CLUSTAL_O)  

        #Add ref sequence to full db.
        #pp.add_seq_to_fasta(args.ref_file, fasta_filename)

        hmmalignment_fasta_filename = pp.hmmer_build_and_align(rep_seq_aligned_filename, args.output, fasta_filename, HMMER_BIN)
        
        pp.build_dataset(hmmalignment_fasta_filename, cluster_filename, args.output)

    except Exception as e:
        print ("An unexpected error occured")
        raise
        

if __name__ == "__main__":
    main()
    

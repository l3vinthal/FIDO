#Program directories
HMMER_BIN = "/home/ubuntu/cookstore/FIDO/hmmer/bin" #Bin directory of hmmer ex. /home/ubuntu/programs/hmmer/bin
CLUSTAL_O = "/home/ubuntu/cookstore/FIDO/clustalo"  #Directory of clustalo executable ex. /home/ubuntu/programs/clustalo

#MSA Processing
MIN_SEQ_LEN = 10
MAX_SEQ_LEN = 1000

MIN_IDEN = 30
MAX_IDEN = 99 #Set greater than 100 to include identical sequences.
QUERY_COVERAGE = 50

## MMSeq2
MIN_IDEN_CLUST = 0.90
from Bio import Blast
import requests
import re
import time
from typing import List, Union, Dict, Optional

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
import os

class BlastP():
    def __init__(self, seq, output_file, hitnum = 500):
        Blast.email = "cookie2004@gmail.com"
        print ('Emailing issues to', Blast.email)
        self.hitnum = hitnum
        self.data = self.protein_blast(seq)

        #Write the stream to an XML file.
        with open("my_blast.xml", "wb") as out_stream:
            out_stream.write(self.data.read())

        file_to_delete = "my_blast.xml"
        #Read the generated XML file.
        blast_record = Blast.read(file_to_delete)

        #Read accession numbers and fetch the full sequences.
        self.seqs = fetch_full_sequences([i.target.id.split('|')[1] for i in blast_record])
        save_sequences_to_fasta(self.seqs, output_file)

        try:
            os.remove(file_to_delete)
            print(f"File '{file_to_delete}' deleted successfully.")
        except FileNotFoundError:
            print(f"Error: File '{file_to_delete}' not found.")
        except Exception as e:
            print(f"An error occurred: {e}")


    def protein_blast(self, seq):
        return Blast.qblast("blastp", "nr", seq, alignments = self.hitnum, hitlist_size=self.hitnum, expect = 0.0001, descriptions=self.hitnum)

def save_sequences_to_fasta(sequences, output_file):
    """
    Save sequences to a FASTA file.
    
    Args:
        sequences: Dictionary mapping accessions to sequence data
        output_file: Path to save the FASTA file
    """
    with open(output_file, 'w') as f:
        for acc, data in sequences.items():
            f.write(f">{data['id']} {data['description']}\n")
            
            # Write sequence in chunks of 80 characters
            seq = data['sequence']
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")
    
    print(f"Saved {len(sequences)} sequences to {output_file}")

def fetch_full_sequences(accessions, batch_size=50, db="protein"):
    """
    Fetch full sequences for a list of accession numbers using Entrez.
    
    Args:
        accessions: List of accession numbers
        batch_size: Number of sequences to fetch in each batch
        db: Entrez database to use (protein or nucleotide)
        
    Returns:
        Dictionary mapping accession numbers to full sequences
    """
    print(f"Fetching {len(accessions)} full sequences from NCBI...")
    sequences = {}
    
    # Process in batches to avoid overloading NCBI servers
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        
        try:
            # Fetch the sequences
            handle = Entrez.efetch(
                db=db,
                id=",".join(batch),
                rettype="fasta",
                retmode="text"
            )
            
            # Parse the sequences
            records = list(SeqIO.parse(handle, "fasta"))
            
            # Store in dictionary
            for record in records:
                # Extract the accession from the ID
                acc = record.id.split("|")[-1] if "|" in record.id else record.id
                sequences[acc] = {
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq)
                }
                
            handle.close()
            
            # Be nice to NCBI servers
            time.sleep(1)
            
            print(f"Fetched batch {i//batch_size + 1}/{(len(accessions)-1)//batch_size + 1}")
            
        except Exception as e:
            print(f"Error fetching batch {i//batch_size + 1}: {str(e)}")
            # Wait a bit longer before retrying
            time.sleep(5)
    
    return sequences

if __name__ == "__main__":
    import argparse
    import sys
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Run BLAST search and retrieve full sequences.')
    parser.add_argument('input', type=str, help='Input sequence string or path')
    parser.add_argument('-o', '--output', type=str, default="blast_full_sequences.fasta", 
                        help='Output FASTA file name (default: blast_full_sequences.fasta)')
    parser.add_argument('-num_aligns', '--number_of_alignments', type=str, default="500", 
                        help='Max number of alignments to generate.')
        
    args = parser.parse_args()

    print(f"Running BLAST search with sequence of length {len(args.input)}...")
    blaster = BlastP(args.input, args.output, hitnum=int(args.number_of_alignments))

    print(f"Successfully retrieved and saved {len(blaster.seqs)} sequences to {args.output}")


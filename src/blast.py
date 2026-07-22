import os
import time
import argparse
from typing import List, Dict, Any

# Correct modern Biopython imports
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO

class BlastP():
    def __init__(self, seq: str, output_file: str, hitnum: int = 500):
        # Entrez needs the email setup directly
        Entrez.email = "myemail@gmail.com"
        print('Emailing NCBI issues to:', Entrez.email)
        
        self.hitnum = hitnum
        print("Submitting BLAST query to NCBI... (This may take a minute)")
        self.data = self.protein_blast(seq)

        # Write the XML stream to a file.
        file_to_delete = "my_blast.xml"
        with open(file_to_delete, "w") as out_stream:
            out_stream.write(self.data.read())
        self.data.close()

        # Read the generated XML file using modern NCBIXML parser
        with open(file_to_delete, "r") as blast_file:
            blast_record = NCBIXML.read(blast_file)

        # Extract accessions cleanly from alignments
        accessions = []
        for alignment in blast_record.alignments:
            # alignment.accession pulls out the clean NCBI accession directly (e.g., 'NP_001011')
            if alignment.accession:
                accessions.append(alignment.accession)

        if not accessions:
            print("No hits found or accessions could not be parsed.")
            self.seqs = {}
            return

        # Fetch full sequences using Entrez and save them
        self.seqs = fetch_full_sequences(accessions[:self.hitnum])
        save_sequences_to_fasta(self.seqs, output_file)

        # Cleanup
        try:
            os.remove(file_to_delete)
            print(f"Temporary file '{file_to_delete}' deleted successfully.")
        except FileNotFoundError:
            print(f"Error: File '{file_to_delete}' not found.")
        except Exception as e:
            print(f"An error occurred during cleanup: {e}")

    def protein_blast(self, seq: str):
        # NCBIWWW.qblast handles the web request directly
        return NCBIWWW.qblast(
            program="blastp", 
            database="nr", 
            sequence=seq, 
            alignments=self.hitnum, 
            hitlist_size=self.hitnum, 
            expect=0.0001
        )

def save_sequences_to_fasta(sequences: Dict[str, Any], output_file: str):
    """Save sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for acc, data in sequences.items():
            f.write(f">{data['id']} {data['description']}\n")
            
            # Write sequence in chunks of 80 characters
            seq = data['sequence']
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")
    
    print(f"Saved {len(sequences)} sequences to {output_file}")

def fetch_full_sequences(accessions: List[str], batch_size: int = 50, db: str = "protein") -> Dict[str, Any]:
    """Fetch full sequences for a list of accession numbers using Entrez."""
    print(f"Fetching {len(accessions)} full sequences from NCBI...")
    sequences = {}
    
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
            
            for record in records:
                acc = record.id.split("|")[-1] if "|" in record.id else record.id
                sequences[acc] = {
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq)
                }
                
            handle.close()
            
            # Respect NCBI servers
            time.sleep(1)
            print(f"Fetched batch {i//batch_size + 1}/{(len(accessions)-1)//batch_size + 1}")
            
        except Exception as e:
            print(f"Error fetching batch {i//batch_size + 1}: {str(e)}")
            time.sleep(5)
    
    return sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run BLAST search and retrieve full sequences.')
    parser.add_argument('input', type=str, help='Input sequence string or path')
    parser.add_argument('-o', '--output', type=str, default="blast_full_sequences.fasta", 
                        help='Output FASTA file name (default: blast_full_sequences.fasta)')
    parser.add_argument('-num_aligns', '--number_of_alignments', type=int, default=500, 
                        help='Max number of alignments to generate.')
        
    args = parser.parse_args()

    print(f"Running BLAST search with sequence of length {len(args.input)}...")
    blaster = BlastP(args.input, args.output, hitnum=args.number_of_alignments)

    print(f"Successfully retrieved and saved {len(blaster.seqs)} sequences to {args.output}")
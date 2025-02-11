from Bio import Entrez , SeqIO
from Bio.Seq import Seq
import os
Entrez.email = "sidharth@random.com"
file_path = "EcoliFasta"
def fetch_fasta(record_id, output_file):
    try:
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        with open(output_file, "w") as file:
            file.write(fasta_data)

        print(f"FASTA for {record_id} saved to {output_file}")

    except Exception as e:
        print(f"Error fetching {record_id}: {e}")

fetch_fasta("NC_000913.3",file_path)

def analyze_fasta(file_path):
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = record.seq.upper()
        seq_length = len(sequence)

        a_count = sequence.count("A")
        t_count = sequence.count("T")

        at_richness = (a_count + t_count) / seq_length * 100 if seq_length > 0 else 0

        print(f"Record ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Length: {seq_length}")
        print(f"AT Richness: {at_richness:.2f}%")
        print("-" * 40)
        return at_richness

genomeAT = analyze_fasta(file_path)

def motifRichness(motif):
    length = len(motif)
    cost = 0
    for c in motif:
        if c=='A' or c=='T':
            cost+=1
    return (cost/length)*100 if length>0 else 0

def motifCheck(motif , file_path,genomeAT):
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            print(f"Currently checking {motif} motif")
            motifAT = motifRichness(motif)
            print(f"AT richness of this motif is {motifAT:.2f}%")
            print(f"It has {(motifAT-genomeAT) :.2f}% higher AT content compared to the EColi K12 Genome")
            occurences = 0
            genome_seq = str(record.seq)
            for i in range(len(genome_seq) - len(motif) + 1):
                if genome_seq[i:i+len(motif)] == motif:
                    print(f"Motif found at position {i+1}")
                    occurences+=1
    print(f"This motif {motif} had {occurences} occurences in the entire Genome")                
    print("-" * 40)

motiflist = ["TTATACACA","TTATTCACA","TTATGCACA","TTATCCACA"]
for motif in motiflist:
    motifCheck(motif,file_path,genomeAT)

print("Successfully Executed")

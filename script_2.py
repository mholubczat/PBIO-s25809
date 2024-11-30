from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import csv

iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
test_seq = next(iterator)


def nucleotide_counts(record):
    return {"A": record.seq.count("A"),
            "C": record.seq.count("C"),
            "G": record.seq.count("G"),
            "T": record.seq.count("T")}

print(nucleotide_counts(test_seq))

def gc_content(record):
    return "{:.2f}%".format(gc_fraction(record) * 100)

print(gc_content(test_record))

def rev_complement(record):
    return str(record.seq.reverse_complement())

print(rev_complement(test_record))



with open('../../../../../Desktop/pythonProject/sequence_analysis.csv', 'w', encoding='utf-8') as output:
    writer = csv.writer(output)
    writer.writerow(["SeqID", "Nucleotide_Counts", "GC_Content", "Motif_Positions", "Reverse_Complement", "Translation_Lengths"])
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        writer.writerow([seq_record.id,
                        nucleotide_counts(seq_record),
                        gc_content(seq_record),

                         ])

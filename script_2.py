from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import csv

iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
test_record = next(iterator)
print(test_record)


def nucleotide_counts(record):
    return {"A": record.seq.count("A"),
            "T": record.seq.count("T"),
            "C": record.seq.count("C"),
            "G": record.seq.count("G")}


print(nucleotide_counts(test_record))


def gc_content(record):
    return "{:.2f}%".format(gc_fraction(record) * 100)


print(gc_content(test_record))


def rev_complement(record):
    return record.seq.reverse_complement()


print(rev_complement(test_record))

motif1 = "ATG"
motif2 = "TATA"
motif3 = "GAATTC"
motifs_to_find = [motif1, motif2, motif3]


def get_motif_positions(record, motifs):
    result = dict()
    for motif in motifs:
        result.update({motif: []})
    for (index, motif) in record.seq.search(motifs):
        result[motif].append(index)
    return result


print(get_motif_positions(test_record, motifs_to_find))


def get_protein_lengths_per_frame(record):
    rev = rev_complement(record)
    return {
        "Frame1": len(record.seq.translate(to_stop=True)),
        "Frame2": len(record.seq[1:].translate(to_stop=True)),
        "Frame3": len(record.seq[2:].translate(to_stop=True)),
        "Frame4": len(rev.translate(to_stop=True)),
        "Frame5": len(rev[1:].translate(to_stop=True)),
        "Frame6": len(rev[2:].translate(to_stop=True))
    }


print(get_protein_lengths_per_frame(test_record))

with open('./sequence_analysis.csv', 'w', encoding='utf-8', newline='') as output:
    writer = csv.writer(output)
    writer.writerow(["SeqID", "Nucleotide_Counts", "GC_Content", "Motif_Positions",
                     "Reverse_Complement", "Translation_Lengths"])
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        writer.writerow([seq_record.id,
                         nucleotide_counts(seq_record),
                         gc_content(seq_record),
                         get_motif_positions(seq_record, motifs_to_find),
                         rev_complement(seq_record),
                         get_protein_lengths_per_frame(seq_record)
                         ])

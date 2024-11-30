from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

seq = Seq("AAGAAATTCCAAGTCCAGGGATACACAAACAGGTGTACAGC"
          "AAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATATT")

countA = seq.count("A")
print("countA =", countA)

countT = seq.count("T")
print("countT =", countT)

countC = seq.count("C")
print("countC =", countC)

countG = seq.count("G")
print("countG =", countG)

GC_fraction = gc_fraction(seq)
print("countGC =", GC_fraction)

seq_rna = seq.transcribe()
print("seqRna =", seq_rna)

seq_protein = seq.translate()
print("seq_protein =", seq_protein)

seq_rev_complement = seq.reverse_complement()
print("seq_rev_complement =", seq_rev_complement)

with open('../../../../../Desktop/pythonProject/sequence_analysis.txt', 'w', encoding='utf-8') as output:
    output.write("Oryginalna sekwencja DNA: " + str(seq) + "\n")
    output.write("Liczba nukleotydów:\n")
    output.write(" A: " + str(countA) + "\n")
    output.write(" T: " + str(countT) + "\n")
    output.write(" C: " + str(countC) + "\n")
    output.write(" G: " + str(countG) + "\n")
    output.write("Zawartość GC: " + str("{:.2f}%".format(GC_fraction * 100)) + "\n")
    output.write("Transkrybowany RNA: " + str(seq_rna) + "\n")
    output.write("Translowane białko: " + str(seq_protein) + "\n")
    output.write("Odwrotne dopełnienie: " + str(seq_rev_complement) + "\n")

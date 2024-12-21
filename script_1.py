from Bio import Align

seq_1 = "AAAACCCCTTGGTTTT"
seq_2 = "AAACACCCTTTGTTTA"

aligner = Align.PairwiseAligner()
alignments = aligner.align(seq_1, seq_2)
score = aligner.score(seq_1, seq_2)

def interpret_alignment(coordinates):
    matches = []
    mismatches = []
    for i in range(len(coordinates[0]) - 1):
        if coordinates[0][i] == coordinates[0][i+1] or coordinates[1][i] == coordinates[1][i+1]:
            mismatches.append({"target": [coordinates[0][i].item(), coordinates[0][i+1].item()], "query": [coordinates[1][i].item(), coordinates[1][i+1].item()]})
        else:
            matches.append({"target": [coordinates[0][i].item(), coordinates[0][i+1].item()], "query": [coordinates[1][i].item(), coordinates[1][i+1].item()]})
    print("matches and mismatches in sequences:")
    print("matches", matches)
    print("mismatches", mismatches)


def count_mismatches(all_alignments):
    count = 0
    for alg in all_alignments:
        coordinates = alg.coordinates
        for i in range(len(coordinates[0]) - 1):
            if coordinates[0][i] == coordinates[0][i + 1] or coordinates[1][i] == coordinates[1][i + 1]:
                count += 1
                break
    return count


print("Alignment score between sequences", seq_1, seq_2, "is", score)
print("*******************************************************")
print("Displaying top 3 alignments")
iteration = 0
for alignment in alignments:
    iteration += 1
    print("-------------------------------------------------------")
    print("Alignment", iteration)
    print(alignment)
    print("coordinates")
    print(alignment.coordinates)
    interpret_alignment(alignment.coordinates)
    if iteration == 3:
        break

print("*******************************************************")
print("There are", count_mismatches(alignments), "alignments containing at least one mismatch")



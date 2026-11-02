# CIF421 Term Project - Progressive Multiple Sequence Alignment (ClustalW style)
# Name / Student no: Nisa ÜNSAL 220602016
# Language: Python 3
# How to run: python msa.py   (Input.txt and Blosum62.txt must be in the same folder)
#             the program asks for the gap penalty, e.g. -8 or -4

# read input sequences from a fasta formatted text file
def read_fasta(filename):
    names = []
    seqs = []
    with open(filename) as f:
        for line in f:
            # strip whitespaces, tabs, and newline characters
            line = line.strip()
            if line == "":
                continue
            # header line contains the sequence identifier
            if line.startswith(">"):
                names.append(line[1:])  # remove the leading '>' character
                seqs.append("")         # allocate a new empty string for incoming sequence
            else:
                # append multi-line sequence chunks to the current sequence
                seqs[-1] += line
    return names, seqs


# parse the substitution matrix (Blosum62) into a nested dictionary
def read_blosum(filename):
    with open(filename) as f:
        lines = f.read().strip().split("\n")
    # skip empty lines and comment lines (starting with '#')
    lines = [ln for ln in lines if ln.strip() != "" and not ln.strip().startswith("#")]
    # the first row gives column amino acid headers
    letters = lines[0].split()
    blosum = {}
    for line in lines[1:]:
        parts = line.split()
        row_letter = parts[0]  # row amino acid
        values = parts[1:]     # corresponding score values
        blosum[row_letter] = {}
        # map each column letter to its integer score
        for col, val in zip(letters, values):
            blosum[row_letter][col] = int(val)
    return blosum


# standard global alignment (Needleman-Wunsch) with specific tie-breaking rule
def needleman_wunsch(a, b, gap, score):
    n = len(a)
    m = len(b)

    # F stores dynamic programming scores, B stores traceback directions
    F = [[0] * (m + 1) for _ in range(n + 1)]
    B = [[""] * (m + 1) for _ in range(n + 1)]

    # initialize first column: aligning sequence 'a' against gaps
    for i in range(1, n + 1):
        F[i][0] = i * gap
        B[i][0] = "U"  # move up: letter in 'a' aligns with gap in 'b'

    # initialize first row: aligning sequence 'b' against gaps
    for j in range(1, m + 1):
        F[0][j] = j * gap
        B[0][j] = "L"  # move left: letter in 'b' aligns with gap in 'a'

    # fill the DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i-1][j-1] + score(a[i-1], b[j-1])  # match or mismatch
            up   = F[i-1][j] + gap                     # gap in second sequence
            left = F[i][j-1] + gap                     # gap in first sequence

            # assignment tie-breaking priority: Diagonal > Up > Left
            if diag >= up and diag >= left:
                F[i][j] = diag
                B[i][j] = "D"
            elif up >= left:
                F[i][j] = up
                B[i][j] = "U"
            else:
                F[i][j] = left
                B[i][j] = "L"

    return F, B


# backtrack from bottom-right (n, m) to top-left (0, 0) using direction pointers
def traceback(a, b, B):
    i = len(a)
    j = len(b)
    top = []
    bottom = []

    # trace moves until reaching the origin
    while i > 0 or j > 0:
        move = B[i][j]
        if move == "D":
            top.append(a[i-1])
            bottom.append(b[j-1])
            i -= 1
            j -= 1
        elif move == "U":
            top.append(a[i-1])
            bottom.append("-")
            i -= 1
        else:  # move == "L"
            top.append("-")
            bottom.append(b[j-1])
            j -= 1

    # reverse accumulated characters because we collected them backwards
    top.reverse()
    bottom.reverse()
    return "".join(top), "".join(bottom)


# calculate similarity (exact matches / alignment length) including gap columns
def similarity(top, bottom):
    length = len(top)
    matches = 0
    for idx in range(length):
        # count matching characters excluding gap columns
        if top[idx] == bottom[idx] and top[idx] != "-":
            matches += 1
    # similarity ratio = matching positions divided by total alignment length
    return matches, length, matches / length


# construct a guide tree and merge order using UPGMA clustering
def upgma(sim_matrix, names):
    k = len(names)
    # create a deep copy of the similarity matrix to avoid modifying the original
    S = [row[:] for row in sim_matrix]
    labels = names[:]
    sizes = [1] * k                 # track number of original sequences in each cluster
    members = [[i] for i in range(k)]  # track original sequence indices for progressive alignment
    merge_order = []

    iteration = 1
    while len(labels) > 1:
        best = -1.0
        bi, bj = -1, -1

        # search for the maximum off-diagonal similarity
        # strict '>' with i < j naturally breaks ties by picking earlier pairs
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if S[i][j] > best:
                    best = S[i][j]
                    bi, bj = i, j

        print(f"Iteration {iteration}: {labels[bi]} + {labels[bj]}  similarity = {best:.4f}")
        merge_order.append((members[bi][:], members[bj][:], best))

        # calculate weighted average similarity to all other clusters
        for c in range(len(labels)):
            if c != bi and c != bj:
                new_sim = (sizes[bi] * S[bi][c] + sizes[bj] * S[bj][c]) / (sizes[bi] + sizes[bj])
                S[bi][c] = new_sim
                S[c][bi] = new_sim
        S[bi][bi] = 1.0  # self-similarity remains 1.0

        # update cluster metadata at index bi
        labels[bi] = f"({labels[bi]},{labels[bj]})"
        sizes[bi] += sizes[bj]
        members[bi] += members[bj]

        # delete cluster bj (bi < bj ensures deleting bj doesn't shift bi's index)
        del labels[bj]
        del sizes[bj]
        del members[bj]

        # remove row and column bj from similarity matrix
        del S[bj]
        for row in S:
            del row[bj]

        # print updated similarity matrix for the lab report
        header = f"{'':25}" + "".join(f"{lbl:>25}" for lbl in labels)
        print(header)
        for r in range(len(labels)):
            row_str = f"{labels[r]:<25}" + "".join(f"{S[r][c]:25.4f}" for c in range(len(labels)))
            print(row_str)
        print()

        iteration += 1

    return labels[0], merge_order


# compute profile column-to-column score as average of all pairwise pairs
def col_score(colA, colB, gap, blosum):
    total = 0
    for x in colA:
        for y in colB:
            if x == "-" and y == "-":
                total += 0       # gap-gap has no penalty
            elif x == "-" or y == "-":
                total += gap     # letter-gap penalty
            else:
                total += blosum[x][y]  # amino acid substitution score
    return total / (len(colA) * len(colB))


# align two profiles (groups of pre-aligned sequences) using profile NW
def align_groups(A, B, gap, blosum):
    n = len(A[0])
    m = len(B[0])

    # extract columns from both alignment groups
    colsA = [[row[c] for row in A] for c in range(n)]
    colsB = [[row[c] for row in B] for c in range(m)]

    # dummy all-gap columns representing gap insertions
    gapA = ["-"] * len(A)
    gapB = ["-"] * len(B)

    # float score matrix and direction matrix
    F = [[0.0] * (m + 1) for _ in range(n + 1)]
    B_matrix = [[""] * (m + 1) for _ in range(n + 1)]

    # boundary initializations
    for i in range(1, n + 1):
        F[i][0] = F[i-1][0] + col_score(colsA[i-1], gapB, gap, blosum)
        B_matrix[i][0] = "U"
    for j in range(1, m + 1):
        F[0][j] = F[0][j-1] + col_score(gapA, colsB[j-1], gap, blosum)
        B_matrix[0][j] = "L"

    # fill profile DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i-1][j-1] + col_score(colsA[i-1], colsB[j-1], gap, blosum)
            up   = F[i-1][j] + col_score(colsA[i-1], gapB, gap, blosum)
            left = F[i][j-1] + col_score(gapA, colsB[j-1], gap, blosum)

            # round to 9 decimals to avoid floating-point precision tie issues
            d_r = round(diag, 9)
            u_r = round(up, 9)
            l_r = round(left, 9)

            # preserve priority: D > U > L
            if d_r >= u_r and d_r >= l_r:
                F[i][j] = diag
                B_matrix[i][j] = "D"
            elif u_r >= l_r:
                F[i][j] = up
                B_matrix[i][j] = "U"
            else:
                F[i][j] = left
                B_matrix[i][j] = "L"

    newA = [[] for _ in range(len(A))]
    newB = [[] for _ in range(len(B))]

    # profile traceback
    i = n
    j = m
    while i > 0 or j > 0:
        move = B_matrix[i][j]
        if move == "D":
            # align column from A with column from B
            for idx in range(len(A)):
                newA[idx].append(colsA[i-1][idx])
            for idx in range(len(B)):
                newB[idx].append(colsB[j-1][idx])
            i -= 1
            j -= 1
        elif move == "U":
            # align column from A with gap inserted into group B
            for idx in range(len(A)):
                newA[idx].append(colsA[i-1][idx])
            for idx in range(len(B)):
                newB[idx].append("-")
            i -= 1
        else:  # move == "L"
            # align column from B with gap inserted into group A
            for idx in range(len(A)):
                newA[idx].append("-")
            for idx in range(len(B)):
                newB[idx].append(colsB[j-1][idx])
            j -= 1

    # reverse each row and join into string
    resA = ["".join(reversed(row)) for row in newA]
    resB = ["".join(reversed(row)) for row in newB]
    return resA + resB


# --- Main Execution Flow ---

# 1. Load FASTA sequences and BLOSUM62 matrix
names, seqs = read_fasta("Input.txt")
blosum = read_blosum("Blosum62.txt")

# diagonal bonus requested by assignment (+2 to all self-scores)
for letter in blosum:
    blosum[letter][letter] += 2

# user input for gap penalty
gap = int(input("Gap penalty (negative integer): "))
score = lambda x, y: blosum[x][y]

# 2. Compute all pairwise alignments (10 pairs)
k = len(seqs)
sim_matrix = [[0.0] * k for _ in range(k)]
for i in range(k):
    sim_matrix[i][i] = 1.0  # diagonal is 1.0

alignments = {}

print("\nAlignment Results:")
print(f"{'Pair':<20} {'Score':<8} {'Matches':<10} {'Length':<10} {'Similarity':<10}")
print("-" * 60)

for i in range(k):
    for j in range(i + 1, k):
        F, B = needleman_wunsch(seqs[i], seqs[j], gap, score)
        top, bottom = traceback(seqs[i], seqs[j], B)
        matches, length, sim = similarity(top, bottom)

        # store alignment strings and similarity matrix values
        alignments[(i, j)] = (top, bottom)
        sim_matrix[i][j] = sim
        sim_matrix[j][i] = sim

        pair_name = f"{names[i]} - {names[j]}"
        print(f"{pair_name:<20} {F[-1][-1]:<8} {matches:<10} {length:<10} {sim:<10.4f}")

# 3. Print complete similarity matrix
print("\nSimilarity Matrix:")
header = f"{'':12}" + "".join(f"{name:>12}" for name in names)
print(header)
for i in range(k):
    row = f"{names[i]:<12}" + "".join(f"{sim_matrix[i][j]:12.4f}" for j in range(k))
    print(row)

# 4. Display the most similar pair alignment in blocks of 50
best_pair = max(alignments.keys(), key=lambda p: sim_matrix[p[0]][p[1]])
top, bottom = alignments[best_pair]

print(f"\nMost Similar Pair: {names[best_pair[0]]} - {names[best_pair[1]]} (Similarity: {sim_matrix[best_pair[0]][best_pair[1]]:.4f})\n")
for p in range(0, len(top), 50):
    print(top[p:p+50])
    print(bottom[p:p+50])
    print()

# 5. Build guide tree using UPGMA
tree, merge_order = upgma(sim_matrix, names)
print("\nGuide tree:", tree)

# 6. Progressive Multiple Sequence Alignment along guide tree order
groups = {}
for i in range(k):
    groups[(i,)] = [seqs[i]]

print("\nProgressive Alignment Steps:")
for left, right, sim in merge_order:
    A = groups[tuple(left)]
    B = groups[tuple(right)]
    merged = align_groups(A, B, gap, blosum)
    new_key = tuple(left + right)
    groups[new_key] = merged

    left_str = "+".join(names[x] for x in left)
    right_str = "+".join(names[x] for x in right)
    print(f"\nMerging ({left_str}) + ({right_str}):")
    for idx, member_idx in enumerate(new_key):
        print(f"{names[member_idx]:<12} {merged[idx]}")

# 7. Reorder final MSA back to the original input sequence order
final_key = tuple(merge_order[-1][0] + merge_order[-1][1])
final_alignment = groups[final_key]

ordered_msa = [None] * k
for idx, member_idx in enumerate(final_key):
    ordered_msa[member_idx] = final_alignment[idx]

# 8. Sanity and consistency checks
msa_len = len(ordered_msa[0])
all_same_length = len(set(len(s) for s in ordered_msa)) == 1
all_match_original = all(ordered_msa[i].replace("-", "") == seqs[i] for i in range(k))

# count conserved positions across all 5 sequences
total_stars = 0
for col_idx in range(msa_len):
    chars = [ordered_msa[r][col_idx] for r in range(k)]
    if len(set(chars)) == 1 and chars[0] != "-":
        total_stars += 1

# 9. Print final MSA in 50-character blocks with consensus row
print("\nFinal Multiple Sequence Alignment:")
for p in range(0, msa_len, 50):
    for i in range(k):
        print(f"{names[i]:<12} {ordered_msa[i][p:p+50]}")
    # place '*' only when all sequences have identical non-gap residue
    block_stars = "".join("*" if len(set(ordered_msa[r][c] for r in range(k))) == 1 and ordered_msa[0][c] != "-" else " " for c in range(p, min(p + 50, msa_len)))
    print(f"{'':<12} {block_stars}\n")

# summary verification report
print(f"Final alignment length: {msa_len}")
print(f"Total fully conserved columns (*): {total_stars}")
print(f"All sequences have identical length: {all_same_length}")
print(f"All sequences match original without gaps: {all_match_original}")

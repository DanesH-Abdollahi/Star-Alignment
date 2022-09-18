import numpy as np


def global_align(x, y, s_match, s_mismatch, s_gap):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))

    for i in range(len(y) + 1):
        A[i][0] = s_gap * i

    for i in range(len(x) + 1):
        A[0][i] = s_gap * i

    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1]
                + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != "-") else 0)
                + (
                    s_mismatch
                    if (y[i - 1] != x[j - 1] and y[i - 1] != "-" and x[j - 1] != "-")
                    else 0
                )
                + (s_gap if (y[i - 1] == "-" or x[j - 1] == "-") else 0),
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        current_score = A[j][i]
        if (
            i > 0
            and j > 0
            and (
                (
                    (x[i - 1] == y[j - 1] and y[j - 1] != "-")
                    and current_score == A[j - 1][i - 1] + s_match
                )
                or (
                    (y[j - 1] != x[i - 1] and y[j - 1] != "-" and x[i - 1] != "-")
                    and current_score == A[j - 1][i - 1] + s_mismatch
                )
                or (
                    (y[j - 1] == "-" or x[i - 1] == "-")
                    and current_score == A[j - 1][i - 1] + s_gap
                )
            )
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1

        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1

    return (align_X, align_Y, A[len(y)][len(x)])


def MSA_Score(N, seqs):
    Score = 0
    for l in range(len(seqs[0])):
        for i in range(N):
            for j in range(i + 1, N):
                if seqs[i][l] == seqs[j][l]:
                    if seqs[i][l] == "-":
                        Score += 0

                    else:
                        Score += 3

                else:
                    if seqs[i][l] != "-" and seqs[j][l] != "-":
                        Score += -1

                    else:
                        Score += -2

    return int(Score)


def Star_Alignment(N, Sequences, s_match, s_mismatch, s_gap):
    Aligned_Sequences = Sequences
    Similarity_Matrix = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i != j:
                Similarity_Matrix[i][j] = global_align(
                    Sequences[i], Sequences[j], s_match, s_mismatch, s_gap
                )[2]

    Max = sum(Similarity_Matrix[0, :])
    Max_i = 0
    for i in range(1, N):
        if sum(Similarity_Matrix[i, :]) > Max:
            Max_i = i

    Star_Sequence = Sequences[Max_i]
    array = Similarity_Matrix[Max_i, :]
    array = np.argsort(-1 * array)
    Alignment_Order = []
    for i in array:
        if i != Max_i:
            Alignment_Order.append(i)
    Seqs_in_Orde_To_Align = []

    for ii in Alignment_Order:
        Seqs_in_Orde_To_Align.append(Sequences[ii])

    counter = 0

    for seq in Seqs_in_Orde_To_Align[0:]:

        Answ = global_align(seq, Star_Sequence, s_match, s_mismatch, s_gap)

        Temp = Answ[1]
        Aligned_Sequences[Alignment_Order[counter]] = Answ[0]

        char_loc = []

        iii = 0
        while iii < len(Star_Sequence):
            if Temp[iii] != Star_Sequence[iii]:
                char_loc.append(iii)
                Temp = Temp[0:iii] + Temp[iii + 1 :]
            else:
                iii += 1

        last_gaps = len(Temp) - len(Star_Sequence)

        char_loc = np.flip(char_loc)

        if counter >= 1:

            for k in range(counter):
                for loc in char_loc:

                    Aligned_Sequences[Alignment_Order[k]] = (
                        Aligned_Sequences[Alignment_Order[k]][0:loc]
                        + "-"
                        + Aligned_Sequences[Alignment_Order[k]][loc:]
                    )

                Aligned_Sequences[Alignment_Order[k]] = (
                    Aligned_Sequences[Alignment_Order[k]] + last_gaps * "-"
                )

        Star_Sequence = Answ[1]
        counter += 1

        Aligned_Sequences[Max_i] = Star_Sequence

    return Aligned_Sequences


def Block(N, Sequences):
    s_match = 3
    s_mismatch = -1
    s_gap = -2

    Where_is_equal = []

    Where_is_equal.append(-1)

    for l in range(len(Sequences[0])):
        counter = 0
        temp_char = Sequences[0][l]
        for i in range(1, N):
            if Sequences[i][l] == temp_char:
                counter += 1

        if counter == N - 1:
            Where_is_equal.append(l)

    # Find Block's
    Where_is_equal.append(len(Sequences[0]))

    Blocked_Seqs = []
    Locations = []
    for i in range(1, len(Where_is_equal)):
        end = Where_is_equal[i]
        first = Where_is_equal[i - 1]
        if (end - first) >= 3:
            Locations.append([first, end])
            Blocked_Seqs.append([])
            for j in range(N):
                Blocked_Seqs[len(Blocked_Seqs) - 1].append(
                    Sequences[j][first + 1 : end]
                )

    # Remove All gaps from Block's Sequences
    for j in range(len(Blocked_Seqs)):
        for i in range(N):
            Blocked_Seqs[j][i] = Blocked_Seqs[j][i].replace("-", "")

    Scores = []
    Scores_Seq = []

    for i in range(len(Blocked_Seqs)):
        Temp_Seq = Sequences.copy()
        New_Block = Star_Alignment(N, Blocked_Seqs[i], s_match, s_mismatch, s_gap)

        for j in range(N):
            Temp_Seq[j] = (
                Temp_Seq[j][0 : Locations[i][0] + 1]
                + New_Block[j]
                + Temp_Seq[j][Locations[i][1] :]
            )

        Scores.append(MSA_Score(N, Temp_Seq))
        Scores_Seq.append(Temp_Seq)

    max_value = max(Scores)
    index = Scores.index(max_value)

    Sequences = Scores_Seq[index]

    return Sequences


s_match = 3
s_mismatch = -1
s_gap = -2


N = int(input())
Sequences = []
for i in range(N):
    Sequences.append(input())

Aligned_Sequences = Star_Alignment(N, Sequences, s_match, s_mismatch, s_gap)

for i in range(50):
    Aligned_Sequences = Block(N, Aligned_Sequences)


print(MSA_Score(N, Aligned_Sequences))

for seq in Aligned_Sequences:
    print(seq)

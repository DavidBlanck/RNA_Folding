test_seq = 'gaggaaaguccgggcUAGCACACACCUUAUGGGUGUGUAGUGUUUGUGCUAAGGGAAAUCAUAACCUUAGGUAUGUUGUAUAAACAUAACGGCAAACUAGUUAUAGCUAAGGUGUUUCACUACGUUAUAACUUAAAUUAAAGUGCCACAGAGACGAAUCUAUUUAGAAAUAAAUAGAGUGAAACGCGGUAAACCCCUCAAGCUAGCAACCCAAAUUAGGUAGGGGCACAUGAUGUGUAGCAAUACAACAUCAUGCAAGAUUUGAAUCUUGAGAUUAAUAGUCACAAAAGAAGAAAUUCUUUacagaacgcggcuua'


"""
This function evaluates an RNA sequence string of nucleotide encodings (A, U, G, or C) and returns a dynamic programming
table and a trace back table storing the results of the evaluation. The last element in the first nested array of the 
score table stores the optimal number of pairings. The trace back table can be used to find the pairings.

Parameters - The RNA sequence to be evaluated

Returns - Score table and trace table in that order. Both are nested arrays. The score table can be used to find the
optimal score of the evaluation. The trace table can be used to determine the pairs in the optimal alignment.
"""


def evaluate_seq(seq):
    seq = seq.upper();
    score_table = []
    for s in range(len(seq)):
        score_table.append([0] * len(seq))

    trace_table = []
    for s in range(len(seq)):
        trace_table.append([0] * len(seq))

    for k in range(5, len(seq)):
        for i in range(0, len(seq) - k):
            j = i + k
            opt = 0
            for t in range(i, j - 4):
                if (seq[j] == 'A' and seq[t] == 'U'
                        or seq[j] == 'U' and seq[t] == 'A'
                        or seq[j] == 'C' and seq[t] == 'G'
                        or seq[j] == 'G' and seq[t] == 'C'):
                    if i == t:
                        temp = 1 + score_table[t + 1][k - 2]
                    else:
                        temp = 1 + score_table[i][t - i - 1] + score_table[t + 1][j - t - 2]
                    if temp > opt:
                        opt = temp
                        optt = t
            if score_table[i][k - 1] >= opt:
                score_table[i][k] = score_table[i][k - 1]
                trace_table[i][k] = (i, k - 1)
            elif optt == i:
                score_table[i][k] = opt
                trace_table[i][k] = (i + 1, k - 2, 'p')
            else:
                score_table[i][k] = opt
                trace_table[i][k] = (i, optt - 1 - i, optt + 1, j - optt - 2)

    return score_table, trace_table


"""
Evaluates a trace table to determine the pairings in an optimal RNA alignment
Parameters -
    trace_table - a nested array storing pointers to the previous node or nodes that each node was derived from.
Returns - 
    A list of the positions for each pairing in the optimal alignment
"""


def traceback(trace_table):
    nodes = [trace_table[0][-1]]
    pairs = []
    while len(nodes) > 0:
        temp = nodes.pop()
        if temp == 0:
            pass
        elif len(temp) == 2:
            nodes.append(trace_table[temp[0]][temp[1]])
        elif len(temp) == 3:
            pairs.append((temp[0] - 1, temp[0] + temp[1] + 1))
            nodes.append(trace_table[temp[0]][temp[1]])
        else:
            pairs.append((temp[2] - 1, temp[3] + temp[2] + 1))
            nodes.append(trace_table[temp[0]][temp[1]])
            nodes.append(trace_table[temp[2]][temp[3]])
    return pairs


st, tt = evaluate_seq(test_seq)

pairs = traceback(tt)

dot_paren = "." * len(test_seq)
for pair in pairs:
    dot_paren = dot_paren[:pair[0]] + '(' + dot_paren[pair[0] + 1:]
    dot_paren = dot_paren[:pair[1]] + ')' + dot_paren[pair[1] + 1:]

print(dot_paren)

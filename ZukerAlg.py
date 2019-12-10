"""
This function evaluates an RNA sequence string of nucleotide encodings (A, U, G, or C) and returns a dynamic programming
table and a trace back table storing the results of the evaluation. The last element in the first nested array of the 
score table stores the optimal number of pairings. The trace back table can be used to find the pairings.

Parameters - The RNA sequence to be evaluated

Returns - Score table and trace table in that order. Both are nested arrays. The score table can be used to find the
optimal score of the evaluation. The trace table can be used to determine the pairs in the optimal alignment.
"""


def evaluate_seq_weighted(seq):
    seq = seq.upper();
    score_table = []
    for s in range(len(seq)):
        score_table.append([0] * len(seq))

    trace_table = []
    for s in range(len(seq)):
        trace_table.append([0] * len(seq))

    # Iterate through range of distances between the two indices
    for k in range(5, len(seq)):
        # Iterate through range of starting positions for the first index value
        for i in range(0, len(seq) - k):
            j = i + k
            opt = 0
            # Now we consider every case. Either j forms a bond with any of the nucleotides
            # between i and j or j does not bond.
            for t in range(i, j - 4):
                if (seq[j] == 'A' and seq[t] == 'U'
                        or seq[j] == 'U' and seq[t] == 'A'
                        or seq[j] == 'C' and seq[t] == 'G'
                        or seq[j] == 'G' and seq[t] == 'C'):
                    if seq[j] == 'C' or seq[j] == 'G':
                        bonus = 5.53
                    else:
                        bonus = 4.42
                    if i == t:
                        temp = bonus + score_table[t + 1][k - 2]
                    else:
                        temp = bonus + score_table[i][t - i - 1] + score_table[t + 1][j - t - 2]
                    if temp > opt:
                        opt = temp
                        optt = t
            # If the score at i, k - 1 is greater than the optimal binding score, then j does not bond. Otherwise,
            # j is bonded to an intermediate position.
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



def dotparen(pairs, length):
    dot_paren = "." * length
    for pair in pairs:
        dot_paren = dot_paren[:pair[0]] + '(' + dot_paren[pair[0] + 1:]
        dot_paren = dot_paren[:pair[1]] + ')' + dot_paren[pair[1] + 1:]
    return dot_paren



def topairs(dotparen):
    opening = []
    pairs = []
    for i in range(len(dotparen)):
        if dotparen[i] == '(':
            opening.append(i)
        elif dotparen[i] == ')' and len(pairs) > 0:
            pairs.append((opening.pop(), i))
        elif dotparen[i] == '>' or dotparen[i] == '}':
            pairs.append('x')
    return pairs


def findcommon(l1, l2):
    s2 = set()
    common = 0
    for e in l2:
        s2.add(e)
    for e in l1:
        if (e in s2):
            common += 1
    return common



test_seq = 'UUUUAUGGAGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAGCGAACGGACGAGAAGCUUGCUUCUCUGAUGUUAGCGGCGGACGGGUGAGUAACACGUGGAUAACCUACCUAUAAGACUGGGAUAACUUCGGGAAACCGGAGCUAAUACCGGAUAAUAUUUUGAACCGCAUGGUUCAAAAGUGAAAGACGGUCUUGCUGUCACUUAUAGAUGGAUCCGCGCUGCAUUAGCUAGUUGGUAAGGUAACGGCUUACCAAGGCAACGAUGCAUAGCCGACCUGAGAGGGUGAUCGGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUAGGGAAUCUUCCGCAAUGGGCGAAAGCCUGACGGAGCAACGCCGCGUGAGUGAUGAAGGUCUUCGGAUCGUAAAACUCUGUUAUUAGGGAAGAACAUAUGUGUAAGUAACUGUGCACAUCUUGACGGUACCUAAUCAGAAAGCCACGGCUAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGUGGCAAGCGUUAUCCGGAAUUAUUGGGCGUAAAGCGCGCGUAGGCGGUUUUUUAAGUCUGAUGUGAAAGCCCACGGCUCAACCGUGGAGGGUCAUUGGAAACUGGAAAACUUGAGUGCAGAAGAGGAAAGUGGAAUUCCAUGUGUAGCGGUGAAAUGCGCAGAGAUAUGGAGGAACACCAGUGGCGAAGGCGACUUUCUGGUCUGUAACUGACGCUGAUGUGCGAAAGCGUGGGGAUCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGAGUGCUAAGUGUUAGGGGGUUUCCGCCCCUUAGUGCUGCAGCUAACGCAUUAAGCACUCCGCCUGGGGAGUACGACCGCAAGGUUGAAACUCAAAGGAAUUGACGGGGACCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAAGCAACGCGAAGAACCUUACCAAAUCUUGACAUCCUUUGACAACUCUAGAGAUAGAGCCUUCCCCUUCGGGGGACAAAGUGACAGGUGGUGCAUGGUUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAAGCUUAGUUGCCAUCAUUAAGUUGGGCACUCUAAGUUGACUGCCGGUGACAAACCGGAGGAAGGUGGGGAUGACAUCAAAUCAUCAUGCCCCUUAUGAUUUGGGCUACACACGUGCUACAAUGGACAAUACAAAGGGCAGCGAAACCGCGAGGUCAAGCAAAUCCCAUAAAGUUGUUCUCAGUUCGGAUUGUAGUCUGCAACUCGACUACAUGAAGCUGGAAUCGCUAGUAAUCGUAGAUCAGCAUGCUACGGUGAAUACGUUCCCGGGUCUUGUACACACCGCCCGUCACACCACGAGAGUUUGUAACACCCGAAGCCGGUGGAGUAACCUUUUAGGAGCCAGCCGUCGAAGGUGGGACAAAUGAUUGGGGUGAAGUCGUAACAAGGUAGCCGUAUCGGAAGGUGCGGCUGGAUCACCUCCUUUCU'
biostruct = '.........(((((...<<<.))))).((((.(((((.(((((((((....(((.(((..(((..((((((....((((((....))))))....)))))))))......(((......((((((((..((...(((((((.((((....(((((((....))))))).....))))......(((((((((....)))))))))......((((((...)))))).)))))))..))))))))))(((..(.(((..((((((((.......))))))))))).....))))..((((((((....))))...))))))).((((((..........)))))).((((....))))...)))))).).....(.(((...(((((....))))).)))).)).))))))..((((......((((....)))).....)))).(..((((((...(.....((((((((......)))))))).....)....))))))..)...((((([[[...(((((.....((.]]])).......)))))))))).))))))))))..........((({{.....((((.(.(((.(((((((.(((((((((((....(((((((.....)))))))..)))))))))..)))))))))...(((((((((..(((((((((..((((((((...(((......)))......))))))))..))....(..((....)))))))))).)))))).)))...))))..))))....((((((...((...((((.........))))...))))))))..........((((((..(((((((((((((....)))))))))))))...((..}})).....)))))).))).(((......((((....))))....)))...>>>..(((((.(((((((.((..(((((((((((((((((....((((........))))........(((((((.....((((((((..(((((((....)))))))..((((....)))).))))).))).((.((((..(((((((((...(((((((((....)))..((((......))))..)))))).....((((.(((((((...((..((......))))....)))))))..((((((((.....)))))))).....))))....)))).)))...))))))))....)))))))...)).))))))))))...(((((((.....(((..((...(((....)))...))....))).....)))))))......(...((((((((........))))))))...).....))))).....((((((((.......))))))))......))...)))))))))).))....((.((...((((((((..((((((((((((....((((((..(((..(((...))).))).))))))...))))))))))))..))))))))....))..))....((((((((((....))))))))))...............'

st1, tt1 = evaluate_seq_weighted(test_seq)

ps1 = traceback(tt1)

#print(ps1)

print(dotparen(ps1, len(test_seq)))

ps2 = topairs(biostruct)

#print(ps2)

#print(biostruct)

print(findcommon(ps1, ps2))
#print(len(ps1))
print(len(ps2))

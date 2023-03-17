import numpy as np

rna_code = {'UUU': 'F', 'UUC': 'F',
            'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
            'AUG': 'M',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'UAU': 'Y', 'UAC': 'Y',
            'CAU': 'H', 'CAC': 'H',
            'CAA': 'Q', 'CAG': 'Q',
            'AAU': 'N', 'AAC': 'N',
            'AAA': 'K', 'AAG': 'K',
            'GAU': 'D', 'GAC': 'D',
            'GAA': 'E', 'GAG': 'E',
            'UGU': 'C', 'UGC': 'C',
            'UGG': 'W',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
            'AGU': 'S', 'AGC': 'S',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            'UAA': '', 'UAG': '', 'UGA': ''}


def Counting_DNA_Nucleotides(x: str):
    A = x.count('A')
    C = x.count('C')
    G = x.count('G')
    T = x.count('T')
    return [A, C, G, T]


def Transcribing_DNA_into_RNA(x: str):
    return x.replace('T', 'U')


def Complementing_a_Strand_of_DNA(x: str):
    x = x[::-1]
    x = x.replace('A', 'Q')
    x = x.replace('G', 'W')
    x = x.replace('C', 'E')
    x = x.replace('T', 'R')
    x = x.replace('Q', 'T')
    x = x.replace('W', 'C')
    x = x.replace('E', 'G')
    x = x.replace('R', 'A')
    return x


def Rabbits_and_Recurrence_Relations(n: int, k: int):
    rabbit = [0] * 50
    rabbit[0] = 1
    rabbit[1] = 1
    if n > 2:
        for i in range(2, n):
            rabbit[i] = rabbit[i - 1] + rabbit[i - 2] * k
    return rabbit[n - 1]


def GC_Content(x: str):
    x = x.replace(' ', '')
    cg = x.count('C') + x.count('G')
    return (cg / len(x)) * 100


def Counting_Point_Mutations(x: str, y: str):
    dh = 0
    x = x.replace(' ', '')
    y = y.replace(' ', '')
    for i, j in zip(x, y):
        if i != j:
            dh += 1
    return dh


def Translating_RNA_into_Protein(x: str):
    lenx = len(x)
    res = ''
    for i in range(1, lenx + 1, 3):
        s = x[i - 1] + x[i] + x[i + 1]
        res += rna_code[s]
    return res


def Finding_a_Motif_in_DNA(s: str, t: str):
    lent, lens = len(t), len(s)
    res = []
    for i in range(1, lens + 1 - lent):
        mid = s[i - 1:i + lent - 1]
        if mid == t:
            res.append(i)
    return res


def Consensus_and_Profile(s):
    res = {'A': [], 'C': [], 'G': [], 'T': []}
    s = np.matrix(s)
    trans = np.transpose(s)
    s = trans.tolist()
    max_s = ''
    acgt = ['A', 'C', 'G', 'T']
    for mid in s:
        if mid.count('A') != 0:
            res['A'].append(mid.count('A'))
        else:
            res['A'].append(0)
        if mid.count('C') != 0:
            res['C'].append(mid.count('C'))
        else:
            res['C'].append(0)
        if mid.count('G') != 0:
            res['G'].append(mid.count('G'))
        else:
            res['G'].append(0)
        if mid.count('T') != 0:
            res['T'].append(mid.count('T'))
        else:
            res['T'].append(0)
    for mid_list in list(zip(res['A'], res['C'], res['G'], res['T'])):
        max_s += acgt[mid_list.index(max(mid_list))]
    return res, max_s


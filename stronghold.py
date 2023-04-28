import numpy as np
import re
from itertools import combinations, permutations

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
            'UAA': ' ', 'UAG': ' ', 'UGA': ' '}

pro_code = {
    'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Y': ['UAU', 'UAC'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['UGU', 'UGC'],
    'W': ['UGG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    ' ': ['UAA', 'UAG', 'UGA']}


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


def Mortal_Fibonacci_Rabbits(n: int, m: int):
    rabbits = [0] * 150
    rabbits[0] = 1
    rabbits[1] = 1
    rabbits[2] = 1
    for i in range(2, m):
        rabbits[i] = rabbits[i - 1] + rabbits[i - 2]
    for i in range(m, n + 1):
        mid = 0
        for j in range(i - m, i - 1):
            mid += rabbits[j]
        rabbits[i] = mid
    print(rabbits)
    return rabbits[n - 1]


def RNA_Splicing(s: str):
    s = s.replace('T', 'U')
    ls = len(s)
    pro = ''
    for i in range(0, ls, 3):
        mid = ''
        mid += s[i] + s[i + 1] + s[i + 2]
        pro += rna_code[mid]
    return pro


def sub_str(input_string):
    length = len(input_string)
    return [input_string[i:j + 1] for i in range(length) for j in range(i, length)]


def common_str(str1, str2):
    ls = sub_str(str1)
    print(ls)
    target = []
    for i in ls:
        if str2.find(i) != -1:
            target.append(i)
    return target


def Finding_a_Shared_Motif(s):
    mid = common_str(s[1], s[0])
    print(mid)
    mid = set(mid)
    mid = list(mid)
    print(mid)
    for dna in s[2:]:
        new_mid = []
        for i in mid:
            if dna.find(i) != -1:
                new_mid.append(i)
        mid = new_mid
    res = ''
    for i in mid:
        if len(i) > len(res):
            res = i
    print(mid)
    return res


def read_txt(path):
    with open(path, "r") as f:
        data = f.readlines()
    dna_list = []
    mid = ''
    for i in data:
        if i[0] == '>':
            dna_list.append(mid)
            mid = ''
        else:
            mid += i[:-1]
    return dna_list[1:]


def RNA_proList(rna: str):
    ls = len(rna)
    pro_list = []
    pro = ''
    for i in range(0, ls, 3):
        mid = ''
        mid += rna[i] + rna[i + 1] + rna[i + 2]
        if rna_code[mid] == ' ':
            pro_list.append(pro)
            pro = ''
            continue
        pro += rna_code[mid]
    pro = ''
    for i in range(1, ls - 2, 3):
        mid = ''
        mid += rna[i] + rna[i + 1] + rna[i + 2]
        if rna_code[mid] == ' ':
            pro_list.append(pro)
            pro = ''
            continue
        pro += rna_code[mid]
    pro = ''
    for i in range(2, ls - 1, 3):
        mid = ''
        mid += rna[i] + rna[i + 1] + rna[i + 2]
        if rna_code[mid] == ' ':
            pro_list.append(pro)
            pro = ''
            continue
        pro += rna_code[mid]
    return pro_list


def Open_Reading_Frames(dna: str):
    dna_s = Complementing_a_Strand_of_DNA(dna)
    rna = dna.replace('T', 'U')
    rna_s = dna_s.replace('T', 'U')
    prolist = RNA_proList(rna)
    prolist += RNA_proList(rna_s)
    res = []
    for pro in prolist:
        i = pro.find('M')
        if i != -1:
            res.append(pro[i:])
        while True:
            i = pro.find('M', i + 1)
            if i == -1:
                break
            res.append(pro[i:])
    res = set(res)
    res = list(res)
    return res


def cut(s: str):
    results = []
    for x in range(len(s)):
        for i in range(len(s) - x):
            results.append(s[i:i + x + 1])
    return results


def Locating_Restriction_Sites(dna: str):
    inp = cut(dna)
    inp = set(inp)
    inp = list(inp)
    res = []
    for s in inp:
        if len(s) > 3 and len(s) < 13:
            rs = Complementing_a_Strand_of_DNA(s)
            if rs == s:
                pos = -1
                while True:
                    pos = dna.find(s, pos + 1)
                    if pos != -1:
                        res.append([pos + 1, len(s)])
                    else:
                        break
    return res


def Inferring_mRNA_from_Protein(pro: str):
    pro += ' '
    res = 1
    for mrna in pro:
        res *= len(pro_code[mrna])
    return res


def Finding_a_Spliced_Motif(dna: str, s: str):
    len_dna = len(dna)
    len_s = len(s)
    j = 0
    res = []
    for i in range(len_dna):
        if j == len_s:
            break
        if dna[i] == s[j]:
            j += 1
            res.append(i + 1)
    return res


def Enumerating_Gene_Orders(n: int):
    nums = []
    for i in range(n):
        nums.append(i + 1)
    result = []
    for mid in permutations(nums, len(nums)):
        result.append(list(mid))
    return result


def read_txt_2(path):
    a = ""
    b = ""
    with open(path, "r") as f:
        data = f.readlines()
    rosa_list = {}
    a = data[0][1:-1]
    for i in data[1:]:
        if i[0] != '>':
            b += i[:-1]
        else:
            rosa_list[a] = b
            b = ""
            a = i[1:-1]
    rosa_list[a] = b
    return rosa_list


def Overlap_Graphs(rosa_list):
    res = []
    for i in rosa_list:
        mid = rosa_list[i]
        for j in rosa_list:
            if rosa_list[j] == mid:
                continue
            if mid[-3:] == rosa_list[j][0:3]:
                res.append((i, j))
    # lth = len(res)
    # iex = []
    # result = []
    # for i in range(lth):
    #     for j in range(i+1, lth):
    #         if res[i][0] == res[j][1] and res[i][1] == res[j][0]:
    #             iex.append(j)
    # print(iex)
    # for i in range(lth):
    #     if i not in iex:
    #         result.append(res[i])
    return res
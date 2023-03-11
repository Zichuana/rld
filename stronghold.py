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
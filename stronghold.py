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

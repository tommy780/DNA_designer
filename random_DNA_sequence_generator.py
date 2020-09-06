"""
Random DNA sequence generator
Author: Chaser
Version: 0.1
Description:
    Generate random DNA sequence without self-complementary domain for idea single DNA strand.
"""




from random import randint
import numpy as np
import itertools
from Thermodynamic import calc_dG, melting_temp


# def randomDNA(num, base):
#     """Random DNA sequence generator from base species"""
#
#     dna = ''
#     j = num // 5
#     i = 0
#
#     for i in range(1, num + 1):
#         k = int(randint(0, len(base) - 1))
#         dna += base[k]
#
#     return dna


# generate DNA sequence sequentially. 0, 1, 2, 3 refer to base[0], base[1], base[2], base[3]
def sequential_seq(num, base, flag):
    dna = ''
    code = len(base)
    for i in range(num-1, -1, -1):
        bit = flag // code**i
        flag %= code**i
        dna += base[bit]
    return dna



# GC_ratio_restrict at ratio
def GC_ratio_restrict(seq, ratio=0.5):
    num = len(seq)
    GC_limit = round(ratio * num)
    seq_GC = seq_base_format(seq, 'S')
    GC_num = 0
    for i in range(0, num):
        if seq_GC[i] == 'S':
            GC_num += 1
    if GC_num != GC_limit:
        return 0
    else:
        return 1


# check the base number in seq among (num//3 to ((num//3)+1))  eg. num = 5  -> [1,2]
'''
def seq_num_check(seq, base):
    num = len(seq)
    if num < 4:
        return 1

    c_num = 0
    for i in range(0, num - 1):
        if seq[i] == base:
            c_num += 1
    if (c_num < (num // 3)) or (c_num > ((num // 3) + 1)):
        return 0
    else:
        return 1
'''


# if the seq contains N consecutive same base then return 0, and vice verse.
def secondary_structure(seq, n, base):
    if len(seq) < 6:
        return 1

    seq_count = [(k, len(list(g))) for k, g in itertools.groupby(seq)]
    for i in range(0, len(seq_count)):
        if seq_count[i][0] == base:
            if seq_count[i][1] >= n:
                return 0
    return 1


# format seq with IUPAC nucleotide code
"""
*******************IUPAC degenerate nucleotide codes for DNA**********************
Code Nucleotides
M       A or C
R       A or G
W       A or T
S       C or G
Y       C or T
K       G or T
V       A, C, or G
H       A, C, or T
D       A, G, or T
B       C, G, or T
N       A, C, G, or T
****************************************
"""


def seq_base_format(seq, code):
    dna_base_f = ''
    if code == 'M':  # A or C
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'C':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'R':  # A or T
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'W':  # A or T
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'S':  # C or G
        for i in range(0, len(seq)):
            if seq[i] == 'C' or seq[i] == 'G':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'Y':  # C or T
        for i in range(0, len(seq)):
            if seq[i] == 'C' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'K':  # G or T
        for i in range(0, len(seq)):
            if seq[i] == 'G' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'V':  # A or C or G
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'C' or seq[i] == 'G':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'H':  # A or C or T
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'C' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'D':  # A or G or T
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'G' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'B':  # C or G or T
        for i in range(0, len(seq)):
            if seq[i] == 'C' or seq[i] == 'G' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    if code == 'N':  # A or C or G or T
        for i in range(0, len(seq)):
            if seq[i] == 'A' or seq[i] == 'C' or seq[i] == 'G' or seq[i] == 'T':
                dna_base_f += code
            else:
                dna_base_f += seq[i]
    return dna_base_f


# seq prevent AAA GGG CCC TTT KKKKKK MMMMMM RRRRRR SSSSSS WWWWWW YYYYYY ,
# if yes, return 0, and vice verse
def seq_prevent(seq):
    """AAA"""
    if secondary_structure(seq, 3, 'A') == 0:
        return 0
    """GGG"""
    if secondary_structure(seq, 3, 'G') == 0:
        return 0
    """CCC"""
    if secondary_structure(seq, 3, 'C') == 0:
        return 0
    """TTT"""
    if secondary_structure(seq, 3, 'T') == 0:
        return 0

    """KKKKKK"""
    seq_f = seq_base_format(seq, 'K')
    if secondary_structure(seq_f, 6, 'K') == 0:
        return 0
    """MMMMMM"""
    seq_f = seq_base_format(seq, 'M')
    if secondary_structure(seq_f, 6, 'M') == 0:
        return 0
    """RRRRRR"""
    seq_f = seq_base_format(seq, 'R')
    if secondary_structure(seq_f, 6, 'R') == 0:
        return 0
    """SSSSSS"""
    seq_f = seq_base_format(seq, 'S')
    if secondary_structure(seq_f, 6, 'S') == 0:
        return 0
    """WWWWWW"""
    seq_f = seq_base_format(seq, 'W')
    if secondary_structure(seq_f, 6, 'W') == 0:
        return 0
    """YYYYYY"""
    seq_f = seq_base_format(seq, 'Y')
    if secondary_structure(seq_f, 6, 'Y') == 0:
        return 0
    return 1


# DNA seq divided by ' ' each 5 bases
def seq_format(seq):
    dna_a = ''
    for i in range(1, len(seq) + 1):
        dna_a += seq[i - 1]
        if i % 5 == 0 and i != len(seq):
            dna_a += ' '
    return dna_a


# delete blank in seq
def seq_format_no_blank(seq):
    dna_a = ''
    for i in range(len(seq)):
        if (seq[i] != 'A') & (seq[i] != 'T') & (seq[i] != 'C') & (seq[i] != 'G'):
            continue
        dna_a += seq[i]
    return dna_a


# calculate the hamming distance of str1 and str2
def hamming_distance(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(a != b for a, b in zip(str1, str2))


# overlapping limit
def overlapping_limit(seq, sub_seq, lapping_num):
    if lapping_num > len(seq) | lapping_num > len(sub_seq):
        raise ValueError("Undefined for sequences of overflow lapping_num")
    else:
        # choose subunit from sub_seq
        for i in range(0, len(sub_seq) - lapping_num + 1):
            temp_sub_seq = ''
            for j in range(0, lapping_num):
                temp_sub_seq += sub_seq[i + j]
            # choose subunit from seq
            for k in range(0, len(seq) - lapping_num + 1):
                temp_seq = ''
                for j in range(0, lapping_num):
                    temp_seq += seq[k + j]
                hamming_d = hamming_distance(temp_seq, temp_sub_seq)
                if hamming_d <= (lapping_num/2):  # judge hamming distance
                    return 0
        return 1


#  hairpin resist
def hairpin_resist(seq, strength=3):
    sub_seq = ''
    for i in range(len(seq)-2*strength):
        for j in range(strength):
            sub_seq += seq[i+j]
        revers_seq = complimentary_sequence(sub_seq)
        for k in range(i+strength+1, len(seq)-strength):
            sub_seq2 = ''
            for m in range(strength):
                sub_seq2 += seq[k+m]
            if revers_seq == sub_seq2:
                return 0
    return 1

def DNA_Generator(num, gc_ratio=0.5, source_str='ATC', lapping_num=3, flag=0):
    #
    while 1:
        # dna = randomDNA(num, source_str)
        if flag >= len(source_str)**num:
            print("All the cases have been traversed, but no sequence was found")
            return 0
        dna = sequential_seq(num, source_str, flag)
        flag += 1
        if GC_ratio_restrict(dna, gc_ratio) == 0:  # 1 for right gc_ratio, 0 for wrong gc_ratio
            continue
        elif seq_prevent(dna) == 0:  # 0 for sequence with secondary structure and vice verse
            continue
        elif overlapping_limit(dna, complimentary_sequence(dna), lapping_num) == 0:  # 0 refers to self complementary sequence
            continue
        elif hairpin_resist(dna, 3) == 0:  # 0 refers to hairpin structure, 3 refers to the least stem length
            continue
        break
    return dna, flag-1


def complimentary_sequence(seq):
    temp = ''
    seq = str(seq)
    i = len(seq)-1
    while i >= 0:
        if seq[i] == 'A':
            temp += 'T'
        elif seq[i] == 'T':
            temp += 'A'
        elif seq[i] == 'C':
            temp += 'G'
        elif seq[i] == 'G':
            temp += 'C'
        elif seq[i] == ' ':
            temp += ' '
        i -= 1
    return temp


#  reverse sequence
def reverse_sequence(seq):
    rever_seq = ''
    for i in range(len(seq)-1, -1, -1):
        rever_seq += seq[i]
    return rever_seq


# total seq num fit the resist
def total_seq_num(num, gc_ratio, source_str, lapping_num):
    total_num = 0
    flag = 0
    while 1:
        if flag >= len(source_str)**num:
            print("All the cases have been traversed")
            break
        dna = sequential_seq(num, source_str, flag)
        flag += 1
        if GC_ratio_restrict(dna, gc_ratio) == 0:  # 1 for right gc_ratio, 0 for wrong gc_ratio
            continue
        elif seq_prevent(dna) == 0:  # 0 for sequence with secondary structure and vice verse
            continue
        elif overlapping_limit(dna, complimentary_sequence(dna), lapping_num) == 0:  # 0 refers to self complementary sequence
            continue
        elif hairpin_resist(dna, 3) == 0:  # 0 refers to hairpin structure, 3 refers to the least stem length
            continue
        total_num += 1
    return total_num

# GC_ratio_restrict test
'''
dna = randomDNA(5, 'ACTG')
while GC_ratio_restrict(dna, 0.37) == 0:
    dna = randomDNA(5, 'ACTG')
print("GC ratio 37%:", dna)
'''

# hamming_distance test
'''
dna1 = "AAAAA"
dna2 = "TTAAA"
print("hamming distance is: ", hamming_distance(dna1, dna2))
'''

# run for test
# num = int(input('Input DNA number:'))
# DNA_Generator(num, gc_ratio=0.5, source_str='ATC', lapping_num=6)
# dna = DNA_Generator(25, 0.5, 'ACGT')
# print("Random DNA sequence 5'->3':" + seq_format(dna))


def main():
    # num = int(input("Input domain length (nt):"))
    # ratio = float(input("Input GC ratio (range(0,1)(eg. 0.37)):"))
    # alphabet = str(input("Input code alphabet ('ATCG') (eg.'ACT'):"))
    # lapping = int(input("Input lapping number (eg. 7):"))

    num = 4  # domain length (nt)
    ratio = 0.5  # GC ratio
    alphabet = 'ACT'  # code alphabet ('ATCG')
    lapping = 3  # lapping number
    flag = 0  # flag number
    conc = 1e-6  # total strands concentration (M)

    # print(f'total seq num is {total_seq_num(num, ratio, alphabet, lapping)}')  # count the total seq num

    domain, flag = DNA_Generator(num, ratio, alphabet, lapping, flag)
    M_temperature = melting_temp(domain, conc)
    print(f'flag = {flag}')
    print("Result sequence:")
    print(f'5\'-{domain}-3\'')
    print(' '*3 + '|'*num)
    print(f'3\'-{reverse_sequence(complimentary_sequence(domain))}-5\'')
    print(f'dG of the duplex is {calc_dG(domain):.2f} kcal/mol')
    print(f'Melting temperature equals to {M_temperature:.2f}°C')


if __name__ == '__main__':
    main()

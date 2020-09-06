import matplotlib.pyplot as plt
import random_DNA_sequence_generator as rDsg

'''seq refers to the sequence to be analysed (eg.'ATCGAAT'), 
sub_range refers to the accuracy of the analysis process'''
def distribution_GC_ratio(seq, sub_range=5):
    temp_seq = ''
    dis_ratio = []
    seq = rDsg.seq_format_no_blank(seq)
    for i in range(0, len(seq)-sub_range+1):
        for j in range(0, sub_range):
            temp_seq += seq[i+j]
        dis_ratio.append(get_GC_num(temp_seq)/len(temp_seq))
        temp_seq = ''
    return dis_ratio
def get_GC_num(seq):
    num = 0
    for i in range(0, len(seq)):
        if (seq[i] == 'G') | (seq[i] == 'C'):
            num += 1
    return num


'''
DNA = 'ACTGAT'
a = distribution_GC_ratio(DNA)
print(a)
'''
if __name__ == '__main__':
    dna_seq = str(input("Input the DNA sequence to analyse (5'->3'):"))
    gc_ratio = distribution_GC_ratio(dna_seq, 5)
    plt.figure(figsize=(8, 6), dpi=100)
    plt.plot(range(len(gc_ratio)), gc_ratio, color='red', linestyle='-', linewidth=3, label='GC_ratio_distribution')
    plt.legend(loc='upper right')
    plt.title('GC Ratio Distribution')
    plt.xlabel('Base index (nt)')
    plt.ylabel('GC Ratio')
    plt.ylim((0, 1))
    plt.xlim((0, len(gc_ratio)-1))
    plt.show()

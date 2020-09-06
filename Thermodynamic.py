"""
Gibbs free energy generator
Created by: Chaser
"""
from complimentary_seq import complimentary_sequence
from math import log

"""
Nearest neighbor thermodynamics
"""

_dH_dict = {'AA': -7.6, 'AT': -7.2, 'TA': -7.2, 'CA': -8.5,
            'GT': -8.4, 'CT': -7.8, 'GA': -8.2, 'CG': -10.6,
            'GC': -9.8, 'GG': -8.0, 'AG': -8.2, 'CC': -8.0,
            'TT': -7.6, 'AC': -8.5, 'TG': -8.4, 'TC': -7.8,
            'ini': 0.2, 'AT_terminal_penalty': 2.2, 'symmetry_correction': 0.0}  # unit: kcal/mol
_dS_dict = {'AA': -21.3, 'AT': -20.4, 'TA': -21.3, 'CA': -22.7,
            'GT': -22.4, 'CT': -21.0, 'GA': -22.2, 'CG': -27.2,
            'GC': -24.4, 'GG': -19.9, 'AG': -22.2, 'CC': -19.9,
            'TT': -21.3, 'AC': -22.7, 'TG': -22.4, 'TC': -21.0,
            'ini': 5.7, 'AT_terminal_penalty': 6.9, 'symmetry_correction': -1.4}  # unit: cal/(K*mol)
_dG_dict = {'AA': -1.00, 'AT': -0.88, 'TA': -0.58, 'CA': -1.45,
            'GT': -1.44, 'CT': -1.28, 'GA': -1.30, 'CG': -2.17,
            'GC': -2.24, 'GG': -1.84, 'AG': -1.30, 'CC': -1.84,
            'TT': -1.00, 'AC': -1.45, 'TG': -1.44, 'TC': -1.28,
            'ini': 1.96, 'AT_terminal_penalty': 0.05, 'symmetry_correction': 0.43}  # unit: kcal/mol
# The symmetry correction is applied only to self-complementary duplexes.
# The terminal AT penalty applies to each end of a duplex that has terminal AT.
# A duplex with both ends closed by AT pairs has a penalty of +1.0 kcal/mol for dG.


# 输入单链序列以及是否自我互补，输出对应互补配对的双链序列自由能
def calc_dG(seq, self_complementary=False):
    dg = 0
    AT_flag = 0
    dg += _dG_dict['ini']
    if self_complementary:
        dg += _dG_dict['symmetry_correction']
    for i in range(len(seq)-1):
        sub_str = seq[i]+seq[i+1]
        if i == 0 or i == len(seq)-2:
            if sub_str == 'AT' or sub_str == 'TA':
                AT_flag += 1
        dg += _dG_dict[sub_str]
    if AT_flag == 2:
        dg += 1.00
    elif AT_flag == 1:
        dg += _dG_dict['AT_terminal_penalty']
    return dg


# 输入单链序列以及是否自我互补，输出对应互补配对的双链序列焓变
def calc_dH(seq, self_complementary=False):
    dh = 0
    dh += _dH_dict['ini']
    if self_complementary:
        dh += _dH_dict['symmetry_correction']
    for i in range(len(seq)-1):
        sub_str = seq[i]+seq[i+1]
        dh += _dH_dict[sub_str]
    return dh


# 输入单链序列以及是否自我互补，输出对应互补配对的双链序列熵变
def calc_dS(seq, self_complementary=False):
    ds = 0
    ds += _dS_dict['ini']
    if self_complementary:
        ds += _dS_dict['symmetry_correction']
    for i in range(len(seq)-1):
        sub_str = seq[i]+seq[i+1]
        ds += _dS_dict[sub_str]
    return ds


#  溶解温度计算
def melting_temp(seq, total_molar_concentration=1e-6, self_complementary=False):  # total_molar_concentration (M)
    r = 1.98718  # 气体常数 R = 1.98718 cal/(K*mol)
    z = 4 if self_complementary is False else 1
    return (calc_dH(seq)*1000/(calc_dS(seq)+r*log(total_molar_concentration/z)))-273.15  # 返回溶解温度 单位为°C


def main():
    seq = 'CGAATGAAGAATGGAACGCA'

    print(f'dG of \'{seq}\' is {calc_dG(seq):.02f} kcal/mol')
    print(f'Melting temperature equals to {melting_temp(seq, 1e-6, False):.02f}°C')


if __name__ == '__main__':
    main()

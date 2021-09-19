#1.create csv file of positions which are different from each other or if any of the file includes mixed bases (Ex. A and T)
#2.create csv file that includes details of the position

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='identify positions which are different from each other')
parser.add_argument('input_initial_control', help='file name for initial control sequence after "/mnt/c/Users/naoki/Downloads/remdesivir/raw_data/"')
parser.add_argument('input_control', help='file name for control sequence after "/mnt/c/Users/naoki/Downloads/remdesivir/raw_data/"')
parser.add_argument('input_Remdesivir', help='file name for virus adapted to Remdesivir after "/mnt/c/Users/naoki/Downloads/remdesivir/raw_data/"')
parser.add_argument('input_high_concentration_Remdesivir', help='file name for virus adapted to high concentrarion Remdesivir after "/mnt/c/Users/naoki/Downloads/remdesivir/raw_data/"')
parser.add_argument('output_strain_number', help='number which represents the number of strain in the bracket "/mnt/c/Users/naoki/Downloads/remdesivir/analysis/differences_compared_with_con-seq/()_sequence.csv"')

args = parser.parse_args()

#read csv to dataframe
df_ini = pd.read_csv("/mnt/c/Users/naoki/Downloads/Remdesivir/raw_data/{}".format(args.input_initial_control), error_bad_lines=False, sep='\t')
df_con = pd.read_csv("/mnt/c/Users/naoki/Downloads/Remdesivir/raw_data/{}".format(args.input_control), error_bad_lines=False, sep='\t')
df_rem = pd.read_csv("/mnt/c/Users/naoki/Downloads/Remdesivir/raw_data/{}".format(args.input_Remdesivir), error_bad_lines=False, sep='\t')
df_hrem = pd.read_csv("/mnt/c/Users/naoki/Downloads/Remdesivir/raw_data/{}".format(args.input_high_concentration_Remdesivir), error_bad_lines=False, sep='\t')

#create initial most frequent sequence
def most_freq_nuc(df_file):
    freq_seq_list = []
    for nuc_A, nuc_C, nuc_G, nuc_T, nuc_gap, nuc_ins in zip(df_file['A'], df_file['C'], df_file['G'], df_file['T'], df_file['Gap'], df_file['Ins']):
        nuc_list = []
        nuc_list.append(nuc_A)
        nuc_list.append(nuc_C)
        nuc_list.append(nuc_G)
        nuc_list.append(nuc_T)
        nuc_list.append(nuc_gap)
        nuc_list.append(nuc_ins)
        if nuc_list[0] == max(nuc_list) and nuc_list[0] >= 90:
            freq_seq_list.append('A')
        elif nuc_list[0] == max(nuc_list) and nuc_list[0] < 90:
            freq_seq_list.append('(A)')
        elif nuc_list[1] == max(nuc_list) and nuc_list[1] >= 90:
            freq_seq_list.append('C')
        elif nuc_list[1] == max(nuc_list) and nuc_list[1] < 90:
            freq_seq_list.append('(C)')
        elif nuc_list[2] == max(nuc_list) and nuc_list[2] >= 90:
            freq_seq_list.append('G')
        elif nuc_list[2] == max(nuc_list) and nuc_list[2] < 90:
            freq_seq_list.append('(G)')
        elif nuc_list[3] == max(nuc_list) and nuc_list[3] >= 90:
            freq_seq_list.append('T')
        elif nuc_list[3] == max(nuc_list) and nuc_list[3] < 90:
            freq_seq_list.append('(T)')
        elif nuc_list[4] == max(nuc_list) and nuc_list[4] >= 90:
            freq_seq_list.append('-')
        elif nuc_list[4] == max(nuc_list) and nuc_list[4] < 90:
            freq_seq_list.append('(-)')
        elif nuc_list[5] == max(nuc_list) and nuc_list[5] >= 90:
            freq_seq_list.append('+')
        elif nuc_list[5] == max(nuc_list) and nuc_list[5] < 90:
            freq_seq_list.append('(+)')

    return freq_seq_list

ini_seq_list = most_freq_nuc(df_ini)
con_seq_list = most_freq_nuc(df_con)
rem_seq_list = most_freq_nuc(df_rem)
hrem_seq_list = most_freq_nuc(df_hrem)

#create dataframe which contains each sequence
df = pd.DataFrame(list(zip(df_ini['Pos'], df_ini['RefN'], ini_seq_list, con_seq_list, rem_seq_list, hrem_seq_list)), columns=['Pos', 'RefN', 'Ini_N', 'con', 'rem', 'high_rem'])
print(df)

#make a list of positions with differences
nuc_list = ['(A)', '(C)', '(G)', '(T)', '(-)', '(+)']
diff_list = []
for pos, ref, ini, con, rem, hrem in zip(df['Pos'], df['RefN'], df['Ini_N'], df['con'], df['rem'], df['high_rem']):
    if ini == con and con == rem and rem == hrem:
        pass
    elif ini != con or con != rem or rem!=hrem:
        diff_list.append(pos)
    elif ini in nuc_list or con in nuc_list or rem in nuc_list or hrem in nuc_list:
        diff_list.append(pos)

print(diff_list)

#convert the list of differences to dataframe
dataframe_list = []
for i in diff_list:
    dataframe_list.append(df.iloc[i-1])

df2 = pd.DataFrame(dataframe_list)

print(df2)

#to csv
df2.to_csv("/mnt/c/Users/naoki/Downloads/remdesivir/analysis/differences_compared_with_con-seq/sequence/{}_sequence_final_version.csv".format(args.output_strain_number), index=False)

#insert 'status' as a column
df_ini['status'] = 'Ini_control'
df_ini = df_ini.iloc[:, [0,22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]
df_con['status'] = 'control'
df_con = df_con.iloc[:, [0,22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]
df_rem['status'] = 'rem'
df_rem = df_rem.iloc[:, [0,22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]
df_hrem['status'] = 'high_rem'
df_hrem = df_hrem.iloc[:, [0,22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]

dataframe_concat_list = []
for i in diff_list:
    dataframe_concat_list.append(df_ini.iloc[i-1])
    dataframe_concat_list.append(df_con.iloc[i-1])
    dataframe_concat_list.append(df_rem.iloc[i-1])
    dataframe_concat_list.append(df_hrem.iloc[i-1])

df_detail = pd.DataFrame(dataframe_concat_list)

#drop column 'Unnamed: 21'
df_detail = df_detail.drop('Unnamed: 21', axis=1)
print(df_detail)

#to csv
df_detail.to_csv("/mnt/c/Users/naoki/Downloads/remdesivir/analysis/differences_compared_with_con-seq/details_of_diff/{}_details_final_version.csv".format(args.output_strain_number), index=False)
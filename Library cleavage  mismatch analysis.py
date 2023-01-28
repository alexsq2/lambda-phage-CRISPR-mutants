#script only works for gene J targets (398) unless "to_find" sequences are switched out
# skeleton of script that checks R1 and R2 files to see if they match




import xlsxwriter
import glob
import statistics

# take two sequences that are supposed to be reverse complement of each other, from stack overflow
def reverse_comp(in_sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = in_sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return  reverse_complement

# takes in a list with file pairs and outputs a dictionary with info on each sequence
def file_pair_analysis(file_pair, random_list, control = False):

    sequence_lines = 0  # this represents the total number of valid sequences found, valid = R1R2 match
    seq_count_lib = {}  # defining this library to add sequences with info to

    R1_file = open(file_pair[0])
    R2_file = open(file_pair[1])

    lines_R1 = R1_file.readlines()
    lines_R2 = R2_file.readlines()

    print(len(lines_R1), 'lines in R1 file')
    print(len(lines_R2), 'lines in R2 file')

    line_count_r1 = 0
    line_count_r2 = 0

    same_length_sequences = 0
    matching_sequences = 0

    for line_num in range(len(lines_R1)):  # lines_R1 and R2 SHOULD have the same length
        temp_sequence_R1 = ''
        temp_sequence_R2 = ''
        if lines_R1[line_num].find(to_find1_R1) != -1 and lines_R1[line_num].find(to_find2_R1) != -1:
            temp_sequence_R1 = lines_R1[line_num][lines_R1[line_num].find(to_find1_R1) + len(to_find1_R1): lines_R1[line_num].find(to_find2_R1)]
            line_count_r1 += 1
            # print(temp_sequence_R1, 'r1 line')
            # print(temp_sequence_R1)

        if lines_R2[line_num].find(to_find1_R2) != -1 and lines_R2[line_num].find(to_find2_R2) != -1:
            temp_sequence_R2 = lines_R2[line_num][lines_R2[line_num].find(to_find1_R2) + len(to_find1_R2): lines_R2[line_num].find(to_find2_R2)]
            # print(temp_sequence_R2, 'r2 line')
            line_count_r2 += 1


        if len(temp_sequence_R1) == len(temp_sequence_R2) and len(temp_sequence_R1) > 20 and len(
                temp_sequence_R2) > 20:  # if they are same length, rev complement R2, match with R1
            same_length_sequences += 1

            if temp_sequence_R2 == reverse_comp(temp_sequence_R1) :
                # print('found R1R2 pair')
                matching_sequences += 1
                sequence_lines += 1
                # print(temp_sequence_R1)
                temp_sequence = temp_sequence_R2 # the final sequence used for analysis and comparison
                # print(temp_sequence)

                if temp_sequence not in seq_count_lib and len(temp_sequence) == len(perfect_target):
                    mm_count_temp = 0
                    mm_pos_1 = 'NA'
                    mm_bin_1 = 'NA'
                    mm_pos_2 = 'NA'
                    mm_bin_2 = 'NA'


                    for x in range(len(temp_sequence)):  # counting total mismatches
                        if temp_sequence[x] != perfect_target[x]:
                            mm_count_temp += 1
                            if mm_count_temp == 1:
                                mm_pos_1 = x
                                for region in target_regions:
                                    if x in target_regions[region]:
                                        mm_bin_1 = region
                            if mm_count_temp == 2:
                                mm_pos_2 = x
                                for region in target_regions:
                                    if x in target_regions[region]:
                                        mm_bin_2 = region
                    if mm_count_temp < 3:
                        seq_count_lib[temp_sequence] = {'count': 1, 'length': len(temp_sequence),
                             'mismatches': mm_count_temp, 'mismatch1': mm_pos_1,'mismatch2' : mm_pos_2,
                                    'mm_bin_1' : mm_bin_1, 'mm_bin_2': mm_bin_2, 'instances': 1}
                        seq_count_lib[temp_sequence]['random'] = temp_sequence in random_list

                    if mm_count_temp > 3 and temp_sequence in random_list:
                        seq_count_lib[temp_sequence] = {'count': 1, 'length': len(temp_sequence),
                                                        'mismatches': mm_count_temp, 'mismatch1': mm_pos_1,
                                                        'mismatch2': mm_pos_2,
                                                        'mm_bin_1': mm_bin_1, 'mm_bin_2': mm_bin_2, 'instances': 1}
                        seq_count_lib[temp_sequence]['random'] = temp_sequence in random_list

                elif len(temp_sequence) == len(perfect_target):
                    seq_count_lib[temp_sequence]['count'] = seq_count_lib[temp_sequence]['count'] + 1

    print(matching_sequences, 'matching sequences')
    print(sequence_lines)
    print(same_length_sequences, 'same length sequences')


    for item in seq_count_lib:
        seq_count_lib[item]['fraction'] = seq_count_lib[item]['count'] / matching_sequences

        if seq_count_lib[item]['mm_bin_1'] != 'NA' and seq_count_lib[item]['mm_bin_2'] != 'NA':
            seq_count_lib[item]['mm_combo'] = seq_count_lib[item]['mm_bin_1'] + ' + ' + seq_count_lib[item]['mm_bin_2']
        elif seq_count_lib[item]['mm_bin_1'] != 'NA' and seq_count_lib[item]['mismatches'] == 1:
            seq_count_lib[item]['mm_combo'] = seq_count_lib[item]['mm_bin_1']

        elif seq_count_lib[item]['mismatches'] == 0:
            seq_count_lib[item]['mm_combo'] = "perfect"


    if control == False:
        for item in seq_count_lib:
            if item in analysis_control_lib_combined:
                seq_count_lib[item]['enrichment'] = seq_count_lib[item]['fraction'] / analysis_control_lib_combined[item]['fraction']
            else:
                seq_count_lib[item]['enrichment'] = 1

        for item in seq_count_lib:
            if item in analysis_control_lib_combined:
                seq_count_lib[item]['z_numerator'] = seq_count_lib[item]['fraction'] - analysis_control_lib_combined[item]['fraction']
            else:
                seq_count_lib[item]['z_numerator'] = 0
        #now making a z-score value using stdev of all z_numberators

        fraction_sums = []
        for item in seq_count_lib:
            fraction_sums.append(seq_count_lib[item]['fraction'])
        standard_dev_exp = statistics.stdev(fraction_sums)

        for item in seq_count_lib:
            seq_count_lib[item]['z-score'] = seq_count_lib[item]['z_numerator'] / standard_dev_exp


    return seq_count_lib


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#start of actual analysis
#data that is specific for the target being analyzed
workbook = xlsxwriter.Workbook('(top) data_L_library.xlsx')

perfect_target = 'TTTGTACTGTCCACGCCGACGGAA' #FORWARD orientation

to_find1_R1 = 'CCGAGCTTGACCA'  #extract the target sequence between these two sequences from the R1 file
to_find2_R1 = 'ATCCTGAAGCAC'

to_find1_R2 = 'AGGAT'  #extract the target sequence between these two sequences from the R2 file
to_find2_R2 = 'TGGTC'
random_sequences_list = ['GCATTCTTATTAATACATTTGAAA','CGCGCCCAACTGACGCTAGGCAAG','TCAGTGCAGGCTCCCGTGTTAGGA','TAAGGGTAAACATACAAGTCGATA','CATCCAGCACTCTACGGTTCCTCC']

complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]
target_regions = {'PAM': [0,1,2,3], 'seed': [4,5,6,7,8,9], 'mid': [10,11,12,13,14,15,16], 'distal':[17,18,19,20,21,22,23] }

# selecting files for control sample, which the experimental samples will be compared with


control_file_list = []
control_file_list += glob.glob('L_Fn_Supercoil_10_Control*')
control_file_list_pairs = []

control_file_list_length_half = int(len(control_file_list) / 2)
for x in range(control_file_list_length_half):
    control_file_list_pairs.append([control_file_list[2 * x]] + [control_file_list[2 * x + 1]])
print(control_file_list_pairs[0])

analysis_control_lib_combined = {}
for pair in control_file_list_pairs:
    print(pair)
    control_lib = file_pair_analysis(pair, random_sequences_list, control = True) #analyzing control files, making library
    if len(control_lib) > 2400:
        for item in control_lib:
            if item in analysis_control_lib_combined:
                analysis_control_lib_combined[item]['fraction'] += analysis_control_lib_combined[item]['fraction']
                analysis_control_lib_combined[item]['instances'] += 1
            else:
                analysis_control_lib_combined[item] = control_lib[item]


for item in analysis_control_lib_combined:
    analysis_control_lib_combined[item]['fraction'] = analysis_control_lib_combined[item]['fraction']/analysis_control_lib_combined[item]['instances']



#collecting experimental samples, arranging so 1minute and 30minute samples are next to each other
exp_file_list_30min = []
exp_file_list_1min = []
exp_file_list = []
# exp_file_list_1min += glob.glob('L_As_Supercoil_1_1*')
# exp_file_list_30min += glob.glob('L_As_Supercoil_1_30*')
# exp_file_list_1min += glob.glob('L_As_Supercoil_2_1*')
# exp_file_list_30min += glob.glob('L_As_Supercoil_2_30*')
# exp_file_list_1min += glob.glob('L_As_Supercoil_5_1*')
# exp_file_list_30min += glob.glob('L_As_Supercoil_5_30*')
# exp_file_list_1min += glob.glob('L_As_Supercoil_10_1*')
# exp_file_list_30min += glob.glob('L_As_Supercoil_10_30*')
exp_file_list_1min += glob.glob('L_Lb_Supercoil_1_1*')
exp_file_list_30min += glob.glob('L_Lb_Supercoil_1_30*')
exp_file_list_1min += glob.glob('L_Lb_Supercoil_2_1*')
exp_file_list_30min += glob.glob('L_Lb_Supercoil_2_30*')
exp_file_list_1min += glob.glob('L_Lb_Supercoil_5_1*')
exp_file_list_30min += glob.glob('L_Lb_Supercoil_5_30*')
exp_file_list_1min += glob.glob('L_Lb_Supercoil_10_1*')
exp_file_list_30min += glob.glob('L_Lb_Supercoil_10_30*')

# exp_file_list_1min += glob.glob('*1min*')
# exp_file_list_30min += glob.glob('*30min*')

exp_file_list_1min = [item for item in exp_file_list_1min if 'Li' not in item]
exp_file_list_1min = [item for item in exp_file_list_1min if 'C1' not in item]
exp_file_list_1min = [item for item in exp_file_list_1min if 'C2' not in item]

exp_file_list_30min = [item for item in exp_file_list_30min if 'Li' not in item]
exp_file_list_30min = [item for item in exp_file_list_30min if 'C1' not in item]
exp_file_list_30min = [item for item in exp_file_list_30min if 'C2' not in item]

print(exp_file_list_1min)

file_list_length_half = int(len(exp_file_list_1min) / 2)
total_file_list_pairs_1min = [] #this list will contain pairs of files ordered by their order in total file list, check the list is in order
total_file_list_pairs_30min = []


for x in range(file_list_length_half):
    total_file_list_pairs_1min.append([exp_file_list_1min[2 * x]] + [exp_file_list_1min[2 * x + 1]])
    total_file_list_pairs_30min.append([exp_file_list_30min[2 * x]] + [exp_file_list_30min[2 * x + 1]])



total_file_list_pairs = []
if len(total_file_list_pairs_1min) == len(total_file_list_pairs_30min):
    for pos in range(len(total_file_list_pairs_1min)):
        total_file_list_pairs.append(total_file_list_pairs_1min[pos])
        total_file_list_pairs.append(total_file_list_pairs_30min[pos])

for pos in range(len(total_file_list_pairs)):
    print(total_file_list_pairs[pos])


all_sample_info = [] # this will be written on the firstsheet containing info for all of the samples processed

custom_format = workbook.add_format({'num_format': ';;;'})

#making worksheets containing info from all of the samples
worksheet_all = workbook.add_worksheet('All info')

worksheet_single_heat = workbook.add_worksheet('Singles')
single_heat_row = 2

worksheet_double_heat = workbook.add_worksheet('Doubles')
double_heat_row = 2

counter_1min_30min = 0 #keeps track of alternating 1min and 30minute samples for heat maps

row_random_seq = 0

for pair in total_file_list_pairs:
    file_name = pair[0][0:22]
    print(file_name)

    worksheet = workbook.add_worksheet(file_name)

    # Add the worksheet data to be plotted.
    dictionary_file_pair_analysis = file_pair_analysis(pair,random_sequences_list) #analyzing experimental files, making library
    print(dictionary_file_pair_analysis)

    data = []
    for item in dictionary_file_pair_analysis:
        data.append(dictionary_file_pair_analysis[item]['fraction'])
    data.sort()
    data.reverse()

    all_sample_info.append(len(dictionary_file_pair_analysis))

    worksheet.write_column('A1', data)

    # Create a new chart object.
    chart = workbook.add_chart({'type': 'column'})

    # Add a series to the chart.
    file_name_string = '=' + file_name + '!$A$1:$A$2500'
    # chart.add_series({'values': '=file_name!$A$1:$A$100'})
    chart.add_series({'values': file_name_string})

    chart.set_x_axis({'name': 'unique sequences','visible': False})
    chart.set_y_axis({'name': 'proportion of reads','min': 0, 'max': 0.004})

    chart.set_title({'name': file_name})
    # Insert the chart into the worksheet.
    worksheet.insert_chart('BL1', chart)

    row = 2
    column = 2
    for item in dictionary_file_pair_analysis:
        worksheet.write(row, column+1, item)
        worksheet.write(row, column+2, dictionary_file_pair_analysis[item]['enrichment'])
        worksheet.write(row, column+3, dictionary_file_pair_analysis[item]['mm_bin_1'])
        worksheet.write(row, column+4, dictionary_file_pair_analysis[item]['mm_bin_2'])
        row += 1

    enrich_dict = {}
    singles_enrich_dict = {}
    perfect_found = False
    for item in dictionary_file_pair_analysis:
        if dictionary_file_pair_analysis[item]['mm_combo'] not in enrich_dict:
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']] = [  dictionary_file_pair_analysis[item]['enrichment'], 1, 0, 0]
            singles_enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']] = [dictionary_file_pair_analysis[item]['enrichment']] # each entry will have all enrichement values for all sequences with that combo
        else:
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][1] += 1
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][2] += dictionary_file_pair_analysis[item]['enrichment']
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][3] = enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][2] / enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][1]
            singles_enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']].append(dictionary_file_pair_analysis[item]['enrichment'])

        if 'perfect' in enrich_dict and perfect_found == False:
            perfect_found = True
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][2] += dictionary_file_pair_analysis[item]['enrichment']
            enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][3] = enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][2] / enrich_dict[dictionary_file_pair_analysis[item]['mm_combo']][1]

    print(enrich_dict)
    print(singles_enrich_dict)

#making 2D heat map for double mismatches
    grid24_list = []

    for x in range(24):
        zero_list = []
        for x in range(24):
            zero_list.append(0)
        grid24_list.append(zero_list)
    print(grid24_list)

    row = 1
    column = 49
    for num in range(-4, 21):
        if num != 0:
            worksheet.write(row, column, num)
            row += 1

    for item in dictionary_file_pair_analysis:
        if dictionary_file_pair_analysis[item]['mismatch1'] != 'NA' and dictionary_file_pair_analysis[item]['mismatch2'] != 'NA':
            grid24_list[dictionary_file_pair_analysis[item]['mismatch1']][dictionary_file_pair_analysis[item]['mismatch2']] += dictionary_file_pair_analysis[item]['enrichment']
    print(grid24_list)


    if counter_1min_30min == 0: # 1 minutes samples
        for row in range(len(grid24_list)):
            for column in range(len(grid24_list)):
                worksheet.write(row + 1, column + 25, grid24_list[row][column], custom_format)
                worksheet_double_heat.write(row + 1 + double_heat_row, column*2, grid24_list[row][column], custom_format)

    if counter_1min_30min == 1:
        for row in range(len(grid24_list)):
            for column in range(len(grid24_list)):
                worksheet.write(row + 1, column + 25, grid24_list[row][column], custom_format)
                worksheet_double_heat.write(row + 1 + double_heat_row, (column * 2) + 1, grid24_list[row][column],
                                            custom_format)

    row = 0
    column = 25
    for num in range(-4, 21):
        if num != 0:
            worksheet.write(double_heat_row -1, column, num)
            column += 1




    for num in range(1000):
        worksheet_double_heat.set_row(num, 45)
    worksheet_double_heat.set_column(0, 100, 3.4)

    row = 0
    column = 0
    for num in range(-4,21):
        if num != 0:
            worksheet_double_heat.write(row + double_heat_row, column * 2, str(num))
            column += 1
    for num in range(-4,21):
        if num != 0:
            worksheet_double_heat.write(row + double_heat_row + 1, 48, str(num))
            row += 1

    worksheet_double_heat.write(double_heat_row - 1, 0, file_name)
    worksheet_double_heat.conditional_format(0, 0, 1000, 1000,
                                             {'type': '2_color_scale', 'min_color': 'white', 'max_color': 'red'})


#making heat map for single mismatches only

    perfect_target_list = list(perfect_target)
    row = 27
    column = 25
    for base in perfect_target_list:
        worksheet.write(row, column, base)
        column += 1
    base_num_dict = {'A': 0, 'T':1, 'C':2, 'G':3 }
    for item in base_num_dict:
        worksheet.write(base_num_dict[item] + 28,24,item)
    row = 28
    column = 25
    for item in dictionary_file_pair_analysis:
        if dictionary_file_pair_analysis[item]['mismatches'] == 1:
            if dictionary_file_pair_analysis[item]['mismatch1'] != 'NA':
                worksheet.write(base_num_dict[   item[dictionary_file_pair_analysis[item]['mismatch1']]   ] + row,
                                dictionary_file_pair_analysis[item]['mismatch1'] + column,
                                dictionary_file_pair_analysis[item]['enrichment'],
                                custom_format)

        if counter_1min_30min == 0:
            if dictionary_file_pair_analysis[item]['mismatches'] == 1:
                if dictionary_file_pair_analysis[item]['mismatch1'] != 'NA':
                    worksheet_single_heat.write(base_num_dict[   item[dictionary_file_pair_analysis[item]['mismatch1']]   ] + single_heat_row + 1,
                                    dictionary_file_pair_analysis[item]['mismatch1'] * 2 ,
                                    dictionary_file_pair_analysis[item]['enrichment'],
                                    custom_format)
        if counter_1min_30min == 1:
            if dictionary_file_pair_analysis[item]['mismatches'] == 1:
                if dictionary_file_pair_analysis[item]['mismatch1'] != 'NA':
                    worksheet_single_heat.write(
                        base_num_dict[item[dictionary_file_pair_analysis[item]['mismatch1']]] + single_heat_row + 1,
                        (dictionary_file_pair_analysis[item]['mismatch1'] *2 ) + 1,
                        dictionary_file_pair_analysis[item]['enrichment'],
                        custom_format)

    worksheet_single_heat.write(single_heat_row -1, 0, file_name)

    for item in base_num_dict:
        worksheet_single_heat.write(base_num_dict[item] + single_heat_row + 1 ,48,item)
    column = 0
    for base in perfect_target_list:
        worksheet_single_heat.write(single_heat_row, column, base)
        column += 1

    for num in range(1000):
        worksheet_single_heat.set_row(num, 45)
    worksheet_single_heat.set_column(0, 100, 3.4)

    worksheet_single_heat.conditional_format(0, 0, 1000, 1000,
                                             {'type': '2_color_scale', 'min_color': 'white', 'max_color': 'red'})

    if counter_1min_30min == 1:
        double_heat_row += 27
        counter_1min_30min = 0
        single_heat_row += 7
    else:
        counter_1min_30min += 1


#writing more info in each separate sheet for each sample

    for num in range(100):
        worksheet.set_row(num, 15)
    worksheet.set_column(0, 100, 2.2)

    worksheet.conditional_format('Z2:AW28', {'type': '2_color_scale', 'min_color': 'white', 'max_color': 'red'})
    worksheet.conditional_format('Z29:AW32', {'type': '2_color_scale', 'min_color': 'white', 'max_color': 'red'})


    enrich_list = []
    for item in enrich_dict:
        enrich_list.append([item, enrich_dict[item][3]])
    enrich_list.sort()
    row = 2
    column = 7
    for item in enrich_list:
        worksheet.write(row, column + 1, item[0])
        worksheet.write(row, column + 3, item[1])
        row += 1

    worksheet.write(0,9, 'average enrichment')

    #writing single enrich dict

    row = 20
    column = 9

    for item in singles_enrich_dict:
        worksheet.write(row,column,item)
        column += 1

    column = 9
    for item in singles_enrich_dict:
        row = 21
        for freq in singles_enrich_dict[item]:
            worksheet.write(row, column, freq)
            row += 1
        column += 1


#writing enrichment of random NT sequences
    column = 1
    worksheet_all.write(row_random_seq, column, file_name)
    column += 1

    for item in dictionary_file_pair_analysis:
        if dictionary_file_pair_analysis[item]['random'] == True:
            print(item,dictionary_file_pair_analysis[item]['enrichment'])
            worksheet_all.write(row_random_seq,column,dictionary_file_pair_analysis[item]['enrichment'])
            column += 1
    row_random_seq += 1


worksheet_all.write_column('A1', all_sample_info)
workbook.close()






import pyhmmer
import subprocess
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
from Bio import SeqIO
from pyhmmer.easel import TextSequence, TextMSA, DigitalSequenceBlock
import pandas as pd
import matplotlib.pyplot as plt
import csv


"""
Train L2_i1.hmm with L2_i1.fasta
"""
# Set the amino acid alphabet for the HMM
alphabet = pyhmmer.easel.Alphabet.amino()

# Path to the FASTA file for L2 sequence training (Salima's dataset)
fasta_file_path = '../../dataset/reproduce_experiment_for_iteration1/train/L2_i1.fasta' # salima

# Load the sequences from the FASTA file
L2_sequences_train_load = SeqIO.parse(fasta_file_path, "fasta")

# Convert the sequences to a format suitable for the HMM (TextSequence objects)
L2_sequences_train = [
    TextSequence(name=bytes(seq.id, "utf-8"), sequence=str(seq.seq))
    for seq in L2_sequences_train_load
]

# Create a Multiple Sequence Alignment (MSA) object for the L2 sequences
msa_L2 = TextMSA(
    name=b"L2_i1", # salima
    sequences=L2_sequences_train
)

# Digitize the MSA using the defined alphabet (amino acids)
digital_msa_L2 = msa_L2.digitize(alphabet)

# Initialize an HMM builder using the alphabet
builder = pyhmmer.plan7.Builder(alphabet)

# Set the background probabilities (default for amino acids)
background = pyhmmer.plan7.Background(alphabet)

# Build the HMM from the MSA and background model
hmm_L2, _, _ = builder.build_msa(digital_msa_L2, background)

# Write the resulting HMM to a file for L2 sequences
with open("../../hmm/reproduce_experiment_for_iteration1/L2_i1.hmm", "wb") as output_file: 
    hmm_L2.write(output_file)

"""
Train L3_i1.hmm with L3_i1.fasta
"""
# Repeat the process for L3 sequences
# Set the amino acid alphabet for the HMM
alphabet = pyhmmer.easel.Alphabet.amino()

# Path to the FASTA file for L3 sequence training (Salima's dataset)
fasta_file_path = '../../dataset/reproduce_experiment_for_iteration1/train/L3_i1.fasta' # salima

# Load the sequences from the FASTA file
L3_sequences_train_load = SeqIO.parse(fasta_file_path, "fasta")

# Convert the sequences to a format suitable for the HMM (TextSequence objects)
L3_sequences_train = [
    TextSequence(name=bytes(seq.id, "utf-8"), sequence=str(seq.seq))
    for seq in L3_sequences_train_load
]

# Create a Multiple Sequence Alignment (MSA) object for the L3 sequences
msa_L3 = TextMSA(
    name=b"L3_i1", # salima
    sequences=L3_sequences_train
)

# Digitize the MSA using the defined alphabet (amino acids)
digital_msa_L3 = msa_L3.digitize(alphabet)

# Initialize an HMM builder using the alphabet
builder = pyhmmer.plan7.Builder(alphabet)

# Set the background probabilities (default for amino acids)
background = pyhmmer.plan7.Background(alphabet)

# Build the HMM from the MSA and background model
hmm_L3, _, _ = builder.build_msa(digital_msa_L3, background)

# Write the resulting HMM to a file for L3 sequences
with open("../../hmm/reproduce_experiment_for_iteration1/L3_i1.hmm", "wb") as output_file: 
    hmm_L3.write(output_file)


"""
Train All_i1.hmm with all_i1.fasta
"""
# Repeat the process for L2+L3 (All) sequences
# Set the amino acid alphabet for the HMM
alphabet = pyhmmer.easel.Alphabet.amino()

# Path to the FASTA file for L3 sequence training (Salima's dataset)
fasta_file_path = '../../dataset/reproduce_experiment_for_iteration1/train/All_i1.fasta' # salima

# Load the sequences from the FASTA file
All_sequences_train_load = SeqIO.parse(fasta_file_path, "fasta")

# Convert the sequences to a format suitable for the HMM (TextSequence objects)
All_sequences_train = [
    TextSequence(name=bytes(seq.id, "utf-8"), sequence=str(seq.seq))
    for seq in All_sequences_train_load
]

# Create a Multiple Sequence Alignment (MSA) object for the L2+L3 (All) sequences
msa_All = TextMSA(
    name=b"All_i1", # salima
    sequences=All_sequences_train
)

# Digitize the MSA using the defined alphabet (amino acids)
digital_msa_All = msa_All.digitize(alphabet)

# Initialize an HMM builder using the alphabet
builder = pyhmmer.plan7.Builder(alphabet)

# Set the background probabilities (default for amino acids)
background = pyhmmer.plan7.Background(alphabet)

# Build the HMM from the MSA and background model
hmm_All, _, _ = builder.build_msa(digital_msa_All, background)

# Write the resulting HMM to a file for L2+L3 (All) sequences
with open("../../hmm/reproduce_experiment_for_iteration1/All_i1.hmm", "wb") as output_file: 
    hmm_All.write(output_file)


"""
Use L2 & L3 hmm to get score of protein sequence
"""
with pyhmmer.plan7.HMMFile("../../hmm/reproduce_experiment_for_iteration1/L2_i1.hmm") as hmm_file:
    hmm_L2 = next(hmm_file)
with pyhmmer.plan7.HMMFile("../../hmm/reproduce_experiment_for_iteration1/L3_i1.hmm") as hmm_file:
    hmm_L3 = next(hmm_file)


# 創建字母表
alphabet = pyhmmer.easel.Alphabet.amino()
fasta_file_path = '../../dataset/reproduce_experiment_for_iteration1/test/interpro_L7.fasta' # salima
sequences = SeqIO.parse(fasta_file_path, "fasta")

search_sequences = []
sequence_dict = {}
for seq in sequences:
    # 打印序列ID
    id_name = bytes(seq.id, "utf-8")
    # 添加 TextSequence 對象到列表
    search_sequences.append(TextSequence(name=id_name, sequence=str(seq.seq)))
    # 同時存儲序列信息到字典
    sequence_dict[id_name] = (len(seq.seq), str(seq.seq))


# 將序列轉換為 DigitalSequence
digital_search_sequences = [seq.digitize(alphabet) for seq in search_sequences]

# 創建搜索管道並搜索 HMM
background = pyhmmer.plan7.Background(alphabet)
pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

sequences = pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_search_sequences)

# 使用 DigitalSequenceList 包裝序列
hits_L2 = pipeline.search_hmm(query=hmm_L2, sequences=sequences)
hits_L3 = pipeline.search_hmm(query=hmm_L3, sequences=sequences)

fixed_length = 33

# 構建名稱和分數字典以便於查找
l2_scores = {hit.domains[0].alignment.target_name.decode(): hit.score for hit in hits_L2}
l3_scores = {hit.domains[0].alignment.target_name.decode(): hit.score for hit in hits_L3}


# 獲取所有名稱以確保完整打印
all_names = list(set(l2_scores.keys()) | set(l3_scores.keys()))
NA_name = []
# 打開一個 CSV 文件以寫入數據
with open('../../result/reproduce_experiment_for_iteration1/i1_interpro_L7.csv', 'w', newline='') as csvfile:
    fieldnames = ['Name', 'L2 Score', 'L3 Score', 'Length', 'Sequence']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # 寫入標題行
    writer.writeheader()

    for name in all_names:
        # formatted_name = name.ljust(fixed_length)
        formatted_name = name
        l2_score = l2_scores.get(name, "N/A")
        l3_score = l3_scores.get(name, "N/A")
        length, sequence = sequence_dict.get(bytes(name, "utf-8"), (0, "N/A"))
        if l2_score == "N/A" or l3_score == "N/A":
            NA_name.append(name)
            continue

                # 寫入數據行
        writer.writerow({
            'Name': formatted_name.split('|')[0],
            'L2 Score': f"{l2_score:.2f}",
            'L3 Score': f"{l3_score:.2f}",
            'Length': length,
            'Sequence': sequence
        })
        # print(f"{formatted_name.split('-')[1]},{l2_score:.2f},{l3_score:.2f},{length},{sequence}") # for decoy
        # print(f"{formatted_name.split('|')[0]},{l2_score:.2f},{l3_score:.2f},{length},{sequence}") 
  

file_path = '../../result/reproduce_experiment_for_iteration1/i1_NA_names.txt'
with open(file_path, 'w') as file:
    for name in NA_name:
        file.write(name + '\n')  
# print(f"Names with N/A scores have been saved to {file_path}.")  


"""
plot result based on the csv
"""
import pandas as pd
import matplotlib.pyplot as plt

# 定义CSV文件路径
files = ['../../result/reproduce_experiment_for_iteration1/i1_interpro_L7.csv']

# 标签列表
labels = ['i1_interpro']


# L2 名稱集合
L2_name = {
    "sp|Q9NQ29|LUC7L_HUMAN", # 1
    "sp|Q9CYI4|LUC7L_MOUSE", # 2
    "tr|A0A3Q2UH51|A0A3Q2UH51_CHICK", # 3
    "tr|Q6GLC1|Q6GLC1_XENTR", # 4
    "tr|H3AZN1|H3AZN1_LATCH", # 5
    "tr|Q6IQD4|Q6IQD4_DANRE", # 6
    "sp|Q9Y383|LC7L2_HUMAN", # 7
    "sp|Q7TNC4|LC7L2_MOUSE", # 8
    "tr|Q5ZLW7|Q5ZLW7_CHICK", # 9
    "tr|Q28EN5|Q28EN5_XENTR", # 10
    "tr|H2ZXS0|H2ZXS0_LATCH", # 11
    "tr|A0A8M1NDA7|A0A8M1NDA7_DANRE", # 12
    "tr|A0A8C4NH22|A0A8C4NH22_EPTBU", # 13
    "tr|F6S0J6|F6S0J6_CIOIN", # 14
    "tr|A0A182H045|A0A182H045_AEDAL", # 15
    "tr|Q9VVI1|Q9VVI1_DROME", # 16
    "tr|B3RNA2|B3RNA2_TRIAD", # 17
    "tr|T2M333|T2M333_HYDVU", # 18
    "sp|Q09217|YP68_CAEEL", # 19
    "tr|A0A1X7V7G6|A0A1X7V7G6_AMPQE", # 20

    "tr|F4IZU3|F4IZU3_ARATH", # 21
    "tr|Q940U9|Q940U9_ARATH", # 22
    "tr|A0A3Q7GLL9|A0A3Q7GLL9_SOLLC", # 23
    "tr|U5D3Z2|U5D3Z2_AMBTC", # 24
    "tr|Q75LD6|Q75LD6_ORYSJ", # 25
    "tr|B4FUS0|B4FUS0_MAIZE", # 26
    "tr|B4FXT5|B4FXT5_MAIZE", # 27
    "tr|Q6K8R4|Q6K8R4_ORYSJ", # 28
    "tr|A0A2R6WMH7|A0A2R6WMH7_MARPO", # 29
    "tr|A0A388LKU4|A0A388LKU4_CHABU", # 30

    "tr|Q7S615|Q7S615_NEUCR", # 31
    "tr|Q5BCF4|Q5BCF4_EMENI", # 32
    "tr|F9X0V8|F9X0V8_ZYMTI", # 33
    "sp|Q07508|LUC7_YEAST", # 34
    "sp|Q9USM4|LUC7_SCHPO", # 35
    "tr|Q5KBY8|Q5KBY8_CRYNJ", # 36
    "tr|E3L016|E3L016_PUCGT", # 37
    "tr|A0A139ATP5|A0A139ATP5_GONPJ" # 38
}


# L3 名稱集合
L3_name = {
    "sp|Q5SUF2|LC7L3_MOUSE", # 1
    "sp|O95232|LC7L3_HUMAN", # 2
    "tr|F1NXR8|F1NXR8_CHICK", # 3
    "tr|Q28G85|Q28G85_XENTR", # 4
    "tr|H3ABC3|H3ABC3_LATCH", # 5
    "tr|Q1ED13|Q1ED13_DANRE", # 6
    "tr|B3RQ49|B3RQ49_TRIAD", # 7
    "tr|A0A8B6XHB2|A0A8B6XHB2_HYDVU", # 8
    "tr|F6Z8S6|F6Z8S6_CIOIN", # 9
    "tr|A0A8C4N567|A0A8C4N567_EPTBU", # 10
    "tr|A0A182GZD3|A0A182GZD3_AEDAL", # 11
    "tr|Q9W3X8|Q9W3X8_DROME", # 12
    "tr|Q8ITY5|Q8ITY5_CAEEL", # 13
    "tr|H2L0S8|H2L0S8_CAEEL", # 14
    "tr|A0A1X7UMY2|A0A1X7UMY2_AMPQE", # 15

    "tr|P94088|P94088_ARATH", # 16
    "tr|A0A3Q7GAS3|A0A3Q7GAS3_SOLLC", # 17
    "tr|W1P1B9|W1P1B9_AMBTC", # 18
    "tr|B4FE46|B4FE46_MAIZE", # 19
    "tr|Q84YS0|Q84YS0_ORYSJ", # 20
    "tr|A0A176W202|A0A176W202_MARPO", # 21
    "tr|A0A388KYX9|A0A388KYX9_CHABU", # 22
    "tr|A0A2R6XJM8|A0A2R6XJM8_MARPO", # 23
    "tr|A0A2K3CY04|A0A2K3CY04_CHLRE", # 24
    "tr|A8HP50|A8HP50_CHLRE", # 25
    "tr|A4SBD0|A4SBD0_OSTLU" # 26
}


classification = {
    "sp|Q9NQ29|LUC7L_HUMAN": 'Animal', # 1
    "sp|Q9CYI4|LUC7L_MOUSE": 'Animal', # 2
    "tr|A0A3Q2UH51|A0A3Q2UH51_CHICK": 'Animal', # 3
    "tr|Q6GLC1|Q6GLC1_XENTR": 'Animal', # 4
    "tr|H3AZN1|H3AZN1_LATCH": 'Animal', # 5
    "tr|Q6IQD4|Q6IQD4_DANRE": 'Animal', # 6
    "sp|Q9Y383|LC7L2_HUMAN": 'Animal', # 7
    "sp|Q7TNC4|LC7L2_MOUSE": 'Animal', # 8
    "tr|Q5ZLW7|Q5ZLW7_CHICK": 'Animal', # 9
    "tr|Q28EN5|Q28EN5_XENTR": 'Animal', # 10
    "tr|H2ZXS0|H2ZXS0_LATCH": 'Animal', # 11
    "tr|A0A8M1NDA7|A0A8M1NDA7_DANRE": 'Animal', # 12
    "tr|A0A8C4NH22|A0A8C4NH22_EPTBU": 'Animal', # 13
    "tr|F6S0J6|F6S0J6_CIOIN": 'Animal', # 14
    "tr|A0A182H045|A0A182H045_AEDAL": 'Animal', # 15
    "tr|Q9VVI1|Q9VVI1_DROME": 'Animal', # 16
    "tr|B3RNA2|B3RNA2_TRIAD": 'Animal', # 17
    "tr|T2M333|T2M333_HYDVU": 'Animal', # 18
    "sp|Q09217|YP68_CAEEL": 'Animal', # 19
    "tr|A0A1X7V7G6|A0A1X7V7G6_AMPQE": 'Animal', # 20

    "tr|F4IZU3|F4IZU3_ARATH": 'Plant', # 21
    "tr|Q940U9|Q940U9_ARATH": 'Plant', # 22
    "tr|A0A3Q7GLL9|A0A3Q7GLL9_SOLLC": 'Plant', # 23
    "tr|U5D3Z2|U5D3Z2_AMBTC": 'Plant', # 24
    "tr|Q75LD6|Q75LD6_ORYSJ": 'Plant', # 25
    "tr|B4FUS0|B4FUS0_MAIZE": 'Plant', # 26
    "tr|B4FXT5|B4FXT5_MAIZE": 'Plant', # 27
    "tr|Q6K8R4|Q6K8R4_ORYSJ": 'Plant', # 28
    "tr|A0A2R6WMH7|A0A2R6WMH7_MARPO": 'Plant', # 29
    "tr|A0A388LKU4|A0A388LKU4_CHABU": 'Plant', # 30

    "tr|Q7S615|Q7S615_NEUCR": 'Fungi', # 31
    "tr|Q5BCF4|Q5BCF4_EMENI": 'Fungi', # 32
    "tr|F9X0V8|F9X0V8_ZYMTI": 'Fungi', # 33
    "sp|Q07508|LUC7_YEAST": 'Fungi', # 34
    "sp|Q9USM4|LUC7_SCHPO": 'Fungi', # 35
    "tr|Q5KBY8|Q5KBY8_CRYNJ": 'Fungi', # 36
    "tr|E3L016|E3L016_PUCGT": 'Fungi', # 37
    "tr|A0A139ATP5|A0A139ATP5_GONPJ": 'Fungi', # 38

    #################################################

    "sp|Q5SUF2|LC7L3_MOUSE": 'Animal', # 1
    "sp|O95232|LC7L3_HUMAN": 'Animal', # 2
    "tr|F1NXR8|F1NXR8_CHICK": 'Animal', # 3
    "tr|Q28G85|Q28G85_XENTR": 'Animal', # 4
    "tr|H3ABC3|H3ABC3_LATCH": 'Animal', # 5
    "tr|Q1ED13|Q1ED13_DANRE": 'Animal', # 6
    "tr|B3RQ49|B3RQ49_TRIAD": 'Animal', # 7
    "tr|A0A8B6XHB2|A0A8B6XHB2_HYDVU": 'Animal', # 8
    "tr|F6Z8S6|F6Z8S6_CIOIN": 'Animal', # 9
    "tr|A0A8C4N567|A0A8C4N567_EPTBU": 'Animal', # 10
    "tr|A0A182GZD3|A0A182GZD3_AEDAL": 'Animal', # 11
    "tr|Q9W3X8|Q9W3X8_DROME": 'Animal', # 12
    "tr|Q8ITY5|Q8ITY5_CAEEL": 'Animal', # 13
    "tr|H2L0S8|H2L0S8_CAEEL": 'Animal', # 14
    "tr|A0A1X7UMY2|A0A1X7UMY2_AMPQE": 'Animal', # 15

    "tr|P94088|P94088_ARATH": 'Plant', # 16
    "tr|A0A3Q7GAS3|A0A3Q7GAS3_SOLLC": 'Plant', # 17
    "tr|W1P1B9|W1P1B9_AMBTC": 'Plant', # 18
    "tr|B4FE46|B4FE46_MAIZE": 'Plant', # 19
    "tr|Q84YS0|Q84YS0_ORYSJ": 'Plant', # 20
    "tr|A0A176W202|A0A176W202_MARPO": 'Plant', # 21
    "tr|A0A388KYX9|A0A388KYX9_CHABU": 'Plant', # 22
    "tr|A0A2R6XJM8|A0A2R6XJM8_MARPO": 'Plant', # 23
    "tr|A0A2K3CY04|A0A2K3CY04_CHLRE": 'Plant', # 24
    "tr|A8HP50|A8HP50_CHLRE": 'Plant', # 25
    "tr|A4SBD0|A4SBD0_OSTLU": 'Plant' # 26

}

# # 添加分類信息
# classification = {
#     'sp|Q9CYI4|LUC7L_MOUSE': 'Animal',
#     'tr|A0A3Q2UH51|A0A3Q2UH51_CHICK': 'Animal',
#     'tr|Q6GLC1|Q6GLC1_XENTR': 'Animal',
#     'tr|H3AZN1|H3AZN1_LATCH': 'Animal',
#     'tr|Q6IQD4|Q6IQD4_DANRE': 'Animal',
#     'sp|Q7TNC4|LC7L2_MOUSE': 'Animal',
#     'tr|Q5ZLW7|Q5ZLW7_CHICK': 'Animal',
#     'tr|Q28EN5|Q28EN5_XENTR': 'Animal',
#     'tr|H2ZXS0|H2ZXS0_LATCH': 'Animal',
#     'tr|A0A8M1NDA7|A0A8M1NDA7_DANRE': 'Animal',
#     'tr|A0A8C4NH22|A0A8C4NH22_EPTBU': 'Animal',
#     'tr|F6S0J6|F6S0J6_CIOIN': 'Animal',
#     'tr|A0A182H045|A0A182H045_AEDAL': 'Animal',
#     'tr|Q9VVI1|Q9VVI1_DROME': 'Animal',
#     'tr|B3RNA2|B3RNA2_TRIAD': 'Animal',
#     'tr|T2M333|T2M333_HYDVU': 'Animal',
#     'sp|Q09217|YP68_CAEEL': 'Animal',
#     'tr|A0A1X7V7G6|A0A1X7V7G6_AMPQE': 'Animal',
#     'tr|A0A3Q7GLL9|A0A3Q7GLL9_SOLLC': 'Plant',
#     'tr|U5D3Z2|U5D3Z2_AMBTC': 'Animal',
#     'tr|Q75LD6|Q75LD6_ORYSJ': 'Plant',
#     'tr|B4FUS0|B4FUS0_MAIZE': 'Plant',
#     'tr|B4FXT5|B4FXT5_MAIZE': 'Plant',
#     'tr|Q6K8R4|Q6K8R4_ORYSJ': 'Plant',
#     'tr|A0A2R6WMH7|A0A2R6WMH7_MARPO': 'Plant',
#     'tr|A0A388LKU4|A0A388LKU4_CHABU': 'Fungi',
#     'tr|Q7S615|Q7S615_NEUCR': 'Fungi',
#     'tr|Q5BCF4|Q5BCF4_EMENI': 'Fungi',
#     'tr|F9X0V8|F9X0V8_ZYMTI': 'Fungi',
#     'tr|Q5KBY8|Q5KBY8_CRYNJ': 'Fungi',
#     'tr|E3L016|E3L016_PUCGT': 'Fungi',
#     'tr|A0A139ATP5|A0A139ATP5_GONPJ': 'Fungi',
#     'tr|Q28G85|Q28G85_XENTR': 'Animal',
#     'tr|F1NXR8|F1NXR8_CHICK': 'Animal',
#     'sp|Q5SUF2|LC7L3_MOUSE': 'Animal',
#     'tr|F6Z8S6|F6Z8S6_CIOIN': 'Animal',
#     'tr|A0A1X7UMY2|A0A1X7UMY2_AMPQE': 'Animal',
#     'tr|Q1ED13|Q1ED13_DANRE': 'Animal',
#     'tr|H3ABC3|H3ABC3_LATCH': 'Animal',
#     'tr|B3RQ49|B3RQ49_TRIAD': 'Animal',
#     'tr|Q84YS0|Q84YS0_ORYSJ': 'Plant',
#     'tr|A0A2R6XJM8|A0A2R6XJM8_MARPO': 'Plant',
#     'tr|A0A3Q7GAS3|A0A3Q7GAS3_SOLLC': 'Plant',
#     'tr|A0A182GZD3|A0A182GZD3_AEDAL': 'Animal',
#     'tr|Q8ITY5|Q8ITY5_CAEEL': 'Animal',
#     'tr|B4FE46|B4FE46_MAIZE': 'Plant',
#     'tr|A0A176W202|A0A176W202_MARPO': 'Plant',
#     'tr|Q9W3X8|Q9W3X8_DROME': 'Animal',
#     'tr|H2L0S8|H2L0S8_CAEEL': 'Animal',
#     'tr|A0A388KYX9|A0A388KYX9_CHABU': 'Fungi',
#     'tr|W1P1B9|W1P1B9_AMBTC': 'Animal',
#     'tr|A8HP50|A8HP50_CHLRE': 'Plant',
#     'tr|A0A8C4N567|A0A8C4N567_EPTBU': 'Animal',
#     'tr|A0A8B6XHB2|A0A8B6XHB2_HYDVU': 'Animal',
#     'tr|A0A2K3CY04|A0A2K3CY04_CHLRE': 'Plant',
#     'tr|A4SBD0|A4SBD0_OSTLU': 'Plant'
# }


# 简化分类信息
simplified_classification = {k.split('|')[1]: v for k, v in classification.items()}

# L2和L3名单
simplified_L2_name = {name.split('|')[1] for name in L2_name}
simplified_L3_name = {name.split('|')[1] for name in L3_name}


# 创建空的DataFrame
all_data = pd.DataFrame()

# 颜色和标记设置
colors = {'L2': 'blue', 'L3': 'red', 'Other': 'grey'}
markers = {'Animal': 'o', 'Plant': 's', 'Fungi': '^', 'Unknown': 'x'}

# markers = {'Animal': 'o', 'Plant': 's', 'Fungi': '^', 'Unknown': 'x'}  # 为未知分类添加标记


# 颜色和标记设置
colors = {'L2': 'blue', 'L3': 'red', 'Other': 'LightGrey'}
markers = {'Animal': 'o', 'Plant': 'x', 'Fungi': '^', 'Unknown': '.'}

# 创建空的DataFrame
all_data = pd.DataFrame()

# 读取并处理每个文件
for file, label in zip(files, labels):
    df = pd.read_csv(file)
    df['Name'] = df['Name']
    df['x'] = df['L2 Score'] / df['Length']
    df['y'] = df['L3 Score'] / df['Length']
    df['Label'] = label

    plt.figure(figsize=(10, 6))
    handles, labels = [], []

    # 先绘制Unknown
    unknown_df = df[df['Name'].apply(lambda x: simplified_classification.get(x, 'Unknown') == 'Unknown')]
    known_df = df[df['Name'].apply(lambda x: simplified_classification.get(x, 'Unknown') != 'Unknown')]

    for df_subset, is_unknown in [(unknown_df, True), (known_df, False)]:
        for index, row in df_subset.iterrows():
            if row['Name'] in simplified_L2_name:
                group = 'L2'
                color = colors['L2']
            elif row['Name'] in simplified_L3_name:
                group = 'L3'
                color = colors['L3']
            else:
                group = 'Other'
                color = colors['Other']

            category = simplified_classification.get(row['Name'], 'Unknown')
            marker = markers[category]
            plot_label = f'{group} {category}' if is_unknown else f'{group} {category}'

            # 绘制并更新图例
            sc = plt.scatter(row['x'], row['y'], color=color, marker=marker, label=plot_label)
            if plot_label not in labels:
                handles.append(sc)
                labels.append(plot_label)

    # 添加对角虚线
    plt.plot([0, 2], [0, 2], 'k--', label='Diagonal line')

    plt.xlabel('L2 Score / Length')
    plt.ylabel('L3 Score / Length')
    plt.xlim(0, 2.0)
    plt.ylim(0, 2.0)
    plt.title(f'{label} Analysis')

    # 排序图例，鲜红色（L3）在蓝色（L2）之前
    handles, labels = zip(*sorted(zip(handles, labels), key=lambda t: t[1], reverse=True))
    plt.legend(handles, labels)

    plt.savefig(f'../../result/reproduce_experiment_for_iteration1/i1_csv.png')






############################################################################################
############################################################################################
############################################################################################
# # 定義CSV文件路徑
# files = [
#     '../../result/reproduce_experiment_for_iteration1/i1_interpro_L7.csv'
#     # './i1_csv/i1_8872_interpro.csv',
#     # './i1_csv/i1_all_L7_10000_decoy.csv',
#     # './i1_csv/i1_L2_10000_decoy.csv',
#     # './i1_csv/i1_L3_10000_decoy.csv'
# ]

# # 顏色列表，為每個文件指定不同顏色
# # colors = ['red', 'green', 'blue', 'purple']
# # labels = ['i1_interpro', 'i1_all_L7_decoy', 'i1_L2_decoy', 'i1_L3_decoy']
# colors = ['blue']
# labels = ['i1_interpro']

# # 創建空的DataFrame
# all_data = pd.DataFrame()

# # 讀取每個文件
# for file, color in zip(files, colors):
#     df = pd.read_csv(file)
#     df['x'] = df['L2 score'] / df['Length']  # 計算L2 score/Length
#     df['y'] = df['L3 score'] / df['Length']  # 計算L3 score/Length
#     df['Color'] = color  # 為數據集添加顏色
#     all_data = pd.concat([all_data, df])

# # 繪製圖表
# i=0
# plt.figure(figsize=(10, 6))
# for color in colors:
#     subset = all_data[all_data['Color'] == color]
#     plt.scatter(subset['x'], subset['y'], color=color, label=f'{labels[i]}')
#     i+=1

# # 添加對角虛線
# plt.plot([0, 2], [0, 2], 'k--', label='Diagonal line')  # 假設x和y的比例範圍最大為2

# plt.xlabel('L2 Score / Length')
# plt.ylabel('L3 Score / Length')
# plt.xlim(0, 2.0)  # 設置x軸範圍最大為2.0
# plt.ylim(0, 2.0)  # 設置y軸範圍最大為2.0
# plt.title('Iterative 1 Normalized Scores Plot')
# plt.legend()
# plt.savefig('../../result/reproduce_experiment_for_iteration1/i1_csv.png')





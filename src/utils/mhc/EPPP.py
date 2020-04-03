#!/usr/bin/env python3

import numpy as np, pandas as pd
import sys, os, re, getopt, pickle
import theano
import keras
from keras.models import load_model
from collections import defaultdict, OrderedDict

AA_dict = {'A': 1,
           'C': 2,
           'D': 3,
           'E': 4,
           'F': 5,
           'G': 6,
           'H': 7,
           'I': 8,
           'K': 9,
           'L': 10,
           'M': 11,
           'N': 12,
           'P': 13,
           'Q': 14,
           'R': 15,
           'S': 16,
           'T': 17,
           'V': 18,
           'W': 19,
           'Y': 20,
           'Z': 0}  # 'Z' is a padding variable

len_AA_dict, len_peptide, hla_per_sample = len(AA_dict), 11, 6
HLAs = pickle.load(open('/home/danxu/python_nbs/Dumps/MS_HLAs.list', 'rb'))
model_path = '/home/danxu/python_nbs/data/DL_NeoAntigen/ms_model/models/ms_model_sigmoid.BNL1.S0.1_1epoch.h5'
Pr_thld = 0.1 # 新生抗原概率阈值
Len_thld = 7 # 新抗原肽段长度阈值

def category_encode(data, categories):
    if isinstance(data, str):
        data = [data]
    if isinstance(data, np.ndarray):
        data = data.tolist()
    encoded = []
    for datum in data:
        if datum not in categories:
            #raise ValueError('Category not found!: %s' % datum)
            encoded.append(0)
        else:
            encoded.append(categories.index(datum) + 1)
    return encoded

def hla_encode(alleles, hla_per_sample=hla_per_sample, HLAs=HLAs):
    if isinstance(alleles, np.ndarray):
        alleles = alleles.tolist()
    type_check = [isinstance(sample, list) for sample in alleles]
    if any(type_check):
        assert all(type_check), \
            'Must provide either a list of alleles or a list of allele lists!'
    else:
        alleles = [alleles]
    onehots = []
    for sample in alleles:
        onehot = category_encode(sample, HLAs)
        onehot = [0] * (hla_per_sample - len(onehot)) + onehot
        onehots.append(onehot)
    return np.array(onehots)

def peptide_encode(peptides, maxlen=None, AA_dict=AA_dict):
    if isinstance(peptides, str):
        peptides = [peptides]
    num_peptides = len(peptides)
    if maxlen is None:
        maxlen = max(map(len, peptides))
    onehot = np.zeros((num_peptides, maxlen), dtype=np.uint8)
    for i, peptide in enumerate(peptides):
        if len(peptide) > maxlen:
            msg = 'Peptide %s has length %d > maxlen = %d.'
            raise ValueError(msg % (peptide, len(peptide), maxlen))
        o = list(map(lambda x: AA_dict[x], peptide))
        k = len(o)
        o = o[:k // 2] + [0] * (maxlen - k) + o[k // 2:]
        if len(o) != maxlen:
            msg = 'Peptide %s has length %d < maxlen = %d, but pad is "none".'
            raise ValueError(msg % (peptide, len(peptide), maxlen))
        onehot[i, :] = o
    return np.array(onehot)

def supplement_hla2six(hla_list):
    if(len(hla_list) != 6):
        new_arr = []
        for i in ['HLA-A', 'HLA-B', 'HLA-C']:
            tmp = [hla for hla in hla_list if i in hla]
            if(len(tmp) == 1):
                new_arr.extend(tmp); new_arr.extend(tmp)
            elif(len(tmp) == 2):
                new_arr.extend(tmp)
            elif(len(tmp) == 0):
                pass
            else:
                raise Exception('hla_list error:', hla_list)
        return new_arr
    else:
        return hla_list

def cut_3Lpeps(P1, P2):
    t1, t2 = [], []
    midx = 0
    if(len(P1) == len(P2)):
        for i in range(len(P1)):
            if(P1[i] != P2[i]): midx = i; break
        for L in [9,10,11]:
            t1.extend([P1[midx-i:midx+L-i] for i in range(L)])
            t2.extend([P2[midx-i:midx+L-i] for i in range(L)])
    else:
        for i in range(min(len(P1), len(P2))):
            if(P1[i] != P2[i]): midx = i; break
        for L in [9,10,11]:
            t1.extend([P1[midx-i:midx+L-i] for i in range(L)])
            t2.extend([P2[midx-i:midx+L-i] for i in range(L)])
            # frameshift extension peptides
            #t1.extend([P1[midx+1+i:midx+1+i+L] for i in range(len(P1)-midx-L)])
            #t2.extend(['' for i in range(len(P1)-midx-L)])
    return(t1, t2)

def trans_arr_col_row(in_arr):
    in_arr = np.array(in_arr); col_num = in_arr.shape[1]
    return([in_arr[:,i] for i in range(col_num)])
        
def usage():
    print("""
    Epitopes presentation probability prediction
    
    Usage: python3 EPPP.py -i <hla1_file>
                           -m <mutationtype file>
                           -w <wildtype file>
                           -p <peptideID file>
                           -s <output prefix sampleName>
                           -o <output path>
                           """)

if len(sys.argv) == 1:
    usage()
    sys.exit(1)
    
config = {
    "hla1_file" : "",
    "mt_file" : "",
    "wt_file" : "",
    "pid_file": "",
    "prefix":"",
    "outpath" : ""
}

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:p:m:w:s:o:")
except getopt.error as v:
    print(v)
    sys.exit(1)
    
for option, value in opts:
    if option == "-m":
        config["mt_file"] = value
    elif option == "-w":
        config["wt_file"] = value
    elif option == "-p":
        config["pid_file"] = value
    elif option == "-i":
        config["hla1_file"] = value
    elif option == "-s":
        config["prefix"] = value
    elif option == "-o":
        config["outpath"] = os.path.abspath(value)

with open(config["mt_file"], 'r') as f:
    C = f.readlines()
    pepid_list = [l.strip().replace('>','') for n,l in enumerate(C) if n%2 == 0]
    mt_list = [l.strip().replace('*','') for n,l in enumerate(C) if n%2 == 1]
with open(config["wt_file"], 'r') as f:
    wt_list = [l.strip().replace('*','') for n,l in enumerate(f.readlines()) if n%2 == 1]

mt_peps, wt_peps = [], []
for P1, P2 in zip(mt_list, wt_list):
    t1, t2 = cut_3Lpeps(P1, P2)
    mt_peps.append(t1); wt_peps.append(t2)
        
with open(config["hla1_file"], 'r') as f:
    hla_onehot = f.readline().strip(); hla_onehot = hla_onehot.split(',')
    hla_onehot = [hla[:5] + '*' + hla[5:] for hla in hla_onehot if '*' not in hla] #convert into model readable format
    hla_onehot = supplement_hla2six(hla_onehot)
    
with open(config["pid_file"], 'r') as f:
    dic_pid = {}
    for l in f:
        l = l.strip().split('\t')
        #Gene Chr Pos  Transcript cDNA Function ProteinChange ProteinChangeAbbr ExtendedPeptide AAchangePosition
        dic_pid['Peptide' + l[0]] = [l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10]]

ms_model = load_model(model_path, custom_objects={'len_AA_dict':len_AA_dict})
header = ['PeptidePostion','HLA','mt_peptides','PeptideID','mt_probs','NeoepitopeType','wt_probs',
          'wt_peptides','Gene','Chr','Position','Transcript','CDSMutation','Function','ProteinChange',
          'ProteinChangeAbbr','ExtendedPeptide','AAchangePosition']
scols = ['PeptidePostion','HLA','mt_peptides','PeptideID','mt_probs','wt_probs','wt_peptides']

df_final_all, df_final = pd.DataFrame(), pd.DataFrame()
for P1_list, P2_list, pid in zip(mt_peps, wt_peps, pepid_list):
    Len = len(P1_list)
    hla_list = [hla_onehot for i in range(Len)]; pepid = [pid for i in range(Len)]
    gene, Chr, Pos, Trans, cDNA, Func, ProteinChange, ProteinChangeAbbr, ExPep, AAPos = trans_arr_col_row([dic_pid[i] for i in pepid])
    mt_input = {'peptide': peptide_encode(P1_list, maxlen=len_peptide),
                'hla_onehot': hla_encode(hla_list)
                }
    wt_input = {'peptide': peptide_encode(P2_list, maxlen=len_peptide),
                'hla_onehot': hla_encode(hla_list)
                }
    mt_probs = ms_model.predict(mt_input); mt_probs = mt_probs[:,0]
    wt_probs = ms_model.predict(wt_input); wt_probs = wt_probs[:,0]
    pepp = [ExPep[i].index(P1_list[i])+1 for i in range(Len)]
    df = pd.DataFrame({'PeptidePostion':pepp, 'HLA':[','.join(i) for i in hla_list], 'mt_peptides':P1_list, 'PeptideID':pepid,
                      'mt_probs':mt_probs, 'NeoepitopeType':['HP' for i in range(Len)], 'wt_probs':wt_probs, 'wt_peptides':P2_list,
                      'Gene':gene, 'Chr':Chr, 'Position':Pos, 'Transcript':Trans, 'CDSMutation':cDNA, 'Function':Func,
                      'ProteinChange':ProteinChange, 'ProteinChangeAbbr':ProteinChangeAbbr, 'ExtendedPeptide':ExPep, 'AAchangePosition':AAPos})
    #if((df.wt_probs > 0.1).mean() != 0): continue
    df.loc[df.wt_peptides == '', 'wt_probs'] = 1e-10
    df_final_all = pd.concat([df_final_all, df[scols][((df.mt_peptides.apply(lambda x:len(x))>Len_thld) & (df.mt_probs>Pr_thld))]])
    df = df[(df.mt_peptides.apply(lambda x:len(x))>Len_thld)]
    df_final = pd.concat([df_final, df.sort_values('mt_probs', ascending=False).head(3)])

prefix = config["prefix"]
#阴性结果
if 'df' not in vars():
    with open(f'{config["outpath"]}/{prefix}.msP.txt', 'w') as fo:
        fo.write('\t'.join(header) + '\n')
    with open(f'{config["outpath"]}/{prefix}.msP.all.txt', 'w') as fo:
        fo.write('')
    print('No mutation Found.'); sys.exit(0)

#去重
df_final_all.drop_duplicates(subset=['mt_peptides','mt_probs'], inplace=True)
df_final.drop_duplicates(subset=['mt_peptides','mt_probs'], inplace=True)
                    
#排序
dic_rank = df_final.groupby(['Chr','Position'], sort=False)['mt_probs'].max().rank(ascending=False).to_dict()
df_final['NeoRank'] = [dic_rank[(c,p)] for c,p in zip(df_final.Chr, df_final.Position)]
df_final = df_final.sort_values(['NeoRank', 'mt_probs'], ascending=[1,0]).drop('NeoRank', axis=1)

#Write Results
df_final_all['mt_probs'] = 1/df_final_all['mt_probs']; df_final_all['wt_probs'] = 1/df_final_all['wt_probs']
df_final_all.insert(5, 'NeoepitopeType', 'HP')
df_final_all.to_csv(f'{config["outpath"]}/{prefix}.msP.all.txt', index=None, header=False, sep='\t')

#df_final['mt_probs'] = 1/df_final['mt_probs']; df_final['wt_probs'] = 1/df_final['wt_probs']
df_final.to_csv(f'{config["outpath"]}/{prefix}.msP.txt', index=None, sep='\t')

                
                
'''
     AUTHOR: D.X, danxu@yucebio.com
    VERSION: 1.0
    CREATED: 2019-6-4
  COPYRIGHT: Copyright belongs to Yucebio company and author, any use of this script must ask for authority from owner,
            otherwise the company and author have right to ban the usage.
'''

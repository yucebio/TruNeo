#!/usr/bin/env python3
import pandas as pd
import sys,os


hlaloh_file = sys.argv[1]
hla_file = sys.argv[2]
output_summary = sys.argv[3]

hla = pd.read_csv(hla_file,sep="\t",dtype=str,header=None)
hlaloh = pd.read_csv(hlaloh_file,sep="\t",dtype=str)

def hlaInDup(cnv,hla_a_b_c):
    hla_dup = False
    if cnv.empty:
        return hla_dup

    for index, row in cnv.iterrows():
        if row['cnv'] != "DUP":
            continue
        if row['chr'] != hla_a_b_c[0]:
            continue
        if int(row['begin']) < hla_a_b_c[2] and int(row['end']) > hla_a_b_c[1]:
            hla_dup = True
    return hla_dup

def AlleleImbalance(hlaloh):
    hlaloh_list = []
    if hlaloh.empty:
        return hlaloh_list
    for index, row in hlaloh.iterrows():
        #if float(row['PVal_unique']) < 0.0002 and (float(row['HLA_type2copyNum_withBAFBin']) < 0.5 or float(row['HLA_type1copyNum_withBAFBin']) < 0.5):
        if float(row['PVal_unique']) < 0.0002:
            hlaloh_list.append(row['LossAllele'])
    return hlaloh_list

hla[1]='阴性'
hla_a_b_c=['chr6',29908247,31241913]
hla.columns=["血液样本HLA分型","肿瘤样本HLA缺失检测结果"]

if len(sys.argv) == 5:
    cns_file = sys.argv[4]
    cnv = pd.read_csv(cns_file,sep="\t",dtype=str)
    if hlaInDup(cnv,hla_a_b_c):
        hla.to_csv(output_summary,sep="\t",index=False)
        exit(0)

hlaloh_list=AlleleImbalance(hlaloh)
hla['肿瘤样本HLA缺失检测结果'] = hla['血液样本HLA分型'].apply(lambda x: '阳性' if x in hlaloh_list else '阴性')
hla.to_csv(output_summary,sep="\t",index=False)



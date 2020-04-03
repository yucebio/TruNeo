#!/usr/bin/python

import pandas as pd
import vcf
import re
import os
import sys
import time
import coloredlogs
import logging
import argparse
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')

class variant:
    def __init__(self,chromosome, position, ref_str, alt_str, maf=0):
        self.chr = chromosome
        self.pos = position
        self.ref = ref_str
        self.alt = alt_str
        self.maf = maf
        
    def same(self,v):
        # if the same variant
        # input: another variant
        if v.chr != self.chr:
            return False
        if v.pos != self.pos and (v.pos+1)!= self.pos and (v.pos-1)!= self.pos:
            return False
        return True

def vcf2variant_list(vcf_file):
    vcf_reader = vcf.Reader(filename=vcf_file)
    variant_list = []
    for record in vcf_reader:
        af = record.samples[0]['AF']
        my_variant = variant(record.CHROM, record.POS, record.REF, record.ALT, af)
        variant_list.append(my_variant)
    return variant_list

def get_af(record, vlist):
    # a dataframe record get af from variant list
    my_variant = variant(record['chr'],record['allelepos'],record['ref'],record['alt'])
    for i in vlist:
        if my_variant.same(i):
            my_variant.maf = p2f(i.maf)
    return my_variant.maf

def p2f(x):
    if isinstance(x,str):
        if "%" in x:
            y = round(float(x.strip('%'))/100,4)
    return y

def main():
    program_license = '''
    
    ~~~~~~~~~~~~~~

    :Description: This module will filter neoantigens by RNA expression and transform to yuce format

    ::Inputs::
    vcf: variant file to extract VAF information
    neo: neoantigen result from NeoPredPipe

    ::Output::
    output: yuce format neoantigen report file
    '''

    parser = ArgumentParser(
        description = program_license,
        formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v","--input_vcf",dest="input_vcf", type=str, required=True,help="Full path to the input vcf file")
    parser.add_argument("-o","--output_neo",dest="output_neo", type=str, required=True,help="Full path to the output neoantigen file")
    parser.add_argument("-i","--input_neo",dest="input_neo", type=str, required=True, help="Full path to the inputput neoantigen file, from NeoPredPipe")
    parser.add_argument("-e","--expression",dest="expression", type=float, required=False, default=1.0, help="threthold for expression (defaul 1.0)")

    args = parser.parse_args()

    df=pd.read_csv(args.input_neo,sep="\t",header=None)
    if len(df.columns) == 24:
        name=["Sample","Line","chr","allelepos","ref","alt","GeneName","Expression","pos","HLA","Peptide","core","Of","Gp","Gl","Ip","Il","Icore","Identity","Score","BA","Rank","Candidate","BindLevel"]
        df.columns=name
        df = df[df['Expression']>args.expression]
    else:
        name=["Sample","Line","chr","allelepos","ref","alt","GeneName","pos","HLA","Peptide","core","Of","Gp","Gl","Ip","Il","Icore","Identity","Score","BA","Rank","Candidate","BindLevel"]
        df.columns=name

    vcf_file = args.input_vcf

    v_list = vcf2variant_list(vcf_file)
    df['af']=df.apply(lambda x:get_af(x,v_list),axis=1)
    df = df.sort_values(by=['BindLevel','af'],ascending=[True,False])
    df['NeoRank'] = range(1,len(df)+1)
    df2=pd.DataFrame()
    df2[['Peptide','NeoRank','HLA', 'NeoepitopeType','Gene','Chr','Position','Psite','DNAVAF']] = df[['Peptide','NeoRank','HLA','BindLevel','GeneName','chr','allelepos','pos','af']]
    #df2['Score'] = round(df['Score'],4)
    df2['NeoepitopeScore'] =df2['Score'] = df['Score']
    df2['WildtypeScore']=df2['WildtypePeptide']=df2['Transcript']=df2['CDSMutation']=df2['ProteinChangeAbbr']=df2['CCLEexpStatus']=df2['CCF']=df2['Clonal']="NA"
    df2=df2[['NeoRank','Score','HLA','Peptide','NeoepitopeScore', 'NeoepitopeType','WildtypeScore','WildtypePeptide','Gene','Chr','Position','Transcript','CDSMutation','ProteinChangeAbbr','Psite','DNAVAF','CCLEexpStatus','CCF','Clonal']]
    df2.to_csv(args.output_neo,index=False,sep="\t")

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logger.info("vardict_filter.py: Elapsed time was %g seconds", totaltime)



import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def Process_HBonds_File(filename):
    df = pd.read_csv(filename,delim_whitespace=True)
    df = df.drop(columns=["AvgDist","Frames","AvgAng"])
    df["#Acceptor"] = df["#Acceptor"].str.split("_").str[1]
    df["#Acceptor"] = df["#Acceptor"].str.split("@").str[0]
    df["DonorH"] = df["DonorH"].str.split("_").str[1]
    df["DonorH"] = df["DonorH"].str.split("@").str[0]
    df["Donor"] = df["Donor"].str.split("_").str[1]
    df["Donor"] = df["Donor"].str.split("@").str[0]
    df["#Acceptor"] = pd.to_numeric(df["#Acceptor"])
    df["Donor"] = pd.to_numeric(df["Donor"])
    df["DonorH"] = pd.to_numeric(df["DonorH"])
    return df

def Clear_Non_Interface_HBonds(df,res_num):
    index_names = df[ (df['#Acceptor'] >= res_num) & (df['Donor'] >= res_num)].index 
    df.drop(index_names, inplace = True) 
    index_names = df[ (df['#Acceptor'] < res_num) & (df['Donor'] < res_num)].index 
    df.drop(index_names, inplace = True) 
    return df

def Group_Multiple_HBonds(df):
    df = df.groupby(["#Acceptor","Donor"])["Frac"].apply(list).reset_index()
    return df

def Mod_ResNums(df,res_num):
    df["Mod_Acceptor"] = df["#Acceptor"]
    df["Mod_Donor"] = df["Donor"]
    df.loc[df["Donor"] >= res_num, "Mod_Donor"] = df["Donor"] - res_num + 5
    df.loc[df["#Acceptor"] >= res_num, "Mod_Acceptor"] = df["#Acceptor"] - res_num + 5
    df.loc[df["Donor"] < res_num, "Mod_Donor"] = df["Donor"] + 4
    df.loc[df["#Acceptor"] < res_num, "Mod_Acceptor"] = df["#Acceptor"] + 4
    return df

def Compare_HBond_DF(excel,sheetlabel,ref,sub,res_num):
    compare = pd.merge(ref,sub,how='outer',left_on=["#Acceptor","Donor"],right_on=["#Acceptor","Donor"])
    compare = compare.fillna(0)
    compare = compare.rename(columns={"Frac_x": "Reference", "Frac_y": "Subject"})
    compare = compare.drop(columns=["Mod_Acceptor_x","Mod_Donor_x","Mod_Acceptor_y","Mod_Donor_y"])
    compare = Mod_ResNums(compare,res_num)
    compare.to_excel(excel,sheet_name = sheetlabel,index=False)


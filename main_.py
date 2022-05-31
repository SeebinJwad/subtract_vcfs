import pandas as pd
import numpy as np
import os
from functools import partial
import jdc

class Mech():
    def __init__(self, chromosomes, sample_id, folder):
        usecols = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Otherinfo10", "Otherinfo11"]
        dtype = {"Chr":"category",
          "Start":np.int32,
          "End":np.int32,
          "Ref":"category",
          "Alt":"category",
          "Func.refGene":"category",
          "Otherinfo10":"category",
          "Otherinfo11":"string"}
        self.chromosomes = chromosomes
        self.sample_id = sample_id
        self.folder = folder
        
        # IF YOU ARE COMMA DELIMITED REMOVE sep='\t' THAT IS FOR TAB DELIMITED 
        map_concat = partial(pd.read_csv, dtype=dtype, usecols=usecols, sep='\t')
        
        # INPUT FOLDER PATH BELOW AND REPLACE \ WITH \\
        # "C:\\Users\\kazij\\Documents\\annovar-full"
        file_list = [f"{self.folder}\\chr{chrom}_{self.sample_id}.csv" for chrom in self.chromosomes]
        self.df = pd.concat(map(map_concat, file_list), ignore_index=True)
        col_names = self.df.columns.values.tolist()
        new_col_names = [''.join(char for char in name if char.isalnum()) for name in col_names]
        self.df.rename(columns=dict(zip(col_names, new_col_names)))
        self.df = self.df.astype(dtype)

    def category_filter(self, **kwargs):
        for key, value in kwargs.items():
            if key == "FuncrefGene":
                key = "Func.refGene"
            cats = self.df[key].cat.categories
            inv_cats = [cat for cat in cats if not cat in value]
            for cat in inv_cats:
                self.df = self.df[self.df[key] != cat]
    def num_filter(self, **kwargs):
        pd.eval("DP = self.df.Otherinfo11.str.split(';',3).str[2].str.split('DP=').str[1].astype('int32')", target=self.df, inplace=True)
        pd.eval("AF = self.df.Otherinfo11.str.split(';',5).str[4].str.split('AF=').str[1].astype('float')", target=self.df, inplace=True)
        for key, value in kwargs.items():
            self.df.loc[:,key].where(self.df.loc[:,key] > value, inplace=True)
        self.df.dropna(inplace=True)

def find_mid(indices):
    copy = indices[:-1].copy()
    copy = np.insert(copy, 0, indices[0] - 1)
    sub = np.subtract(indices, copy)
    # [1,2,3,1,2,3] - [0,1,2,3,1,2]
    return sub.argmin()

def inv_arr(max_val, idx_arr):
    arr = set(range(0,max_val))
    diff = arr.difference(idx_arr)
    return list(diff)

class Uqid():
    def __init__(self, series, len1, len2):
        self.series = series
        self.len1 = len1
        self.len2 = len2
        
    def create_idx_col(self):
        self.idxs = np.append(np.arange(0, self.len1), np.arange(0, self.len2))
        self.df = pd.DataFrame({'idx':self.idxs}, index=self.series)
        self.df.reset_index(inplace=True)
        
    def filter_uqs(self):
        self.df = self.df.sort_index()
        self.df = self.df.groupby('index')
        self.df = self.df.filter(lambda x: len(x) == 1)

def combine_uqids(mech1_,mech2_):
    uqid1 = mech1_.df[mech1_.df.columns[1:5]].apply(lambda x: ''.join(x.astype(str)),axis=1).squeeze()
    uqid2 = mech2_.df[mech2_.df.columns[1:5]].apply(lambda x: ''.join(x.astype(str)),axis=1).squeeze()    
    combined_uqids = pd.concat([uqid1, uqid2], ignore_index=True)
    return Uqid(combined_uqids, len(uqid1.index), len(uqid2.index))

def subtract(mech1, mech2):
    combine_mechs = combine_uqids(mech1,mech2)
    combine_mechs.create_idx_col()
    combine_mechs.filter_uqs()
    idx = combine_mechs.df['idx'].to_numpy()
    mid_val = find_mid(idx)
    drop_vals = inv_arr(combine_mechs.len1, idx[:mid_val])
    uq_df = mech1.df.drop(mech1.df.index[[drop_vals]])
    return uq_df

def chrom_range(up_to_chrom_num, include_x=False, include_y=False):
    range_ = [*range(1, up_to_chrom_num+1)]
    range_ = list(map(str, range_))
    if include_x:
        range_ += ['X']
    if include_y:
        range_ += ['Y']
    return range_

def main():
   # CHANGE THESE VALUES BELOW
    ids = {'22241':['22334']}
    path = "C:\\Users\\kazij\\Documents\\annovar-full"
    chromosomes = chrom_range(up_to_chrom_num=1, include_x=False, include_y=False)
    min_cat_params = dict(Otherinfo10=['PASS'])
    min_num_params = dict(DP=15, AF=.15)
    sub_cat_params = dict()
    sub_num_params = dict()
    # CHANGE THESE VALUES ABOVE
    
    for key, value in ids.items():
        min_parent = Mech(chromosomes, key, path)
        min_parent.category_filter(**min_cat_params)
        min_parent.num_filter(**min_num_params)
        sub_parent = Mech(chromosomes, key, path)
        sub_parent.category_filter(**sub_cat_params)
        sub_parent.num_filter(**sub_num_params)
        for val in value:
            min_variant = Mech(chromosomes, val, path)
            min_variant.category_filter(**min_cat_params)
            min_variant.num_filter(**min_num_params)
            sub_variant = Mech(chromosomes, val, path)
            sub_variant.category_filter(**sub_cat_params)
            sub_variant.num_filter(**sub_num_params)
            
            min_parent_sub_variant = subtract(min_parent, sub_variant)
            min_variant_sub_parent = subtract(min_variant, sub_parent)

            if not os.path.exists("results"):
                os.makedirs("results")
            min_parent_sub_variant.to_csv(f"results\\{key}-{val}.csv", index=False)
            min_variant_sub_parent.to_csv(f"results\\{val}-{key}.csv", index=False)

main()

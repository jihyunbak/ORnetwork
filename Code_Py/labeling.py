# labeling.py
# perceptual labeling of receptor-code groups
# python code for demonstration

# Copyright 2019 Ji Hyun Bak


# === initialize

import numpy as np
import pandas as pd

# JHB custom
from Code_Py import aux, bim
Tray = aux.Tray

# /

# === data file loaders

# --- set directories

indir = 'Data/database/' # web-scraped percept data (GoodScents)
passdir = 'Data/database/' # intermediate data being passed
corrdir = 'Data/database/' # cross-list correspondence (Mainland2015 and GoodScents)

# --- functions

def load_mols_with_descriptors_goodscent(add_fixup=False, filter_corr=False):
    ''' load the list of molecules with odor descriptors '''
    # mfilename = indir + 'goodscents_out_merged.csv'
    mfilename = indir + 'odorants_with_descriptors.csv'
    mtab = aux.csv_to_df(mfilename)
    # optionally, add fix-up table
    if add_fixup:
        fixtab = fixup_mols_with_descriptors(mtab.columns, mtab.shape[0])
        mtab = mtab.append(fixtab, ignore_index=True)
    # then optionally, filter odorants in Mainland2015 dataset (keep Ltab order)
    if filter_corr:
        mtab_full = mtab
        Ltab = load_odorants_Mainland15_with_goodscent_corr()
        mtab = aux.collect_df_rows_by_index(mtab_full, Ltab.mtab_idx)
    return mtab

def fixup_mols_with_descriptors(columns, ibase):
    ''' 
    add odor information for test molecules
    that are not covered by GoodScents CAS-numbered list 
    '''
    # add information (TODO: should replace by a text file)
    addlist = [\
     ['benzene', 'aromatic', 'aromatic;sweet;gasoline-like'],
     ['5,5-Dimethyl-1,3-cyclohexanedione', '', 'odorless;weak'],
     ['nonanedioic acid', '', 'fatty'],
     ['thioglycolic acid', '', 'unpleasant;pungent;rotten'],
     ['1-octanethiol', 'sulfurous', 'sulfurous'],
     ['TMT', '', 'fox-odor'],
     ['(+)-menthol', 'mentholic', 'mentholic;cooling;minty'],
     ['androstenone', 'sweaty', 'sweaty;urinous;woody;floral'],
     ['banana', 'fruity', 'banana;fruity;creamy;tropical'],
     ['androstadienone', '', 'sweaty'],
     ['1-formylpiperidine', '', 'odorless'],
     ['3-methyl-2-hexenoic acid', '', 'sweaty'],
     ['butyric acid', '', 'unpleasant;rancid;penetrating;obnoxious']]
    
    # match format to mtab
    def match_format_to_mtab(fixlist, ibase, columns):
        fixlist_ext = []
        for i in range(len(fixlist)):
            nr = ibase + i
            front_matter = [nr, 0, 0, 'N/A', '']
            fixlist_ext.append(front_matter + fixlist[i])
        fix_df = pd.DataFrame(fixlist_ext, columns=columns)
        return fix_df
    
    # "fixed" dataset with added information
    fix_df = match_format_to_mtab(addlist, ibase, columns)
    return fix_df

def load_mols_with_descriptor_indices_goodscent(add_fixup=False, filter_corr=False):
    # load nonzero indices
    fixtag = 'fixup_' if add_fixup else ''
    # ifilename = passdir + 'goodscents_' + fixtag + 'out_merged2_basis.csv'
    ifilename = passdir + 'odorants_with_descriptor_indices.csv'
    ilist, iheader = aux.mycsvread(ifilename, 1)
    itab = pd.DataFrame(data=ilist, columns=iheader[0])
    # optionally, filter odorants in Mainland2015 dataset (keep Ltab order)
    if filter_corr:
        Ltab = load_odorants_Mainland15_with_goodscent_corr()
        itab = aux.collect_df_rows_by_index(itab, Ltab.mtab_idx)
        ilist = [ilist[midx] for midx in Ltab.mtab_idx.astype('int')]
    return itab, ilist

def load_descriptor_basis_goodscent(add_fixup=False):
    # load basis words
    fixtag = 'fixup_' if add_fixup else ''
    # bfilename = passdir + 'goodscents_' + fixtag + 'descriptor_basis.tsv' # 5/29/2019 updated
    bfilename = passdir + 'descriptor_basis.tsv'
    bdat = aux.mycsvread(bfilename, delimiter='\t')
    btab = pd.DataFrame(data=bdat, columns=['dim', 'word'])
    return btab

def load_odorants_Mainland15_with_goodscent_corr():
    # load Mainland2015 odorant list (with goodScents correspondence)
    # lfilename = corrdir + 'Ltab_corr3.csv' # with out-of-GoodScents input
    lfilename = corrdir + 'odorants_corr.csv' # with out-of-GoodScents input
    Ltab = aux.csv_to_df(lfilename)
    return Ltab

def load_descriptors_with_Castro13_category():
    # load manually annotated Castro2013 grouping
    # filename = passdir + 'descriptor_with_rcode_c13grp_manual.csv'
    filename = passdir + 'descriptor_proxy_category_map.csv'
    ctab = aux.csv_to_df(filename, delimiter=',')
    return ctab

# --- wrapper ---

def load_descriptor_space(add_fixup=True, filter_corr=False, dataset='goodscents'):
    if dataset == 'goodscents':
        mtab = load_mols_with_descriptors_goodscent(add_fixup, filter_corr) # mols with descriptors
        itab, ilist = load_mols_with_descriptor_indices_goodscent(add_fixup, filter_corr) # nonzero idx
        btab = load_descriptor_basis_goodscent(add_fixup) # basis words
    # else:
        # todo
    return mtab, itab, ilist, btab

# /


# === perceptual odor groups

def split_collect_idx(mystr, delimiter=';'):
    idxlist = []
    if (mystr is not None) and (len(mystr) > 0):
        nb_words_list = mystr.split(delimiter)
        for word in nb_words_list:
            idx_word = int(word)
            idxlist.append(idx_word)
    return idxlist

def unpack_itab_to_type_all(itab):
    ilist2 = []
    for j in range(len(itab)):
        il_type = split_collect_idx(itab.iloc[j].OdorTypeN)
        il_words = split_collect_idx(itab.iloc[j].OdorDescriptorN)
        il_all = list(set(il_type + il_words)) # union
        mylist = [il_type, il_all]
        ilist2.append(mylist)
    return ilist2

def assign_wgroup_to_odorants(itab, ctab, Klist=None):
    ''' assign to each odorant a perceptual odor category index '''
    
    # unpack index data: [[typeN], [typeN and descriptorN merged]]
    ilist2 = unpack_itab_to_type_all(itab)
    
    # unpack category indices
    u = ctab.widx.astype('int') # index to Wtab
    cgrp = ctab.c13grp_manual.astype('int') # manually assigned Castro2013 groups
    
    # detect Klist if not provided
    if not Klist:
        Klist_full = np.sort(np.unique(cgrp)).tolist()
        Klist = [k for k in Klist_full if k>0]
    
    # make a list of descriptor categories
    cglist2 = list_descriptor_categories(ilist2, cgrp, u)
    
    # fractional assignment
    odorant_cgrp = fractional_assignment(cglist2, Klist)
    
    # add more fields to output
    out = odorant_cgrp
    out.list = cglist2
    out.columns = Klist
    out.rows = ilist2
    return out

def list_descriptor_categories(ilist2, cgrp, u):
    ''' collect the odor category labels for the descriptors'''
    
    def collect_list(sublist):
        mycg = []
        for wi in sublist:
            wi1 = (wi + 1) # make 1-based indices
            v2 = aux.findby(cgrp, u == wi1)
            mycg = mycg + v2
            if not v2: # indicates error
                print(wi1) 
        return mycg
    
    # loop over odorants
    cglist2 = []
    for j in range(len(ilist2)):
        # list all associated descriptor categories
        mycg = []
        for sublist in ilist2[j]:
            mycg_sub = collect_list(sublist)
            mycg.append(mycg_sub)
        cglist2.append(mycg)
    return cglist2

def fractional_assignment(cglist2, Klist):
    ''' fractional assignment '''
    
    # local function
    def count_hits(clist):
        mymat = np.equal(np.matrix(clist).transpose(), np.array(Klist)).astype('int')
        myvec = np.sum(mymat, axis=0)
        return myvec
    
    # loop
    vecs_cnt = []
    vecs_norm = []
    inds_best = []
    for j in range(len(cglist2)):
        # odor type
        myvec_type = count_hits(cglist2[j][0])
        # odor descriptors
        myvec_words = count_hits(cglist2[j][1])
        # vote
        myvec = myvec_words + 0.5 * myvec_type # let OdorType break ties
        myvecsum = np.sum(myvec)
        if myvecsum > 0:
            myvec_norm = list(myvec / myvecsum)
            idx_max = np.argmax(myvec) + 1 # make 1-based
        else:
            myvec_norm = myvec # all zeros
            idx_max = 0 # keep 0
        vecs_cnt.append(myvec)
        vecs_norm.append(myvec_norm)
        inds_best.append(idx_max)
    odorant_cgrp = Tray()
    odorant_cgrp.vecs_cnt = vecs_cnt
    odorant_cgrp.vecs_norm = vecs_norm
    odorant_cgrp.best = inds_best
    return odorant_cgrp


# /

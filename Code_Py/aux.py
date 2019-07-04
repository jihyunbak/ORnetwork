# aux.py
# auxiliary functions

# Copyright 2019 Ji Hyun Bak

import numpy as np
import pandas as pd

# for stat
from scipy.sparse import coo_matrix
from scipy import stats

# for io
import csv

# for plot
import matplotlib as mpl
import matplotlib.pyplot as plt

# === ds: custom data structure 

class Tray:
    ''' empty class, to emulate Matlab's struct '''
    def __init__(self):
        pass
    def get_attr_keys(self):
        dkey = self.__dict__.keys()
        return dkey

# /

# === dm: data manipulation

# --- pandas DataFrame specific

def collect_df_rows_by_index(df, idx_input, drop=True):
    # should extend for the bad-index case (NaN)
    idx = idx_input.astype('int')
    df_new = df.iloc[idx].reset_index(drop=drop)
    return df_new

def convert_data_types(df, fields, type):
    for myfield in fields:
        myvalue = getattr(df, myfield).astype(type)
        setattr(df, myfield, myvalue)
    return df

def sort_and_reset_index(intab, columns, drop=True):
    ''' sort by columns and reset index '''
    sorttab = intab.sort_values(columns)
    outtab = sorttab.reset_index(drop=drop)
    return outtab

# --- other

def find_equal(listlike, targ):
    idx_hit = []
    for m in range(len(listlike)):
        if targ == listlike[m]:
            idx_hit.append(m)
    return idx_hit

def find_idx(testlist_bool):
    # https://stackoverflow.com/questions/364621/how-to-get-items-position-in-a-list
    myidx = [i for i,x in enumerate(testlist_bool) if x == 1]
    return myidx

def findby(vlist, testlist_bool):
    myidx_list = find_idx(testlist_bool)
    val = [vlist[i] for i in myidx_list]
    return val

def isin_lists(list, testlist):
    y_array = np.isin(np.array(list), np.array(testlist))
    y = y_array.tolist()
    return y

def normalize_by(mat, axis):
    mysum = np.sum(mat, axis=axis)
    newmat = np.true_divide(mat, mysum)
    return newmat

def center_by(mat, axis):
    mymean = np.mean(mat, axis=axis)
    newmat = mat - mymean
    return newmat

# /

# === stat: reusable statistics

# --- counting & probability estimation

def count_with_weight(vec, wgt=None, *args):
    # v_uniq, v_cnt = np.unique(vec, return_counts=True)
    if wgt is None:
        wgt = np.ones(np.size(vec))
    v_uniq = np.unique(vec).tolist()
    v_wgtcnt = []
    for vu in v_uniq:
        myidx = find_idx(isin_lists(vec, vu))
        mywgtcnt = sum([wgt[i] for i in myidx])
        v_wgtcnt.append(mywgtcnt)
    return v_uniq, v_wgtcnt

def samp_prob1(vec, wgt=None, normalize=True):
    ''' sampled probability for one variable with discrete values '''
    v_uniq, v_cnt = count_with_weight(vec, wgt)
    cnt_mat = np.matrix(v_cnt).transpose()
    if normalize:
        cnt_mat = normalize_by(cnt_mat, axis=None) # single dimension
    return cnt_mat, v_uniq

def samp_joint_prob(v1, v2, wgt=None, normalize=True):
    ''' sampled joint probability for two variables v1 and v2 '''
    
    if not wgt:
        wgt = np.ones(np.size(v1))
    
    # use COO matrix
    v1uniq, v1iinv = np.unique(v1, return_inverse=True) # renumber
    v2uniq, v2iinv = np.unique(v2, return_inverse=True)
    mat_shape = (len(v1uniq), len(v2uniq))
    cnt_mat_sparse = coo_matrix((wgt, (v1iinv, v2iinv)), shape=mat_shape)
    cnt_mat = cnt_mat_sparse.todense()
    if normalize:
        cnt_mat = cnt_mat / np.sum(cnt_mat) # normalize by all-entries sum
    
    return cnt_mat, v1uniq, v2uniq

def get_joint_counts(vars, wgt, names=('v1', 'v2')):
    ''' 
    given simultaneous samples of two variables v1 and v2, 
    compute the joint counts and probabilities and return DataFrame objects.
    each row is a distinct value of v1 (first input); 
    each column is a distinct value of v2 (second input).
    
    INPUT: vars = (v1, v2) and names = (v1name, v2name) are tuples.
    OUTPUT: (cnts, probs) with Tray objects cnts and probs.
    '''
    
    # unpack input
    (h2, b2) = vars
    (v1name, v2name) = names
    
    # -- count matrices
    # receptor code groups (marginal counts)
    p_h, h2_uniq1 = samp_prob1(h2, wgt=wgt, normalize=True)
    cnt_h, _ = samp_prob1(h2, wgt=wgt, normalize=False)
    dat_h = np.concatenate((cnt_h.astype('int'), p_h), axis=1)
    # perceptual odor categories (marginal counts)
    p_b, b2_uniq1 = samp_prob1(b2, wgt=wgt, normalize=True)
    cnt_b, _ = samp_prob1(b2, wgt=wgt, normalize=False)
    dat_b = np.concatenate((cnt_b.astype('int'), p_b), axis=1)
    # joint statistics
    p_hb, _, _ = samp_joint_prob(h2, b2, wgt=wgt, normalize=True)
    cnt_hb, _, _ = samp_joint_prob(h2, b2, wgt=wgt, normalize=False)
    # expected joint distribution (product of marginals)
    dat_p_exp = np.multiply(np.matrix(p_h), np.matrix(p_b).transpose())
    
    # -- make DataFrame objects
    names_h = [v1name + '=' + str(h) for h in h2_uniq1]
    names_b = [v2name + '=' + str(b) for b in b2_uniq1]
    cnt_h_df = pd.DataFrame(data=dat_h, index=names_h, columns=['cnt', 'p'])
    cnt_b_df = pd.DataFrame(data=dat_b, index=names_b, columns=['cnt', 'p'])
    cnt_hb_df = pd.DataFrame(data=cnt_hb.astype('int'), index=names_h, columns=names_b)
    p_hb_df = pd.DataFrame(data=p_hb, index=names_h, columns=names_b)
    p_exp_df = pd.DataFrame(data=dat_p_exp, index=names_h, columns=names_b)
    
    # -- pack output and return
    # raw counts
    cnts = Tray()
    setattr(cnts, v1name, cnt_h_df)
    setattr(cnts, v2name ,cnt_b_df)
    cnts.joint = cnt_hb_df
    # joint probabilities
    probs = Tray()
    probs.obs = p_hb_df
    probs.exp = p_exp_df
    
    return cnts, probs

# --- statistical test

def chisq_to_pvals(chisq, dof):
    pval_lo = stats.chi2.cdf(chisq, dof)
    pval_up = 1 - stats.chi2.cdf(chisq, dof)
    return (pval_lo, pval_up)

# /

# === io: file input/output

def csv_to_df(filename, delimiter=','):
    '''
    assuming a single header line,
    read a csv file and return a pandas DataFrame
    '''
    dat, header = mycsvread(filename, 1, delimiter=delimiter)
    df = pd.DataFrame(dat, columns=header[0])
    return df

def mycsvread(filename, nheader=0, row_filter=None, \
    encoding='utf-8', delimiter=','):
    '''
    reads from a csv file and returns a list (or two lists)
    optionally reads the first n lines seperately as header (default is 0)
    optinally specify the encoding (default is utf-8)
    '''
    
    # -- default is to read each row as-is
    if not row_filter:
        row_filter = lambda row: row # dummy function to just return the input
    
    # -- read the file content
    mylist = []
    myheader = []
    cnt = 0
    with open(filename, 'r', newline='', encoding=encoding) as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            # read row as header
            if(cnt < nheader):
                myheader.append(row)
                cnt = cnt + 1
                continue
            # read row as body
            myrow = row_filter(row)
            mylist.append(myrow)
    if nheader>0:
        return mylist, myheader
    else:
        return mylist

# /

# === plot: reusable plots

# --- 2D heatmap ---

def draw_heatmap(data, row_labels, col_labels, filename=None, extend='neither', **kwargs):
    fig, ax = plt.subplots()
    im = ax.imshow(data, **kwargs)

    # tick labels
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # colorbar
    cbar = ax.figure.colorbar(im, ax=ax, extend=extend)
    
    # TODO: also annotate values in each cell?
    
    if not filename:
        pass
    else:
        plt.savefig(filename)
        print('figure saved to:' + filename)
    plt.show()

def draw_heatmap_df(mydf, **kwargs):
    draw_heatmap(mydf.get_values(), mydf.index, mydf.columns, **kwargs)

# /

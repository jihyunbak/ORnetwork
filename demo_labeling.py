# demo_labeling.py

# Demo script for the perceptual labeling of 
# descriptors, odorants, and the receptor-code-based groups of odorants.

# Copyright 2019 Ji Hyun Bak


# === initialize 

import pandas as pd
import numpy as np

from Code_Py import aux
import Code_Py.labeling as lbl

Tray = aux.Tray

# /

# === load data

# NOTE: all input data files are in Data/database/

# load descriptor space
_, itab, _, _ = lbl.load_descriptor_space(filter_corr=True)

# load Mainland2015 odorant list (with goodScents correspondence)
Ltab = lbl.load_odorants_Mainland15_with_goodscent_corr()

# load manually annotated Castro2013 grouping
ctab = lbl.load_descriptors_with_Castro13_category()

# /

# === perceptual odor labeling

# assign a K-dim vector to each odorant, averaging over its descriptors
odorcat_names = ['citrus', 'chemical', 'minty', 'floral', 'fruity', 'green', 'spicy', 'animal', 'unsorted']
Klist = [n for n in range(1,len(odorcat_names))]
odorant_cgrp = lbl.assign_wgroup_to_odorants(itab, ctab, Klist)  # function is in module 'actions'

# odorant group indices (either Rgroup or Wgroup) and odorant degrees
h2 = Ltab.h2.astype('int').tolist()
b2 = odorant_cgrp.best
b2_0to9 = [(b if b>0 else 9) for b in b2] # replace 0's with 9's (for unsorted)
degL = Ltab.deg.astype('int').tolist()

# take the 6 largest groups + all others (7 groups total)
ngrpcut = 6
h2_trimmed = [(s if s<=ngrpcut else ngrpcut+1) for s in h2]

# reorder odor categories, for plotting
def reorder_odor_cats(b2):
    newlist = []
    bmap = [9, 8, 4, 3, 2, 7, 1, 6, 5]
    newlabels = ['green', 'floral', 'minty', 'chemical', 'animal', 'spicy', 'fruity', 'citrus', 'unsorted']
    for bfull in b2:
        bnew = bmap[bfull]
        newlist.append(bnew)
    return (newlist, newlabels)

# regroup odor categories (reduced labels)
def regroup_odor_cats(b2):
    newlist = []
    bmap = [5, 4, 1, 1, 1, 3, 1, 3, 2]
    # newlabels = ['green/floral/minty', 'animal', 'spicy/fruity', 'citrus', 'chemical/unsorted']
    newlabels = ['plant-like', 'animal-like', 'culinary', 'citrus', 'other']
    for bfull in b2:
        bnew = bmap[bfull]
        newlist.append(bnew)
    return (newlist, newlabels)

b2_reordered, odorcat_names_reordered = reorder_odor_cats(b2)
b2_regrouped, odorcat_re_names = regroup_odor_cats(b2)

# /


# === Rgroup vs Wgroup: joint distribution and observed/expected ratio

# --- count group frequencies ---

# get count tables
grp_names = ('grp', 'cat')
cnt1t, _ = aux.get_joint_counts((h2_trimmed, b2_0to9), None, grp_names) # unweighted counts
cntwt, _ = aux.get_joint_counts((h2_trimmed, b2_0to9), degL, grp_names) # weighted by degree

# view counts (unweighted and weighted by degree)
cnt1t.joint
cntwt.joint

# ---- observed/expected ratio ----

def compute_diff_cnt(h2, b2, *args, ngrpcut=None):

    # count h2 (raw counts)
    cnt_h, h2_uniq = aux.samp_prob1(h2, *args, normalize=False) # raw counts
    cnt_mat0 = np.concatenate((np.matrix([h2_uniq]).transpose(), cnt_h), axis=1)
    
    # marginal probabilities (h2 and b2 separately)
    p_h, h2_uniq0 = aux.samp_prob1(h2, *args) # normalized
    p_b, b2_uniq = aux.samp_prob1(b2, *args)
    
    # joint probability
    p_hb, h2_uniq2, b2_uniq2 = aux.samp_joint_prob(h2, b2, *args)
    
    # observed/expected ratio
    cutidx = slice(0,ngrpcut) # slice(0,None) gives entire range
    p_h_cut = aux.normalize_by(p_h[cutidx], axis=None) # re-normalize
    mat_exp = np.multiply(p_h_cut, p_b.transpose())
    mat_obs = aux.normalize_by(p_hb[cutidx,:], axis=None) # re-normalize
    cdiff_mat = np.true_divide(mat_obs, mat_exp)
    
    
    # make DataFrame objects
    # - probs
    cnt_mat = np.concatenate((cnt_mat0, p_hb), axis=1)
    cnt_header = ['h2_idx', 'h2_cnt'] + ['cond_k' + str(k) for k in b2_uniq2]
    cnt_df = pd.DataFrame(data=cnt_mat, columns=cnt_header)
    # - prob ratios
    cdiff_header = ['cdiff_' + str(k) for k in b2_uniq2]
    cdiff_df = pd.DataFrame(data=cdiff_mat, columns=cdiff_header)
    
    # pass additional data
    more_out = Tray()
    more_out.p_h = p_h
    more_out.p_h_cut = p_h_cut
    more_out.p_b = p_b
    # more_out.p_hb = p_hb
    more_out.mat_obs_cut = mat_obs
    more_out.mat_exp_cut = mat_exp
    more_out.mat_oe = cdiff_mat
    more_out.cnt_df = cnt_df
    more_out.cdiff_df = cdiff_df
    
    return p_hb, more_out

def proc_cdiff_evaluate_and_plot(h2, b2, b2_labels=None, use_degree=True, figname_header=None, ngrpcut=None, **kwargs):

    # optionally weighted by odorant degree (# receptors activated)
    vec_ones = [1 for j in range(len(degL))]
    wgt_degree = degL if use_degree else vec_ones
    
    # count
    _, out = compute_diff_cnt(h2, b2, wgt_degree, ngrpcut=ngrpcut)
    mat_obs = out.mat_obs_cut # p(h,b)
    mat_exp = out.mat_exp_cut # p(h)p(b)
    mat_oe = out.mat_oe[0:ngrpcut,:] # p(h,b)/p(h)p(b)
    cfdata = np.log(mat_oe)

    # prepare for plot
    if b2_labels is None:
        b2_labels = ['cat ' + str(b) for b in range(1, np.add(1, np.size(cfdata,1)))]
    h2_labels = ['group ' + str(h) for h in range(1, np.add(1, np.size(cfdata,0)))]
    
    if not figname_header:
        figname = None
    else:
        wgttag = '_wgtdeg' if use_degree else '_wgt1'
        gcuttag = '_gcut' + str(ngrpcut) if (ngrpcut is not None) else ''
        figname = figname_header + wgttag + gcuttag + condtag + '.pdf'
    
    # make plot
    aux.draw_heatmap(cfdata, h2_labels, b2_labels, filename=figname, **kwargs)
    
# evaluate and plot
np.seterr(divide='ignore') # let's ignore the divide-by-zero warning
kwargs_colorbar = {'extend':'both', 'cmap':'bwr', 'vmin':-1.5, 'vmax':1.5}
proc_cdiff_evaluate_and_plot(h2_trimmed, b2_0to9, odorcat_names, **kwargs_colorbar)

# /

# === KL divergence test

def kld_density_to_pval_up(kld_summand, num_samples):
    
    kld_all = np.nansum(kld_summand, axis=None)
    Gvalue = 2*num_samples*kld_all
    (nr, nc) = kld_summand.shape
    dof = (nr - 1) * (nc - 1)
    _, pval_up = aux.chisq_to_pvals(Gvalue, dof)
    return pval_up, Gvalue

def get_kldiv(h2, b2, wgt, ngrpcut=None):

    _, out = compute_diff_cnt(h2, b2, wgt, ngrpcut=ngrpcut)

    mat_obs = out.mat_obs_cut # p(h,b)
    mat_exp = out.mat_exp_cut # p(h)p(b)

    mat_oe = out.mat_oe[0:ngrpcut,:] # p(h,b)/p(h)p(b)
    mat_oe_b = np.true_divide(aux.normalize_by(mat_obs, axis=1), np.sum(mat_exp, axis=1)) # p(h|b)/p(h), at *fixed* b
    mat_oe_h = np.true_divide(aux.normalize_by(mat_obs, axis=0), np.sum(mat_exp, axis=0)) # p(b|h)/p(b), at *fixed* h

    # use natural log (unit in nats)
    kld_summand_all = np.multiply(mat_obs, np.log(mat_oe))  # np.log: returns -Inf for input 0
    kld_summand_b = np.multiply(aux.normalize_by(mat_obs, axis=1), np.log(mat_oe_b)) # at *fixed* b
    kld_summand_h = np.multiply(aux.normalize_by(mat_obs, axis=0), np.log(mat_oe_h)) # at *fixed* h

    # np.nansum: treats NaNs as zero
    kld_all = np.nansum(kld_summand_all, axis=None) # KLD[p(h,b) || p(h)p(b)]
    kld_b = np.nansum(kld_summand_b, axis=0) # KLD[p(h|b) || p(h)] at fixed b
    kld_h = np.nansum(kld_summand_h, axis=1) # KLD[p(b|h) || p(b)] at fixed h

    # make dataframes
    _, h2_uniq1 = aux.samp_prob1(h2, wgt=wgt, normalize=True)
    names_h = ['h=' + str(h) for h in h2_uniq1[0:ngrpcut]]
    _, b2_uniq1 = aux.samp_prob1(b2, wgt=wgt, normalize=True)
    names_b = ['b=' + str(b) for b in b2_uniq1]
    
    kld_summand_all_df = pd.DataFrame(kld_summand_all, columns=names_b, index=names_h)
    kld_summand_b_df = pd.DataFrame(kld_summand_b, columns=names_b, index=names_h)
    kld_summand_h_df = pd.DataFrame(kld_summand_h, columns=names_b, index=names_h)
    kld_b_df = pd.DataFrame(data=kld_b, columns=names_b)
    kld_h_df = pd.DataFrame(data=kld_h, index=names_h)
    
    # p-value estimate
    num_samples = len(b2)
    pval_up, Gvalue = kld_density_to_pval_up(kld_summand_all, num_samples)
    
    # pack 
    kld = Tray()
    kld.all = kld_all
    kld.b = kld_b_df
    kld.h = kld_h_df
    
    summand = Tray()
    summand.all = kld_summand_all_df
    summand.b = kld_summand_b_df
    summand.h = kld_summand_h_df
    
    kl_out = Tray()
    kl_out.kld = kld
    kl_out.summand = summand
    kl_out.pval = pval_up
    kl_out.G = Gvalue
    
    return kl_out

def proc_kldiv_evaluate_and_show(wgt, regroup_categories, figname=None, **kwargs):
    b2_input = b2_regrouped if regroup_categories else b2_0to9
    b2_names = odorcat_re_names if regroup_categories else odorcat_names
    out = get_kldiv(h2_trimmed, b2_input, wgt)
    # -- print values
    print('KLD = ' + str(out.kld.all)) # in nats
    print('G = ' + str(out.G))
    print('p = ' + str(out.pval))
    # -- test plot
    myplotdf = out.summand.all.copy()
    myplotdf.columns = b2_names
    aux.draw_heatmap_df(myplotdf, cmap='OrRd', filename=figname, **kwargs)
    return out


# evaluate KL divergence and print values

np.seterr(invalid='ignore') # let's also ignore the invalid value warning
Gtest0 = proc_kldiv_evaluate_and_show(degL, False)

# just set a lower bound to colorbar
Gtest0 = proc_kldiv_evaluate_and_show(degL, False, extend='min', vmin=0)

# with reduced labels

Gtest1 = proc_kldiv_evaluate_and_show(degL, True)

# just set a lower bound to colorbar
Gtest1 = proc_kldiv_evaluate_and_show(degL, True, extend='min', vmin=0)

# /

# === chi-square test ===

def chisq_density_to_pval_up(chisq_density_mat):
    (nr, nc) = chisq_density_mat.shape
    dof = (nr - 1) * (nc - 1)
    chisq_sum = np.sum(chisq_density_mat, axis=None)
    _, pval_up = aux.chisq_to_pvals(chisq_sum, dof)
    return pval_up

def get_chisq(h2, b2, wgt, ngrpcut=None):

    _, out = compute_diff_cnt(h2, b2, wgt, ngrpcut=ngrpcut)

    num_samples = len(h2) # ad hoc un-normalizing
    mat_obs = num_samples * out.mat_obs_cut # p(h,b)
    mat_exp = num_samples * out.mat_exp_cut # p(h)p(b)

    mat_oediffsq = np.power(mat_obs - mat_exp, 2)
    mat_chisq_each = np.true_divide(mat_oediffsq, mat_exp)

    # np.nansum: treats NaNs as zero
    chisq_all = np.nansum(mat_chisq_each, axis=None)

    # make dataframes
    _, h2_uniq1 = aux.samp_prob1(h2, wgt=wgt, normalize=True)
    names_h = ['h=' + str(h) for h in h2_uniq1[0:ngrpcut]]
    _, b2_uniq1 = aux.samp_prob1(b2, wgt=wgt, normalize=True)
    names_b = ['b=' + str(b) for b in b2_uniq1]
    
    chisq_each_df = pd.DataFrame(mat_chisq_each, columns=names_b, index=names_h)
    
    pval_up = chisq_density_to_pval_up(mat_chisq_each)
    
    # pack 
    chisq = Tray()
    chisq.sum = chisq_all
    chisq.each = chisq_each_df
    chisq.numer = pd.DataFrame(mat_oediffsq, columns=names_b, index=names_h)
    chisq.denom = pd.DataFrame(mat_exp, columns=names_b, index=names_h)
    chisq.pval = pval_up
    
    return chisq


# --- chi-sq test: evaluate and plot

def proc_chisq_evaluate_and_show(wgt, regroup_categories, figname=None, cmap='OrRd', **kwargs):
    b2_input = b2_regrouped if regroup_categories else b2_0to9
    b2_names = odorcat_re_names if regroup_categories else odorcat_names
    chisq = get_chisq(h2_trimmed, b2_input, wgt)
    # -- print values
    print('chisq = ' + str(chisq.sum))
    print('p = ' + str(chisq.pval))
    # -- test plot
    myplotdf = chisq.each.copy()
    myplotdf.columns = b2_names
    aux.draw_heatmap_df(myplotdf, cmap=cmap, filename=figname, **kwargs)
    return chisq

_ = proc_chisq_evaluate_and_show(degL, False, cmap='Purples', vmin=0)


_ = proc_chisq_evaluate_and_show(degL, True, cmap='Purples', vmin=0)


# /

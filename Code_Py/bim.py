# bim.py
# a class for pairwise matrices representing bipartite interactions

# Copyright 2019 Ji Hyun Bak

import numpy as np
from scipy.sparse import dok_matrix, coo_matrix

class PairMat:
    def __init__(self, dat, rows=None, cols=None):
        self.set_data(dat)
        self.set_rows(rows)
        self.set_columns(cols)
        
    def set_data(self, dat):
        self.rawdat = dat
        if type(dat) is np.matrix:
            self.dat = dat
        elif type(dat) is dok_matrix:
            # allows a DOK (dictionary of keys) sparse matrix
            self.dat = dat
        elif (type(dat) is list):
            self.dat = dok_matrix_from_rows(dat)
        elif type(dat) is np.ndarray:
            self.dat = np.matrix(dat)
        elif type(dat) is tuple:
            # ijv type
            self.dat = coo_matrix_from_ijk(dat).todok()
        else:
            print('unknown data format.')
            self.dat = None
        
    def set_rows(self, rows):
        num_rows = np.size(self.dat, axis=0)
        if (type(rows) is list) and (len(rows) == num_rows):
            myrows = rows
        else:
            myrows = ['row' + str(i) for i in range(num_rows)]
        self.row_names = myrows 
        self.num_rows = num_rows
        
    def set_columns(self, columns):
        num_cols = np.size(self.dat, axis=1)
        if (type(columns) is list) and (len(columns) == num_cols):
            mycols = columns
        else:
            mycols = ['col' + str(i) for i in range(num_cols)]
        self.col_names = mycols
        self.num_cols = num_cols
        
    # ----- non-constructor methods -----
        
    def shape(self):
        return get_matrix_shape(self.dat)
        
    def degree(self, axis=1):
        return get_degree(self.dat, axis)
        
    def projection(self, axis=1):
        '''
        projection onto one node type.
        axis specifies the (other) node type that is summed over: 0 or 1
        '''
        mymat = self.dat
        if axis==0:
            mycomat = mymat.T.dot(mymat) # [num_cols num_cols] matrix
        elif axis==1:
            mycomat = mymat.dot(mymat.T) # [num_rows num_rows] matrix
        else:
            mycomat = None # error
        return mycomat
        
# -----------------------------------------------------------------------------

def coo_matrix_from_ijk(ijk_tuple, **kwargs):
    ''' 
    build a COOrdinate sparse matrix using ijv format
    each of i,j,k can be a np array (what other formats?)
    '''
    (row, col, data) = ijk_tuple    # unpack
    # if not shape:
    #     num_rows = max(row) + 1     # 0-based
    #     num_cols = max(col) + 1
    #     shape = (num_rows, num_cols)
    mymat = coo_matrix((data, (row, col)), **kwargs)
    return mymat

def dok_matrix_from_rows(list_rows, num_cols=None):
    '''
    builds a DOK (dictionary of keys) sparse pairwise matrix 
    from a list of nonzero indices in each row.
    for now binary matrix only.
    returns a dok_matrix object.
    '''
    
    # -- detect size
    num_rows = len(list_rows)
    if num_cols is None:
        num_cols = 0
        for row in list_rows:
            if not row:
                continue
            num_cols = max([num_cols, max(row) + 1])
    
    # -- make a binary matrix
    Wmat = dok_matrix((num_rows, num_cols))
    for irow in range(0, num_rows):
        myidxlist = list_rows[irow]
        myidxlist = list(set(myidxlist)) # remove repeated values
        for icol in myidxlist:
            Wmat[irow,icol] = 1 # for now binary only 
    return Wmat

def get_matrix_shape(M):
    '''
    M: a numpy matrix (or sparse matrix)
    returns a list [#rows, #cols]
    '''
    num_rows = np.size(M, axis=0)
    num_cols = np.size(M, axis=1)
    myshape = [num_rows, num_cols]
    return myshape

def get_degree(M, axis=1):
    '''
    M: a numpy matrix (or sparse matrix) represeenting a bipartite graph.
    returns the degree of one node type;
    axis specifies the (other) node type that is summed over: 0 or 1
    '''
    mydeg = np.sum(M, axis=axis)
    return mydeg

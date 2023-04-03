#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 21:00:20 2019

@author: ben
"""

import numpy as np
from ATL11.RDE import RDE
from scipy.linalg import lstsq

def my_lstsq(A, b):
    return np.linalg.solve(A.T.dot(A), A.T.dot(b))
    

class poly_ref_surf(object):
    def __init__(self, degree_xy=None, exp_xy=None, xy0=[0,0], skip_constant=False, xy_scale=1.0):
        self.x0=xy0[0]
        self.y0=xy0[1]
        if degree_xy is not None:
            poly_exp_x, poly_exp_y=np.meshgrid(np.arange(0, degree_xy[0]+1), np.arange(0, degree_xy[1]+1))
            temp=np.asarray(list(set(zip(poly_exp_x.ravel(), poly_exp_y.ravel()))))  # makes one array out of two
            # Remove exponents with exp_x+exp_y larger that max(exp_x, exp_y)
            temp=temp[np.sum(temp, axis=1)<=np.array(degree_xy).max(),:];
            # if the skip_constant option is chosen, eliminate the constant term
            if skip_constant:
                temp=temp[np.all(temp>0, axis=1)]
            # sort the exponents first by x, then by y
            temp=temp[(temp[:,0]+temp[:,1]/(temp.shape[0]+1.)).argsort()]  # orders the array: 0,0  0,1  1,0
            self.exp_x=temp[:,0]  # [0,0,1]
            self.exp_y=temp[:,1]  # [0,1,0]
        if exp_xy is not None:
            self.exp_x=np.array(exp_xy[0]).ravel()
            self.exp_y=np.array(exp_xy[1]).ravel()
        # establish arrays for fitting
        self.poly_vals=np.NaN+np.zeros(self.exp_x.shape)
        self.model_cov_matrix=None
        self.xy_scale=xy_scale
        self.skip_constant=skip_constant
        self.valid=None
        self.fit_matrix=None

    def build_fit_matrix(self, x, y):
        self.fit_matrix=np.zeros([x.size, self.exp_x.size])  # size is the ravel of shape, here G is len(x) * 3
        for col, ee in enumerate(zip(self.exp_x, self.exp_y)):
            self.fit_matrix[:,col]=((x.ravel()-self.x0)/self.xy_scale)**ee[0] * ((y.ravel()-self.y0)/self.xy_scale)**ee[1]
        return self

    def z(self, x=None, y=None, mask=None):
        # evaluate the poltnomial at [x0, y0]
        if mask is None:
            mask=np.ones_like(self.poly_vals)
        if x is not None and y is not None:
            self.build_fit_matrix(x, y)
        z=np.dot(self.fit_matrix, self.poly_vals*mask)
        if x is not None:
            z.shape=x.shape
        return z
    
    def fit(self, zd, xd=None, yd=None, sigma_d=None, mask=None, max_iterations=1, min_sigma=0):

        # asign poly_vals and cov_matrix with a linear fit to zd at points xd, yd
        # build the design matrix:
        if self.fit_matrix is None:
            self.build_fit_matrix(xd, yd)
        G=self.fit_matrix
        #initialize outputs
        m=np.zeros(G.shape[1])+np.NaN
        residual=np.zeros_like(xd)+np.NaN
        rows=np.ones_like(xd, dtype=bool)
        # build a sparse covariance matrix
        if sigma_d is None:
            sigma_d=np.ones_like(zd.ravel())
        chi2r=len(zd)*2
        if mask is None:
            mask=np.ones_like(zd.ravel(), dtype=bool)
        sigma_inv=1./sigma_d.ravel()
        sigma_inv.shape=[len(sigma_inv), 1]
        
        for k_it in np.arange(max_iterations):
            rows=mask
            if rows.sum()==0:
                chi2r=np.NaN
                break
            Gsub=G[rows,:]
            cols=(np.amax(Gsub,axis=0)-np.amin(Gsub,axis=0))>0
            if self.skip_constant is False:
                cols[0]=True
            m=np.zeros([Gsub.shape[1],1])
            # compute the LS coefficients
            # note: under broadcasting rules, sigma_inv*G = diags(sigma_inv).dot(G)
            # and  sigma_inv.ravel()*zd.ravel() = diags(sigma_inv)*zd.ravel()
            #msub = my_lstsq(sigma_inv[rows]*(Gsub[:,cols]), sigma_inv[rows].ravel()*(zd.ravel()[rows]))
            msub = lstsq(sigma_inv[rows]*(Gsub[:,cols]), sigma_inv[rows].ravel()*(zd.ravel()[rows]))[0]

            msub.shape=(len(msub), 1)
            m[np.where(cols)] = msub  # only takes first three coefficients?
            residual=zd.ravel()-G.dot(m).ravel()
            rs=residual/sigma_d.ravel()
            chi2r_last=chi2r
            deg_of_freedom=(rows.sum()-cols.sum())
            if deg_of_freedom  > 0:
                chi2r=sum(rs**2)/(deg_of_freedom)
            else:
                # the inversion is even-determined or worse, no further improvement is expected
                break
            #print('chi2r is ',chi2r)
            if np.abs(chi2r_last-chi2r)<0.01 or chi2r<1:
                break
            sigma=RDE(rs[rows])
            if not np.isfinite(sigma):
                sigma=0
            #print('sigma from RDE ',sigma)
            threshold=3.*np.max([sigma, min_sigma])
            mask=np.abs(rs)<threshold
            #print("\tsigma=%3.2f, f=%d/%d" % (sigma, np.sum(mask), len(mask)))
        self.poly_vals=m

        return m, residual, chi2r, rows

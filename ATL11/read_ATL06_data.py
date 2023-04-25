#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:27:09 2019

@author: ben
"""

import ATL11
import numpy as np
import re
import pointCollection as pc
from ATL11.check_ATL06_hold_list import check_ATL06_hold_list, read_hold_files


def get_seg_index(seg_range, D6):
    '''
    Find the indexes in an ATL06 to read to match seg_range

    required arguments:
         seg_range:  first and last segments to read
         D6:         ATL06 data structure, containing at least a segment_id field
    returns:
         list of indices needed to span seg_range
    '''
    index_range=np.searchsorted(D6.segment_id[:,0], seg_range)
    if index_range[0] > 0:
        index_range[0] -= 1
    return index_range

def get_ATL06_release(D6_list):
    """Get the release number for ATL06 data."""
    ATL06_re=re.compile('ATL06_\d+_\d+_(\d{3})_\d{2}.h5')
    for D6i in D6_list:
        release=int(ATL06_re.search(D6i.filename).group(1))
        D6i.assign({'release':np.zeros_like(D6i.delta_time)+release})

def read_ATL06_data(ATL06_files, beam_pair=2, cycles=[1, 12], use_hold_list=False, minimal=False, ATL06_dict=None, seg_range=None):
    '''
    Read ATL06 data from a list of files for a specific beam pair

    required arguments:
        ATL06_files: a list of ATL06 files
        beam_pair: pair number to read from the files
        cycles: first and last cycles to include
    '''
    params_11=ATL11.defaults()
    ATL06_re=re.compile('ATL06_\d+_\d\d\d\d(\d\d)\d\d_')

    if use_hold_list:
        hold_list = read_hold_files()
    else:
        hold_list = None

    # check the files against cycle bounds
    for filename in ATL06_files.copy():
        try:
            m=ATL06_re.search(filename)
            if (int(m.group(1)) < cycles[0]) or (int(m.group(1)) > cycles[1]) :
                ATL06_files.remove(filename)
                continue
        except Exception:
            pass
        # check ATL06 filename against hold list
        if check_ATL06_hold_list([filename], hold_list=hold_list)[0]:
            ATL06_files.remove(filename)
    if len(ATL06_files)==0:
        print("edited D6 list has no files")
        return None

    if minimal:
        # only read in the segment_id, latitude, and longitude fields, return in a dictionary
        minimal_field_dict={None:['segment_id','latitude','longitude'],'ground_track':['x_atc', 'y_atc']}
        D6_dict={}
        for filename in ATL06_files:
            temp=pc.ATL06.data(field_dict=minimal_field_dict, beam_pair=beam_pair).from_h5(filename)
            D6_dict[filename]=temp
        return D6_dict

    # read in the ATL06 data from all the repeats
    D6_list=[]

    for filename in ATL06_files:
        try:
            if seg_range is not None:
                this_index_range = get_seg_index(seg_range, ATL06_dict[filename])
            else:
                this_index_range=None
            D6_list.append(pc.ATL06.data(field_dict=params_11.ATL06_field_dict, beam_pair=beam_pair).from_h5(filename, index_range=this_index_range))
        except KeyError:
            pass

    get_ATL06_release(D6_list)

    D6=pc.ATL06.data(beam_pair=beam_pair).from_list(D6_list)

    # index D6 to remove segments with large y_atc offsets:
    D6.index(np.any(np.abs(D6.y_atc-params_11.pair_yatc_ctr[beam_pair])< params_11.pair_yatc_ctr_tol, axis=1))

    if D6.size == 0:
        return None
    # reorder data rows from D6 by cycle
    D6.index(np.argsort(D6.cycle_number[:,0],axis=0))

    # choose the hemisphere and project the data to polar stereographic
    if np.max(D6.latitude) < 0:
        D6.get_xy(None, EPSG=3031)
    else:
        D6.get_xy(None, EPSG=3413)

    return D6

def select_ATL06_data(D6, lonlat_bounds=None, first_ref_pt=None, last_ref_pt=None, num_ref_pts=None):
    """
    Select a subset of input ATL06 data, and find the associated reference-point numbers
    Required input arguments:
        D6: An ATL06 data structure (PointDatabase.ATL06_data)

    Optional input agruments
        first_ref_pt: first reference point to attempt to fit
        num_ref_pts: number of reference points to include in the fit
        last_ref_pt: last reference point to include in the fit

    """
    params_11=ATL11.defaults()

    if lonlat_bounds is not None:
        keep = (D6.longitude >= lonlat_bounds[0])
        keep &= (D6.latitude >= lonlat_bounds[1])
        keep &= (D6.longitude <= lonlat_bounds[2])
        keep &= (D6.latitude <= lonlat_bounds[3])
        keep = np.any(keep, axis=1)
        if not np.any(keep):
            return None, None, None
        D6.index(keep)

    # get list of reference points
    uId, iId=np.unique(D6.segment_id.ravel(), return_index=True)
    ctrSegs=np.mod(uId, params_11.seg_number_skip)==0
    ref_pt_numbers=uId[ctrSegs]
    ref_pt_x=D6.x_atc.ravel()[iId[ctrSegs]]

    # apply input arguments to the input reference points
    if first_ref_pt is not None:
        these=ref_pt_numbers>=first_ref_pt
        ref_pt_numbers=ref_pt_numbers[these]
        ref_pt_x=ref_pt_x[these]

    if last_ref_pt is not None and len(ref_pt_numbers)>0:
        these=ref_pt_numbers<=last_ref_pt
        ref_pt_numbers=ref_pt_numbers[these]
        ref_pt_x=ref_pt_x[these]

    if num_ref_pts is not None and len(ref_pt_numbers)>0:
        ref_pt_numbers=ref_pt_numbers[0:int(num_ref_pts)]
        ref_pt_x=ref_pt_x[0:int(num_ref_pts)]
    if len(ref_pt_numbers) > 0:
    # subset D6 to match the reference point numbers
        keep = D6.segment_id >= np.min(ref_pt_numbers) - params_11.seg_number_skip
        keep &= D6.segment_id <= np.max(ref_pt_numbers) + params_11.seg_number_skip
        D6.index(np.any(keep, axis=1))

        return D6, ref_pt_numbers, ref_pt_x
    else:
        return None, None, None

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:07:54 2020

@author: ben
"""

import numpy as np
import pointCollection as pc
import os
import h5py
from ATL11.h5util import create_attribute


class data(pc.data):
    def to_h5(self, fileOut=None,
          h5f_out=None,
          replace=True,
          group='/', extensible=True,
          dimension_fields = None,
          compression='gzip',
          chunks=True,
          meta_dict = None, DEBUG = False):
        """
        Write the contents of self to an hdf5 file.

        Saves the current data object to a file.  Optional parameters specify the
        group within the output file, and specify whether the file will be
        overwritten or if the data in self will be written to new fields and
        groups within the file.

        The meta_dict parameter allows specification of dataset attributes.  Special
        values include:
            'compression': overwrites the default gzip compression
            'source_field': alternative fieldname within self from which data
                will be copied
            'precision':
                value passed to the 'scaleoffset' keyword in h5py
        Other values in meta_dict are assigned as attributes to each dataset.

        Parameters
        ----------
        fileOut : str, optional
            Filename to which the data are written.
        h5f_out : h5py.File, optional
            File handle to which the data are written (overrides fileOut)
        replace : bool, optional
            If true, any existing file is overwritten. The default is True.
        group : str, optional
            Group into which the data are written. The default is '/'.
        extensible : bool, optional
            If true, datasets are created as extensivble. The default is True.
        meta_dict : dict optional
            Dictionary containing attributes of output datasets. The default is None.
        DEBUG : bool, optional
            If True, report error messages from dataset creation

        Returns
        -------
        None.

        """
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'
        if h5f_out is None:
            h5f_out=h5py.File(fileOut, mode)
            close_file = True
        else:
            close_file = False

        if group is not None:
            if not group in h5f_out:
                h5f_out.create_group(group.encode('ASCII'))
        field_dict={}
        if meta_dict is None:
            field_dict = {field:field for field in self.fields}
        else:
            field_dict = {}
            for out_field, this_md in meta_dict.items():
                if 'source_field' in this_md and this_md['source_field'] is not None and this_md['source_field'] in self.fields:
                    field_dict[out_field] = this_md['source_field']
                elif out_field in self.fields:
                    field_dict[out_field] = out_field

        if meta_dict is None:
            meta_dict = {out_field:{} for out_field in field_dict.keys()}

        # establish the coordinate fields first
        dimension_fields=[]
        non_dimension_fields=[]
        for out_field in field_dict.keys():
            if 'dimension' in meta_dict[out_field] and meta_dict[out_field]['dimension']:#(out_field==meta_dict[out_field]['dimensions'] or out_field in meta_dict[out_field]['dimensions']):
                dimension_fields += [out_field]
            else:
                non_dimension_fields += [out_field]

        for out_field in dimension_fields + non_dimension_fields:
            this_data=getattr(self, field_dict[out_field])
            maxshape=this_data.shape
            if extensible:
                maxshape=list(maxshape)
                maxshape[0]=None
            # try prepending the 'group' entry of meta_dict to the output field.
            # catch the exception thrown if meta_dict is None or if the 'group'
            # entry is not there
            out_field_name=out_field
            try:
                out_field_name = meta_dict[out_field]['group'] + out_field
            except (TypeError, KeyError) as e:
                if DEBUG:
                    print(e)
            out_field_name = (group + '/' +out_field).encode('ASCII')
            kwargs = dict( compression=compression,
                          maxshape=tuple(maxshape))
            if meta_dict is not None and 'precision' in meta_dict[out_field] and meta_dict[out_field]['precision'] is not None:
                kwargs['scaleoffset']=int(meta_dict[out_field]['precision'])
            if 'datatype' in meta_dict[out_field]:
                dtype = meta_dict[out_field]['datatype'].lower()
                kwargs['dtype']=dtype
                if 'int' in dtype:
                    kwargs['fillvalue'] = np.iinfo(np.dtype(dtype)).max
                else:
                    kwargs['fillvalue'] = np.finfo(np.dtype(dtype)).max
            # use an explicit '_FillValue' for preference over an overall fillvalue
            if '_FillValue' in meta_dict[out_field]:
                this_data = np.nan_to_num(this_data,nan=meta_dict[out_field]['_FillValue'])
            elif 'fillvalue' in kwargs:
                this_data = np.nan_to_num(this_data,nan=kwargs['fillvalue'])
            if 'dtype' in kwargs:
                this_data = this_data.astype(kwargs['dtype'])

            # Create the dataset
            dset = h5f_out.create_dataset(out_field_name,
                                data=this_data,  **kwargs)
            if 'dimension' in meta_dict[out_field] and meta_dict[out_field]['dimension']:
                dset.make_scale()
            if out_field in meta_dict:
                for key, val in meta_dict[out_field].items():
                    if key.lower() not in ['group','source_field','precision','dimensions']:
                        if isinstance(val, str):
                            create_attribute(h5f_out[out_field_name].id,
                                             key, [], str(val))
                        else:
                            h5f_out[out_field_name].attrs[key] = val
            if 'dimensions' in meta_dict[out_field]:
                dims = meta_dict[out_field]['dimensions']
                if isinstance(dims, str):
                    dims = dims.split(',')
                for ind, dim in enumerate(dims):
                    dset.dims[ind].label=dim
                    if '../' in dim:
                        group_path = group.split('/')
                        while '../' in dim:
                            dim=dim.lstrip('../')
                            group_path=group_path[:-1]
                        try:
                            dset.dims[ind].attach_scale("/"+h5f_out['/'.join(group_path+[dim])])
                        except Exception as e:
                            print("-----")
                            print([group,out_field_name])
                            print('/'.join(group_path+[dim]))
                            print("------")
                            raise e
                    else:
                        dset.dims[ind].attach_scale(h5f_out[dim])
            if 'fillvalue' in kwargs and '_FillValue' not in meta_dict[out_field]:
                dset.attrs['_FillValue'.encode('ASCII')] = kwargs['fillvalue']

        for key, val in self.attrs.items():
            if val is not None:
                h5f_out[group].attrs[key]=val

        for key in ['EPSG','SRS_proj4']:
            val=getattr(self, key)
            if val is not None:
                h5f_out[group].attrs[key] = val
        if close_file:
            h5f_out.close()

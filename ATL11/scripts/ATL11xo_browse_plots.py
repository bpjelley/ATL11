#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
import os, h5py, glob
import shutil
import pointCollection as pc
#from PointDatabase.mapData import mapData

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from importlib import resources
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.io.img_tiles as cimgt

import imageio
from io import BytesIO
import datetime as dt
import re

def ATL11xo_browse_plots(ATL11xo_file,
                       mosaic,
                       args=None,
                       region='AA',
                       out_path=None,
                       pdf=False,
                       nolog=False):
    print('File to plot',os.path.basename(ATL11xo_file))
    # establish output files
    ATL11xo_file_str = os.path.basename(ATL11xo_file).split('.')[0]
    if out_path is None:
        out_path = os.path.dirname(ATL11xo_file)
    if not nolog:
        log_file = '{}/ATL11_BrowsePlots_{}.log'.format(out_path, dt.datetime.now().date())
        fhlog = open(log_file,'a')
    D0 = pc.data().from_h5(args.ATL11xo_file, fields=['x','y','dem_h','fit_quality'])
    Dx = pc.data().from_h5(args.ATL11xo_file, group='crossing_track',fields=['h_corr','delta_time','atl06_quality_summary'])
    Dd = pc.data().from_h5(args.ATL11xo_file, group='datum_track',fields=['h_corr','delta_time'])


    tile_re = re.compile(r'ATL11XO_(..)_E(.*)_N(.*)_c(\d\d)_(\d\d\d)_(\d\d).h5')
    region, x0, y0, cycle, release, version = \
        tile_re.search(os.path.basename(args.ATL11xo_file)).groups()

    x0 = int(x0)*1000
    y0 = int(y0)*1000
    bounds=[x0+np.array([-args.tile_width/2, args.tile_width/2]), y0+np.array([-args.tile_width/2, args.tile_width/2])]

    DEM = pc.grid.data().from_geotif(args.mosaic, bounds=bounds)
    if 'z' in DEM.fields:
        DEM.calc_gradient()

    # make plots,
    png_buffers = []
    fig1, ax1 = plt.subplots(1,3,sharex=True,sharey=True) #, subplot_kw=dict(projection=projection))

    for ax in ax1:
        if 'z_x' in DEM.fields and np.any(np.isfinite(DEM.z_x)):
            DEM.show(ax=ax, field='z_x', xy_scale=1/1000, cmap='gray', \
                     clim = np.percentile(DEM.z_x[~np.isnan(DEM.z_x)], [5,95]),
                     interpolation='nearest', aspect='equal')

    plt.sca(ax1[0])
    good = (D0.fit_quality==0) & (Dx.atl06_quality_summary==0)
    plt.colorbar(
        plt.scatter(D0.x[good]/1000, D0.y[good]/1000, 2, c=Dx.h_corr[good], cmap='terrain'),
        label='height, m', orientation='horizontal', extend='both')

    dh = Dx.h_corr-D0.dem_h
    this_clim = np.percentile(dh[np.isfinite(dh)], [5, 95])
    this_clim = np.array([-1, 1])*np.max(np.abs(this_clim))
    plt.sca(ax1[1])
    plt.colorbar(plt.scatter(D0.x[good]/1000, D0.y[good]/1000, 2,
                             c=dh[good], cmap='RdBu', clim=this_clim),
                 orientation='horizontal', extend='both', label='height - DEM, m')

    def pt_count(D, ii):
        return len(ii)
    if np.any(good):
        n_pts = pc.apply_bin_fn(D0[good], 1000, fn=pt_count, fields='N')
    else:
                n_pts = pc.data().from_dict(
                    {'x':np.array([np.mean(D0.x)]),
                     'y':np.array([np.mean(D0.y)]),
                    'N':np.array([0])})

    plt.sca(ax1[2])
    ind=np.argsort(n_pts.N)
    plt.colorbar(plt.scatter(n_pts.x[ind]/1000, n_pts.y[ind]/1000, 2, c=n_pts.N[ind],
                             clim=[0,np.percentile(n_pts.N, 95)]),
                 orientation='horizontal',
                label='points per square km', extend='max')

    for ax in ax1:
        ax.set_xlim(bounds[0]/1000)
        ax.set_ylim(bounds[1]/1000)
        ax.set_aspect(1)

    if region=='AA':
        coordinate_description='south-polar Antarctic stereographic projection (EPSG=3031)'
        DEM_name = 'REMA DEM'
    elif region=='AR':
        coordinate_description='north-polar Arctic Sea Ice stereographic projection (EPSG=3413)'
        DEM_name= 'ArcticDEM'
    plt.figtext(0.1,0.01,f'Figure 1. Crossing-track heights from cycle {cycle} (1st panel). Crossing-track height difference from ICESat-2 DEM (2nd panel). Number of valid crossing-track heights per square kilometer (3rd panel). All overlaid on gradient of DEM. Coordinates are in {coordinate_description}',wrap=True)
    ax1[0].set_ylabel('y [km]', fontdict={'fontsize':10})
    ax1[0].set_title('Crossing-track height', fontdict={'fontsize':10});
    ax1[1].set_title('Height difference\nfrom DEM', fontdict={'fontsize':10});
    ax1[2].set_title('Point density', fontdict={'fontsize':10});
    fig1.suptitle('{}'.format(os.path.basename(args.ATL11xo_file)))
    plt.subplots_adjust(bottom=0.23, top=0.9)
    buffer = BytesIO()
    fig1.savefig(buffer, format='png')
    buffer.seek(0)
    png_buffers.append(buffer.getvalue())
#    fig1.savefig('{0}/{1}_Figure1_h_corr_NumValids_dHdtOverDEM.png'.format(out_path,ATL11xo_file_str),format='png', bbox_inches='tight')
#    fig1.savefig('{0}/{1}_BRW_default1.png'.format(out_path,ATL11xo_file_str),format='png')

    fig2,ax2 = plt.subplots(1, 3, sharey=True, figsize=[8, 6],
                            gridspec_kw = {'bottom':0.5,'top':0.8,'left':0.1,'right':0.9})
    bins = np.arange(-0.1, 3.1, 0.2)
    xl=[-0.2, 1.2]
    plt.sca(ax2[0])
    plt.hist(D0.fit_quality.astype(float), bins)
    ax2[0].set_xticks([0, 1, 2, 3],
                      ['valid',r'poly coeff \n $\sigma$ > 10','|slope| > 0.2', 'both'],
                      rotation='vertical')
    ax2[0].set_xlim([-0.2, 3.2])
    ax2[0].set_title('reference surface fit')

    bins=np.arange(-0.1, 1.1, 0.2)
    plt.sca(ax2[1])
    plt.hist(np.isnan(Dd.h_corr).astype(float), bins)
    plt.gca().set_xticks([0, 1])
    ax2[1].set_xlim(xl)
    plt.gca().set_xticklabels(['valid','invalid'])
    plt.gca().set_title(f'datum-track height\ncycle {cycle}')
    plt.gca().set_xticks([0, 1], ['valid','invalid'] , rotation='vertical')

    plt.sca(ax2[2])
    plt.hist(Dx.atl06_quality_summary.astype(float), bins)
    ax2[2].set_xlim(xl)
    plt.gca().set_xticks([0, 1], ['valid','potential\nproblems'], rotation='vertical')
    plt.gca().set_title(f'crossing-track height\ncycle {cycle}')

    ax2[0].set_ylabel('count')

    fig2.suptitle('{}'.format(os.path.basename(ATL11xo_file)))
    plt.figtext(0.1,0.01,'Figure 2. Histograms reference surface fit quality (left panel), valid datum-track heights (middle panel) and valid crossing-track heights (right panel)',wrap=True)
    plt.subplots_adjust(bottom=0.2, top=0.9)
    buffer = BytesIO()
    fig2.savefig(buffer, format='png')
    buffer.seek(0)
    png_buffers.append(buffer.getvalue())
#    fig2.savefig('{0}/{1}_Figure2_flags_hist.png'.format(out_path,ATL11xo_file_str),format='png')
#    fig2.savefig('{0}/{1}_BRW_default2.png'.format(out_path,ATL11xo_file_str),format='png')

    if pdf:    #save all to one .pdf file
        figs = list(map(plt.figure, plt.get_fignums()))
        with PdfPages('{0}/{1}.pdf'.format(out_path,ATL11xo_file_str)) as pdff:
            for fig in figs:
                pdff.savefig(fig)

    # put images into browse file
    ATL11xo_file_brw='{}/{}_BRW.h5'.format(out_path,ATL11xo_file_str)
    if os.path.isfile(ATL11xo_file_brw):
        os.remove(ATL11xo_file_brw)
    brw_template = os.path.join( str(resources.files('ATL11')),
                                "package_data",
                                "BRW_template.h5")
    shutil.copyfile(brw_template,ATL11xo_file_brw)

    with h5py.File(ATL11xo_file_brw,'r+') as hf:
#        for ii, name in enumerate(sorted(glob.glob('{0}/{1}_BRW_def*.png'.format(out_path,ATL11xo_file_str)))):
        for ii in range(len(png_buffers)):
            hf.require_group('/default')
#            img = imageio.imread(name, pilmode='RGB')
#            img = png_buffers[ii]
            img = imageio.imread(png_buffers[ii], pilmode='RGB')

#            namestr = os.path.splitext(name)[0]
#            namestr = os.path.basename(namestr).split('BRW_')[-1]

            namestr = 'default'+str(ii+1)
            dset = hf.create_dataset('default/'+namestr, img.shape, data=img.data, \
                                     chunks=img.shape, compression='gzip',compression_opts=6)
            dset.attrs['CLASS'] = np.bytes_('IMAGE')
            dset.attrs['IMAGE_VERSION'] = np.bytes_('1.2')
            dset.attrs['IMAGE_SUBCLASS'] = np.bytes_('IMAGE_TRUECOLOR')
            dset.attrs['INTERLACE_MODE'] = np.bytes_('INTERLACE_PIXEL')
#        for ii, name in enumerate(sorted(glob.glob('{0}/{1}_Figure*.png'.format(out_path,ATL11xo_file_str)))):
#            if 'Figure1' not in name and 'Figure3' not in name:
#                img = imageio.imread(name, pilmode='RGB')
#
#                namestr = os.path.splitext(name)[0]
#                namestr = os.path.basename(namestr).split('Figure')[-1]
#                dset = hf.create_dataset(namestr[2:], img.shape, data=img.data, \
#                                         chunks=img.shape, compression='gzip',compression_opts=6)
#                dset.attrs['CLASS'] = np.bytes_('IMAGE')
#                dset.attrs['IMAGE_VERSION'] = np.bytes_('1.2')
#                dset.attrs['IMAGE_SUBCLASS'] = np.bytes_('IMAGE_TRUECOLOR')
#                dset.attrs['INTERLACE_MODE'] = np.bytes_('INTERLACE_PIXEL')
        del hf['ancillary_data']
        with h5py.File(ATL11xo_file,'r') as g:
            g.copy('ancillary_data',hf)

    # remove individual png files
#    for name in sorted(glob.glob('{0}/{1}_Figure*.png'.format(out_path,ATL11xo_file_str))):
#        if os.path.isfile(name): os.remove(name)
    fhlog.close()

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('ATL11xo_file', type=str)
    parser.add_argument('mosaic', type=str)
    parser.add_argument('--tile_width', type=float, default=200e3, help='tile width, defaults to 200 km ')
    parser.add_argument('--out_path', '-o', type=str, help='default is ATL11xo_file path')
    parser.add_argument('--pdf', action='store_true', default=False, help='write images to .pdf file')
    parser.add_argument('--nolog', action='store_true', default=False, help='no writing errors to .log file')
    args=parser.parse_known_args()[0]
    ATL11xo_browse_plots(args.ATL11xo_file, args.mosaic, args=args, out_path=args.out_path, pdf=args.pdf)

if __name__=="__main__":
    main()

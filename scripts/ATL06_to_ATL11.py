#! /usr/bin/env python3

'''
Executable script to generate ATL11 files based on ATL06 data.
'''

import os
os.environ['MKL_NUM_THREADS']="1"
os.environ['NUMEXPR_NUM_THREADS']="1"
os.environ['OMP_NUM_THREADS']="1"
os.environ['OPENBLAS_NUM_THREADS']="1"

import numpy as np
import ATL11
#import write_METADATA
import glob
import sys, h5py
import matplotlib.pyplot as plt


#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -b -101. -76. -90. -74.5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5"
#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5" 

def get_proj4(hemisphere):
    if hemisphere==-1:
        return'+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' 
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

def main(argv):
    # account for a bug in argparse that misinterprets negative agruents
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    # command-line interface: run ATL06_to_ATL11 on a list of ATL06 files
    import argparse
    parser=argparse.ArgumentParser(description='generate an ATL11 file from a collection of ATL06 files.')
    parser.add_argument('rgt', type=int, help="reference ground track number")
    parser.add_argument('subproduct', type=int, help="ICESat-2 subproduct number (latltude band)")
    parser.add_argument('--directory','-d', default=os.getcwd(), help="directory in which to search for ATL06 files")
    parser.add_argument('--pair','-p', type=int, default=None, help="pair number to process (default is all three)")
    parser.add_argument('--Release','-R', type=int, default=2, help="Release number")
    parser.add_argument('--Version','-V', type=int, default=1, help="Version number")
    parser.add_argument('--cycles', '-c', type=int, nargs=2, default=[3, 4], help="first and last cycles")
    parser.add_argument('--GI_file_glob','-G', type=str, default=None, help="Glob (wildcard) string used to match geoindex files for crossing tracks")
    parser.add_argument('--out_dir','-o', default=None, required=True, help="Output directory")
    parser.add_argument('--first_point','-f', type=int, default=None, help="First reference point")
    parser.add_argument('--last_point','-l', type=int, default=None, help="Last reference point")
    parser.add_argument('--num_points','-N', type=int, default=None, help="Number of reference points to process")
    parser.add_argument('--Hemisphere','-H', type=int, default=-1)
    parser.add_argument('--bounds', '-b', type=float, nargs=4, default=None, help="latlon bounds: west, south, east, north")
    parser.add_argument('--test_plot', action='store_true', help="plots locations, elevations, and elevation differences between cycles")
    parser.add_argument('--Blacklist','-B', action='store_true')
    parser.add_argument('--verbose','-v', action='store_true')
    args=parser.parse_args()

    # output file format is ATL11_RgtSubprod_c1c2_rel_vVer.h5
    out_file="%s/ATL11_%04d%02d_%02d%02d_%03d_%02d.h5" %( \
            args.out_dir,args.rgt, args.subproduct, args.cycles[0], \
            args.cycles[1], args.Release, args.Version)
    if os.path.isfile(out_file):
        os.remove(out_file)

    if args.verbose:
        print('ATL11 output filename',out_file)
    glob_str='%s/*ATL06*_*_%04d??%02d_*.h5' % (args.directory, args.rgt, args.subproduct)
    files=glob.glob(glob_str)

    print("found ATL06 files:" + str(files))

    if args.pair is None:
        pairs=[1, 2, 3]
    else:
        pairs=[args.pair]
    
    if args.GI_file_glob is not None:
        GI_files=glob.glob(args.GI_file_glob)
    else:
        GI_files=None   
    print("found GI files:"+str(GI_files))
    
    for pair in pairs:
        D6 = ATL11.read_ATL06_data(files, beam_pair=pair, cycles=args.cycles, use_blacklist=args.Blacklist)
        if D6 is None:
            continue
        D6, ref_pt_numbers, ref_pt_x = ATL11.select_ATL06_data(D6, \
                            first_ref_pt=args.first_point,\
                            last_ref_pt=args.last_point, \
                            lonlat_bounds=args.bounds, 
                            num_ref_pts=args.num_points)

        if D6 is None or len(ref_pt_numbers)==0: 
            continue
        D11=ATL11.data().from_ATL06(D6, ref_pt_numbers=ref_pt_numbers, ref_pt_x=ref_pt_x,\
                      cycles=args.cycles, beam_pair=pair, verbose=args.verbose, \
                      GI_files=GI_files, hemisphere=args.Hemisphere) # defined in ATL06_to_ATL11
        if D11 is None:
            print("ATL06_to_ATL11: Not enough good data to calculate an ATL11, nothing written")
            return()
        # fill cycle_number list in cycle_stats and ROOT
        setattr(D11.cycle_stats,'cycle_number',list(range(args.cycles[0],args.cycles[1]+1)))
        setattr(D11.ROOT,'cycle_number',list(range(args.cycles[0],args.cycles[1]+1)))
        # add dimensions to D11
        D11.N_pts, D11.N_cycles = D11.ROOT.h_corr.shape
        
        if isinstance(D11.crossing_track_data.h_corr, np.ndarray):
            D11.Nxo = D11.crossing_track_data.h_corr.shape[0]
        
        if D11 is not None:
            D11.write_to_file(out_file)

    out_file = ATL11.write_METADATA.write_METADATA(out_file,files)

    print("ATL06_to_ATL11: done with "+out_file)
        
    if args.test_plot:
        ATL11.ATL11_test_plot.ATL11_test_plot(out_file)
#        ATL11.ATL11_browse_plots.ATL11_browse_plots(out_file,args.Hemispher,mosaic=mosaic)

if __name__=="__main__":
    main(sys.argv)

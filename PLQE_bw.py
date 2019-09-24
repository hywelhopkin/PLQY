#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# written by bernard.wenger@heliodm.com

# Using the procedure of DeMello et al https://doi.org/10.1002/adma.19970090308

# History
# 6/8/2019, v1.0, converted from Matlab
# 9/8/2019, v1.1, first working version (no added functions)
# 22/9/2019, v1.2, added QE pro spectrometer

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.io import loadmat

# define calibration files
cal_dict = {
        "Maya red" : 'cal_Maya_red.txt',
        "Maya steel" : 'cal_Maya_steel.txt',
        "QE red 25"  : 'cal_QEpro_red_25um.txt',
        "QE red 200" : 'cal_QEpro_red_200um.txt',
        "QE steel 25" : 'cal_QEpro_steel_25um.txt',
        "QE steel 200" : 'cal_QEpro_steel_200um.txt',
    }

def get_args():
    """Get arguments and options"""
    parser = argparse.ArgumentParser(description='Calculation of PL quantum efficiency.')
    
    parser.add_argument('-sn', '--short_name', type=str, required=True, help="Prefix of the file names (e.g 'filename_') to be followed by 'in.txt', 'out.txt', 'bckg.txt' or 'empty.txt'")
    
    opt = parser.add_argument_group('optional arguments')
    opt.add_argument('-d', '--directory', default='.', type=str, help="folder in which the raw files are stored (e.g. 'folder_name/sub_folder/'")    
    opt.add_argument('-cfg', '--config', default='Maya red', type=str, help="Fiber and spectrometer configurations. Choices: 'Maya red', 'Maya steel', 'QE red 25', 'QE red 200', 'QE steel 25' or 'QE steel 200', Default = 'Maya red'")    
    opt.add_argument('-lr', '--laser_range', nargs = 2, default= [395, 410], type=int, help="Laser wavelength range in nm, Default = [397, 407]")
    opt.add_argument('-plr', '--pl_range', nargs = 2, default= [550, 850], type=int, help="PL detection range in nm, Default = [550, 850]")
    opt.add_argument('-st', '--short_time', default= 1, type=int, help="Integration time for short measurement in ms")
    opt.add_argument('-ln', '--long_name', default= '', type=str, help="Prefix of the long file name (e.g 'filename_long_') to be followed by 'in.txt', 'out_.txt', 'bckg.txt' or 'empty.txt'")
    opt.add_argument('-lt', '--long_time', default= 0, type=int, help="Integration time for long measurement in ms")
    opt.add_argument('-c', '--common', default= True, type=bool, help="Indicates that common background and empty files are used. Default = True")
    opt.add_argument('-cb', '--common_bckg', default= 'bckg.txt', type=str, help=" Name of the common background file. Default = 'bckg.txt'")
    opt.add_argument('-ce', '--common_empty', default= 'empty.txt', type=str, help=" Name of the common empty file. Default = 'empty.txt'")
    opt.add_argument('-clb', '--common_long_bckg', default= 'long_bckg.txt', type=str, help=" Name of the common long background file. Default = 'long_bckg.txt'")
    opt.add_argument('-cle', '--common_long_empty', default= 'long_empty.txt', type=str, help=" Name of the common long empty file. Default = 'long_empty.txt'")
   
    args = parser.parse_args()
    
    return args


def PLQE(args):
    # Function running all calculations
    data = loadit(args) # Load data and calibration file
    cal = np.loadtxt('cal/' + cal_dict.get(args.config))

    def trim(data, vr):
        # function used to select the valid range (could also be used for removing hot pixels)
        return data[vr[0]:vr[1], :]

    if 'QE' in args.config:
        vr = [4, -5] # valid range limits
        data = trim(data, vr)
    elif 'Maya' in args.config:
        vr = [5, -6] # valid range limits
        data = trim(data, vr)
    
    cal = np.interp(data[:, 0], cal[:, 0], cal[:, 1])

    # Continue here so that it works for every calibration file

    # Unpack data
    short_in = data[:, 0:2]
    short_out = data[:, 2:4]
    short_bckg = data[:, 4:6]
    short_empty = data[:, 6:8]
    wl = short_in[:, 0]

    # Process
    short_in[:, 1] -= short_bckg[:, 1]
    short_out[:, 1] -= short_bckg[:, 1]
    short_empty[:, 1] -= short_bckg[:, 1]

    def scale(d, time):
        # divide by integration time
        return (d[:, 1] - np.mean(d[20:50, 1])) / time

    short_in[:, 1] = scale(short_in, args.short_time)
    short_out[:, 1] = scale(short_out, args.short_time)
    short_empty[:, 1] = scale(short_empty, args.short_time)
    
    def repl(short_, long_, wl, laser_range):
        # function to combine long and short wl
        low = np.max(np.argwhere(wl < laser_range[0]))
        high = np.min(np.argwhere(wl > laser_range[1]))

        combined = np.append(long_[0:low+1, 1], short_[(low+1):(high+1), 1])
        combined = np.append(combined, long_[(high+1):, 1])
        return combined

    if args.long_name != '':
        # process long files if existing
        print('Combining short and long measurements...')
        long_in = data[:, 8:10]
        long_out = data[:, 10:12]
        long_bckg = data[:, 12:14]
        long_empty = data[:, 14:16]
        
        long_in[:, 1] -= long_bckg[:, 1]
        long_out[:, 1] -= long_bckg[:, 1]
        long_empty[:, 1] -= long_bckg[:, 1]

        long_in[:, 1] = scale(long_in, args.long_time)
        long_out[:, 1] = scale(long_out, args.long_time)
        long_empty[:, 1] = scale(long_empty, args.long_time)

        _in = repl(short_in, long_in, wl, args.laser_range)
        _out = repl(short_out, long_out, wl, args.laser_range)
        _empty = repl(short_empty, long_empty, wl, args.laser_range)
    else:
        _in  = short_in[:, 1]
        _out = short_out[:, 1]
        _empty = short_empty[:, 1]
    
    # Apply calibration
    _in  = _in * cal
    _out = _out * cal
    _empty = _empty * cal

    spectrum = np.c_[wl, (_in - _out)]

    def inte(d, x, _range):
        # select data over a  within given wavelength range
        d_int = np.trapz(d[np.where((x > _range[0]) & (x < _range[1]))], x=x[np.where((x > _range[0]) & (x < _range[1]))])
        return d_int

    # The calibration gives a scaled spectrum with a response proportional to the number of photons, not the power 

    absorbed_full = 1 - inte(_in, wl, args.laser_range) / inte(_out, wl, args.laser_range)
    PL_full  = inte(_in, wl, args.pl_range) - inte(_empty, wl, args.pl_range) - (1 - absorbed_full) * (inte(_out, wl, args.pl_range) - inte(_empty, wl, args.pl_range))
    QE_full = PL_full / (inte(_empty, wl, args.laser_range) * absorbed_full)

    # Print results
    print('')
    print('RESULTS')
    print('-------')
    print('PLQY = {:.3f} %'.format(QE_full*100))
    print('Absorbtance = {:.2f} %'.format(absorbed_full*100))
    print('Absorbance = {:.2f} '.format(-np.log10(1-absorbed_full)))

    # Plot results
    fig = plt.figure(figsize=(11,8))
    gs = gridspec.GridSpec(2,2)

    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[:, 1])

    for ax in [ax1, ax2]:
        ax.semilogy(wl, _in, label='in')
        ax.semilogy(wl, _out, label='out')
        ax.semilogy(wl, _empty, label='empty')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Counts [a.u.]')
        ax.legend()

    ax1.set_xlim(args.laser_range[0]-15, args.laser_range[1] + 15)
    ax1.set_ylim(bottom = 1e-2)
    ax1.axvline(args.laser_range[0], linestyle='--', color='k')
    ax1.axvline(args.laser_range[1], linestyle='--', color='k')
    ax2.set_xlim(args.pl_range[0]-25, args.pl_range[1] + 25)
    ax2.axvline(args.pl_range[0], linestyle='--', color='k')
    ax2.axvline(args.pl_range[1], linestyle='--', color='k')

    ax3.plot(wl, _in - _empty, label='in')
    ax3.plot(wl, _out - _empty, label='out')
    ax3.set_xlabel('Wavelength (nm)')
    ax3.set_ylabel('Counts [a.u.]')
    ax3.set_xlim(args.pl_range[0]-25, args.pl_range[1] + 25)
    ax3.set_ylim(0, 1.2* np.max(_in[np.where((wl > args.pl_range[0]) & (wl < args.pl_range[1]))] - _empty[np.where((wl > args.pl_range[0]) & (wl < args.pl_range[1]))]))

    ax3.legend()
    ax3.annotate('PLQY = {0:.3f} % \nAbsorbance = {1:.2f}'.format(QE_full*100, -np.log10(1-absorbed_full)), xy=(0.05, 0.92), xycoords='axes fraction')

    return fig, spectrum

    
def loadit(args):
    # Function to load files
    short_in = np.loadtxt(args.directory + args.short_name + 'in.txt', delimiter='\t')
    short_out = np.loadtxt(args.directory + args.short_name + 'out.txt', delimiter='\t')
    
    if args.common == True:
        short_bckg = np.loadtxt(args.directory + args.common_bckg, delimiter='\t')
        short_empty = np.loadtxt(args.directory + args.common_empty, delimiter='\t')
    else:
        short_bckg = np.loadtxt(args.directory + args.short_name + 'bckg.txt', delimiter='\t')
        short_empty = np.loadtxt(args.directory + args.short_name + 'empty.txt', delimiter='\t')
    
    data = np.c_[short_in, short_out, short_bckg, short_empty]

    if args.long_name != '':
        long_in = np.loadtxt(args.directory + args.long_name + 'in.txt', delimiter='\t')
        long_out = np.loadtxt(args.directory + args.long_name + 'out.txt', delimiter='\t')
        
        if args.common == True:
            long_bckg = np.loadtxt(args.directory + args.common_long_bckg, delimiter='\t')
            long_empty = np.loadtxt(args.directory + args.common_long_empty, delimiter='\t')
        else:
            long_bckg = np.loadtxt(args.directory + args.long_name + 'bckg.txt', delimiter='\t')
            long_empty = np.loadtxt(args.directory + args.long_name + 'empty.txt', delimiter='\t')

        data = np.c_[data, long_in, long_out, long_bckg, long_empty]

    return data

def save_res(fig, spectrum, args):
    # saving results
    plt.savefig(args.directory + args.short_name + 'fig.pdf', format='pdf')
    np.savetxt(args.directory + args.short_name + 'spectrum.txt', spectrum, delimiter='\t', fmt='%.5f')



# Run functions
if __name__ == "__main__":
    args = get_args()
    fig, spectrum = PLQE(args)
    save_res(fig, spectrum, args)

plt.show(block=True)

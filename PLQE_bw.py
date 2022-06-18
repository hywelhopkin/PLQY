#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# written by bernard.wenger@heliodm.com

# Using the procedure of DeMello et al https://doi.org/10.1002/adma.19970090308

# History
# 06/08/2019, converted from Matlab
# 09/08/2019, first working version (no added functions)
# 22/09/2019, added QE pro spectrometer
# 27/09/2019, added subtraction of stray light
# 08/10/2019, added Gooey (see Github for further history)
# 16/4/2021, created a helio branch specific for use in Begbroke

from pathlib import Path, PurePath
from os import chdir
from gooey import Gooey, GooeyParser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lmfit.models import VoigtModel, ConstantModel
from scipy.constants import Planck, speed_of_light
import csv

# define calibration files
cal_dict = {
        "4 pi" : 'cal_PLQY_helio.txt',
        "fQY tube" : 'cal_fQY_tube_helio.txt',
        "fQY sphere" : 'cal_fQY_sphere_helio.txt',
         }

cfgs = ['4 pi', 'fQY tube', 'fQY sphere']

# Use flag --ignore-gooey if you want to use the command line
@Gooey(advanced=True,          # toggle whether to show advanced config or not 
       default_size=(800, 600),   # starting size of the GUI
       show_success_modal = False,
       return_to_config = True,
       navigation = "TABBED",
       tabbed_groups=True,
)

def get_args():
    """Get arguments and options"""
    parser = GooeyParser(description='Calculation of PL quantum efficiency.')

    req = parser.add_argument_group('Select directory and file', gooey_options={'columns': 1})
    req.add_argument('-sp', '--short_path', type=str, widget="FileChooser", help="Path to the '_in.txt' file (e.g. 'C:/folder_name/sub_folder/short_name_in.txt'", gooey_options={'wildcard':"'in' files (*_in.txt)|*_in.txt|" "All files (*.*)|*.*"})    
    req.add_argument('-fd', '--folder', type=Path , widget="DirChooser", help="Select folder containing files")
    req.add_argument('-st', '--short_time', default= 10, type=int, help="Integration time for short measurement in ms")
    req.add_argument('-c', '--common', action='store_true', help="Indicates that common background and empty files are used.")
    req.add_argument('-f', '--fQY', action='store_true', help="Check this for fQY mode.")
    req.add_argument('-nf', '--no_fig', action='store_true', help="Don't show figure")
    

    opt = parser.add_argument_group('optional arguments', gooey_options={'columns': 2})
    opt.add_argument('-lr', '--laser_range', nargs = 2, default = "440 460", type=int, help="Laser wavelength range in nm, Default = 440 460")
    opt.add_argument('-plr', '--pl_range', nargs = 2, default = "550 850", type=int, help="PL detection range in nm, Default = 550 850")
    opt.add_argument('-cfg', '--config', default='4 pi', widget="Dropdown", choices=cfgs, type=str, help="Fiber and spectrometer configurations. Choices: '4 pi', '2 pi', 'fQY' 'Maya red', 'Maya steel', 'QE red 25', 'QE red 200', 'QE steel 25' or 'QE steel 200', Default = '4 pi'")    
    opt.add_argument('-sl', '--stray_light', action='store_true', help="Removes stray light background. Default = 'False'")
    opt.add_argument('-cb', '--common_bckg', default= 'bckg.txt', type=str, help=" Name of the common background file. Default = 'bckg.txt'")
    opt.add_argument('-ce', '--common_empty', default= 'empty.txt', type=str, help=" Name of the common empty file. Default = 'empty.txt'")

    long_group = parser.add_argument_group("Using long integration time", gooey_options={'columns': 3})
    long_group.add_argument('-lt', '--long_time', default= 100, type=int, help="Integration time for long measurement in ms")
    long_group.add_argument('-lp', '--long_path', type=str, default='', widget="FileChooser", help="Path to the long '_in.txt' file (e.g. 'C:/folder_name/sub_folder/long_name_in.txt'", gooey_options={'wildcard':"'in' files (*_in.txt)|*_in.txt|" "All files (*.*)|*.*"})    
    long_group.add_argument('-cl', '--common_long', action='store_true', help="Indicates that common background and empty files are used. Default = True")
    long_group.add_argument('-clb', '--common_long_bckg', default= 'long_bckg.txt', type=str, help=" Name of the common long background file. Default = 'long_bckg.txt'")
    long_group.add_argument('-cle', '--common_long_empty', default= 'long_empty.txt', type=str, help=" Name of the common long empty file. Default = 'long_empty.txt'")

    args = parser.parse_args()

    # Need to find a way to treat these even if we choose the folder option.
    # To test if a string is empty one can use "if not myString:"

    # args.directory = Path(args.short_path).resolve().parent
    # args.short_name = Path(args.short_path).name
    # if (args.long_path != '') : args.long_name = Path(args.long_path).name
    # args.cwd = Path(__file__).resolve().parent
    return args

def PLQE(args):
    # Function running all calculations
    data = loadit(args) # Load data and calibration file
    cal = np.loadtxt(str(PurePath(args.cwd).joinpath('cal', cal_dict.get(args.config))))

    def trim(data, vr):
        # function used to select the valid range (could also be used for removing hot pixels)
        return data[vr[0]:vr[1], :]

    if 'QE' in args.config:
        vr = [4, -5] # valid range limits
        data = trim(data, vr)
    elif 'Maya' in args.config: 
        vr = [15, -6] # valid range limits
        data = trim(data, vr)
    else: 
        vr = [15, -6] # valid range limits
        data = trim(data, vr)
    
    cal = np.interp(data[:, 0], cal[:, 0], cal[:, 1])

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

    if args.long_path != '':
        # process long files if existing
        print('\nCombining short and long measurements...')
        if args.long_time == 0:
            print("\n!!! Don't forget to enter integration times !!!\n")
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
    
    # Apply calibration (generates spectra in µW/nm)
    _in  = _in * cal
    _out = _out * cal
    _empty = _empty * cal



    def inte(d, x, _range):
        # select data over a  within given wavelength range
        d_int = np.trapz(d[np.where((x > _range[0]) & (x < _range[1]))], x=x[np.where((x > _range[0]) & (x < _range[1]))])
        return d_int

    def remove_stray(_in, _out, _empty, wl, wl_range, pl_range):
        # remove background due to stray light
        def get_avg(x, wl, wl_range):
            idx = [np.argmin(abs(wl - wl_range[0])), np.argmin(abs(wl - wl_range[1]))]
            avg = np.trapz(x[idx[0]:idx[1]], wl[idx[0]:idx[1]])
            return avg
        
        # first take the average in a region below the laser peak
        avg_empty = get_avg(_empty, wl, wl_range)
        avg_out = get_avg(_out, wl, wl_range)
        avg_in = get_avg(_in, wl, wl_range)

        # Apply correction only in PL range
        def corr_pl_range(d, x, _empty, avg_data, avg_empty, _range):
            d_pl_range = d[np.where((x > _range[0]) & (x < _range[1]))]
            _empty_pl_range = _empty[np.where((x > _range[0]) & (x < _range[1]))]
            d_corr = (d_pl_range * (avg_empty / avg_data) - _empty_pl_range) * (avg_data / avg_empty)
            d[np.where((x > _range[0]) & (x < _range[1]))] = d_corr
            return d
        
        _in_corr = corr_pl_range(_in, wl, _empty, avg_in, avg_empty, pl_range)
        _out_corr = corr_pl_range(_out, wl, _empty, avg_out, avg_empty, pl_range)
        _empty_corr = corr_pl_range(_empty, wl, _empty, avg_empty, avg_empty, pl_range)
        
        return _in_corr, _out_corr, _empty_corr
    
    if args.stray_light == True:
        # Removes stray light background
        print('With stray light correction')
        wl_range = [370, 390]
        _in, _out, _empty = remove_stray(_in, _out, _empty, wl, wl_range, args.pl_range)
        
    spectra = np.c_[wl, _empty, _in, _out, (_in - _out)] # in µW/nm

    p2e = Planck * speed_of_light / (1e-9 * wl)  # photon to energy conversion J/photon
    
    # Applying the correction 'cal' gives a spectrum in W/nm, we have to convert into photons using p2e
    absorbed_full = 1 - inte(_in / p2e , wl, args.laser_range) / inte(_out / p2e, wl, args.laser_range)
    PL_full  = inte(_in / p2e, wl, args.pl_range) - inte(_empty / p2e, wl, args.pl_range) - (1 - absorbed_full) * (inte(_out / p2e, wl, args.pl_range) - inte(_empty / p2e, wl, args.pl_range))
    QE_full = PL_full / (inte(_empty / p2e, wl, args.laser_range) * absorbed_full)

    # Values in power scale (µW/nm)
    absorbed_W = 1 - inte(_in, wl, args.laser_range) / inte(_out, wl, args.laser_range) 
    PL_W  = inte(_in, wl, args.pl_range) - inte(_empty, wl, args.pl_range) - (1 - absorbed_W) * (inte(_out, wl, args.pl_range) - inte(_empty, wl, args.pl_range))
    PE_full = PL_W / (inte(_empty, wl, args.laser_range) * absorbed_W)


    OD = -np.log10(inte(_in, wl, args.laser_range) / inte(_out, wl, args.laser_range)) # this is in energy scale

    laser_power = np.trapz(_empty[(wl > args.laser_range[0]) & (wl < args.laser_range[1])], x=wl[(wl > args.laser_range[0]) & (wl < args.laser_range[1])])/1000 # in mW
    
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
        ax.set_xlabel('Wavelength [nm]')
        ax.set_ylabel('Intensity [µW/nm]')
        ax.legend()

    ax1.set_xlim(args.laser_range[0]-15, args.laser_range[1] + 15)
    ax1.set_ylim(bottom = 1e-4)
    ax1.axvline(args.laser_range[0], linestyle='--', color='k')
    ax1.axvline(args.laser_range[1], linestyle='--', color='k')
    ax1.annotate(f'Laser power: {laser_power:.2f} mW', xy=(0.1, 0.89), xycoords='axes fraction')
    ax2.set_xlim(args.pl_range[0]-25, args.pl_range[1] + 25)
    ax2.axvline(args.pl_range[0], linestyle='--', color='k')
    ax2.axvline(args.pl_range[1], linestyle='--', color='k')


    ax3.plot(wl, _in - _empty, label='in')
    ax3.plot(wl, _out - _empty, label='out')
    ax3.set_xlabel('Wavelength [nm]')
    ax3.set_ylabel('Intensity [µW/nm]')
    ax3.set_xlim(args.pl_range[0]-25, args.pl_range[1] + 25)
    ax3.set_ylim(0, 1.2* np.max(_in[np.where((wl > args.pl_range[0]) & (wl < args.pl_range[1]))] - _empty[np.where((wl > args.pl_range[0]) & (wl < args.pl_range[1]))]))

    ax3.legend()
    center, fwhm = fit_voigt(ax3, spectra, args)

    ax3.annotate(
        f'PLQY = {QE_full*100:.2f} % \n'
        f'OD = {-np.log10(1-absorbed_full):.2f} \n'
        f'Peak center = {center:.2f} nm \n'
        f'FWHM = {fwhm:.2f} nm',
        xy=(0.05, 0.89), xycoords='axes fraction')

    # Print results
    string_pc = f'<html><table><tr><td>{QE_full*100:.2f}</td><td>{OD:.2f}</td><td>{center:.2f}</td><td>{fwhm:.2f}</td><td>{laser_power:.2f}</td></tr></table></html>'

    print(
        f'\nRESULTS ({args.short_name})\n'
        f'-------\n'
        f'PLQY = {QE_full*100:.2f} %\n'
        f'PE = {PE_full*100:.2f} %\n'
        f'OD = {OD:.2f}\n'
        f'Peak WL = {center:.2f} nm\n'
        f'FWHM = {fwhm:.2f} nm\n'
        f'Laser power = {laser_power:.2f} mW\n'
        # f'\n/// string to copy-paste for Excel ////\n'
        # f'{string_pc}'
    )
    # pc.copy(string_pc) # string to Excel is copied in the clipboard

    to_log = [args.short_name, f'{QE_full*100:.3f}', f'{PE_full*100:.3f}', f'{OD:.3f}', f'{center:.2f}', f'{fwhm:.2f}', f'{laser_power:.2f}']

    return fig, spectra, to_log

def loadit(args):
    # Function to load files
    chdir(args.directory)
    if args.common == True:
        try:
            short_bckg = np.loadtxt(args.common_bckg, delimiter='\t')
            short_empty = np.loadtxt(args.common_empty, delimiter='\t')
        except:
            short_bckg = np.loadtxt(args.common_bckg, delimiter='\t', skiprows=14)
            short_empty = np.loadtxt(args.common_empty, delimiter='\t', skiprows=14)
    else:
        short_bckg = np.loadtxt(str(args.short_name).replace('in.txt', 'bckg.txt'), delimiter='\t')
        short_empty = np.loadtxt(str(args.short_name).replace('in.txt', 'empty.txt'), delimiter='\t')

    try:
        short_in = np.loadtxt(args.short_name, delimiter='\t')
        if args.fQY == False: 
            short_out = np.loadtxt(str(args.short_name).replace('in.txt', 'out.txt'), delimiter='\t')
        else:
            short_out = short_empty.copy()
    except:
        short_in = np.loadtxt(args.short_name, delimiter='\t', skiprows=14)
        if args.fQY == False: 
            short_out = np.loadtxt(str(args.short_name).replace('in.txt', 'out.txt'), delimiter='\t', skiprows=14)
        else:
            short_out = short_empty.copy()

    data = np.c_[short_in, short_out, short_bckg, short_empty]

    if args.long_path != '':
        if args.common_long == True:
            try:
                long_bckg = np.loadtxt(args.common_long_bckg, delimiter='\t')
                long_empty = np.loadtxt(args.common_long_empty, delimiter='\t')
            except:
                long_bckg = np.loadtxt(args.common_long_bckg, delimiter='\t', skiprows=14)
                long_empty = np.loadtxt(args.common_long_empty, delimiter='\t', skiprows=14)
        else:
            long_bckg = np.loadtxt(str(args.long_name).replace('in.txt', 'bckg.txt'), delimiter='\t')
            long_empty = np.loadtxt(str(args.long_name).replace('in.txt', 'empty.txt'), delimiter='\t')

        try:
            long_in = np.loadtxt(args.long_name, delimiter='\t')
            if args.fQY == False: 
                long_out = np.loadtxt(str(args.long_name).replace('in.txt', 'out.txt'), delimiter='\t')
            else:
                long_out = long_empty.copy()
        except:
            long_in = np.loadtxt(args.long_name, delimiter='\t', skiprows=14)
            if args.fQY == False: 
                long_out = np.loadtxt(str(args.long_name).replace('in.txt', 'out.txt'), delimiter='\t', skiprows=14)
            else:
                long_out = long_empty.copy()

        data = np.c_[data, long_in, long_out, long_bckg, long_empty]

    return data

def fit_voigt(ax, spectra, args):
    fit_range = args.pl_range
    wl = spectra[:, 0]
    sp = spectra[:, 2]

    wl_fit = wl[(wl>fit_range[0]) & (wl<fit_range[1])]
    sp_fit = sp[(wl>fit_range[0]) & (wl<fit_range[1])]


    mod = VoigtModel() + ConstantModel()
    pars = mod.make_params(amplitude=np.max(sp_fit), center=wl_fit[np.argmax(sp_fit)], 
                            sigma=10, gamma=10, c=0)

    out = mod.fit(sp_fit, pars, x=wl_fit)

    ax.plot(wl_fit, out.best_fit, 'k--', alpha=0.8)
    
    return out.params['center'].value, out.params['fwhm'].value, 

def save_res(fig, spectra, args):
    # saving results
    plt.savefig(str(PurePath(args.directory).joinpath(str(args.short_name).replace('in.txt', 'fig.pdf'))), format='pdf')
    spectra_header = 'Wavelength\t empty\t in\t out\t proc'
    np.savetxt(str(PurePath(args.directory).joinpath(str(args.short_name).replace('in.txt', 'spectra.txt'))), spectra, delimiter='\t', fmt='%.5e', header=spectra_header, comments='')

def save_log(args, to_log):
    logfile = Path(args.directory).joinpath('log_plqy.csv')
    header = ['Name', 'PLQY (%)', 'PE (%)', 'OD (-)', 'Peak WL (nm)', 'FWHM (nm)', 'Laser power (mW)']

    if not logfile.is_file():
        with open(logfile, 'w', newline='') as f:
            f_write = csv.writer(f, dialect='excel')
            f_write.writerow(header)

    with open(logfile, 'a', newline='') as f:
        f_write = csv.writer(f, dialect='excel')
        f_write.writerow(to_log)

# Run functions
if __name__ == "__main__":
    args = get_args()
    # # if (args.folder) check if string is empty. If it is continue as normal with args.short_name. If not run a loop where args.short_name is changed at every iteration
    for in_file in Path(args.folder).glob('*in  .txt'):
        print(in_file)
    # # do something with "txt_file
    # fig, spectra, res_to_log = PLQE(args)
    # save_res(fig, spectra, args)
    # save_log(args, res_to_log)

if args.no_fig == False:
    plt.show()

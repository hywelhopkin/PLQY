#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# written by bernard.wenger@heliodm.com
# History
# 6/8/2019, v1.0, converted from Matlab

import sys
import numpy as np
import matplotlib.pyplot as plt

# Only change the settings in these top two sections
def PLQE_bw(short_name, folder='.', fiber='red', laser_range=[397, 407], pl_range=[450, 850], short_time=1, long_name='', long_time=0, common=True, common_bckg='bckg.txt', common_empty='empty.txt'):
    """
    Calculation of PL quantum efficiency.

    Parameters
    ----------
    folder : str, optional
        folder in which the raw files are stored (e.g. 'folder_name/sub_folder/'
        Default = '.'
    fiber : str, optional
        Fiber used. 'red' or 'silver'
        Default = 'red'
    laser_range : list of integer, optional
        Wavelength range for the excitation laser (e.g [397, 407]) 
        Default = [397, 407]
    pl_range : list of integer, optional
        Wavelength range used to calculate the emission (e.g [550, 750])
        Default = [450, 850]
    short_name : str, required
        Prefix of the file names (e.g 'filename_') to be followed by 'in.txt', 'out_.txt', 'bckg.txt' or 'empty.txt'
    short_time : int, optional
        Integration time for short measurement in ms
        Default = 1
    long_name : str, optional
        Prefix of the long file name (e.g 'filename_long_') to be followed by 'in.txt', 'out_.txt', 'bckg.txt' or 'empty.txt'
        Default = ''
    long_time : int, optional
        Integration time for long measurement in ms
        Default = 0
    common : bool, optional
        Indicates that common background and empty files are used
        Default = True
    common_bckg : str, optional
        Name of the common background file
        Default = 'bckg.txt'
    common_empty : str, optional
        Name of the common empty file
        Default = 'empty.txt'

    Returns
    -------
    tbdout : ndarray
        Data read from the text file.

    Examples
    --------
    >>> write down an example
    >>> 
    """

    print('There you go')
    print(short_name)

#loadit = @(fn)(importdata(['perovTest_100ms_',fn])');
# loadit_long = @(fn)(dlmread([data_folder,long_name,fn],'\t', 0, 0)');
# loadit_short = @(fn)(dlmread([data_folder,short_name,fn],'\t', 0, 0)');
# loadit_generic = @(fn)(dlmread([data_folder,fn],'\t', 0, 0)');

# if strcmp(common_bckg,'') == 0
#     short_bckg  = loadit_generic(common_bckg);
# else
#     short_bckg    = loadit_short('bckg.txt');
# end
# short_in    = loadit_short('in.txt');
# short_out   = loadit_short('out.txt');
# if strcmp(common_empty,'') == 0
#     short_empty = loadit_generic(common_empty);
# else
#     short_empty = loadit_short('empty.txt');
# end

# if nargin < 4
#     long_bckg = short_bckg;
#     long_empty = short_empty;
#     long_in = short_in;
#     long_out = short_out;
    
# else
# #     long_bckg   = loadit_generic('long_bckg.txt');
#     long_bckg   = loadit_long('bckg.txt');
#     long_in     = loadit_long('in.txt');
#     long_out    = loadit_long('out.txt');
# #     long_empty  = loadit_generic('long_empty.txt');
#     long_empty  = loadit_long('empty.txt'); 
# end


# ## Select Fiber calibration file

# if fiber == 1 # red fiber
#     cal_file ='SpecResp_MayaPro_200umRedFiber_WavelengthCorrected_Updated_16052014';
# elseif fiber == 2 # silver fiber
#     cal_file = 'SpecResp_MayaPro_50umSteelFiber_WavelengthCorrected_Updated_16052014';
# else
#     str = 'Select fiber: 1 (red) or (2) silver'
#     return
# end

# ## BackgroundSubtract

# short_in(2,:) = short_in(2,:)-short_bckg(2,:);
# short_out(2,:) = short_out(2,:)-short_bckg(2,:);
# short_empty(2,:) = short_empty(2,:)-short_bckg(2,:);

# long_in(2,:) = long_in(2,:)-long_bckg(2,:);
# long_out(2,:) = long_out(2,:)-long_bckg(2,:);
# long_empty(2,:) = long_empty(2,:)-long_bckg(2,:);




# ## Process
# proc   = @(d,time)((d(2,:)-mean(d(2,1:24)))./time);

# x = short_in(1,:);

# short_in       = proc(short_in,short_time);
# short_out      = proc(short_out,short_time);
# short_empty    = proc(short_empty,short_time);

# long_in       = proc(long_in,long_time);
# long_out      = proc(long_out,long_time);
# long_empty    = proc(long_empty,long_time);

# ## Combine
# repl   = @(x,high,low,range)([low(1:search(x,range(1))),high(search(x,range(1))+1:search(x,range(2))),low(search(x,range(2))+1:end)]);

# in    = repl(x,short_in,long_in,laser_range);
# out   = repl(x,short_out,long_out,laser_range);
# empty = repl(x,short_empty,long_empty,laser_range);

# ## Calibrate
# #load('IntSphere_calib_532Laser_nofilt','cal');
# load(cal_file,'cal');
# in    = in.*cal;
# out   = out.*cal;
# empty = empty.*cal;
# # trick to remove background
# # out(538:end) = mean(in(1200:1250));
# # empty(538:end) = mean(in(1200:1250));

# ## Calculate quantum efficiency (naive approach, no reabsorbtion)

# inte = @(x,y,n)(sum(y(search(x,n(1)):search(x,n(2))))); 
# #load('PowerCalibration','PowCal');
# #Power = inte(x,short_empty.*cal,laser_range).*PowCal(p)/0.16;
# #PowerText = [num2str(Power,3),'mW/cm^2']

# absorbed = [0,1-(inte(x,in,laser_range)./(inte(x,out,laser_range)))]; # Fraction of photons absorbed
# PL       = [0,inte(x,in,pl_range)-inte(x,empty,pl_range)];   #PL (not removing PL from reabsorption)

# QE = PL(2)./(inte(x,empty,laser_range).*absorbed(2));
# QEText = [num2str(QE*100,3),'#'];

# # From Adv Mater, 9, 230 (1997) (including reabsorption)

# absorbed_full = [0,1-(inte(x,in,laser_range)./(inte(x,out,laser_range)))];   # Fraction of photons absorbed on first pass
# PL_full       = [0,inte(x,in,pl_range)-inte(x,empty,pl_range)-(1-absorbed_full(2)).*(inte(x,out,pl_range)-inte(x,empty,pl_range))];  # PL from first pass

# QE_full = PL_full(2)./(inte(x,empty,laser_range).*absorbed_full(2));
# QEText_full = [num2str(QE_full*100,3),'#'];

# # print results
# fprintf('\nRESULTS\n')
# fprintf('-------\n')
# fprintf('PLQY = #.3f ## (#.3f ## w/o reabs) \n', QE_full*100, QE*100)
# fprintf('Absorbtance = #.2f ## \n', absorbed_full(2)*100)
# fprintf('Absorbance = #.2f \n', - log10(1 - absorbed_full(2)))

# ##Giles' export function
# fileID = fopen(['PLQE_log.txt'], 'a');
# fwrite(fileID,[13 10],'char');
# fprintf(fileID,[short_name char(9)]);
# fprintf(fileID,QEText_full);
# fclose(fileID);

# #save spectra
# dlmwrite([data_folder, short_name, 'spectrum.txt'], [x.' (in-out).'], '\t');

# ## Plotting
# clf;
# subplot(221)
# title('All data');
# warning('off','MATLAB:Axes:NegativeDataInLogAxis')
# semilogy(x,in,x,out,x,empty);
# axis([380 1000 5e-2 2*max(empty)]);
# legend('In','Out','Laser');
# xlabel('Wavelength (nm)');
# ylabel('PL (counts/sec)');
# title('All PL (calibrated)');
# subplot(223);
# plot(x,in-empty,x,out-empty);
# axis([380 900 0 max(in(548:1582)-empty(548:1582))]);
# title('Just PL (empty baseline is subtracted)');
# legend('In','Out', 'Location', 'NorthWest');
# xlabel('Wavelength (nm)');
# ylabel('PL (counts/sec) - corrected');
# subplot(122);
# plot(absorbed,PL,'d:');
# xlabel('Absorbed photons');
# ylabel('PL photons');
# legend(['  PLQE = ' QEText_full],'location','best');
# #legend(['Power =', PowerText, '  PLQE = ' QEText_full],'location','best');
# hgexport(gcf, [data_folder long_name '.jpeg'],hgexport('factorystyle'),'Format', 'jpeg')

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    PLQE_bw(*sys.argv[1:])
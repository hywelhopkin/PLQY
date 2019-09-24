import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from scipy.io import loadmat


cal_dict = {
        "Maya red" : 'SpecResp_MayaPro_200umRedFiber_WavelengthCorrected_Updated_16052014',
        "Maya steel" : 'SpecResp_MayaPro_50umSteelFiber_WavelengthCorrected_Updated_16052014',
        "Maya red BW" : 'cal_HL3_red_sphere.txt',
        "QE red 25"  : 'cal_QEpro_red_25um.txt',
        "QE red 200" : 'cal_QEpro_red_200um.txt',
        "QE steel 25" : 'cal_QEpro_steel_25um.txt',
        "QE steel 200" : 'cal_QEpro_steel_200um.txt',
    }

def get_args():
    """Get arguments and options"""
    parser = argparse.ArgumentParser(description='Calculation of PL quantum efficiency.')
    
    parser.add_argument('-cfg', '--config', default = 'Maya Pro', type=str, required=True, help="Prefix of the file names (e.g 'filename_') to be followed by 'in.txt', 'out.txt', 'bckg.txt' or 'empty.txt'")
    
    args = parser.parse_args()
    
    return args

def test(args):
    print(args.config)
    cal = np.loadtxt(cal_dict.get(args.config))
    print(cal.shape)

# Run functions
if __name__ == "__main__":
    args = get_args()
    test(args)


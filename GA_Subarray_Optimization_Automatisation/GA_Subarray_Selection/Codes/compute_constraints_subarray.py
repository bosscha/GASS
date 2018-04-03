#! /usr/bin/env python2
# -*- coding: utf-8 -*-
__author__="Ignacio Toledo"
__version__="0.1"

'''
    File name: compute_constraints_subarray.py
    Author: Ignacio Toledo 
    Modified by : Roxane Lassis 
    Description : Functions to compute the real constraint values of a subarray
    Context : ALMA internship
    Date created: 06/2017
    Date last modified: 08/2017
    Python Version: 2.7
'''

import sys
sys.path.insert(0,'/users/sleon/python/astropy-2.0')
import astropy

import numpy as np
import pandas as pd
import time

from matplotlib import pyplot as plt
from math import *
from astropy import modeling
from astropy.stats import gaussian_sigma_to_fwhm


#==================================================================================================================================
# Useful functions to compute constraints values of a subarray :   
#==================================================================================================================================
   


def calc_baselines(array,num_subarray, lat=radians(-23.0262015)):

    """
    calculate baseline coordinates in geocentrical coordinates, starting from a
    subarray and given the
    site latitude (by default Alma's)
    
    cfg_file : file with the configuration (position 3d of the pads, type of antenna, name of the pads)
    lat : site latitude
    """
    coordinates=[[],[],[]]
    for j in range(0,array.configuration_Manager.get_Number_Pads_Per_Subarray()[num_subarray]):
            pad=array.get_Pad(num_subarray,j)
            coordinates[0].append(pad.get_E())
            coordinates[1].append(pad.get_N())
            coordinates[2].append(pad.get_U())
                

    rot = np.array(
        [[0, -1 * sin(lat), cos(lat)],
         [1, 0, 0],
         [0, cos(lat), sin(lat)]])

    xx, yy, zz = np.dot(
        rot, coordinates)
        
    # Initialize N(N-1)/2 baselines with N antennas
    N=array.configuration_Manager.get_Number_Pads_Per_Subarray()[num_subarray]
    b = np.zeros([N * (N-1) / 2, 3]) 
    c = 0
    for i in range(0, N):
        for j in range(i + 1, N):
            b[c][0] = xx[i] - xx[j]
            b[c][1] = yy[i] - yy[j]
            b[c][2] = zz[i] - zz[j]
            c += 1

    return b


def calc_uv_coverage(ha, dec, b, t=1):
    """
    Calculate the uv coverage using the B matrix created with calc_b,
    given a source hour angle in hours, dec in degrees and a time duration.
    :param ha : hour angle 
    :param dec : declination 
    :param b : baselines computed with calc_b
    :param t : integration time
    :return : u and v coordinates of the sample uv plane
    """
    ha *= 15
    dec = radians(dec)
    UU = np.zeros(0)
    VV = np.zeros(0)

    for ht in np.radians(np.arange(ha, ha+(t*15.)+1.5, 1.5)):
        rot = np.array(
            [[sin(ht), cos(ht), 0],
             [-1 * sin(dec) * cos(ht), sin(dec) * sin(ht), cos(dec)],
             [cos(dec) * cos(ht), -1 * cos(dec) * sin(ht), sin(dec)]])
        uu, vv, ww = np.dot(rot, b.transpose())
        UU = np.concatenate([UU, uu])
        VV = np.concatenate([VV, vv])

    return UU, VV
    

#@numba.jit
def calc_constraints(ha, dec, baselines,t=1, doplot=False, display_values=False, px=512):
    """
    Calculate the beam angular resolution expected from an array with baselines
    described by a B array (created with calc_b), for a source at a given
    dec in degrees and hour angle in hours at the beginning, assuming a
    observing time of one hour by default.

    A briggs weighting with factor 0.5 is used to grid the uv plane.

    :param ha : hour angle
    :param dec : declination of the source
    :param b : baselines computed with calc_b
    :param t : integration time
    :param doplot : boolean to plot or no
    :param px : 
    :return : ar, the resolution
    """

    uv = calc_uv_coverage(ha, dec, baselines, t=t)
    l_ax = 2. * max(np.abs(uv[0]).max(), np.abs(uv[1]).max())
    dx = l_ax * 2. / px
    bin_ax = np.arange(-2 * l_ax, 2 * l_ax, dx) 
    ds = dx / (299792458.0 / (100. * 1e9))
    npix = len(bin_ax)
    H, xd, yd = np.histogram2d(uv[0], uv[1], bins=(bin_ax, bin_ax))

    # Apply briggs weight with robust factor 0.5
    f2 = (5*10**0.5)**2. * (((1/H[H>0])**2.).sum() / H.sum())
    H = H / (1. + H * f2)

    # Apply fourier transform
    h = np.fft.fft2(H)
    hshift = np.abs(np.fft.fftshift(h).real)
    dr = degrees(1. / (npix * ds)) * 3600.
    cent = npix / 2
    size = 50
    coord_a = cent - size / 2
    coord_b = cent + size / 2
    x, y = np.mgrid[:size, :size]
    hshiftcut = hshift.transpose()[coord_a:coord_b, coord_a:coord_b]
    p_init = modeling.models.Gaussian2D(
        amplitude=hshiftcut.max(), x_mean=size/2, y_mean=size/2,
        x_stddev=1, y_stddev=1)
    fit_p = modeling.fitting.LevMarLSQFitter()
    p = fit_p(p_init, x, y, hshiftcut, maxiter=100)
    
    #Compute the sidelobes
    sidelobe_max_value=hshiftcut[(hshiftcut-p(x,y))==np.max(hshiftcut-p(x,y))][0]
    beam_max_value=np.max(np.abs(hshiftcut))
    sidelobe_percentage=sidelobe_max_value*100/beam_max_value
    
    if display_values: 
        print 'sidelobe', sidelobe_max_value
        print 'beam', beam_max_value

    if doplot:
        # Figure of the computed beam 
        plt.figure()
        plt.imshow(hshiftcut)
        CS = plt.contour(hshiftcut)#p(x, y))
        plt.clabel(CS)
        plt.colorbar()
        # Figure of the theoric beam
        plt.figure()
        plt.imshow(p(x,y))
        CS = plt.contour(p(x,y))
        plt.colorbar()
        # Figure of the difference
        plt.figure()
        plt.imshow(hshiftcut - p(x,y))
        plt.colorbar()
        plt.show()

    # Axis of the beam
    a=p.y_stddev.value*(astropy.stats.gaussian_sigma_to_fwhm * dr)
    b=p.x_stddev.value*(astropy.stats.gaussian_sigma_to_fwhm * dr)
    aa=max(a,b)
    bb=min(a,b)
    # Compute the elongation
    elongation=aa/bb
    # Compute the Maximal Recoverable Scale for a frequency of 100Ghertz (ie. wave length=3mm)
    distances=np.zeros(baselines.shape[0])
    for j in range(0,baselines.shape[0]):
        distances[j]=baselines[j][0]**2+baselines[j][1]**2+baselines[j][2]**2
    Lmin=np.sqrt(min(distances))
    #mrs=(0.6*(3*10**(-3))/Lmin)*(3600*180/(np.pi))
    mrs=37100./(Lmin*100)
    # Compute the resolution
    res = (np.sqrt(np.abs(p.x_stddev.value * p.y_stddev.value)) *
          astropy.stats.gaussian_sigma_to_fwhm * dr)
          
    return res,mrs,elongation,sidelobe_percentage

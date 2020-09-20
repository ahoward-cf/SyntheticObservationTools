#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 13:28:02 2020

@author: c1649794
"""

import numpy as np
from astropy.io import fits
from reproject import reproject_exact
from astropy.convolution import convolve_fft

class fitsConvolve(object):
    '''
    A class to read in, regrid, and convolve fits images with beam profiles
    '''
    
    def __init__(self, filename, dirname = './', beampath = './Beams/'):
        self.dirname = dirname
        if type(filename) == list:
            self.filenameList = filename
            self.fitsList = []
            for f in filename:
                current_fits = fits.open(dirname + f)[0]
                self.fitsList.append(current_fits)
        elif type(filename) == str:
            self.filenameSingle = filename
            self.fitsSingle = fits.open(dirname + f)[0]
        else:
            raise TypeError('Type of filename not recognised. filename must be list or str type')
        self.beampath = beampath
            
    def get_wavelength(self, wavelength = None):
        if hasattr(self, 'fitsList'):
            self.wvlnList = []
            for f in self.fitsList:
                self.wvlnList.append(f.header['WVLENGTH'])
        elif hasattr(self, 'fitsSingle'):
            self.wvlnSingle = self.fitsSingle.header['WVLENGTH']
        else:
            raise AttributeError('No images detected. Unable to generate wavelength list.')
            
    def generate_new_header(self, oldHead, beampixel):
        newHead = oldHead.copy()
        newHead['CDELT1'] = -beampixel / 3600.
        newHead['CDELT2'] = beampixel / 3600.
        newHead['NAXIS1'] = int(oldHead['NAXIS1'] * (oldHead['CDELT1'] / newHead['CDELT1']))
        newHead['NAXIS2'] = int(oldHead['NAXIS2'] * (oldHead['CDELT2'] / newHead['CDELT2']))
        newHead['CRPIX1'] = newHead['NAXIS1'] / 2.
        newHead['CRPIX2'] = newHead['NAXIS2'] / 2.
        return newHead
            
    def regrid_to_beam(self, beampixel = 1):
        if hasattr(self, 'fitsList'):
            self.beamReprojectedList = []
            for f in self.fitsList:
                oldHead = f.header
                oldData = f.data
                newHead = self.generate_new_header(oldHead, beampixel)
                newData = reproject_exact((oldData, oldHead), newHead, parallel = False)[0]
                HDUList = fits.HDUList(fits.ImageHDU(newData, newHead))[0]
                self.beamReprojectedList.append(HDUList)
        elif hasattr(self, 'fitsSingle'):
            oldHead = self.fitsSingle.header
            oldData = self.fitsSingle.data
            newHead = self.generate_new_header(oldHead, beampixel)
            newData = reproject_exact((oldData, oldHead), newHead, parallel = False)[0]
            HDUList = fits.HDUList(fits.ImageHDU(newData, newHead))[0]
            self.beamReprojectedImage = HDUList
        else:
            raise AttributeError('No images detected. Unable to reproject to beam grid.')
            
    def convolve_with_beams(self, pad = True):
        if not hasattr(self, 'wvlnList') or hasattr(self, 'wvlnSingle'):
            self.get_wavelength()
        if hasattr(self, 'beamReprojectedList'):
            self.beamConvolvedList = []
            for f in range(len(self.beamReprojectedList)):
                currentWvln = self.wvlnList[f]
                currentBeam = fits.getdata(self.beampath + 'psf_{:04d}.fits'.format(currentWvln))
                currentData = self.beamReprojectedList[f].data
                if pad:
                    currentData = np.pad(currentData, 50, 'edge')
                currentHead = self.beamReprojectedList[f].header
                convolvedImage = convolve_fft(currentData, currentBeam, normalize_kernel = True, allow_huge = True)
                if pad:
                    convolvedImage = convolvedImage[50:-50,50:-50]
                HDUList = fits.HDUList(fits.ImageHDU(convolvedImage, currentHead))[0]
                self.beamConvolvedList.append(HDUList)
        elif hasattr(self, 'beamReprojectedImage'):
            currentWvln = self.wvlnSingle
            currentBeam = fits.getdata(self.beampath + 'psf_{:04d}.fits'.format(currentWvln))
            currentData = self.beamReprojectedImage.data
            if pad:
                currentData = np.pad(currentData, 50, 'edge')
            currentHead = self.beamReprojectedImage.header
            convolvedImage = convolve_fft(currentData, currentBeam, normalize_kernel = True, allow_huge = True)
            if pad:
                convolvedImage = convolvedImage[50:-50,50:-50]
            HDUList = fits.HDUList(fits.ImageHDU(convolvedImage, currentHead))[0]
            self.beamConvolvedImage = HDUList
        else:
            raise AttributeError('No reprojected images detected. Unable to convolve with beams.')            
            
    def reproject_convolved_to_original(self):
        if hasattr(self, 'beamConvolvedList'):
            self.convolvedOriginalList = []
            for f in range(len(self.beamConvolvedList)):
                currentHead = self.beamConvolvedList[f].header
                currentData = self.beamConvolvedList[f].data
                newHead = self.fitsList[f].header
                newData = reproject_exact((currentData, currentHead), newHead, parallel = False)[0]
                HDUList = fits.HDUList(fits.ImageHDU(newData, newHead))[0]
                self.convolvedOriginalList.append(HDUList)
        elif hasattr(self, 'beamReprojectedImage'):
            currentHead = self.beamConvolvedImage.header
            currentData = self.beamConvolvedImage.data
            newHead = self.fitsSingle.header
            newData = reproject_exact((currentData, currentHead), newHead, parallel = False)[0]
            HDUList = fits.HDUList(fits.ImageHDU(newData, newHead))[0]
            self.convolvedOriginalSingle.append(HDUList)
        else:
            raise AttributeError('No convolved images detected. Unable to reproject.')
            
    def save_convolved_to_fits(self):
        if hasattr(self, 'convolvedOriginalList'):
            for f in range(len(self.convolvedOriginalList)):
                HDU = self.convolvedOriginalList[f]
                outName = self.dirname + 'conv_' + self.filenameList[f]
                fits.writeto(outName, HDU.data, header = HDU.header, overwrite = True)
        elif hasattr(self, 'convolvedOriginalSingle'):
            HDU = self.convolvedOriginalSingle
            outName = self.dirname + 'conv_' + self.filenameSingle
            fits.writeto(outName, HDU.data, header = HDU.header, overwrite = True)
        else:
            raise AttributeError('No convolved images detected. Unable to save out files.')

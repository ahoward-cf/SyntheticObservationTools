#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 13:01:29 2020

@author: c1649794
"""

import numpy as np
from astropy.io import fits

class FitsAddNoise(object):
    '''
    A class to add noise to Fits images
    '''
    
    def __init__(self, filename, dirname = './'):
        self.dirname = dirname
        if type(filename) == list:
            self.filenameList = filename
            self.fitsList = []
            for f in filename:
                current_fits = fits.open(dirname + f)[0]
                self.fitsList.append(current_fits)
        elif type(filename) == str:
            self.filenameSingle = filename
            self.fitsSingle = fits.open(dirname + filename)[0]
        else:
            raise TypeError('Type of filename not recognised. filename must be list or str type')
            
    def add_noise (self, SNR = 300):
        self.SNR = SNR
        if (hasattr(self, 'fitsList') and
            ((type(self.SNR) == int or
            type(self.SNR) == float) or
            (type(self.SNR) == list and
            len(self.fitsList) == len(self.SNR)))):
            self.noisyList = []
            for f in range(len(self.fitsList)):
                currentHead = self.fitsList[f].header.copy()
                currentData = self.fitsList[f].data.copy()
                if type(self.SNR) == list:
                    currentSNR = self.SNR[f]
                else:
                    currentSNR = self.SNR
                noise = np.nanmax(currentData) / currentSNR
                noisyHead = currentHead.copy()
                noisyHead['NOISE'] = noise
                currentNoise = np.random.normal(0, noise, currentData.shape)
                noisyData = currentData + currentNoise
                HDUList = fits.HDUList(fits.ImageHDU(noisyData, noisyHead))[0]
                self.noisyList.append(HDUList)
        elif (hasattr(self, 'fitsSingle') and 
            (type(self.SNR) == int or
             type(self.SNR) == float)):
            currentHead = self.fitsSingle.header.copy()
            currentData = self.fitsSingle.data.copy()
            currentSNR = self.SNR
            noise = np.nanmax(currentData) / currentSNR
            noisyHead = currentHead.copy()
            noisyHead['NOISE'] = noise
            currentNoise = np.random.normal(0, noise, currentData.shape)
            noisyData = currentData + currentNoise
            HDUList = fits.HDUList(fits.ImageHDU(noisyData, noisyHead))[0]
            self.noisySingle = HDUList
        else:
            raise ValueError('Length of SNR and image list do not match. SNR must be list of length equal to image list, or single value.')
    
    def save_noisy_to_fits(self):
        if hasattr(self, 'noisyList'):
            for f in range(len(self.noisyList)):
                HDU = self.noisyList[f]
                outName = self.dirname + 'noisy_' + self.filenameList[f]
                fits.writeto(outName, HDU.data, header = HDU.header, overwrite = True)
        elif hasattr(self, 'noisySingle'):
            HDU = self.noisySingle
            outName = self.dirname + 'noisy_' + self.filenameSingle
            fits.writeto(outName, HDU.data, header = HDU.header, overwrite = True)
        else:
            raise AttributeError('No noisy images detected. Unable to save out files.')
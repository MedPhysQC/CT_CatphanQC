#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This code is developed to be used as an analysis module within the WADQC software
# It is basically a wrapper to use the Catphan module developed within the pylinac package which is a requirement for this wrapper to work
# 
#
#
# The WAD Software can be found on https://bitbucket.org/MedPhysNL/wadqc/wiki/Home
#  
#
# Changelog:
#
#


__version__='20181116'
#__author__ = 'DD, tdw'



import sys,os

try:
    import pydicom as dicom
    from pydicom import tag

except ImportError:
    import dicom
    from dicom import tag


import getopt
import numpy as np
from numpy import random as rnd
try:
    # this will fail unless wad_qc is already installed
    from wad_qc.module import pyWADinput
except ImportError:
    import sys
    # add root folder of WAD_QC to search path for modules
    _modpath = os.path.dirname(os.path.abspath(__file__))
    while(not os.path.basename(_modpath) == 'Modules'):
        _new_modpath = os.path.dirname(_modpath)
        if _new_modpath == _modpath:
            raise
        _modpath = _new_modpath
    sys.path.append(os.path.dirname(_modpath))
    from wad_qc.module import pyWADinput

from wad_qc.modulelibs import wadwrapper_lib
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.


import pylinac
import os


def _addBoundary(pathtodicom,enlargepixels=25):
    'This routine manipulates the Catphan series by adding a number of background rows and columns to make the  detection of the phantom easier for pylinac'
    flist = []
    ew = enlargepixels
    for root, subFolders, files in os.walk(pathtodicom):
        for tfile in files:
            tmpfile = os.path.join(root,tfile)
            
            tmpdcm =  dicom.read_file(tmpfile)
        
            tmparray = tmpdcm.pixel_array
            arrshape = np.shape(tmparray)
        
            newarray = np.zeros((arrshape[0] + 2*ew,arrshape[1] + 2*ew))
            newarray[ew:ew+arrshape[0],ew:ew+arrshape[1]] = tmparray

            if newarray.dtype != np.uint16:
                newarray = newarray.astype(np.uint16)
                tmpdcm.PixelData = newarray.tostring()

        
            tag = 'Rows'
            if tag in tmpdcm:
                tmpdcm.data_element(tag).value  = arrshape[0]+2*ew

            tag = 'Columns'
            if tag in tmpdcm:
                tmpdcm.data_element(tag).value  = arrshape[1]+2*ew

        
            tmpdcm.save_as(tmpfile)
    return
    


def Catphan_Analysis(data, results,actions):
    dcmInfile = os.path.dirname(os.path.abspath(data.series_filelist[0][0]))
    
    try:
        params = action['params']
    except:
        params = {}

        
    version = params['version']
    if params['classifier'] == 'true':
       classifier = True
    else:
       classifier = False
    hut = int(params['hu_tolerance'])

    addboundary = params['add_boundary']
    if addboundary == "True":
        _addBoundary(dcmInfile)

    
    if not version in ["503","504","603","604"]:
        print ('Sorry, Catphan version not supported! has to be 503,504,600 or 604')
        sys.exit()
    else:
        if version == "503":
           tmpcat = pylinac.CatPhan503(dcmInfile)
        elif version == "504":
           tmpcat = pylinac.CatPhan504(dcmInfile)

        elif version == "600":
           tmpcat = pylinac.CatPhan600(dcmInfile)
        elif version == "604":
           tmpcat = pylinac.CatPhan604(dcmInfile)


    tmpcat.analyze(hu_tolerance=hut)
    #print (tmpcat.return_results())
    'HU Linearity ROI'
    tmphu =  tmpcat.ctp404.hu_roi_vals
    for key in tmphu.keys():
        results.addFloat('HU_'+key,float(tmphu[key]))

    
    results.addBool('HU Passed',bool(tmpcat.ctp404.passed_hu))
    results.addFloat('Uniformity index', float(tmpcat.ctp486.uniformity_index))
    results.addFloat('Integral non-uniformity',float(tmpcat.ctp486.integral_non_uniformity))
    results.addBool('Uniformity Passed',bool(tmpcat.ctp486.overall_passed))
    results.addFloat('Low contrast visibility',tmpcat.ctp404.lcv)
    results.addFloat('MTF 50 (lp/mm)', tmpcat.ctp528.mtf(50))
    results.addFloat('Geometric Line Average (mm)',tmpcat.ctp404.avg_line_length)
    results.addBool('Geometry Passed', bool(tmpcat.ctp404.passed_geometry))
    results.addFloat('Slice Thickness (mm)', tmpcat.ctp404.meas_slice_thickness)
    results.addBool('Slice Thickness Passed',bool(tmpcat.ctp404.passed_thickness))


    objects = ['linearity','rmtf','hu','un','sp','prof']

    for obj in objects:
        tmpcat.save_analyzed_subimage('%s.jpg'%obj,subimage=obj)
        results.addObject(obj,'%s.jpg'%obj)


def acqdatetime_series(data, results, action):
    """
    Read acqdatetime from dicomheaders and write to IQC database

    Workflow:
        1. Read only headers
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)

    results.addDateTime('AcquisitionDateTime', dt)



from wad_qc.module import pyWADinput
if __name__ == "__main__":
    data, results, config = pyWADinput()

    print(config)
    for name,action in config['actions'].items():

        if name == 'acqdatetime':
            acqdatetime_series(data, results, action)

        elif name == 'Catphan_Analysis':
            Catphan_Analysis(data, results, action)


    results.write()



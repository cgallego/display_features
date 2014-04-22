# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 12:18:04 2014

This script will load volumes, Load a lesion Segmentation (VOI), Visualize volumes
and then extract Dynamic, Morphology and Texture features from the VOI.

Arguments:
============
sys.argv[1] = input text file with one case per line in the following format:
StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID


@ author (C) Cristina Gallego, University of Toronto, 2014
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string

from inputs_init import *
from feature_maps import *

from features_dynamic import *
from features_morphology import *
from features_texture import *


# Open filename list
file_ids = open(sys.argv[1],"r")
for fileline in file_ids:
    # Get the line: StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID    image_pos_pat    image_ori_pat
    fileline = fileline.split()
    cond = fileline[0] 
    StudyID = fileline[1]  
    DicomExamNumber = fileline[2]
    Lesions_id = fileline[3]
    dateID = fileline[4]
    SeriesID = fileline[5] # corresponds to dynamic sequence;
    
    ###### Loading 
    print "Start by loading volumes..."
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = 'Z:\\Cristina\\MassNonmass\\mass'
    load = Inputs_init()
    [series_path, phases_series, lesionID_path] = load.readVolumes(path_rootFolder, StudyID, DicomExamNumber, SeriesID, Lesions_id)
    print "Path to series location: %s" % series_path 
    print "List of pre and post contrast volume names: %s" % phases_series
    print "Path to lesion segmentation: %s" % lesionID_path    
    
    print "\n Load Segmentation..."
    lesion3D = load.loadSegmentation(lesionID_path)
    print "Data Structure: %s" % lesion3D.GetClassName()
    print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
    print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
    
    print "\n Visualize in sub volumes..."
    loadDisplay = Display()
    lesion3D_mesh = loadDisplay.addSegment(lesion3D)
    loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=1, interact=True)

    ###### Extract Dynamic features
    print "\n Visualize certain features maps..."
        
    VOIbounds = load.lesion3D_mesh.GetBounds()
    print "VOI boundas: %d, %d, %d, %d, %d, %d" % VOIbounds
    fmaps = Maps() 
    img_feature = fmaps.featuremap_beta(load.DICOMImages, series_path, phases_series, load.image_pos_pat, load.image_ori_pat, VOIbounds, features=['beta'])
    
    
    
    
    
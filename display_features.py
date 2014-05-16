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
    caseLabel=cond+fileline[9]+'_'+fileline[10]
    
    ###### Loading 
    print "Start by loading volumes..."
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = 'Z:\\Cristina\\MassNonmass\\mass'
    path_outputFolder ="Z:\Cristina\MassNonmass\codeProject\codeBase\display_features"         
    
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
    
    ###### Extract Dynamic features
    print "\n Visualize certain features maps... Extract VOI data first... "
    fmaps = Maps() 
    [img_features, time_points] = fmaps.getVOIdata(load.DICOMImages, series_path, phases_series, load.image_pos_pat, load.image_ori_pat, load.spacing, lesion3D, path_outputFolder, caseLabel)
    
    # extract feature maps and Visualize each one
    fmaps.featureMap(load.DICOMImages, img_features, time_points, ['beta','Tpeak','Kpeak'], caseLabel,  path_outputFolder)
    fmaps.renderer1.Render()
    fmaps.renWin1.Finalize()

    
    
    
    
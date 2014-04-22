# -*- coding: utf-8 -*-
"""
Create visualization with standard vtk actors, renders, windowsn, interactors

USAGE:
=============
from feature_maps import *
featuremaps = Maps()  
featuremaps.dicomTransform(image, image_pos_pat, image_ori_pat)
featuremaps.addSegment(lesion3D)
featuremaps.subImage(Images2Sub, timep)                  
featuremaps.visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Methods:
=============
dicomTransform(image, image_pos_pat, image_ori_pat)
addSegment(lesion3D)
subImage(Images2Sub, timep)                  
visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Instance Attributes:
===============
'origin': (-167.0, -69.0, -145.0)
'spacing': (0.44920000433921814, 0.44920000433921814, 3.0)
'dims': (512, 512, 96), 

VTK Instance objects:
=============
'xImagePlaneWidget': (vtkImagePlaneWidget)
'yImagePlaneWidget': (vtkImagePlaneWidget)
'zImagePlaneWidget': (vtkImagePlaneWidget)
'picker': (vtkCellPicker)
'iren1': (vtkWin32RenderWindowInteractor)
'camera': (vtkOpenGLCamera)
'mapper_mesh': (vtkPainterPolyDataMapper)
'actor_mesh': (vtkOpenGLActor)
'renWin1': (vtkWin32OpenGLRenderWindow)
'renderer1': (vtkOpenGLRenderer)

Created on Mon Apr 21 12:59:02 2014

@ author (C) Cristina Gallego, University of Toronto, 2014
----------------------------------------------------------------------
"""


import os, os.path
import sys
import string
import datetime
from numpy import *
import dicom
from vtk.util.numpy_support import vtk_to_numpy
import vtk

from lmfit import minimize, Parameters, Parameter, report_errors, Minimizer
from lmfit.printfuncs import *
from pandas import DataFrame
import pandas as pd
import pylab

from inputs_init import *
from display import *

#!/usr/bin/env python
class Maps(object):
    """
    USAGE:
    =============
    loadDisplay = Display()
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.005) 
        
        # Create 3 orthogonal view using the ImagePlaneWidget
        self.xImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.yImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.zImagePlaneWidget = vtk.vtkImagePlaneWidget()
        
        #  The 3 image plane widgets
        self.xImagePlaneWidget.DisplayTextOn();
        self.xImagePlaneWidget.SetPicker(self.picker);
        self.xImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.xImagePlaneWidget.SetKeyPressActivationValue('x');
        self.xImagePlaneWidget.GetPlaneProperty().SetColor(1, 0, 0);
        self.xImagePlaneWidget.SetResliceInterpolateToNearestNeighbour();
        
        self.yImagePlaneWidget.DisplayTextOn();
        self.yImagePlaneWidget.SetPicker(self.picker);
        self.yImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.yImagePlaneWidget.SetKeyPressActivationValue('y');
        self.yImagePlaneWidget.GetPlaneProperty().SetColor(0, 1, 0);
        self.yImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        
        self.zImagePlaneWidget.DisplayTextOn();
        self.zImagePlaneWidget.SetPicker(self.picker);
        self.zImagePlaneWidget.SetKeyPressActivationValue('z');
        self.zImagePlaneWidget.GetPlaneProperty().SetColor(0, 0, 1);
        self.zImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        self.zImagePlaneWidget.SetRightButtonAutoModifier(1);
        
        # Create a renderer, render window, and render window interactor to
        # display the results.
        self.renderer1 = vtk.vtkRenderer()
        self.renWin1 = vtk.vtkRenderWindow()
        self.iren1 = vtk.vtkRenderWindowInteractor()
        
        self.renWin1.SetSize(1000, 800);
        self.renWin1.AddRenderer(self.renderer1);
        self.iren1.SetRenderWindow(self.renWin1);
        
        self.xImagePlaneWidget.SetInteractor( self.iren1 )
        self.yImagePlaneWidget.SetInteractor( self.iren1 )
        self.zImagePlaneWidget.SetInteractor( self.iren1 )
        
        # Set Up Camera view
        self.camera = self.renderer1.GetActiveCamera()
        self.renderer1.SetBackground(0.0, 0.0, 0.0)
        self.iren1.SetPicker(self.picker)
                
        self.origin=[]
        
    def __call__(self):       
        """ Turn Class into a callable object """
        Maps()
        
    def visualize_map(self, VOIclip):
        
        # get info from image before visualization
        VOIclip.UpdateInformation()
        self.dims = VOIclip.GetDimensions()
        print "Image Dimensions"
        print self.dims
        (xMin, xMax, yMin, yMax, zMin, zMax) = VOIclip.GetWholeExtent()
        print "Image Extension"
        print xMin, xMax, yMin, yMax, zMin, zMax
        self.spacing = VOIclip.GetSpacing()
        print "Image Spacing"
        print self.spacing
            
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( VOIclip )
        self.xImagePlaneWidget.SetPlaneOrientationToXAxes()
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( VOIclip )
        self.yImagePlaneWidget.SetPlaneOrientationToYAxes()
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( VOIclip )
        self.zImagePlaneWidget.SetPlaneOrientationToZAxes()
        self.zImagePlaneWidget.SetSliceIndex(0)
            
        self.xImagePlaneWidget.On()
        self.yImagePlaneWidget.On()
        self.zImagePlaneWidget.On()
        
        # set up cube actor with Orientation(A-P, S-I, L-R) using transform_cube
        # Set up to ALS (+X=A, +Y=S, +Z=L) source:
        cube = vtk.vtkAnnotatedCubeActor()
        cube.SetXPlusFaceText( "L" );
        cube.SetXMinusFaceText( "R" );
        cube.SetYPlusFaceText( "A" );
        cube.SetYMinusFaceText( "P" );
        cube.SetZPlusFaceText( "S" );
        cube.SetZMinusFaceText( "I" );
        cube.SetFaceTextScale( 0.5 );
            
        # Set UP the axes
        axes2 = vtk.vtkAxesActor()
        axes2.SetShaftTypeToCylinder();
        #axes2.SetUserTransform( transform_cube );         
        axes2.SetTotalLength( 1.5, 1.5, 1.5 );
        axes2.SetCylinderRadius( 0.500 * axes2.GetCylinderRadius() );
        axes2.SetConeRadius( 1.025 * axes2.GetConeRadius() );
        axes2.SetSphereRadius( 1.500 * axes2.GetSphereRadius() );
    
        tprop2 = axes2.GetXAxisCaptionActor2D()
        tprop2.GetCaptionTextProperty();
    
        assembly = vtk.vtkPropAssembly();
        assembly.AddPart( axes2 );
        assembly.AddPart( cube );
        
        widget = vtk.vtkOrientationMarkerWidget();
        widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 );
        widget.SetOrientationMarker( assembly );
        widget.SetInteractor( self.iren1 );
        widget.SetViewport( 0.0, 0.0, 0.4, 0.4 );
        widget.SetEnabled( 1 );
        widget.InteractiveOff();
                    
        # Create a text property for both cube axes
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(1, 1, 1)
        tprop.ShadowOff()
        
        # Create a vtkCubeAxesActor2D.  Use the outer edges of the bounding box to
        # draw the axes.  Add the actor to the renderer.
        axes = vtk.vtkCubeAxesActor2D()
        axes.SetInput(VOIclip)
        axes.SetCamera(self.renderer1.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(1.2)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)      
        self.renderer1.AddViewProp(axes)
        
        ############
        # bounds and initialize camera
        bounds = VOIclip.GetBounds()
        self.renderer1.ResetCamera(bounds)    
        self.renderer1.ResetCameraClippingRange()
        self.camera.SetViewUp(0.0,-1.0,0.0)
        self.camera.Azimuth(315)
        
        # Initizalize
        self.renWin1.Render()
        self.renderer1.Render()
        self.iren1.Start()  
                            
        return
        
        
    def featuremap_beta(self, DICOMImages, series_path, phases_series, image_pos_pat, image_ori_pat, VOIbounds, features):
        """ featuremap_beta: Creates some feature maps from image
        
        INPUTS:
        =======        
        image: (vtkImageData)    Input image to Transform
        image_pos_pat: (list(dicomInfo[0x0020,0x0032].value)))  Image position patient Dicom Tag
        image_ori_pat: (list(dicomInfo[0x0020,0x0037].value))   Image oreintation patient Dicom Tag
        
        OUTPUTS:
        =======
        transformed_image (vtkImageData)    Transformed imaged mapped to dicom coords frame
        transform (vtkTransform)            Transform used
        
        """ 
        if 'beta' in features:
            pixVals = []
            deltaS = {}
            timepoints = [] 
            
            # necessary to read point coords
            VOIPnt = [0,0,0]
            ijk = [0,0,0]
            pco = [0,0,0]
            
            for i in range(len(DICOMImages)):
                abspath_PhaseID = series_path+os.sep+str(phases_series[i]) 
                print phases_series[i]
                
                # Get total number of files
                load = Inputs_init()
                [len_listSeries_files, FileNms_slices_sorted_stack] = load.ReadDicomfiles(abspath_PhaseID)
                mostleft_slice = FileNms_slices_sorted_stack.slices[0]
                
                # Get dicom header, retrieve
                dicomInfo_series = dicom.read_file(abspath_PhaseID+os.sep+str(mostleft_slice)) 
                
                # (0008,0031) AT S Series Time            # hh.mm.ss.frac
                seriesTime = str(dicomInfo_series[0x0008,0x0031].value) 
                # (0008,0033) AT S Image Time             # hh.mm.ss.frac
                imageTime = str(dicomInfo_series[0x0008,0x0033].value)
                
                # (0008,0032) AT S Acquisition Time       # hh.mm.ss.frac
                ti = str(dicomInfo_series[0x0008,0x0032].value) 
                
                acquisitionTimepoint = datetime.time(hour=int(ti[0:2]), minute=int(ti[2:4]), second=int(ti[4:6]))
                timepoints.append( datetime.datetime.combine(datetime.date.today(), acquisitionTimepoint) )
                
                # find mapping to Dicom space  
                [transformed_image, transform_cube] = Display().dicomTransform(DICOMImages[i], image_pos_pat, image_ori_pat)
                
                #################### HERE GET IT AND MASK IT OUT       
#                clipper=vtk.vtkImageClip()
#                clipper.ClipDataOn()
#                clipper.SetOutputWholeExtent(VOIbounds) #where 0 gives all the cells outside of the box, and 1 gives all the cells inside the box.
#                clipper.SetInput( transformed_image )
#                clipper.Update()
                
                extract = vtk.vtkExtractVOI()
                extract.SetVOI(VOIbounds)
                extract.SetSampleRate(1, 1, 1)
                extract.SetInput(transformed_image)
                extract.Update()
                
                VOIimage = vtk.vtkImageData()
                VOIimage.SetWholeExtent(VOIbounds)
                VOIimage.SetOrigin(VOIbounds[0], VOIbounds[2], VOIbounds[4])
                VOIimage.DeepCopy(extract.GetOutput())
                VOIimage.Update()
                
                self.visualize_map(VOIimage)
                
            
            print timepoints
            
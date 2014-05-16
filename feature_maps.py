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
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtk.util import vtkImageImportFromArray
import vtk

from lmfit import minimize, Parameters, Parameter, report_errors, Minimizer
from lmfit.printfuncs import *
from pandas import DataFrame
import pandas as pd
import pylab
from scipy import stats

from inputs_init import *
from display import *

#!/usr/bin/env python
class Maps(object):
    """
    USAGE:
    =============
    maps = Maps()
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.timepoints = [] 
        self.time_points = []
        self.deltaS = {}
        self.model = []
        self.ext = []
        self.ix = []
        self.iy = []
        self.iz = []
        self.maxCr = 0 
        
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
    
    def addSegment(self, lesion3D):        
        # Set the planes based on seg bounds
        bounds = lesion3D.GetBounds()
        print "\n Mesh DICOM bounds: "
        print "xmin, xmax= [%d, %d]" % (bounds[0], bounds[1])
        print "yin, ymax= [%d, %d]" %  (bounds[2], bounds[3]) 
        print "zmin, zmax= [%d, %d]" % (bounds[4], bounds[5])
        
        # Add ICPinit_mesh.vtk to the render
        self.mapper_mesh = vtk.vtkPolyDataMapper()
        self.mapper_mesh.SetInput( lesion3D )
        self.mapper_mesh.ScalarVisibilityOff()
        
        self.actor_mesh = vtk.vtkActor()
        self.actor_mesh.SetMapper(self.mapper_mesh)
        self.actor_mesh.GetProperty().SetColor(0, 1, 0)    #R,G,B
        self.actor_mesh.GetProperty().SetOpacity(0.6)
        self.actor_mesh.GetProperty().SetPointSize(5.0)
        self.actor_mesh.GetProperty().SetRepresentationToWireframe()
        
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetSliceIndex( 0 )
        
        self.renderer1.AddActor(self.actor_mesh)
        self.renWin1.Render()
                
        return
        
    def convertArray2vtkImage(self, nparray, t_ImagedataVTK, npImagesandMask): 
        """ Takes a numpy.ndarray and converts it to a vtkimageData. require npImagesandMask to pass on image info """
        # Create vtk object
        size_array = npImagesandMask['dims'][0]*npImagesandMask['dims'][1]*npImagesandMask['dims'][2]
        flatim = nparray.transpose(2,1,0)
        flatim = flatim.flatten()
        
        # create vtk image
        vtk_image = vtk.vtkImageData()
        vtk_image.DeepCopy(t_ImagedataVTK)
        vtk_image.SetNumberOfScalarComponents(1)
        vtk_image.SetScalarTypeToDouble()
        vtk_image.AllocateScalars()
        
        # Get scalars from numpy
        image_array = vtk.vtkDoubleArray() 
        image_array.SetNumberOfValues(size_array)
        image_array.SetNumberOfComponents(1) 
        
        # not too efficient convertion of np.array to vtk. Far from ideal
        for k in range(size_array):
            image_array.InsertTuple1(k,flatim[k])
            
        vtk_image.GetPointData().SetScalars(image_array) 
        vtk_image.Update()
          
        return vtk_image   
        
    def visualize_map(self, VOIclip):
        # get info from image before visualization
        VOIclip.UpdateInformation()
        self.dims = VOIclip.GetDimensions()
        (xMin, xMax, yMin, yMax, zMin, zMax) = VOIclip.GetWholeExtent()
        self.spacing = VOIclip.GetSpacing()

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
        self.camera.Azimuth(360)
        
        # Initizalize
        #self.renWin1.Render()
        self.renderer1.Render()
        self.iren1.Start()  
                            
        return
        
    def createMaskfromMesh(self, VOI_mesh, im):
        """ Takes an image and a VOI_mesh and returns a boolean image with only 1s inside the VOI_mesh """                
        # Create an Image 
        white_image = vtk.vtkImageData()
        white_image.DeepCopy(im) 
        
        # extract VOI bounds in dicom space
        self.VOIbounds = VOI_mesh.GetBounds()
        print "self.VOIbounds:"
        print self.VOIbounds
        
        self.VOIbounds_expand = []
        self.VOIbounds_expand.append( self.VOIbounds[0] )
        self.VOIbounds_expand.append( self.VOIbounds[1]  )
        self.VOIbounds_expand.append( self.VOIbounds[2]  )
        self.VOIbounds_expand.append( self.VOIbounds[3] )
        self.VOIbounds_expand.append(  self.VOIbounds[4] )
        self.VOIbounds_expand.append( self.VOIbounds[5] )
            
        self.VOIbounds = self.VOIbounds_expand
            
        roiStencil = vtk.vtkCubeSource()
        roiStencil.SetBounds(self.VOIbounds)
        roiStencil.Update()
        print "\nGetXLength roiStencil: %d " % roiStencil.GetXLength()
        print "GetYLength roiStencil: %d " % roiStencil.GetYLength()
        print "GetZLength roiStencil: %d " % roiStencil.GetZLength()
  
        # polygonal data --> image stencil:
        pol2stenc = vtk.vtkPolyDataToImageStencil()
        pol2stenc.SetInput(roiStencil.GetOutput())
        pol2stenc.SetOutputOrigin(im.GetOrigin())
        pol2stenc.SetOutputSpacing(im.GetSpacing())
        pol2stenc.SetOutputWholeExtent(im.GetWholeExtent())
        pol2stenc.Update()
         
        # cut the corresponding white image and set the background:
        imgstenc = vtk.vtkImageStencil()
        imgstenc.SetInput(im)
        imgstenc.SetStencil(pol2stenc.GetOutput())
        imgstenc.ReverseStencilOff()
        imgstenc.SetBackgroundValue(0.0)
        imgstenc.Update()
        
        self.VOIdims = imgstenc.GetOutput().GetDimensions()
                  
        return imgstenc.GetOutput()
        
                
    def convertfeatureMap2vtkImage(self, nparray, imageStencil, factor): 
        """ Takes a numpy.ndarray from a featureMap and writes at the ijk positions corresponding to imageStencil """
        # Create vtk object
        nparray_dims = nparray.shape
        #size_array = nparray_dims[0]*nparray_dims[1]*nparray_dims[2]
                    
        # get non zero elements
        VOI_scalars = imageStencil.GetPointData().GetScalars()
        numpy_VOI_imageStencil = vtk_to_numpy(VOI_scalars)     
        
        numpy_VOI_imageStencil = numpy_VOI_imageStencil.reshape(self.VOIdims[2], self.VOIdims[1], self.VOIdims[0]) 
        numpy_VOI_imageStencil = numpy_VOI_imageStencil.transpose(2,1,0)
        
        # iterate point-by-point to extract feature map       
        for i in range(nparray_dims[0]):
            for j in range(nparray_dims[1]):
                for k in range(nparray_dims[2]):
                    newpixValx = nparray[i,j,k]*factor
                    if newpixValx > 255: 
                        newpixValx=0                        
                    numpy_VOI_imageStencil[ self.ix+i, self.iy+j, self.iz+k] = newpixValx
                    imageStencil.SetScalarComponentFromFloat( self.ix+i, self.iy+j, self.iz+k, 0, newpixValx)
                    imageStencil.UpdateData()
   
        return imageStencil
        
    
    def getVOIdata(self, DICOMImages, series_path, phases_series, image_pos_pat, image_ori_pat, VOIspacing, VOI_mesh, path_outputFolder, caseLabeloutput):
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

            # (0008,0032) AT S Acquisition Time       # hh.mm.ss.frac
            ti = str(dicomInfo_series[0x0008,0x0032].value) 
            
            acquisitionTimepoint = datetime.time(hour=int(ti[0:2]), minute=int(ti[2:4]), second=int(ti[4:6]))
            self.timepoints.append( datetime.datetime.combine(datetime.date.today(), acquisitionTimepoint) )
            
            # find mapping to Dicom space  
            displayFuncs = Display()
            [self.transformed_image, transform_cube] = displayFuncs.dicomTransform(DICOMImages[i], image_pos_pat, image_ori_pat)
             
            #################### HERE GET IT AND MASK IT OUT                       
            ### Get inside of VOI            
            self.imageStencil = self.createMaskfromMesh(VOI_mesh, self.transformed_image)
            if i==0:
                print path_outputFolder
                print caseLabeloutput
                # ## save mask as metafile image
                os.chdir(path_outputFolder)
                vtkmask_w = vtk.vtkMetaImageWriter()
                vtkmask_w.SetInput( self.imageStencil )
                vtkmask_w.SetFileName( 'pre_con.mhd' )
                vtkmask_w.Write()
                vtkmask_w.Update()
                
            # get non zero elements
            VOI_scalars = self.imageStencil.GetPointData().GetScalars()
            numpy_VOI_scalars = vtk_to_numpy(VOI_scalars)     
            
            numpy_VOI_scalars = numpy_VOI_scalars.reshape(self.VOIdims[2], self.VOIdims[1], self.VOIdims[0]) 
            numpy_VOI_scalars = numpy_VOI_scalars.transpose(2,1,0)
            
            print "\n VOIbounds"
            print self.VOIbounds
            
            # compute origin
            IMorigin = displayFuncs.origin

            if i==0:               
                print "\n Compute VOI extent:"
                self.ext.append( round((self.VOIbounds[1]-self.VOIbounds[0])/VOIspacing[0]) )
                self.ext.append( round((self.VOIbounds[3]-self.VOIbounds[2])/VOIspacing[1]) )
                self.ext.append( round((self.VOIbounds[5]-self.VOIbounds[4])/VOIspacing[2]) )
                print self.ext                
                            
            # get non zero elements
            image_scalars = self.transformed_image.GetPointData().GetScalars()
            numpy_VOI_imagedata = vtk_to_numpy(image_scalars)     
            
            numpy_VOI_imagedata = numpy_VOI_imagedata.reshape(self.VOIdims[2], self.VOIdims[1], self.VOIdims[0]) 
            numpy_VOI_imagedata = numpy_VOI_imagedata.transpose(2,1,0)
            
            ## compute ijk extent
            self.ix = round((self.VOIbounds[0]-IMorigin[0])/VOIspacing[0])
            self.iy = round((self.VOIbounds[2]-IMorigin[1])/VOIspacing[1])
            self.iz = round((self.VOIbounds[4]-IMorigin[2])/VOIspacing[2])
            
            numpy_VOI_imagedataext = numpy_VOI_imagedata[ self.ix:self.ix+self.ext[0],  self.iy:self.iy+self.ext[1],  self.iz:self.iz+self.ext[2] ]
            
            print "Shape of numpy_VOI_imagedataext: "
            print size(numpy_VOI_imagedataext)
            
            #################### HERE GET IT AND MASK IT OUT
            # Now collect pixVals
            print "Saving %s" % 'VOIimage'+str(i)
            self.deltaS['VOI'+str(i)] = numpy_VOI_imagedataext
            numpy_VOI_imagedataext = []
            
            # Visualize
            if i==0: 
                self.addSegment(VOI_mesh)
                self.renderer1.RemoveActor(self.actor_mesh) 
                self.visualize_map(self.imageStencil)
                
        
               
        # Collecting timepoints in proper format
        t_delta = []
        t_delta.append(0)
        total_time = 0
        for i in range(len(DICOMImages)-1):
            current_time = self.timepoints[i+1]
            previous_time = self.timepoints[i]
            difference_time = current_time - previous_time
            timestop = divmod(difference_time.total_seconds(), 60)
            t_delta.append( t_delta[i] + timestop[0]+timestop[1]*(1./60))
            total_time = total_time+timestop[0]+timestop[1]*(1./60)
            
        # finally print t_delta
        print "\n time_points:"
        print t_delta
        print "total_time"
        print total_time
        self.time_points = array(t_delta)
                        
        return self.deltaS, self.time_points  
        
        
    def init_features(self, img_features, featuresKeys):
        """ Initializes feature map objects based on request from vector of keywords featuresKeys"""
        self.R_square_map = zeros(shape(img_features['VOI0']))
        
        if 'amp' in featuresKeys:
            self.amp_map = zeros(shape(img_features['VOI0']))
        
        if 'beta' in featuresKeys:
            self.beta_map = zeros(shape(img_features['VOI0']))
            
        if 'alpha' in featuresKeys:
            self.alpha_map = zeros(shape(img_features['VOI0']))
            
        if 'iAUC1' in featuresKeys:
            self.iAUC1_map = zeros(shape(img_features['VOI0']))

        if 'Slope_ini' in featuresKeys:
            self.Slope_ini_map = zeros(shape(img_features['VOI0']))

        if 'Tpeak' in featuresKeys:
            self.Tpeak_map = zeros(shape(img_features['VOI0']))

        if 'Kpeak' in featuresKeys:
            self.Kpeak_map = zeros(shape(img_features['VOI0']))
            
        if 'SER' in featuresKeys:
            self.SER_map = zeros(shape(img_features['VOI0']))
            
        if 'maxCr' in featuresKeys:
            self.maxCr_map = zeros(shape(img_features['VOI0']))
        
        if 'peakCr' in featuresKeys:
            self.peakCr_map = zeros(shape(img_features['VOI0']))
        
        if 'UptakeRate' in featuresKeys:
            self.UptakeRate_map = zeros(shape(img_features['VOI0']))
        
        if 'washoutRate' in featuresKeys:
            self.washoutRate_map = zeros(shape(img_features['VOI0']))
        
        if 'var_F_r_i' in featuresKeys:
            self.allvar_F_r_i_map = zeros(shape(img_features['VOI0']))

        return

    def fcn2min(self, params, t, data):
        """ model EMM for Bilateral DCE-MRI, subtract data"""
        # define objective function: returns the array to be minimized
        # unpack parameters:
        #  extract .value attribute for each parameter
        amp = params['amp'].value    # Upper limit of self.deltaS
        alpha = params['alpha'].value    # rate of signal increase min-1
        beta = params['beta'].value        # rate of signal decrease min-1
                
        self.model = amp * (1- exp(-alpha*t)) * exp(-beta*t)
    
        return self.model - data
        
        
    def featureMap(self, DICOMImages, img_features, time_points, featuresKeys, caseLabeloutput, path_outputFolder):
        """Extracts feature maps per pixel based on request from vector of keywords featuresKeys """
        ## Retrive image data
        VOIshape = img_features['VOI0'].shape
        print VOIshape
        self.init_features(img_features, featuresKeys)
        data_deltaS=[]  
        self.allvar_F_r_i=[]
        
        # append So and to
        data_deltaS.append( 0 )  
        
        # Based on the course of signal intensity within the lesion
        So = array(img_features['VOI0']).astype(float)
        Crk = {'Cr0': mean(So)}  
        C = {}
        Carray = []
        
        # iterate point-by-point to extract feature map       
        for i in range(VOIshape[0]):
            for j in range(VOIshape[1]):
                for k in range(VOIshape[2]):
                    for timep in range(1, len(DICOMImages)):
                        pix_deltaS = (img_features['VOI'+str(timep)][i,j,k].astype(float) - img_features['VOI0'][i,j,k].astype(float))/img_features['VOI0'][i,j,k].astype(float)
                        if pix_deltaS<0: pix_deltaS=0 
                        data_deltaS.append( pix_deltaS )
                        
                        F_r_i =  array(img_features['VOI'+str(timep)]).astype(float)
                        n_F_r_i, min_max_F_r_i, mean_F_r_i, var_F_r_i, skew_F_r_i, kurt_F_r_i = stats.describe(F_r_i)
                        self.allvar_F_r_i.append(var_F_r_i)
                    
                    data = array(data_deltaS)
                    #print data                        
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('amp',   value= 10,  min=0)
                    params.add('alpha', value= 1, min=0) 
                    params.add('beta', value= 0.05, min=0.0001, max=0.9)
                    
                    # do fit, here with leastsq self.model
                    myfit = Minimizer(self.fcn2min,  params, fcn_args=(time_points,), fcn_kws={'data':data})
                    myfit.prepare_fit()
                    myfit.leastsq()
 
                    ####################################                              
                    # Calculate R-square: R_square = sum( y_fitted - y_mean)/ sum(y_data - y_mean)
                    R_square = sum( (self.model - mean(data))**2 )/ sum( (data - mean(data))**2 )
                    #print "R^2:"
                    #print R_square
                    self.R_square_map[i,j,k] = R_square
                    
                    if 'amp' in featuresKeys:
                        amp = params['amp'].value
                        print "amp:"
                        print amp
                        self.amp_map[i,j,k] = amp
                    
                    if 'beta' in featuresKeys:
                        beta = params['beta'].value
                        self.beta_map[i,j,k] = beta
                        
                    if 'alpha' in featuresKeys:
                        alpha = params['alpha'].value
                        print "alpha:"
                        print alpha
                        self.alpha_map[i,j,k] = alpha
                        
                    if 'iAUC1' in featuresKeys:
                        iAUC1 = params['amp'].value *( ((1-exp(-params['beta'].value*t[1]))/params['beta'].value) + (exp((-params['alpha'].value+params['beta'].value)*t[1])-1)/(params['alpha'].value+params['beta'].value) )
                        print "iAUC1"
                        print iAUC1
                        self.iAUC1_map[i,j,k] = iAUC1
                        
                    if 'Slope_ini' in featuresKeys:
                        Slope_ini = params['amp'].value*params['alpha'].value
                        print "Slope_ini"
                        print Slope_ini
                        self.Slope_ini_map[i,j,k] = Slope_ini
                    
                    if 'Tpeak' in featuresKeys:
                        Tpeak = (1/params['alpha'].value)*log(1+(params['alpha'].value/params['beta'].value))
                        self.Tpeak_map[i,j,k] = Tpeak
                    
                    if 'Kpeak' in featuresKeys:
                        Kpeak = -params['amp'].value * params['alpha'].value * params['beta'].value
                        self.Kpeak_map[i,j,k] = Kpeak
                    
                    if 'SER' in featuresKeys:
                        SER = exp( (t[4]-t[1])*params['beta'].value) * ( (1-exp(-params['alpha'].value*t[1]))/(1-exp(-params['alpha'].value*t[4])) )
                        print "SER"
                        print SER
                        self.SER_map[i,j,k] = SER
                        
                    if 'maxCr' in featuresKeys:
                        print "Maximum Upate (Fii_1) = %d " %  self.maxCr
                        self.maxC_map[i,j,k] = self.maxCr
                        
                    if 'peakCr' in featuresKeys:
                        print "Peak Cr (Fii_2) = %d " %  self.peakCr
                        self.peakCr_map[i,j,k] = self.peakCr    
                        
                    if 'UptakeRate' in featuresKeys:
                        self.UptakeRate = float(self.maxCr/self.peakCr)    
                        print "Uptake rate (Fii_3) "
                        print self.UptakeRate
                        self.UptakeRate_map[i,j,k] = self.UptakeRate
                        
                    if 'washoutRate' in featuresKeys:          
                        if( self.peakCr == 4):
                            self.washoutRate = 0
                        else:
                            self.washoutRate = float( (self.maxCr - array(Crk['Cr'+str(4)]).astype(float))/(4-self.peakCr) )
                        print "WashOut rate (Fii_4) "
                        print self.washoutRate
                        self.washoutRate_map[i,j,k] = self.washoutRate
                        
                    if 'var_F_r_i' in featuresKeys:  
                        print("Variance F_r_i: {0:8.6f}".format( mean(self.allvar_F_r_i) ))
                        self.allvar_F_r_i_map[i,j,k] = mean(self.allvar_F_r_i)
                    
                    data_deltaS=[]
                    data_deltaS.append( 0 )
                    
        # convert feature maps to image
        if 'beta' in featuresKeys:       
            beta_map_stencil = self.convertfeatureMap2vtkImage(self.beta_map, self.imageStencil, 1000) 
            print path_outputFolder
            print caseLabeloutput
            # ## save mask as metafile image
            os.chdir(path_outputFolder)
            vtkmask_w = vtk.vtkMetaImageWriter()
            vtkmask_w.SetInput(beta_map_stencil )
            vtkmask_w.SetFileName( 'beta_'+os.sep+caseLabeloutput+'.mhd' )
            vtkmask_w.Write()
            vtkmask_w.Update()
            
            self.xImagePlaneWidget.SetWindowLevel(640,75)
            self.yImagePlaneWidget.SetWindowLevel(640,75)
            self.zImagePlaneWidget.SetWindowLevel(640,75)
            self.renderer1.Render()            
            self.visualize_map(beta_map_stencil)    
            
            
        if 'Tpeak' in featuresKeys:
            Tpeak_map_stencil = self.convertfeatureMap2vtkImage(self.Tpeak_map, self.imageStencil, 1) 
            print path_outputFolder
            print caseLabeloutput
            # ## save mask as metafile image
            os.chdir(path_outputFolder)
            vtkmask_w = vtk.vtkMetaImageWriter()
            vtkmask_w.SetInput( Tpeak_map_stencil )
            vtkmask_w.SetFileName( 'Tpeak_'+os.sep+caseLabeloutput+'.mhd' )
            vtkmask_w.Write()
            vtkmask_w.Update()
            
            self.xImagePlaneWidget.SetWindowLevel(240,35)
            self.yImagePlaneWidget.SetWindowLevel(240,35)
            self.zImagePlaneWidget.SetWindowLevel(240,35)
            self.renderer1.Render()
            self.visualize_map(Tpeak_map_stencil)    
        
        if 'Kpeak' in featuresKeys:
            Kpeak_map_stencil = self.convertfeatureMap2vtkImage(self.Kpeak_map, self.imageStencil, 1) 
            print path_outputFolder
            print caseLabeloutput
            # ## save mask as metafile image
            os.chdir(path_outputFolder)
            vtkmask_w = vtk.vtkMetaImageWriter()
            vtkmask_w.SetInput( Kpeak_map_stencil )
            vtkmask_w.SetFileName( 'Kpeak_'+os.sep+caseLabeloutput+'.mhd' )
            vtkmask_w.Write()
            vtkmask_w.Update()
            
            self.xImagePlaneWidget.SetWindowLevel(118,15)
            self.yImagePlaneWidget.SetWindowLevel(118,15)
            self.zImagePlaneWidget.SetWindowLevel(118,15)
            self.renderer1.Render()           
            self.visualize_map(Kpeak_map_stencil) 
            
               
        
        return
                     
        
        
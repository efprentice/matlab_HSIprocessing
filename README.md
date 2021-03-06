# matlab_HSIprocessing

## These files are needed to run the scripts:

### HSI_wlens.mat  
this is the result of a wavelength calibration (camera specific),
should be done for each camera before each test   

### find_band_num.m  
converts wavelength to approx. band number in image  

## These scripts are for testing HSI images (input is a folder of .png images)

### on_the_fly_testing.m  
used for collecting images in a field scan (terrestrial), output is:  
 - a selection of raw spectral plots per frame,
 used for checking for over-/under-saturation  
  - a spatial image of the scan scene,
  used for verifying spatial targets  
  
 <p align="center">
 <img width="600" height="300" src="sofly.png">
 </p>
  
### select_class_spec.m  
used for getting spectrum of a specific class e.g. water pixels (ui), output is:  
 - spectral plots per class (can use freehand, polygon, of ellipse tool for selection)  
 
  <p align="center">
 <img width="600" height="300" src="classes.png">
 </p>
 
### vert_horz_lines.m 
used for getting spectra of horizontal and vertical lines across a spatial image (ui), output is:  
 - plots of three consecutive vertical lines (or horizontal frames) showing intensity as a function of spatial extent  
 
 <p align="center">
 <img width="600" height="300" src="lines.png">
 </p>
  
### intensity_3d.m  
used for visualizing intensity of a spatial region in the scan image (ui), output is:  
  - 3D plot of intensity values as a function of spected spatial limits  
  
 <p align="center">
 <img width="500" height="300" src="int3d.png">
 </p>

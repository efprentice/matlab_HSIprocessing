# matlab_HSIprocessing

## These files are needed to run the scripts:

### HSI_wlens.mat  
this is the result of a wavelength calibration (camera specific),
should be done for each camera before each test   

### find_band_num.m  
converts wavelength to approx. band number in image  

## These scripts are for testing HSI .png raw outputs (input is a folder of .png images)

### on_the_fly_testing.m  
used for collecting images in a field scan (terrestrial), output is:  
 - a selection of raw spectral plots per frame,
 used for checking for over-/under-saturation  
  - a spatial image of the scan scene,
  used for verifying spatial targets  
  
### select_class_spec.m  
used for getting spectrum of a specific class e.g. water pixels (ui), output is:  
 - spectral plots per class (using muliple vertical pixels and frames)  
 
### vert_horz_lines.m 
used for getting spectra of horizontal and vertical lines across a spatial image (ui), output is  
 - plots of three consecutive lines (frames) showing intensity as a function of spatial extent  
 


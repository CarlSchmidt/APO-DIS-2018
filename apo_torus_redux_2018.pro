;******************************************************************************************
; NAME: 
;      APO_TORUS_REDUCE
;AUTHORS: 
;       Carl Schmidt (LATMOS)
;       
;PURPOSE: Reduce the APO data of the IO Torus from the Dual Imaging Spectrograph 
;

function Scattered_Light, x, p 
  return, p[0]*x^p[1] + p[2]    ; Power law equation 
end 

function Jovian_Centroid_Blue, x, p 
  Common Tilt_Coeffs, Tilt_Coeffs_r, Tilt_Coeffs_b
  ;return, p[0]+p[1]*x+p[2]*(x^2.)  ; Polynomial fiting for jupiter's center
  ;return, p[0] - 0.0024*x - 0. * (x^2.)
  return, p[0] + Tilt_Coeffs_b[1]*x; + Tilt_Coeffs_b[2] * (x^2.)
end 
function Jovian_Centroid_red, x, p 
  Common Tilt_Coeffs, Tilt_Coeffs_r, Tilt_Coeffs_b
  ;return, p[0] + p[1]*x + p[2]*(x^2.)  ; Polynomial fiting for jupiter's center
  ;return, p[0] + P[1]*x 
  return, p[0] + Tilt_Coeffs_r[1]*x; + Tilt_Coeffs_r[2] * (x^2.)
end 
function Gaussian_Ribbon, p, x=x, y=y, err=err, fit=fit
  z = (x - p[1])/p[2]
  fit = p[0]*exp(-z^2/2d) + p[3] + p[4]*x
  return, (y - fit)/err
  ;RETURN, P[0] + GAUSS1(X, P[1:3])
end

function Fit_blue_triplet, p, x=x, y=y, err=err, fit=fit
  Common Dispersions, dispers_b, dispers_r, lines_b, lines_r
  ;P[0] = Gaussian amplitudes of the 3722 SIII line
  ;P[1] = Gaussian amplitudes of the 3726 OII line
  ;P[2] = Gaussian amplitudes of the 3729 OII line
  ;P[3] = Index of the middle 3726 peak. 
  ;P[4] = Gaussian Sigma 
  z_3722 = (x - (p[3] + (lines_b[0] - lines_b[1]) / dispers_b)) / p[4]
  z_3726 = (x - (p[3] + (lines_b[1] - lines_b[1]) / dispers_b)) / p[4]
  z_3729 = (x - (p[3] + (lines_b[2] - lines_b[1]) / dispers_b)) / p[4]
  fit = p[0]*exp(-z_3722^2/2.) + p[1]*exp(-z_3726^2/2.) + p[2]*exp(-z_3729^2/2.); + p[5]
  return, (y - fit)/err
end
function blue_triplet, p, x
  Common Dispersions, dispers_b, dispers_r, lines_b, lines_r
  z_3722 = (x - (p[3] + (lines_b[0] - lines_b[1]) / dispers_b)) / p[4]
  z_3726 = (x - (p[3] + (lines_b[1] - lines_b[1]) / dispers_b)) / p[4]
  z_3729 = (x - (p[3] + (lines_b[2] - lines_b[1]) / dispers_b)) / p[4]
  return, p[0]*exp(-z_3722^2/2d) + p[1]*exp(-z_3726^2/2d) + p[2]*exp(-z_3729^2/2d); + p[5]
end

function Fit_blue_doublet, p, x=x, y=y, err=err, fit=fit
  Common Dispersions, dispers_b, dispers_r, lines_b, lines_r
  ;P[0] = Gaussian amplitudes of the 4069 SIII line
  ;P[1] = Gaussian amplitudes of the 4076 OII line
  ;P[2] = Index of the middle 4069 peak. 
  ;P[3] = Gaussian Sigma 
  z_4069 = (x - (p[2] + (lines_b[3] - lines_b[3]) / dispers_b)) / p[3]
  z_4076 = (x - (p[2] + (lines_b[4] - lines_b[3]) / dispers_b)) / p[3]
  fit = p[0]*exp(-z_4069^2/2d) + p[1]*exp(-z_4076^2/2d); + p[4]
  return, (y - fit)/err
end
function blue_doublet, p, x
  Common Dispersions, dispers_b, dispers_r, lines_b, lines_r
  z_4069 = (x - (p[2] + (lines_b[3] - lines_b[3]) / dispers_b)) / p[3]
  z_4076 = (x - (p[2] + (lines_b[4] - lines_b[3]) / dispers_b)) / p[3]
  return, p[0]*exp(-z_4069^2/2d) + p[1]*exp(-z_4076^2/2d); + p[4]
end
function chunk_subtract, shift
  Common Chunks, torus_chunk, sky_chunk
  t = total(torus_chunk, 2)
  s = smart_shift( total(sky_chunk, 2), shift[0], 0, /interp) ;Hack??? Should probably add , Missing = !values.f_nan here!!! but that'll need testing!
  histogauss, t/s, a, /NOPlot
  return, a[2]
end

;-----------------------------------------------------------------------------------------------------------------------------------------------
;************************************************************Perform the sky subtraction*******************************************************-
;-----------------------------------------------------------------------------------------------------------------------------------------------
pro sky_sub, Torus_image, sky_image, channel, flat, err_flat, debug=debug

  Common sky, Sky_frame, Torus_frame, which_channel, debug_set
  Common smooth, x_smooth, y_smooth
  Common Apertures, aperture_r, aperture_b
  Common Chunks, torus_chunk, sky_chunk
  
  if keyword_set(debug) then debug_set = 1 else debug_set = 0
  which_channel = channel
  
  print,'*********************************************'
  print,' Starting '+which_channel+' Sky Subtraction '
  print,'*********************************************'
  
  ;*************************************************
  ;       Info for Background Subtraction Smoothing
  ;*************************************************
  x_smooth = 200.0d           ;200 from original program
  y_smooth = 0.0d             ;no smoothing in y
  
  ;*******************************************************
  ;                       Images
  ;********************************************************
  
  Torus_Frame = MRDFITS(Torus_image, 0, torus_header, /Dscale, /silent )
  Sky_Frame = MRDFITS(sky_image, 0, sky_header, /Dscale, /silent )
  err_Torus_Frame = MRDFITS(Torus_image, 1, /Dscale, /silent )
  err_Sky_Frame = MRDFITS(sky_image, 1, /Dscale, /silent )

  ;**********************************************************************************
  ;   Use AMOEBA to find the shift between Torus and Sky frames. CANNOT BE RELIABLY PREDICTED
  ;**********************************************************************************

  N_chunks = 60
  chunk_shifts = make_array(N_chunks, value = !values.F_Nan)
  chunk_size = floor(1078 / N_chunks)
  for i = 0, N_chunks - 1 do begin
    T = smooth(Torus_frame, [x_smooth, y_smooth], /EDGE_MIRROR, /NAN)
    S = smooth(Sky_frame, [x_smooth, y_smooth], /EDGE_MIRROR, /NAN)
    s_frame = (T / S)*sky_frame
    
    ;don't seach for optimimum shifts outside the data aperture
    if which_channel eq 'red' then ind = array_indices(T, aperture_r); ind = cgsetDifference(chunk_ind, aperture_r)
    if which_channel eq 'blue' then ind = array_indices(T, aperture_b); ind = cgsetDifference(chunk_ind, aperture_b)
    Y_sampled = mean( [i*chunk_size, (i+1)*chunk_size])
    junk = where(ind[1,*] eq fix(Y_sampled), /NULL, Count)
    if count eq 0 then continue 
    
    torus_chunk = torus_frame[*, i*chunk_size : (i+1)*chunk_size]
    sky_chunk = s_frame[*, i*chunk_size : (i+1)*chunk_size]
    chunk_shifts[i] = AMOEBAX(1.e-3, 1.e-4, function_name='chunk_subtract', SCALE = 2.5, P0 = [-0.2], FUNCTION_VALUE = fval)  
  endfor
  final_shift = median(chunk_shifts) ;taking the median of the shifts helps with moon contaiminations

  ;*******************************************************
  ;                          Shift sky
  ;*******************************************************

  if abs(final_shift) gt 2.5 then final_shift = 0. ;catch crazies  
  shift_x = final_shift
  print, 'Shifted Jupiter frame by ', shift_x
  Sky_frame_shifted = smart_shift(Sky_frame, shift_x, 0, /interp, Missing = !values.f_nan)
  err_Sky_frame_shifted = smart_shift(err_Sky_frame, shift_x, 0., /interp, Missing = !values.f_nan)

  ;********************************************************
  ;                       1st Sky subtraction
  ;********************************************************

  Torus_smooth = smooth(torus_frame,[x_smooth,y_smooth],/EDGE_MIRROR,/NAN)
  Sky_smooth = smooth(Sky_frame_shifted,[x_smooth,y_smooth],/EDGE_MIRROR,/NAN)  
  scale_factor = Torus_smooth / Sky_smooth
  TorusSkySub_Frame = Torus_Frame - scale_factor * Sky_frame_shifted
  ;*******************************************************
  ;                       2nd Sky subtraction
  ;*******************************************************

  ; Re-perform step #1 on the residual, the simple subtraction of the two images
  TorusSky_Frame = (Torus_Frame - TorusSkySub_Frame)
  Torus_smooth = smooth(torussky_frame, [x_smooth,y_smooth], /EDGE_MIRROR, /NAN)
  Sky_smooth = smooth(Sky_frame_shifted, [x_smooth,y_smooth], /EDGE_MIRROR, /NAN)
  scale_factor = Torus_smooth / Sky_smooth
  TorusSkySub_FrameFinal = Torus_Frame - scale_factor * Sky_frame_shifted
  err_TorusSkySub_FrameFinal = sqrt(err_Torus_Frame^2. + (scale_factor * err_Sky_frame_shifted)^2.)

  ;***************************************************************
  ;            Flat Fielding Post Subtraction has best results, FLAT FIELD NOW
  ;***************************************************************

  TorusSkySub_FrameFinal_FF = TorusSkySub_FrameFinal / Flat
  err_TorusSkySub_FrameFinal_FF = abs(TorusSkySub_FrameFinal_FF) * sqrt( (err_TorusSkySub_FrameFinal  / TorusSkySub_FrameFinal)^2. + (err_flat / flat)^2.)
  TorusSkySub_FrameFinal = temporary(TorusSkySub_FrameFinal_FF) & err_TorusSkySub_FrameFinal = temporary(err_TorusSkySub_FrameFinal_FF)

    if keyword_set(debug_set) then begin
      if which_channel eq 'red' then ind = array_indices(TorusSkySub_FrameFinal, aperture_r)
      if which_channel eq 'blue' then ind = array_indices(TorusSkySub_FrameFinal, aperture_b)
        window, 2, xs = max(ind[0,*]) - min(ind[0,*]), ys = max(ind[1,*]) - min(ind[1,*])
        tv, bytscl(TorusSkySub_FrameFinal[ min(ind[0,*]):max(ind[0,*]), min(ind[1,*]):max(ind[1,*]) ], 0, 100)
      ;if which_channel eq 'blue' then stop
      ;if which_channel eq 'red' then stop
    endif

  sxaddpar, header1, 'HISTORY','Background Subtracted: Twice'
  Sky_Subtracted_Filename = STRMID(Torus_image, 0, strpos(Torus_image, '.fits' )) + 'SS.fits' 
  MWRFITS, TorusSkySub_FrameFinal, Sky_Subtracted_Filename, header1, /create, /silent
  err_header = ["COMMENT Extension 1 is the error post sky subtraction", "Extension = 1"] 
  MWRFITS, err_TorusSkySub_FrameFinal, Sky_Subtracted_Filename, err_header, /silent ;Append the fits file with the error -> extension 1

  print,'*******************************************'
  print,'Done with '+which_channel+' Sky Subtraction'
  print,'*******************************************'
end ; sky_sub

;******************************************************************************************
;******************************************************************************************
;;*****************************************************************************************
;******************************************************************************************
;                       APO_Torus_Reduce: Main Program
;******************************************************************************************
;******************************************************************************************
;;*****************************************************************************************
;******************************************************************************************

pro APO_Torus_redux_new, directory, part, debug = debug

Common Apertures, aperture_r, aperture_b
Common Tilt_Coeffs, Tilt_Coeffs_r, Tilt_Coeffs_b
Common rotator_shifts, telrot_array_b, final_shift_array_b, telrot_array_r, final_shift_array_r
Common Dispersions, dispers_b, dispers_r, lines_b, lines_r

;*********************************INPUTS !!!!!!!!************************************
plate_scale_r = mean([0.3898,0.3969,0.3951,0.3965]) ;arcsec / pixel measured from Io, Europa and Callisto separations in frames UT141101/Jupiter.0020-21r.fits
plate_scale_b = mean([0.4195,0.4166,0.4173,0.4221]) ;arcsec / pixel measured from Io, Europa and Callisto separation on UT141101/Jupiter.0020-21b.fits
dispersion_b = .617  ;rough dispserion in angstroms / pixel
dispersion_r = .581  ;rough dispersion in angstroms / pixel
smooth_bright_by = .1      ;Smooth by a running average of (R_J) for the bright SII doublet 
smooth_by = .2             ;Smooth by a running average of (R_J) 
lines_b = [3721.69, 3726.04, 3728.80, 4068.60, 4076.35] 
;lines_b = [3721.69, 3726.04, 3728.80, 4068.60, 3735.] ;offband check
lines_r = [6312.060, 6716.440, 6730.815] ;OII at 7319.92 & 7330.19 is extremely faint --- No clear detection
;lines_r = [6312.060, 6716.440, 6730.815, 7319.92, 7330.19]
Red_colors = ['orange red', 'red', 'dark red']  
Blue_colors = ['Violet', 'Blue Violet', 'Dodger Blue', 'Green', 'Dark Green']  

;Define the directory specific information, like frames used for calibration of the mask throughput, which frames to co-add, etc.
Load_directory: 
CASE 1 OF
  directory eq 'D:\DATA\Apache Point Data\UT131107\': begin
      Jupiter_On_Mask = 'Jupiter1.0010'
      Jupiter_Off_Mask = 'Jupiter1.0011'
      Torus = [5,6,8,9] 
      Sky = [7,10]
      CoAdd_frames = [[0,1,2,3]] ;frame 3 bad
      ;CoAdd_frames = [[0,1,2]] ;frame 3 bad
      key = 1 ;Image to key off of and shift everything else to
      dispersion_b = .617 
      dispersion_r = .58 
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1090.      
   end
   directory eq 'D:\DATA\Apache Point Data\UT131109\': begin
      Jupiter_On_Mask = 'Torus.0009'
      Jupiter_Off_Mask = 'Jupiter.0008'
      Torus = [6,7,8,9,10] 
      Sky = [11,12]
      key = 1 ;Image to key off of and shift everything else to
      CoAdd_frames = [[0,1,2,3,4]]
      ;CoAdd_frames = [[0,1,2,4]] ;frame 3 bad
      dispersion_b = .617 
      dispersion_r = .58
      Center_WL_pixel_b = 1010.
      ;Center_WL_pixel_r = 1082. 
      Center_WL_pixel_r = 1090.      
   end
   directory eq 'D:\DATA\Apache Point Data\UT131217\': begin
      Jupiter_On_Mask = 'Jupiter.0001'
      Jupiter_Off_Mask = 'Jupiter.0002'
      Torus = [6,7,8,9,10,11,12] 
      Sky = [2,3,4]      
      key = 1 ;Image to key off of and shift everything else to
      CoAdd_frames = [[indgen(7)]]
      dispersion_b = .617 
      dispersion_r = .58 
      Center_WL_pixel_b = 1024.
      Center_WL_pixel_r = 1056.      
   end
   directory eq 'D:\DATA\Apache Point Data\UT131224\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT131107\'   
      Torus = [3,4,5,6,7,8,9,10,11] 
      Sky = [0,1,2]
      key = 1 ;Image to key off of and shift everything else to
      CoAdd_frames = [[indgen(9)]]
      dispersion_b = .617 
      dispersion_r = .58 
      Center_WL_pixel_b = 1024.
      Center_WL_pixel_r = 1060.      
   end
   directory eq 'D:\DATA\Apache Point Data\UT140107\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT131107\'
      Torus = [0,2,3,6,8,10,11,13,14,16,17,19,20,22,23] 
      Sky = [1,5,9,12,15,18,21]
      CoAdd_frames = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
      ;CoAdd_frames = [[6,7,8,9,10,11]] 
      key = 12 ;Image to key off of and shift everything else to
      ;CoAdd_frames = [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]
      dispersion_b = .620
      dispersion_r = .58 
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1030.
   end
   directory eq 'D:\DATA\Apache Point Data\UT140209\': begin
      Jupiter_On_Mask = 'Jupiter_on_mask.0020' 
      Jupiter_Off_Mask = 'Jupiter_off_mask.0021'
      key = 1 ;Image to key off of and shift everything else to
      Torus = [5,6,7] 
      Sky = [0,1]
      CoAdd_frames = [[0,1,2]]
      dispersion_b = .617
      dispersion_r = .58 
      Center_WL_pixel_b = 1012.
      Center_WL_pixel_r = 1060.
   end
   directory eq 'D:\DATA\Apache Point Data\UT140214\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT140209\'
      Torus = [4,5,6,7,8,9,10,11] 
      Sky = [0,1,2,3]
      key = 1 ;Image to key off of and shift everything else to
      CoAdd_frames = [[0,1,2,3,4,5,6,7]]
      dispersion_b = .617
      dispersion_r = .58 
      Center_WL_pixel_b = 1012
      Center_WL_pixel_r = 1080
   end
   directory eq 'D:\DATA\Apache Point Data\UT140216\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT140209\'
      Torus = [2,3,4,5] 
      Sky = [0,1]
      key = 1 ;Image to key off of and shift everything else to
      CoAdd_frames = [[2,3]]
      dispersion_b = .617
      dispersion_r = .58 
      Center_WL_pixel_b = 1012
      Center_WL_pixel_r = 1080
   end
  directory eq 'D:\DATA\Apache Point Data\UT140219\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT140209\'
      Torus = [7,8,9,10,11,12,13,14,15,16,17] 
      Sky = [0,1,2,3,4,5]
      key = 10 ;Image to key off of and shift everything else to
      ;CoAdd_frames = [[1,2,3]]  ;best frames as seen in the stack plot
      CoAdd_frames = [[4,5,6,7,8,9,10]]  ;best frames as seen in the stack plot
      dispersion_b = .617
      dispersion_r = .58 
      Center_WL_pixel_b = 1012
      Center_WL_pixel_r = 1065
   end
   directory eq 'D:\DATA\Apache Point Data\UT141030\': begin
      Jupiter_On_Mask = 'Jupiter.0037' 
      Jupiter_Off_Mask = 'Jupiter_Slit_Calibration.0038'  
      Torus = indgen(11) 
      Sky = [11,12,13,14,15,17]
      CoAdd_frames = [[indgen(9)]] 
      key = 0 ;Image to key off of and shift everything else to
      ;CoAdd_frames = [[0,1],[2,3],[4,5],[6,7],[8,9]]
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1060.
   end
   directory eq 'D:\DATA\Apache Point Data\UT141101\': begin 
      ;CoAdd_frames = [[2,3],[4,5],[6,7],[8,9],[10,11]] ;frame
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT141030\'
      Torus = indgen(12) 
      Sky = [12,13,14,15,16,17]
      CoAdd_frames = [[0,1,2,3,4,5,6,8,9,10,11]] ;frames 0,1,7 bad
      key = 7 ;Image to key off of and shift everything else to 
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1060.
   end
   directory eq 'D:\DATA\Apache Point Data\UT150327\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT160406\' ;Actual Mask Unknown!
      Torus = [4,5,6,7,8,9,10] ;ingnore 5 min exposures unless really needed
      Sky = [12,13,14]
      ;CoAdd_frames = [[0,1,2,3,4,5,6]]
      CoAdd_frames = [[5,6,7,8,9,10]]
      key = 0 ;Image to key off of and shift everything else to
      dispersion_b = .617 
      dispersion_r = .58 
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1060.
   end
   directory eq 'D:\DATA\Apache Point Data\UT160406\': begin
      Jupiter_On_Mask = 'Torus.0009' 
      Jupiter_Off_Mask = 'Jupiter_Off_Mask.0006'  
      Torus = indgen(17)+9 
      Sky = [0,1,2,3,4,5,6]
      ;CoAdd_frames = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]] ;Moonless Frames
      CoAdd_frames = [indgen(17)]
      key = 8
      Center_WL_pixel_b = 1013.
      Center_WL_pixel_r = 1090.
   end
   directory eq 'D:\DATA\Apache Point Data\UT160414\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT160406\'
      Torus = indgen(12) + 3
      Sky = indgen(3)
      ;CoAdd_frames = [[0,1,2,3],[4,5,6,7],[8,9,10,11]] ;Moonless Frames
      CoAdd_frames = [[0,1,2,3,4,5,6,7,8,9,10,11]] ;Moonless Frames
      key = 8
      Center_WL_pixel_b = 1014.
      Center_WL_pixel_r = 1085.
   end
   directory eq 'D:\DATA\Apache Point Data\UT160415\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT160406\'
      Torus = indgen(8) + 2
      Sky = indgen(2)
      ;CoAdd_frames = [[0,1,2,3],[4,5,6,7]] ;Moonless Frames
      CoAdd_frames = [[0,1,2,3,4,5,6,7]] ;Moonless Frames
      key = 2
      Center_WL_pixel_b = 1014.
      Center_WL_pixel_r = 1090.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170310\': begin
      Jupiter_On_Mask = 'Torus.0019' 
      Jupiter_Off_Mask = 'Off-Mask.0003'  
      Torus = indgen(18)+5
      Sky = [0,1,2,3]
      CoAdd_frames = [[indgen(18)]]
      key = 8
      Center_WL_pixel_b = 978.
      Center_WL_pixel_r = 1092.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170325\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT170310\' ;for now, there are mask cals this day that should be used though!
      Torus = indgen(19)+5
      Sky = [0,1,2,3,4]
      CoAdd_frames = [[indgen(19)]]
      key = 8
      Center_WL_pixel_b = 1010.
      Center_WL_pixel_r = 1092.
   end
  directory eq 'D:\DATA\Apache Point Data\UT170331\': begin
      Jupiter_On_Mask = 'Torus.0014' 
      Jupiter_Off_Mask = 'Jupiter_Off_Mask_Calibration.0001'  
      Torus = indgen(13)+7 ; there are 14 in fact, but the last has a cosmic ray right on top of the ribbon
      Sky = [0,1,2]
      CoAdd_frames = [[indgen(13)]] 
      key = 8
      Center_WL_pixel_b = 1004.
      Center_WL_pixel_r = 1099.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170511\': begin
      Jupiter_On_Mask = 'Torus.0020' 
      Jupiter_Off_Mask = 'Off_mask.0004'  
      Torus = indgen(21)+9 
      Sky = [0,1,2,3,4,5]
      CoAdd_frames = [[indgen(21)]] 
      key = 8
      Center_WL_pixel_b = 1028.
      Center_WL_pixel_r = 1050.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170519\': begin
      Jupiter_On_Mask = 'Torus.0009' 
      Jupiter_Off_Mask = 'Off_Mask_Calibration.0001'  
      Torus = indgen(8)+2 
      Sky = [0,1]
      CoAdd_frames = [[indgen(8)]] 
      key = 4
      Center_WL_pixel_b = 1022.
      Center_WL_pixel_r = 1080.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170704\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT170519\'
      Torus = indgen(5)+1 
      Sky = [0]
      CoAdd_frames = [[indgen(5)]] 
      key = 4
      Center_WL_pixel_b = 1022.
      Center_WL_pixel_r = 1050.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170707\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT170519\'
      Torus = indgen(12)+3 
      Sky = [0,1,2]
      CoAdd_frames = [[indgen(12)]] 
      key = 4
      Center_WL_pixel_b = 1022.
      Center_WL_pixel_r = 1050.
   end
   directory eq 'D:\DATA\Apache Point Data\UT170714\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT170519\'
      Torus = indgen(7)+2 
      Sky = [0,1]
      CoAdd_frames = [[indgen(7)]] 
      key = 4
      Center_WL_pixel_b = 1026.
      Center_WL_pixel_r = 1074.
   end
   directory eq 'D:\DATA\Apache Point Data\UT180309\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT180322\'
      Torus = indgen(20)+6
      Sky = indgen(6) 
      CoAdd_frames = [[indgen(7)]] 
      key = 4
      Center_WL_pixel_b = 1026.
      Center_WL_pixel_r = 1074.
   end
   directory eq 'D:\DATA\Apache Point Data\UT180320\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT180322\'
      Torus = indgen(12)+4 
      Sky = [0,1,2,3]
      CoAdd_frames = [[indgen(7)]] 
      key = 4
      Center_WL_pixel_b = 988.
      Center_WL_pixel_r = 1118.
   end
   directory eq 'D:\DATA\Apache Point Data\UT180322\': begin
      Jupiter_On_Mask = 'Torus.0001' 
      Jupiter_Off_Mask = 'Off-Mask.0004'  
      Torus = indgen(22)+5
      Sky = indgen(5)
      CoAdd_frames = [[indgen(8)]] 
      key = 4
      Center_WL_pixel_b = 1022.
      Center_WL_pixel_r = 1074.
   end
   directory eq 'D:\DATA\Apache Point Data\UT180330\': begin
      Mask_Cal_Directory = 'D:\DATA\Apache Point Data\UT180322\'
      Torus = indgen(22)+6 
      Sky = indgen(6)
      CoAdd_frames = [[indgen(7)]] 
      key = 4
      Center_WL_pixel_b = 1068.
      Center_WL_pixel_r = 1032.
   end       
ENDCASE

;************************************************************************************

;setup P3D for use with the rectification step
defsysv,'!p3d_path','C:\IDL\Io\Apache_Point_Programs\p3d-2.2.6\'
setenv, 'p3d_path=C:\IDL\Io\Apache_Point_Programs\p3d-2.2.6\'    
    
;Load generic SPICE kernels

        ; Clean any lingering kernels out of memory here:
          cspice_ktotal, 'all', count
          Print, 'Deleting ', strtrim(string(count),2), ' old SPICE kernels from memory'
          i=0
          while i lt count do begin
            cspice_kdata, 0, 'all', file, type, source, handle, found
            cspice_unload, file
            i=i+1
          endwhile
        
        ; Load New Kernels
          CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\lsk\naif0010.tls')         ; leap seconds kernel
          CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\pck\pck00010.tpc')         ; Planet rotational states
          CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Jupiter_System\jup309.bsp')   
          ;CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\pck\de-403-masses.tpc')    ; Planet Masses
          CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\planets\de421.bsp')    ; SPK (ephemeris kernel) for planets
          CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\satellites\sat319.bsp'); SPK (ephemeris kernel) for satellites 
          cspice_ktotal, 'all', count
          Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'

If part eq 1 then begin ;bias / flat / CR
    
    ;*****************************build master bias frames***********************************************
    blue_biases = FILE_SEARCH(strcompress(Directory + 'bias*b.fits'), count = n_bias_b)
    print, n_bias_b, ' blue bias frames found, building master bias. . .'  
    bigarray = fltarr(n_bias_b,2098,1078)
    for i = 0, n_bias_b - 1 do begin
      bigarray[i,*,*] = MRDFITS(blue_biases[i], 0, header, /Dscale, /silent )
    endfor  
    ;Take the pixel-by-pixel median of the bias images in the list
    Blue_bias = MEDIAN(bigarray, DIMENSION=1, /double, /even)
    MWRFITS, Blue_bias, directory+'Processed\Master_blue_bias.fits', header, /create, /silent
  
    red_biases = FILE_SEARCH(strcompress(Directory + 'bias*r.fits'), count = n_bias_r)
    print, n_bias_r, ' red bias frames found, building master bias. . .'  
        bigarray = fltarr(n_bias_r,2098,1078)
    for i = 0, n_bias_r - 1 do begin
      bigarray[i,*,*] = MRDFITS(red_biases[i], 0, header, /Dscale, /silent )
    endfor  
    ;Take the pixel-by-pixel median of the bias images in the list
    red_bias = MEDIAN(bigarray, DIMENSION=1, /double, /even)
    MWRFITS, red_bias, directory+'Processed\Master_red_bias.fits', header, /create, /silent

    ;*****************************build master arcs frames***********************************************
    blue_arcs = FILE_SEARCH(strcompress(Directory + 'Arcs.*b.fits'), count = n_images)                                ;regular blue arc frames
    blue_arcs_long = FILE_SEARCH(strcompress(Directory + 'Arcs_long.*b.fits'), count = n_images_long)                  ;long exposure time blue arc frames
    blue_arcs_long_wide = FILE_SEARCH(strcompress(Directory + 'Arcs_Long_Wide*b.fits'), count = n_images_long_wide)   ;long exposure & Narrow slit mask blue arc frames (for Mercury observation dates only)
    print, n_images, ' blue arcs frames found, building master arcs. . .'  
    print, n_images_long, ' longer expsoure blue arcs frames found, building master arcs. . .'  
    print, n_images_long_wide, ' longer expsoure blue arcs frames found with a narrow slit mask, building master arcs. . .'  
    bigarray = fltarr(n_images,2098,1078)
    for i = 0, N_images - 1 do begin
      bigarray[i,*,*] = MRDFITS(blue_arcs[i], 0, header, /Dscale, /silent )
    endfor  
    ;Take the pixel-by-pixel median of the arcs images in the list
    Blue_arcs = MEDIAN(bigarray, DIMENSION=1, /double, /even) - blue_bias
    MWRFITS, blue_arcs, directory+'Processed\Master_blue_arcs.fits', header, /create, /silent   
        if n_images_long gt 0 then begin
          bigarray = fltarr(n_images_long,2098,1078)
          for i = 0, n_images_long - 1 do begin
            bigarray[i,*,*] = MRDFITS(blue_arcs_long[i], 0, header, /Dscale, /silent )
          endfor  
          ;Take the pixel-by-pixel median of the arcs images in the list
          Blue_arcs_long = MEDIAN(bigarray, DIMENSION=1, /double, /even) - blue_bias
          MWRFITS, blue_arcs_long, directory+'Processed\Master_blue_arcs_long.fits', header, /create, /silent
        endif  
        if n_images_long_wide gt 0 then begin
          bigarray = fltarr(n_images_long,2098,1078)
          for i = 0, n_images_long_wide - 1 do begin
            bigarray[i,*,*] = MRDFITS(blue_arcs_long_wide[i], 0, header, /Dscale, /silent )
          endfor  
          ;Take the pixel-by-pixel median of the arcs images in the list
          Blue_arcs_long_wide = MEDIAN(bigarray, DIMENSION=1, /double, /even) - blue_bias
          MWRFITS, blue_arcs_long_wide, directory+'Processed\Master_blue_arcs_long.fits', header, /create, /silent
        endif  

    red_arcs = FILE_SEARCH(strcompress(Directory + 'Arcs*r.fits'), count = n_images) 
    Potassium_arcs = FILE_SEARCH(strcompress(Directory + 'Arcs*7500*r.fits'), count = n_Potassium_images)
    if n_Potassium_images gt 0 then begin
      match, Potassium_arcs, red_arcs, reg_type, K_type
      keep = cgsetdifference(indgen(n_images), fix(K_type))
      n_images = n_elements(keep) & red_arcs = red_arcs[keep] ;only us the proper arc frames, not potassium for Mercury
    endif
    if directory eq 'D:\DATA\Apache Point Data\UT141030\' then $ ;the wide slit mask for Jupiter, narrow slit was for Mercury
      red_arcs = FILE_SEARCH(strcompress(Directory + 'Arcs_Wide*r.fits'), count = n_images)
    print, n_images, ' red arcs frames found, building master arcs. . .'  
    bigarray = fltarr(n_images, 2098, 1078)
    for i = 0, N_images - 1 do begin
      bigarray[i,*,*] = MRDFITS(red_arcs[i], 0, header, /Dscale, /silent )
    endfor  
    ;Take the pixel-by-pixel median of the arc images in the list
    red_arcs = MEDIAN(bigarray, DIMENSION=1, /double, /even) - red_bias
    MWRFITS, red_arcs, directory+'Processed\Master_red_arcs.fits', header, /create, /silent

    ;*****************************build master flat frames***********************************************
    blue_flats = FILE_SEARCH(strcompress(Directory + 'flat*b.fits'), count = n_flat_b)
    imaging_flats = FILE_SEARCH(strcompress(Directory + 'flat*6731*b.fits'), count = n_imaging_flats)
    if n_imaging_flats gt 0 then begin
      match, imaging_flats, blue_flats, junk, imaging_type
      keep = cgsetdifference(indgen(n_flat_b), fix(imaging_type))
      n_flat_b = n_elements(keep) & blue_flats = blue_flats[keep] ;only us the proper arc frames, not potassium for Mercury
    endif
    if directory eq 'D:\DATA\Apache Point Data\UT141030\' then $ ;the wide slit mask for Jupiter, narrow slit was for Mercury
      blue_flats = FILE_SEARCH(strcompress(Directory + 'Flats_Wide*b.fits'), count = n_flat_b)
    print, n_flat_b, ' blue flat frames found, building master flat. . .'  
    bigarray = fltarr(n_flat_b,2098,1078)
    for i = 0, n_flat_b - 1 do begin
      bigarray[i,*,*] = MRDFITS(blue_flats[i], 0, header, /Dscale, /silent )
    endfor  
    ;Take the pixel-by-pixel median of the flat images in the list
    Blue_flat = MEDIAN(bigarray, DIMENSION=1, /double, /even) - blue_bias
    Gain = float(sxpar(header, 'GAIN'))
    RDnoise = float(sxpar(header, 'RDNOISE'))
    err_blue_flat = sqrt( (!pi*blue_flat / (2.*gain*n_flat_b)) + RDnoise^2. ) ;median's increase noise by sqrt(pi/2)
    MWRFITS, Blue_flat, directory+'Processed\Blue_Master_flat.fits', header, /create, /silent
    MWRFITS, err_blue_flat, directory+'Processed\Blue_Master_flat.fits', /silent ;Append the fits file with the error -> extension 1
  
    red_flats = FILE_SEARCH(strcompress(Directory + 'flat*r.fits'), count = n_flat_r)
    Potassium_flats = FILE_SEARCH(strcompress(Directory + 'flat*7500*r.fits'), count = n_Potassium_flats)
    imaging_flats = FILE_SEARCH(strcompress(Directory + 'flat*6731*r.fits'), count = n_imaging_flats)
    if n_imaging_flats + n_Potassium_flats gt 0 then begin
      match, [imaging_flats, Potassium_flats] , red_flats, junk, imaging_type
      keep = cgsetdifference(indgen(n_flat_r), fix(imaging_type)) 
      n_flat_r = n_elements(keep) & red_flats = red_flats[keep] ;only us the proper arc frames, not potassium for Mercury
    endif
    if directory eq 'D:\DATA\Apache Point Data\UT141030\' then $ ;the wide slit mask for Jupiter, narrow slit was for Mercury
      red_flats = FILE_SEARCH(strcompress(Directory + 'Flats_Wide*r.fits'), count = n_flat_r)
    if directory eq 'D:\DATA\Apache Point Data\UT160414\' then $ ;the wide slit mask for Jupiter, narrow slit was for Mercury
      red_flats = FILE_SEARCH(strcompress(Directory + 'Flat.*r.fits'), count = n_flat_r)
    print, n_flat_r, ' red flat frames found, building master flat. . .'  
    bigarray = fltarr(n_flat_r, 2098, 1078)
    for i = 0, n_flat_r - 1 do begin
      bigarray[i,*,*] = MRDFITS(red_flats[i], 0, header, /Dscale, /silent)
    endfor  
    ;Take the pixel-by-pixel median of the flat images in the list
    red_flat = MEDIAN(bigarray, DIMENSION=1, /double, /even) - red_bias
    Gain = float(sxpar(header, 'GAIN'))
    RDnoise = float(sxpar(header, 'RDNOISE')) 
    err_red_flat = sqrt( (!pi*red_flat / (2.*gain*n_flat_r)) + RDnoise^2. )   ;median's increase noise by sqrt(pi/2)
    MWRFITS, red_flat, directory+'Processed\Red_Master_flat.fits', header, /create, /silent
    MWRFITS, err_red_flat, directory+'Processed\Red_Master_flat.fits', /silent ;Append the fits file with the error -> extension 1

    ;If necessary, bias and CR correct the frames with Jupiter off mask used for mask density calibration 
    if keyword_set(Jupiter_Off_Mask) then begin
      filename = directory+Jupiter_Off_Mask+'b.fits' 
      BS = MRDFITS(filename, 0, header, /Dscale, /silent ) - blue_bias 
      Gain = float(sxpar(header, 'GAIN'))
      RDnoise = float(sxpar(header, 'RDNOISE'))
      new_filename = STRMID(filename, strlen(directory))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      SXADDPAR, Header, 'PIXSCAL2', plate_scale_b, 'Plate scale estimate from moons (arcsec/pixel)' ;write over the default (inaccurate) platescale
      MWRFITS, BS, directory+'Processed\' + new_filename + '.BS.fits', header, /create, /silent
      la_cosmic, directory+'Processed\' + new_filename + '.BS.fits', outsuff = "CR", sigclip = 6., $
        statsec = sxpar( header, 'DETSIZE'), readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.    
      BSCR = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
      err_BSCR = sqrt( (BSCR / gain) + RDnoise^2. )  
      MWRFITS, err_BSCR, directory+'Processed\' + new_filename + '.BSCR.fits', /silent ;Append the fits file with the error -> extension 1
      filename = directory+Jupiter_Off_Mask+'r.fits' 
      BS = MRDFITS(filename, 0, header, /Dscale, /silent ) - blue_bias 
      Gain = float(sxpar(header, 'GAIN'))
      RDnoise = float(sxpar(header, 'RDNOISE'))
      new_filename = STRMID(filename, strlen(directory))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      SXADDPAR, Header, 'PIXSCAL2', plate_scale_b, 'Plate scale estimate from moons (arcsec/pixel)' ;write over the default (inaccurate) platescale
      MWRFITS, BS, directory+'Processed\' + new_filename + '.BS.fits', header, /create, /silent
      la_cosmic, directory+'Processed\' + new_filename + '.BS.fits', outsuff = "CR", sigclip = 6., $
        statsec = sxpar( header, 'DETSIZE'), readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.    
      BSCR = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
      err_BSCR = sqrt( (BSCR / gain) + RDnoise^2. )  
      MWRFITS, err_BSCR, directory+'Processed\' + new_filename + '.BSCR.fits', /silent ;Append the fits file with the error -> extension 1
    endif

    ;*************************cosmic ray correct and bias subtract the data******************************
    Images_b = FILE_SEARCH(strcompress(Directory + '{Torus,Jupiter}*b.fits'), count = n_images) 
    for i = 0, N_images - 1 do begin
      filename = Images_b[i]
      BS = MRDFITS(filename, 0, header, /Dscale, /silent ) - blue_bias  
      Gain = float(sxpar(header, 'GAIN'))
      RDnoise = float(sxpar(header, 'RDNOISE'))
      new_filename = STRMID(filename, strlen(directory))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      SXADDPAR, Header, 'PIXSCAL2', plate_scale_b, 'Plate scale estimate from moons (arcsec/pixel)' ;write over the default (inaccurate) platescale
      MWRFITS, BS, directory+'Processed\' + new_filename + '.BS.fits', header, /create, /silent
      la_cosmic, directory+'Processed\' + new_filename + '.BS.fits', outsuff = "CR", sigclip = 6., $
        statsec = sxpar( header, 'DETSIZE'), readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.    
      BSCR = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
      err_BSCR = sqrt( (BSCR / gain) + RDnoise^2. )  
      MWRFITS, err_BSCR, directory+'Processed\' + new_filename + '.BSCR.fits', /silent ;Append the fits file with the error -> extension 1                          
      if keyword_set(debug) then begin
        window, 0, xs = 2098, ys = 1000
        image = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
        err_image = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 1, /Dscale, /silent )
        tv, bytscl(image, 0, 150)  
      endif ;debug      
    endfor  
    Images_r = FILE_SEARCH(strcompress(Directory + '{Torus,Jupiter}*r.fits'), count = n_images) 
    for i = 0, N_images - 1 do begin
      filename = Images_r[i]
      BS = MRDFITS(filename, 0, header, /Dscale, /silent ) - red_bias
      Gain = float(sxpar(header, 'GAIN'))
      RDnoise = float(sxpar(header, 'RDNOISE'))   
      new_filename = STRMID(filename, strlen(directory))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      SXADDPAR, Header, 'PIXSCAL2', plate_scale_r, 'Plate scale estimate from moons (arcsec/pixel)' ;write over the default (inaccurate) platescale
      MWRFITS, BS, directory+'Processed\' + new_filename + '.BS.fits', header, /create, /silent
      la_cosmic, directory+'Processed\' + new_filename + '.BS.fits', outsuff = "CR", sigclip = 6., $       
        statsec = sxpar( header, 'DETSIZE'), readn=RDnoise, gain=gain  ;best parameters are tested, slow n good.
      BSCR = MRDFITS(directory+'Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
      err_BSCR = sqrt( (BS / gain) + RDnoise^2.)  
      MWRFITS, err_BSCR, directory+'Processed\' + new_filename + '.BSCR.fits', /silent ;Append the fits file with the error -> extension 1
      if keyword_set(debug) then begin
        window, 1, xs = 2098, ys = 1078
        image = MRDFITS(directory+'Processed/' + new_filename + '.BS-mask.fits', 0, header, /Dscale, /silent )
        tv, bytscl(image)  
      endif  
    endfor  
endif ;part 1  

if part eq 2 then begin ;find torus vs sky frames / sky subtract
  Images_r = FILE_SEARCH(strcompress(Directory + '{Torus,Jupiter}*r.fits'), count = n_images) 
  Images_b = FILE_SEARCH(strcompress(Directory + '{Torus,Jupiter}*b.fits'), count = n_images) 
  
  ;figure out which frames are aligned with the torus and which are for Jovian scattered light subtraction     
  for i= 0, n_elements(Images_r) - 1 do begin
    test = MRDFITS(Images_r[i], 0, header, /Dscale, /silent )
    print, i, sxpar(header, 'OBJANGLE'), '    ', sxpar(header, 'DATE-OBS'), sxpar(header, 'EXPTIME'), '   ', Images_r[i], '  ----  Orig Filename = ', sxpar(header, 'FILENAME')
  endfor
  if keyword_set(debug) then stop ;stop here to find indices of torus, and of jupiter pivoted off torus frames

  ;find the mask edges with a flat field
      mask_coords_b = intarr(2) & mask_coords_r = intarr(2)
      blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 0, flat_header_b, /Dscale, /silent )
      red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 0, flat_header_r, /Dscale, /silent )
      err_blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 1, err_flat_header_b, /Dscale, /silent )
      err_red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 1, err_flat_header_r, /Dscale, /silent )   
      aperture_r = where(red_flat gt mean(red_flat)) 
        dummy = red_flat       
        y_aperture_ind = where( total(red_flat,1) gt mean(total(red_flat,1)), /NULL )
        ind = array_indices(red_flat, aperture_r) 
        dummy[0:max(ind[0,*]),y_aperture_ind] = -100
        aperture_r = where(dummy eq -100, /NULL, complement = anti_aperture_red) ;a much better aperture.
      aperture_b = where(blue_flat gt mean(blue_flat)) 
        dummy = blue_flat       
        y_aperture_ind = where( total(blue_flat,1) gt mean(total(blue_flat,1)), /NULL )
        ind = array_indices(blue_flat, aperture_b) 
        dummy[0:max(ind[0,*]),y_aperture_ind] = -100
        aperture_b = where(dummy eq -100, /NULL, complement = anti_aperture_blue) ;a much better aperture.
      profile_r = total(red_flat, 1)
      profile_b = total(blue_flat, 1)
      mask_r = where(profile_r lt mean(profile_r), /NULL, complement = anti_mask_red)
      mask_r = mask_r[where((mask_r gt 350) and (mask_r lt 800))]
      mask_coords_r = [min(mask_r), max(mask_r)]
      mask_b = where(profile_b lt mean(profile_b), /NULL, complement = anti_mask_blue)
      mask_b = mask_b[where((mask_b gt 200) and (mask_b lt 700))]
      mask_coords_b = [min(mask_b), max(mask_b)]
      aperture_coords_b = ARRAY_INDICES(blue_flat, aperture_b)  
      lower_section_lower_row_b = min(aperture_coords_b[1,*]) 
      lower_section_upper_row_b = min(mask_b) - 1
      upper_section_lower_row_b = max(mask_b) + 1 
      upper_section_upper_row_b = max(aperture_coords_b[1,*]) 
      aperture_r = where(red_flat gt mean(red_flat))
      aperture_coords_r = ARRAY_INDICES(red_flat, aperture_r)  
      lower_section_lower_row_r = min(aperture_coords_r[1,*]) 
      lower_section_upper_row_r = min(mask_r) - 1
      upper_section_lower_row_r = max(mask_r) + 1 
      upper_section_upper_row_r = max(aperture_coords_r[1,*])   
      
  ;get Jupiter's instantaneous radius in pixels (rough numbers)
      cspice_utc2ET, sxpar(flat_header_r, 'DATE-OBS'), ET    
      cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET, 'IAU_Jupiter', 'LT+S', 'Earth', Sub_Earth, trgepc, srfvec ;get sub observer point in catesian IAU Jupiter coords
      cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
      cspice_spkpos, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', ptarg, ltime
      R_j = 206264.806* atan(radii[0] / norm(ptarg))         
      R_J_per_pixel_b = float(sxpar(flat_header_b, 'PIXSCAL2')) / R_j     
      R_J_per_pixel_r = float(sxpar(flat_header_r, 'PIXSCAL2')) / R_j  
      
  ;Make bad pixel masks from the flats, fix the bad pixel before normalizing the flats
    dummy = red_flat
    dummy[anti_aperture_red] = max(dummy) 
    dummy[where(dummy lt 100., /Null)] = !values.f_NaN ;Call DN less then 100 a bad pixel
    good = where(finite(dummy), complement = bad_pixels_r, /Null)
    red_flat[bad_pixels_r] = 0.;!values.f_NaN
    fixpix, red_flat, finite(dummy), red_flat_fixed, NPIX = 80
    ;tv, bytscl(test, 0, 20000) 
    red_flat = temporary(red_flat_fixed)
  
  ;NORMALIZE THE FLAT FIELDS to a value of 1. (NO FLAT CORRECTION UNDER THE MASK)
    red_flat[ anti_aperture_red ] = !values.F_NaN
    mean_spec = mean(red_flat, dimension = 2, /NAN)
    red_flat[ anti_aperture_red ] = -1.
    p3d_tool_flatfield, red_flat, red_flat_norm, err_red_flat_norm, din = err_red_flat, /NORMALIZE, spec = mean_spec, /ndisp
    red_flat_norm[ anti_aperture_red ] = 1. & err_red_flat_norm[ anti_aperture_red ] = err_red_flat[ anti_aperture_red ]
    red_flat = temporary(red_flat_norm) & err_red_flat = temporary(err_red_flat_norm)
    
    blue_flat[ anti_aperture_blue ] = !values.F_NaN
    mean_spec = mean(blue_flat, dimension = 2, /NAN)
    blue_flat[ anti_aperture_blue ] = -1.
    p3d_tool_flatfield, blue_flat, blue_flat_norm, err_blue_flat_norm, din = err_blue_flat, /NORMALIZE, spec = mean_spec, /ndisp
    blue_flat_norm[ anti_aperture_blue ] = 1. & err_blue_flat_norm[ anti_aperture_blue ] = err_blue_flat[ anti_aperture_blue ]
    blue_flat = temporary(blue_flat_norm) & err_blue_flat = temporary(err_blue_flat_norm)

  ;-------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ;Finetune the flat fields --- On average, the Jupiter scattered light profiles should have a smooth power law falloff with distance from the centroid. Make that so. 

    ;find the blue polynomical fit to Jupiter's center.
    Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '{Torus,Jupiter}*b.BSCR.fits')) 
    n_images = n_elements([sky, torus])  ;just use sky and torus frames to get this correction
    Images_b = Images_b[[sky, torus]] ;just use sky and torus frames to get this correction
    flat_correction_array = fltarr(n_images, 2098, 1078)
    for i= 0, n_images - 1 do begin  
      
      image = MRDFITS(Images_b[i], 0, header, /Dscale, /silent )
      err_image = MRDFITS(Images_b[i], 1, /Dscale, /silent )
      dummy = image & err_dummy = err_image

      ;get a signal to noise estimate under the mask for every spectral bin
      under_mask = [lower_section_upper_row_b, upper_section_lower_row_b] ; be sure your only searching regions well inside the aperture    
      y = total(dummy[*, under_mask[0]:under_mask[1]], 2, /NAN)
      dy = sqrt(total(err_dummy[*, under_mask[0]:under_mask[1]]^2., 2, /NAN))
      x = findgen(2098)
  
      ;Determine rebinning needed to have a uniform signal to noise, and get the pixel fits of Jupiter's centroid
      new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
      Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
      for j = 0, N_elements(new_y) - 1 do begin 
        x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
        compress = where(BIN_map eq j, /null)
        y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
        junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
        ;fit = jjgaussfit(x, y, a, LIMITamp = [1.e1, 2.^16.], LIMITcen = [min(x)+10., max(x)-10.], FIXsl = 0., guess = [max(y), peak, (2./R_J_per_pixel_b)/3., 0., 0.])
        fit = jjgaussfit(x, y, a, LIMITamp = [0., 2.^16.], LIMITcen = [min(x)+10.,max(x)-10.], FIXsl = 0., guess = [max(y), mean(x), (2./R_J_per_pixel_b)/3., 0., 0.])
                if (total(finite(a)) eq 5) then begin
                  disk = where(y gt .1*a[0], /Null)               
                  if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the ceter of a region 10% the peak differ by >20 pixels
                  disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel_b) and (disk lt 1.2/R_J_per_pixel_b + a[1])), /NULL)]
                  if (disk eq !NULL) then continue ;CS Added, 4/7/2018
                  if (abs(median(x[disk]) - a[1]) gt 20.) then continue
                  Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]
                endif
      endfor 
      reject = where(new_x lt 1400, /NULL) ; not enough light for a centroid
      Jupiter_centroid[reject] = !values.F_nan ; not enough light for a centroid
      keep = where(finite(Jupiter_centroid)) 
  
      ;1D fits
      coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
      center = poly(findgen(2098), coeffs) ; pixel location of Jupiter's center
      if keyword_set(debug) then begin
        window, 2, title = 'Fitting Jupiter''s Centroid'
        plot, new_x[keep], Jupiter_centroid[keep], psym = 4, /ynozero
        oplot, findgen(2098), center
      endif
      
      ;2D unwarping: dispersion axis isn't exactly aligned in x, correct it so that it is.
      Xo = [new_x[keep], new_x[keep], new_x[keep]]
      Yo = [make_array(n_elements(keep), value = mean(center) - 300, /float), $
            make_array(n_elements(keep), value = mean(center), /float), $
            make_array(n_elements(keep), value = mean(center) + 300, /float)] 
      Xi = [new_x[keep], new_x[keep], new_x[keep]]
      Yi = [Jupiter_centroid[keep] - 300., Jupiter_centroid[keep], Jupiter_centroid[keep] + 300.]
      POLYWARP, xi, yi, Xo, yo, 2, P_fwd, Q_fwd ; P and Q are the coefficient matrices to warp the tilted data to make Y exactly spatial and X exactly in the dispersion axis
      POLYWARP, Xo, yo, xi, yi, 2, P_inv, Q_inv ;invert the unwarp, takes, aligned image and make it tilted like the original
 
      ;Perform an image unwarping based on P and Q:
      dummy = dummy / blue_flat & err_dummy = err_dummy / blue_flat
      unwarped = POLY_2D(dummy, P_fwd, Q_fwd, 1) & err_unwarped = POLY_2D(err_dummy, P_fwd, Q_fwd, 1)

      x = findgen(1078)
      y = total(unwarped, 1, /Nan) ;Collapse the wavelength dimension
      err_y = sqrt(total(err_unwarped^2., 1, /Nan))
      
      trim = 15      
      initial_guess = [1.e3 * max(y), -1., 0.]     
      upper_mask_edge = where(y eq max(y[ upper_section_lower_row_b : upper_section_upper_row_b - trim ]), /NULL)
      fit_x = x[upper_mask_edge : upper_section_upper_row_b - trim] 
      fit_y = y[upper_mask_edge : upper_section_upper_row_b - trim]
      err_fit_y = err_y[upper_mask_edge : upper_section_upper_row_b - trim]
      p = mpfitfun('Scattered_Light', fit_x - mean(center), fit_y, err_fit_y, initial_guess, /NAN, /quiet) 
      fit_y_upper = Scattered_Light(fit_x - mean(center), p)      
      lower_mask_edge = where(y eq max(y[ lower_section_lower_row_b  + trim : lower_section_upper_row_b ]), /NULL)
      fit_x = x[lower_section_lower_row_b + trim : lower_mask_edge] 
      fit_y = y[lower_section_lower_row_b + trim : lower_mask_edge]
      err_fit_y = err_y[lower_section_lower_row_b + trim : lower_mask_edge]
      p = mpfitfun('Scattered_Light', reverse(mean(center) - fit_x), reverse(fit_y), reverse(err_fit_y), initial_guess, /NAN, /quiet) 
      fit_y_lower = Scattered_Light(reverse(mean(center) - fit_x), p)

      ;correct the flat field
      model = y
      model[lower_section_lower_row_b + trim : lower_mask_edge] = reverse(fit_y_lower) & model[upper_mask_edge : upper_section_upper_row_b - trim] = fit_y_upper  
      flat_correction = rebin(transpose(y / model), 2098, 1078) ;ratio of the flat corrected data to a fit using an ideal power law fall off for Jovian scatteblue light.
      flat_correction = POLY_2D(flat_correction, P_inv, Q_inv, 1) ;warp it back to the tilted images
      unwarped_post_correction = POLY_2D(image / (blue_flat * flat_correction), P_fwd, Q_fwd, 1) ;for inspection only, unwarp it
      y_corrected = total(unwarped_post_correction, 1) ;for inspection
      flat_correction_array[i, *, *] = flat_correction ; log this correction and average all the corrections ---> all images will still use the same flat. 
      
      ;inspect
      if keyword_set(debug) then begin
        Unwarped[indgen(2098), mean(center)] = 0 
        window, 0, xs = 2098, ys = 1000
        tv, bytscl(unwarped, 0., 1000) ;inspect the unward image, be sure the black line center of Jupiter looks centered 
        Unwarped[indgen(2098), mean(center)] = 1 
        window, 1, xs = 800, ys = 600
        cgplot, x, y, psym = 10, /ynozero, xr = [lower_section_lower_row_b, upper_section_upper_row_b]
        cgplot, x, model, /overplot, color = 'red'
        cgplot, x, y_corrected, /overplot, color = 'green'
        ;if i eq 3 then stop
      endif
  endfor
  blue_flat_correction = median(flat_correction_array, dimension = 1, /even)
  window, 0, xs = 2098, ys = 1000
  tv, bytscl(blue_flat_correction, 0.9, 1.1) ;inspect the flat field correction
  
    ;find the RED polynomical fit to Jupiter's center.
    Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '{Torus,Jupiter}*r.BSCR.fits')) 
    n_images = n_elements([sky, torus])  ;just use sky and torus frames to get this correction
    Images_r = Images_r[[sky, torus]] ;just use sky and torus frames to get this correction
    flat_correction_array = fltarr(n_images, 2098, 1078)
    for i= 0, n_images - 1 do begin
      image = MRDFITS(Images_r[i], 0, header, /Dscale, /silent )
      err_image = MRDFITS(Images_r[i], 1, /Dscale, /silent )
      dummy = image & err_dummy = err_image
      
      ;get a signal to noise estimate under the mask for every spectral bin
      under_mask = [lower_section_upper_row_r, upper_section_lower_row_r] ; be sure your only searching regions well inside the aperture    
      y = total(image[*, under_mask[0]:under_mask[1]], 2, /NAN)
      dy = sqrt(total(err_image[*, under_mask[0]:under_mask[1]]^2., 2, /NAN))
      x = findgen(2098)
  
      ;Determine rebinning needed to have a uniform signal to noise, and get the pixel fits of Jupiter's centroid
      new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
      Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
      for j = 0, N_elements(new_y) - 1 do begin 
        x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
        compress = where(BIN_map eq j, /null)
        y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
        junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
        fit = jjgaussfit(x, y, a, LIMITamp = [1.e3, 2.^16.], LIMITcen = [min(x)+10., max(x)-10.], FIXsl = 0., guess = [max(y), peak, (2./R_J_per_pixel_r)/3., 0., 0.])
        if (total(finite(a)) eq 5) then begin
          disk = where(y gt .1*a[0], /Null)               
          if (disk eq !NULL) then continue 
          disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel_r) and (disk lt 1.2/R_J_per_pixel_r + a[1])), /NULL)]
          if (abs(median(x[disk]) - a[1]) gt 20.) then continue ;reject of the Gaussian centroid and the pixels at the center of a region 10% the peak differ by >20 pixels
          Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]       
        endif
      endfor 
      keep = where(finite(Jupiter_centroid), /NULL) 
     
      ;1D fits
      coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
      center = poly(findgen(2098), coeffs) ; pixel location of Jupiter's center
      if keyword_set(debug) then begin
        window, 2, title = 'Fitting Jupiter''s Centroid'
        plot, new_x[keep], Jupiter_centroid[keep], psym = 4, /ynozero
        oplot, findgen(2098), center
      endif
      
      ;2D unwarping: dispersion axis isn't exactly aligned in x, correct it so that it is.
      Xo = [new_x[keep], new_x[keep], new_x[keep]]
      Yo = [make_array(n_elements(keep), value = mean(center) - 300, /float), $
            make_array(n_elements(keep), value = mean(center), /float), $
            make_array(n_elements(keep), value = mean(center) + 300, /float)] 
      Xi = [new_x[keep], new_x[keep], new_x[keep]]
      Yi = [Jupiter_centroid[keep] - 300., Jupiter_centroid[keep], Jupiter_centroid[keep] + 300.]
      POLYWARP, xi, yi, Xo, yo, 2, P_fwd, Q_fwd ; P and Q are the coefficient matrices to warp the tilted data to make Y exactly spatial and X exactly in the dispersion axis
      POLYWARP, Xo, yo, xi, yi, 2, P_inv, Q_inv ;invert the unwarp, takes, aligned image and make it tilted like the original
 
      ;Perform an image unwarping based on P and Q:
      dummy = dummy / red_flat & err_dummy = err_dummy / red_flat
      unwarped = POLY_2D(dummy, P_fwd, Q_fwd, 1) & err_unwarped = POLY_2D(err_dummy, P_fwd, Q_fwd, 1)
      Unwarped[indgen(2098), mean(center)] = 0 
      window, 0, xs = 2098, ys = 1000
      tv, bytscl(unwarped, 0., 5000) ;inspect the center of Jupiter looks centered
      Unwarped[indgen(2098), mean(center)] = 1 

      x = findgen(1078)
      y = total(unwarped, 1, /Nan) ;Collapse the wavelength dimension
      err_y = sqrt(total(err_unwarped^2., 1, /Nan))
      
      trim = 20      
      initial_guess = [1.e3 * max(y), -1., 0.]     
      upper_mask_edge = where(y eq max(y[upper_section_lower_row_r : upper_section_upper_row_r ]), /NULL)
      fit_x = x[upper_mask_edge : upper_section_upper_row_r - trim] 
      fit_y = y[upper_mask_edge : upper_section_upper_row_r - trim]
      err_fit_y = err_y[upper_mask_edge : upper_section_upper_row_r - trim]
      p = mpfitfun('Scattered_Light', fit_x - mean(center), fit_y, err_fit_y, initial_guess, /NAN, /quiet) 
      fit_y_upper = Scattered_Light(fit_x - mean(center), p)      
      lower_mask_edge = where(y eq max(y[lower_section_lower_row_r : lower_section_upper_row_r ]), /NULL)
      fit_x = x[lower_section_lower_row_r + trim : lower_mask_edge] 
      fit_y = y[lower_section_lower_row_r + trim : lower_mask_edge]
      err_fit_y = err_y[lower_section_lower_row_r + trim : lower_mask_edge]
      p = mpfitfun('Scattered_Light', reverse(mean(center) - fit_x), reverse(fit_y), reverse(err_fit_y), initial_guess, /NAN, /quiet) 
      fit_y_lower = Scattered_Light(reverse(mean(center) - fit_x), p)

      ;correct the flat field
      model = y
      model[lower_section_lower_row_r + trim : lower_mask_edge] = reverse(fit_y_lower) & model[upper_mask_edge : upper_section_upper_row_r - trim] = fit_y_upper  
      flat_correction = rebin(transpose(y / model), 2098, 1078) ;ratio of the flat corrected data to a fit using an ideal power law fall off for Jovian scattered light.
      flat_correction = POLY_2D(flat_correction, P_inv, Q_inv, 1) ;warp it back to the tilted images
      unwarped = POLY_2D(image / (red_flat * flat_correction), P_fwd, Q_fwd, 1) ;for inspection only, unwarp it
      y_corrected = total(unwarped, 1) ;for inspection
      flat_correction_array[i, *, *] = flat_correction ; log this correction and average all the corrections ---> all images will still use the same flat. 
      
      ;inspect
      if keyword_set(debug) then begin
        window, 1, xs = 800, ys = 600
        cgplot, x, y, psym = 10, /ynozero, xr = [lower_section_lower_row_r, upper_section_upper_row_r]
        cgplot, x, model, /overplot, color = 'red'
        cgplot, x, y_corrected, /overplot, color = 'green'
      endif
  endfor
  red_flat_correction = median(flat_correction_array, dimension = 1, /even)
  window, 0, xs = 2098, ys = 1000
  tv, bytscl(red_flat_correction, 0.9, 1.1) ;inspect the flat field correction
  
  ;correct the flat fields
  blue_flat = blue_flat * blue_flat_correction 
  red_flat = red_flat * red_flat_correction 
  ;--------------------------------------------------------------------------------------------------------------------------------------------

  ;Pull the timestamp of the torus frames. Use the sky frame closest in time
    Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '{Torus,Jupiter}*b.BSCR.fits'), count = n_images) 
    Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '{Torus,Jupiter}*r.BSCR.fits'), count = n_images) 
    torus_times = dblarr(N_elements(Torus))
    sky_times = dblarr(N_elements(Sky)) 
    for i = 0, N_elements(Torus) - 1 do begin
          test = MRDFITS(Images_r[Torus[i]], 0, header, /Dscale, /silent )
          cspice_utc2ET, sxpar(header, 'DATE-OBS'), ET
          torus_times[i] = ET + float(sxpar( header, 'EXPTIME'))/2. ;seconds past J2000
    endfor
    for i = 0, N_elements(sky) - 1 do begin
          test = MRDFITS(Images_r[Sky[i]], 0, header, /Dscale, /silent ) 
          cspice_utc2ET, sxpar(header, 'DATE-OBS'), ET
          sky_times[i] = ET + float(sxpar( header, 'EXPTIME'))/2. ;seconds past J2000   
    endfor
    if keyword_set(debug) then print, torus_times, sky_times

  ;SKY SUBTRACT and FLAT FIELD using the index of the appropriate sky frame 
  telrot_array_b = [] & final_shift_array_b = []
  telrot_array_r = [] & final_shift_array_r = []
  for i = 0, N_elements(Torus) - 1 do begin
        Time_difference = min(abs(torus_times[i] - sky_times), sky_location) ;seconds between sky and torus frames 
        sky_sub, Images_b[Torus[i]], Images_b[sky[sky_location]], 'blue', blue_flat, err_blue_flat, debug = debug ;write a fits file of the sky subtracted torus frame with an SS suffix              
        sky_sub, Images_r[Torus[i]], Images_r[sky[sky_location]], 'red' ,  red_flat,  err_red_flat, debug = debug  ;write a fits file of the sky subtracted torus frame with an SS suffix 
        print, Images_r[Torus[i]],'   ', Images_r[sky[sky_location]]
        Sky_Subtracted_Filename_r = STRMID(Images_r[Torus[i]], 0, strpos(Images_r[Torus[i]], '.fits' )) + 'SS.fits' 
        print, Images_b[Torus[i]],'   ', Images_b[sky[sky_location]]
        Sky_Subtracted_Filename_b = STRMID(Images_b[Torus[i]], 0, strpos(Images_b[Torus[i]], '.fits' )) + 'SS.fits' 
        Ssub_torus_image_r = MRDFITS(Sky_Subtracted_Filename_r, 0, header, /Dscale, /silent) / red_flat    ;FLAT FIELD
        Ssub_torus_image_b = MRDFITS(Sky_Subtracted_Filename_b, 0, header, /Dscale, /silent) / blue_flat   ;FLAT FIELD
        ;window, 0, xs = 2098, ys = 1000
        ;tv, bytscl(Ssub_torus_image_r, -10, 100)
        ;window, 1, xs = 2098, ys = 1000        
        ;tv, bytscl(Ssub_torus_image_b, -10, 40)
  endfor  
  ;window, 0
  ;plot, telrot_array_b, final_shift_array_b, psym = 5
  ;window, 1
  ;plot, telrot_array_r, final_shift_array_r, psym = 5
;stop
  ;1) Normalize into counts per second for variable exposure times
  ;2) Put Jupiter back into the sky subtracted frame
  for i = 0, N_elements(Torus) - 1 do begin                 
        orig_torus_image_r = MRDFITS(Images_r[Torus[i]], 0, torus_header, /Dscale, /silent)
        err_orig_torus_image_r = MRDFITS(Images_r[Torus[i]], 1, err_torus_header, /Dscale, /silent)
        Sky_Subtracted_Filename_r = STRMID(Images_r[Torus[i]], 0, strpos(Images_r[Torus[i]], '.fits' )) + 'SS.fits' 
        Ssub_torus_image_r = MRDFITS(Sky_Subtracted_Filename_r, 0, /Dscale, /silent)    
        err_Ssub_torus_image_r = MRDFITS(Sky_Subtracted_Filename_r, 1, /Dscale, /silent)
        Ssub_torus_image_r[*, mask_coords_r[0]:mask_coords_r[1]] = orig_torus_image_r[*, mask_coords_r[0]:mask_coords_r[1]] / 100. ;replace mask subtraction with with 100th the original Jupiter data
        err_Ssub_torus_image_r[*, mask_coords_r[0]:mask_coords_r[1]] = err_orig_torus_image_r[*, mask_coords_r[0]:mask_coords_r[1]] / 100. ;replace mask subtraction with with 100th the original Jupiter data
        Ssub_torus_image_r = Ssub_torus_image_r  / float(sxpar(torus_header, 'EXPTIME')) & err_Ssub_torus_image_r = err_Ssub_torus_image_r  / float(sxpar(torus_header, 'EXPTIME'))
        sxaddpar, torus_header, 'HISTORY', 'Original Jupiter Mask Data / 100'  
        Jupiter_replaced_Filename = STRMID(Images_r[Torus[i]], 0, strpos(Images_r[Torus[i]], '.fits' )) + 'JR.fits' 
        if ((directory eq 'D:\DATA\Apache Point Data\UT170310\') and (i le 8)) then begin
          Ssub_torus_image_r = shift(Ssub_torus_image_r, -34, 0); hack: these were shifted WRT arcs when switching from DIS imaging mode      
          err_Ssub_torus_image_r = shift(err_Ssub_torus_image_r, -34, 0); hack: these were shifted WRT arcs when switching from DIS imaging mode      
        endif
        MWRFITS, Ssub_torus_image_r, Jupiter_replaced_Filename, torus_header, /create, /silent
        MWRFITS, err_Ssub_torus_image_r, Jupiter_replaced_Filename, /silent ;Append the fits file with the error -> extension 1
           
        orig_torus_image_b = MRDFITS(Images_b[Torus[i]], 0, torus_header, /Dscale, /silent)
        err_orig_torus_image_b = MRDFITS(Images_b[Torus[i]], 1, err_torus_header, /Dscale, /silent)
        Sky_Subtracted_Filename_b = STRMID(Images_b[Torus[i]], 0, strpos(Images_b[Torus[i]], '.fits' )) + 'SS.fits' 
        Ssub_torus_image_b = MRDFITS(Sky_Subtracted_Filename_b, 0, /Dscale, /silent)    
        err_Ssub_torus_image_b = MRDFITS(Sky_Subtracted_Filename_b, 1, /Dscale, /silent)
        Ssub_torus_image_b[*, mask_coords_b[0]:mask_coords_b[1]] = orig_torus_image_b[*, mask_coords_b[0]:mask_coords_b[1]] / 100. ;replace mask subtraction with with 100th the original Jupiter data
        err_Ssub_torus_image_b[*, mask_coords_b[0]:mask_coords_b[1]] = err_orig_torus_image_b[*, mask_coords_b[0]:mask_coords_b[1]] / 100. ;replace mask subtraction with with 100th the original Jupiter data
        Ssub_torus_image_b = Ssub_torus_image_b  / float(sxpar(torus_header, 'EXPTIME')) & err_Ssub_torus_image_b = err_Ssub_torus_image_b  / float(sxpar(torus_header, 'EXPTIME'))
        sxaddpar, torus_header, 'HISTORY', 'Original Jupiter Mask Data / 100'  
        Jupiter_replaced_Filename = STRMID(Images_b[Torus[i]], 0, strpos(Images_b[Torus[i]], '.fits' )) + 'JR.fits' 
        if ((directory eq 'D:\DATA\Apache Point Data\UT170310\') and (i le 8)) then begin
          Ssub_torus_image_b = shift(Ssub_torus_image_b, -34, 0); hack: these were shifted WRT arcs when switching from DIS imaging mode      
          err_Ssub_torus_image_b = shift(err_Ssub_torus_image_b, -34, 0); hack: these were shifted WRT arcs when switching from DIS imaging mode      
        endif
        MWRFITS, Ssub_torus_image_b, Jupiter_replaced_Filename, torus_header, /create, /silent
        MWRFITS, err_Ssub_torus_image_b, Jupiter_replaced_Filename, /silent ;Append the fits file with the error -> extension 1
  endfor
endif  

if part eq 3 then begin ;Rectification & Wavelength Calibration (straighten the spectral line curvature). 
  
  ;find the mask edges with a flat field
      mask_coords_b = intarr(2) & mask_coords_r = intarr(2)
      blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 0, torus_header, /Dscale, /silent )
      red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 0, torus_header, /Dscale, /silent )
      profile_r = total(red_flat, 1)
      profile_b = total(blue_flat, 1)
      mask_r = where(profile_r lt mean(profile_r))
      mask_r = mask_r[where((mask_r gt 350) and (mask_r lt 800))]
      mask_coords_r = [min(mask_r), max(mask_r)]
      mask_b = where(profile_b lt mean(profile_b))
      mask_b = mask_b[where((mask_b gt 200) and (mask_b lt 700))]
      mask_coords_b = [min(mask_b), max(mask_b)]
  
  ;find lines that should be in the arc spectrum:
      READCOL,'D:\DATA\Apache Point Data\HeNeAr_Linelist.dat', F='F,A', Wavelength, Type, STRINGSKIP = '#'
      Additional_Lines = [3819.61, 4120.82, 4200.68, 4387.93] ;lines actually verified measurable, but not in HeNeAr_Linelist.dat.
      Wavelength = [Wavelength, Additional_Lines]
      Wavelength = Wavelength[sort(Wavelength)]
      never_measured = [Wavelength[0:7], 3850.57, 4072.20, 4131.73, 4277.55] ;lines that are not actually verified measurable, but present in the APO supplied HeNeAr_Linelist.dat line list
      for i = 0, N_elements(never_measured)-1 do remove, where(wavelength eq never_measured[i], /NULL), wavelength

  ;=====================================blue channel=======================================================  
    blue_arcs_long = FILE_SEARCH(strcompress(Directory + 'Arcs_Long.*b.fits'), count = n_images_long)  
    if n_images_long gt 0 then blue_arcs = MRDFITS(directory+'Processed\Master_blue_arcs_long.fits', 0, header, /Dscale, /silent ) else $
    blue_arcs = MRDFITS(directory+'Processed\Master_blue_arcs.fits', 0, header, /Dscale, /silent )
    if keyword_set(add_blue_arcs) then blue_arcs = blue_arcs + MRDFITS(add_blue_arcs, 0, header, /Dscale, /silent ) ;for when arc frames are too faint, ALWAYS DOUBLE CHECK ALIGNMENT
    aperture_b = where(blue_flat gt mean(blue_flat))
      dummy = blue_flat       
      y_aperture_ind = where( total(blue_flat,1) gt mean(total(blue_flat,1)), /NULL )
      ind = array_indices(blue_flat, aperture_b) 
      dummy[0:max(ind[0,*]),y_aperture_ind] = -100
      aperture_b = where(dummy eq -100, /NULL, complement = anti_aperture_blue) ;a much better aperture.
    center_wl = sxpar(header, 'DISPWC') ;fits header default
    dispersion = sxpar(header, 'DISPDW') ;fits header default (inaccurate)
    search = 20 ;pixel distance from the expected position based of the fits header to search for an arc line 
    dispersion = dispersion_b ;override -- manual inputs at program header     
    center_wl_pixel = center_wl_pixel_b ;override -- manual inputs at program header 
    dummy = blue_arcs ;inspection
    free_spectral_range = [center_wl - dispersion*2098./2., center_wl + dispersion*2098./2.] ;wavelength coverage in Angstroms
    lines_in_range = wavelength[where((wavelength gt free_spectral_range[0]) and (wavelength lt free_spectral_range[1]))]
    expected_pixel = center_wl_pixel + ((lines_in_range - center_wl) / dispersion)       
    keep = where((expected_pixel gt 21.) and (expected_pixel lt 2077.), /NULL, count) ;crop to lines that are fully on chip
    expected_pixel = expected_pixel[keep]
    lines_in_range = lines_in_range[keep]
    aperture_coords = ARRAY_INDICES(blue_arcs, aperture_b)  
    lower_section_lower_row = min(aperture_coords[1,*])
    lower_section_upper_row = min(mask_b) - 1
    upper_section_lower_row = max(mask_b) + 1 
    upper_section_upper_row = max(aperture_coords[1,*])   

    ;----------------------------blue rectification------------------------------------------------
        s = size(blue_arcs[*, lower_section_lower_row:upper_section_upper_row], /dimensions)
        both = blue_arcs[*, lower_section_lower_row:upper_section_upper_row]
        dummy = both
        blank_dummy = fltarr(s[0], s[1])
    ;--------------------find all the lines actually measured--------------------------------------------------------------------   
              identwave = [] & identpixl = []
              collapsed_arcs = total(blue_arcs[*, upper_section_lower_row:upper_section_upper_row], 2, /Nan)
              ;collapsed_arcs = total(blue_arcs[*, lower_section_lower_row:lower_section_upper_row], 2, /Nan)
              collapsed_arcs_WL = [center_wl - reverse(dispersion*findgen(2098./2.)), center_wl + dispersion*findgen(2098./2.)]
              
              ;inspection for initial disperions and center pixel hacks if the fits header is too far off---> adjust manual inputs at program header. 
                window, 0
                plot, collapsed_arcs, /ylog, yr = [50,2.e6], ystyle = 1
                oplot, expected_pixel, make_array(N_elements(expected_pixel), /FLOAT, value = 1.e3), psym = 4

                parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.,0.]}, 4)
                parinfo[2].fixed = 1  
                parinfo[2].value = 1.71  
              window, 1, xs = 600, ys = 400
              for i = 0, n_elements(lines_in_range) - 1 do begin
                result = mpfitpeak(findgen(search*2. + 1), float(collapsed_arcs[expected_pixel[i]-search:expected_pixel[i]+search]), a, $
                                   /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE =2) 
                plot, findgen(search*2. + 1), float(collapsed_arcs[expected_pixel[i]-search:expected_pixel[i]+search])                  
                if not finite(total(a)) then continue
                if ((status gt 0) and (a[0] gt 7.*robust_sigma(collapsed_arcs)+a[3])) then begin ;NEED A MORE ROBUST THRESHOLD
                  oplot, findgen(search*2. + 1), result
                  identwave = [identwave,lines_in_range[i]]
                  identpixl = [identpixl,a[1] - search + expected_pixel[i]]
                endif
              endfor

              ;Inspect predicted vs found lines
                window, 2
                cgplot, identpixl, identwave, psym = 4, yrange = [free_spectral_range[0],free_spectral_range[1]], color = 'blue', $
                  xtitle = 'Pixel #', ytitle = 'Angstroms', Title = 'BLUE = Measured Peaks / RED = Fits Header & Line List', ystyle = 1, xrange = [0,2098], xstyle = 1  
                cgplot, center_wl_pixel + ((lines_in_range - center_wl) / dispersion), lines_in_range, psym = 5, color = 'red', /overplot  

              ;fit dispersion of found lines, throw out high residuals (wrong lines)
                coeff = ROBUST_POLY_FIT(identpixl, identwave, 1, yfit, SIG) 
                cgplot, identpixl, yfit, psym = 5, color = 'green', /overplot  
                residual = abs(yfit - identwave)
                keep = where(residual lt 4.*SIG, /NULL) 

              lines_in_range = identwave[keep] ;replace values with results of the lines actually found.
              expected_pixel = identpixl[keep] ;replace values with results of the lines actually found.

        ;--------------------Trace their measured line shapes--------------------------------------------------------------------   
        print, 'Finding/fitting He Ne Ar emission line locations in the blue arc frame (slow)...'
        line_pixel_positions = fltarr(N_elements(lines_in_range), upper_section_upper_row - lower_section_lower_row + 1)  
        key_search_range = 10. ;range to seach for a line maximum over in pixels
        center = round(N_elements(both[0,*])/2.)
        for i = 0, N_elements(lines_in_range) - 1 do begin 
            
            ;Get an accurate starting position for tracing the line from the lower section
            Fit = mpfitpeak(findgen(2.*key_search_range+1.),both[expected_pixel[i] - key_search_range:expected_pixel[i] + key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row], a)
            location = a[1] ;use the fit peak instead of the max peak
            key = [location + expected_pixel[i] - key_search_range, mean([lower_section_lower_row,lower_section_upper_row]) - lower_section_lower_row]

            ;For each line, get the pixel locations along the cross-dispersion axis
            p3d_wavecal_calculate_mcurv, both, key[0], Key[1], linepos_low, linewidth = 10

            ;if the fit peaks / maxima are far from a measured peak then set them to !nan
            for j = 0, n_elements(linepos_low) - 1 do begin
              check = round(linepos_low[j]) - expected_pixel[i]
              if min(abs(check)) gt 10. then linepos_low[j]= !values.f_nan ;reject the calibration line if there's no measured line within 10 pixels       
            endfor 
            
            ;Get an accurate starting position for tracing the line from the upper section
            Fit = mpfitpeak(findgen(2.*key_search_range+1.),both[expected_pixel[i] - key_search_range:expected_pixel[i] + key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row], a)
            location = a[1] ;use the fit peak instead of the max peak
            key = [location + expected_pixel[i] - key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row]   
            
            ;For each line, get the pixel locations along the cross-dispersion axis
            p3d_wavecal_calculate_mcurv, both, key[0], Key[1], linepos_high, linewidth = 10
            
            ;if the fit peaks / maxima are far from a measured peak then set them to !nan
            for j = 0, n_elements(linepos_high) - 1 do begin
              check = round(linepos_high[j]) - expected_pixel[i]
              if min(abs(check)) gt 10. then linepos_high[j]= !values.f_nan ;reject the calibration line if there's no measured line within 20 pixels       
            endfor 

            line_pixel_positions[i, *] = [linepos_low[0:center], linepos_high[center+1:*]]
        endfor  

        ;correct these line positions using Gaussian fitting at every row.
          refrow = round(mean([lower_section_lower_row,lower_section_upper_row]) - lower_section_lower_row)
          i = 1
          REPEAT BEGIN ;grow the reference row until it's finite and gets all the lines, that might be a lot though!
             lowerpos = line_pixel_positions[*, refrow-i:refrow+i]
             i++
          ENDREP UNTIL finite(total(median(lowerpos, dimension=2)))
          p3d_wavecal_correct_maskpos, both, median(lowerpos, dimension=2), 10, lower_out,  inlines, FWHM = 2.7, method = 'Gaussian', refrow = refrow  
          refrow = round(mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row)
          i = 1
          REPEAT BEGIN ;grow the reference row until it's finite and gets all the lines, that might be a lot though!
             upperpos = line_pixel_positions[*, refrow-i:refrow+i]
             i++
          ENDREP UNTIL finite(total(median(upperpos, dimension=2)))
          p3d_wavecal_correct_maskpos, both, median(upperpos, dimension=2), 10, Upper_out, inlines, FWHM = 2.7, method = 'Gaussian', refrow = refrow  
          line_pixel_positions = [[lower_out[*, 0:center]], [upper_out[*, center+1:*]]]
          
        ;for each arc lamp emission line, fit the line positions in y with a simple polynomial, replace with the fit line pixel position.
        for i = 0, N_elements(lines_in_range) - 1 do begin
          x = findgen(n_elements(line_pixel_positions[i, *]))
          y = line_pixel_positions[i, *]
          reject = where((x gt MASK_COORDS_B[0] - lower_section_lower_row) and (x lt MASK_COORDS_B[1] - lower_section_lower_row), /null, COMPLEMENT = keep) ;reject points under the mask
          x[reject] = !values.F_nan & y[reject] = !values.F_nan
          reject = where(abs(y - expected_pixel[i]) gt 10., /null, COMPLEMENT = keep) ;reject points more than 10 pixels from home
          x[reject] = !values.F_nan & y[reject] = !values.F_nan 
          coeff = ROBUST_POLY_FIT(x[where(finite(y) eq 1, /null)], y[where(finite(y) eq 1, /null)], 2, yfit, SIG)            
          x = findgen(n_elements(line_pixel_positions[i, *]))
          yfit = POLY(X, Coeff)
          if sig lt 1 then line_pixel_positions[i, *] = yfit else line_pixel_positions[i, *] = !values.F_NaN  ;reject fits with high residuals
          blank_dummy[y, indgen(upper_section_upper_row - lower_section_lower_row + 1)] = 10000. ;mark the territory
          dummy[round(line_pixel_positions[i, *]), indgen(upper_section_upper_row - lower_section_lower_row + 1)] = 10000. ;mark the territory         
        endfor
        
        ;truncate to only arc lines with a decent polynomial fit
        keep = where(finite(total(line_pixel_positions, 2)), /NULL)
        lines_in_range = lines_in_range[keep]
        line_pixel_positions = line_pixel_positions[keep,*]
           
        window, 4, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
        tv, bytscl(dummy, 0, 10)
        window, 5, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
        tv, bytscl(blank_dummy, 0, 10)    

        Print, 'Rectifying blue channel using', N_elements(lines_in_range), ' lines... '
        p3d_wavecal_fit_maskwd, both, N_elements(lines_in_range), line_pixel_positions, lines_in_range, 1, 3, -1, warp_coefficients, debug = debug, chisq = chisq, optlow = 1       
        template_both_b = dblarr(2098, N_elements(warp_coefficients.m[0,*]))
        for k = 0L, N_elements(warp_coefficients.m[0,*]) - 1L do begin
          xinv = (2d0*dindgen(2098) - total(warp_coefficients.m[*,k])) / (warp_coefficients.m[1L,k] - warp_coefficients.m[0L,k])
          Template_both_b[*,k] = poly(xinv, warp_coefficients.c[k,*])
        endfor

        ;Rectify & make a dispersion mask
        P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_b, crval, cdelt, both_out, dout
        dispersion_mask_b = both
        for k = 0L, n_elements(both_out[0,*])-1 do begin
          dispersion_mask_b[0:N_elements(both_out[*,k])-1, k] = findgen(N_elements(both_out[*,k]))*cdelt + crval
        endfor
        full_frame_dispersion_mask = make_array(2098,1078,/float,value = !values.f_nan)
        full_frame_dispersion_mask[*, lower_section_lower_row:upper_section_upper_row] = dispersion_mask_b
        dispersion_mask_b = temporary(full_frame_dispersion_mask)   
      
        window, 0, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
        tv, bytscl(both, 0, 10)
        offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges          
        blue_arcs[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = both_out
        should_be = fltarr(n_elements(lines_in_range))
        for i = 0, n_elements(lines_in_range)-1 do begin
          junk = min(abs(dispersion_mask_b[*,500] - lines_in_range[i]), location)
          should_be[i] = location
        endfor
        blue_arcs[should_be,*] = 100
        window, 1, xs = 2098, ys = 800
        tv, bytscl(blue_arcs[*,100:999], 0, 10) 
        ;rdpix, dispersion_mask_b[*,100:999]
        
        ;Find the sigma line width, and it's first order slope with wavelength 
        s = size(both_out, /dimensions)
        tv, bytscl(both_out)
        wid = []
        for i = 0, n_elements(lines_in_range)-1 do begin
          junk = min(abs(dispersion_mask_b[*,500] - lines_in_range[i]), location)
          save_widths = []
          for j = 0, s[1]-1 do begin
            if (j gt lower_section_upper_row-lower_section_lower_row) and (j gt upper_section_lower_row-lower_section_lower_row) then continue ;skip the mask
            Fit = mpfitpeak(findgen(2.*key_search_range+1.), both_out[location - key_search_range:location + key_search_range, j], a)  
            save_widths = [save_widths, a[2]]
          endfor 
          wid = [wid, median( save_widths )]
        endfor
        Lwidth_coeffs_b = ROBUST_POLY_FIT( lines_in_range, wid, 1, yfit)
        if keyword_set(debug) then begin
          window, 2, xs = 600, ys = 400
          plot, lines_in_range, wid, /ynozero, psym=4, xtitle = 'Angstroms', ytitle = 'Linewidth (sigma / pixels)' 
          oplot, lines_in_range, yfit
        endif

        ;unwarp and wavelength calibrate the science frames:  
        Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*b.BSCRJR.fits'), count = n_images) 
        for i = 0, n_elements(Images_b)-1 do begin
          orig_torus_image_b = MRDFITS(Images_b[i], 0, torus_header, /Dscale, /silent)  
          err_orig_torus_image_b = MRDFITS(Images_b[i], 1, /Dscale, /silent)  
          both = orig_torus_image_b[*, lower_section_lower_row:upper_section_upper_row]
          err_both = err_orig_torus_image_b[*, lower_section_lower_row:upper_section_upper_row]
          P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_b, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both    
          offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
          orig_torus_image_b[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = both_out
          err_orig_torus_image_b[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = err_both_out
          window, 1, xs = 2098, ys = upper_section_upper_row -  lower_section_lower_row + 1
          tv, bytscl(orig_torus_image_b[*, lower_section_lower_row:upper_section_upper_row], 0, .05)
          ;rdpix, dispersion_mask_b[*, lower_section_lower_row:upper_section_upper_row]

          ;write the dispersion parameters to the fits header. 
          sxaddpar, torus_header, 'DISPERS', string(cdelt) , 'Calculated Dispersion in A / pix'
          sxaddpar, torus_header, 'WAVELE1', string(crval) , 'Wavelength 1st pix disp mask (A)'
          sxaddpar, torus_header, 'W0_COEF', string(Lwidth_coeffs_b[0]), 'Sigma line width coeff 0'
          sxaddpar, torus_header, 'W1_COEF', string(Lwidth_coeffs_b[1]), 'Sigma line width coeff 1'
          Filename = STRMID(Images_b[i], 0, strpos(Images_b[i], '.fits' )) + 'Rec.fits' 

          MWRFITS, orig_torus_image_b, Filename, torus_header, /CREATE, /silent ;/create overwrites
          MWRFITS, err_orig_torus_image_b, Filename, /silent   ;Append the fits file -> extension 1
          MWRFITS, dispersion_mask_b, Filename, /silent ;Append the fits file -> extension 2
        endfor
        save, Template_both_b, filename = strcompress(Directory + 'Processed/Template_both_b.sav')

  ;=====================================red channel=======================================================  
    red_arcs = MRDFITS(directory+'Processed\Master_red_arcs.fits', 0, header, /Dscale, /silent )
    red_flat = MRDFITS(directory+'Processed\red_Master_flat.fits', 0, flat_header, /Dscale, /silent )
    red_arcs = rotate(red_arcs, 5) ;flip everything so that wavelength increase with pixel number
    red_flat = rotate(red_flat, 5) ;flip everything so that wavelength increase with pixel number
    aperture_r = where(red_flat gt mean(red_flat))
    center_wl = sxpar(header, 'DISPWC')
    dispersion = sxpar(header, 'DISPDW')
    search = 20 ;pixel distance from the expected position based of the fits header to search for an arc line 
    center_wl_pixel = center_wl_pixel_r ;manually input the center wavelength pixel at program header, based on first window 0 results below
    dispersion = dispersion_r ;manually input the Dispersion at program header, based on first window 0 results below
    free_spectral_range = [center_wl - dispersion*2098./2., center_wl + dispersion*2098./2.] ;wavelength coverage in Angstroms
    lines_in_range = wavelength[where((wavelength gt free_spectral_range[0]) and (wavelength lt free_spectral_range[1]))]
    expected_pixel = center_wl_pixel + ((lines_in_range - center_wl) / dispersion)       
    keep = where(((expected_pixel gt search) and (expected_pixel lt 2047.-search)), /NULL)
    expected_pixel = expected_pixel[keep] & lines_in_range = lines_in_range[keep]
    aperture_coords = ARRAY_INDICES(red_arcs, aperture_r)  
    lower_section_lower_row = min(aperture_coords[1,*])
    lower_section_upper_row = min(mask_r) - 1
    upper_section_lower_row = max(mask_r) + 1 
    upper_section_upper_row = max(aperture_coords[1,*])   

    ;----------------------------red rectification------------------------------------------------
        s = size(red_arcs[*, lower_section_lower_row:upper_section_upper_row], /dimensions)
        both = red_arcs[*, lower_section_lower_row:upper_section_upper_row]
        dummy = both
        blank_dummy = fltarr(s[0], s[1])
    
    ;--------------------find all the lines actually measured--------------------------------------------------------------------   
              identwave = [] & identpixl = []
              collapsed_arcs = total(red_arcs[*, upper_section_lower_row:upper_section_upper_row], 2, /Nan)
              collapsed_arcs_WL = [center_wl - reverse(dispersion*findgen(2098./2.)), center_wl + dispersion*findgen(2098./2.)]

              ;inspection for initial disperions and center pixel hacks if the fits header is too far off---> adjust manual inputs at program header. 
                window, 0
                plot, collapsed_arcs, /ylog, yr = [100,1.e7]
                oplot, expected_pixel, make_array(N_elements(expected_pixel), /FLOAT, value = 5.e6), psym = 4
                
                parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.,0.]}, 4)
                parinfo[2].fixed = 1  
                parinfo[2].value = 1.71  
              ;window, 0, xs = 600, ys = 400
              for i = 0, n_elements(lines_in_range) - 1 do begin
                result = mpfitpeak(findgen(search*2. + 1), float(collapsed_arcs[expected_pixel[i]-search:expected_pixel[i]+search]), a, $
                                   /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE =2) 
                ;plot, findgen(search*2. + 1), float(collapsed_arcs[expected_pixel[i]-search:expected_pixel[i]+search])                  
                if not finite(total(a)) then continue
                if ((status gt 0) and (a[0] gt 7.*robust_sigma(collapsed_arcs)+a[3])) then begin ;NEED A MORE ROBUST THRESHOLD
                  oplot, findgen(search*2. + 1), result
                  identwave = [identwave,lines_in_range[i]]
                  identpixl = [identpixl,a[1] - search + expected_pixel[i]]
                endif
              endfor

              ;Inspect predicted vs found lines
                window, 2
                cgplot, identpixl, identwave, psym = 4, yrange = [free_spectral_range[0],free_spectral_range[1]], color = 'blue', $
                  xtitle = 'Pixel #', ytitle = 'Angstroms', Title = 'BLUE = Measured Peaks / RED = Fits Header & Line List', ystyle = 1, xrange = [0,2098], xstyle = 1  
                cgplot, center_wl_pixel + ((lines_in_range - center_wl) / dispersion), lines_in_range, psym = 5, color = 'red', /overplot  

              ;fit dispersion of found lines, throw out high residuals (wrong lines)
                coeff = ROBUST_POLY_FIT(identpixl, identwave, 1, yfit, SIG) 
                cgplot, identpixl, yfit, psym = 5, color = 'green', /overplot  
                residual = abs(yfit - identwave)
                keep = where(residual lt 4.*SIG, /NULL) 

              lines_in_range = identwave[keep] ;replace values with results of the lines actually found.
              expected_pixel = identpixl[keep] ;replace values with results of the lines actually found.

        ;--------------------Trace their measured line shapes--------------------------------------------------------------------   
        print, 'Finding/fitting He Ne Ar emission line locations in the red arc frame (slow)...'
        line_pixel_positions = fltarr(N_elements(lines_in_range), upper_section_upper_row - lower_section_lower_row + 1)  
        key_search_range = 10. ;range to seach for a line maximum over in pixels
        center = round(N_elements(both[0,*])/2.)
        for i = 0, N_elements(lines_in_range) - 1 do begin 
            
            ;Get an accurate starting position for tracing the line from the lower section
            Fit = mpfitpeak(findgen(2.*key_search_range+1.),both[expected_pixel[i] - key_search_range:expected_pixel[i] + key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row], a)
            location = a[1] ;use the fit peak instead of the max peak
            key = [location + expected_pixel[i] - key_search_range, mean([lower_section_lower_row,lower_section_upper_row]) - lower_section_lower_row]

            ;For each line, get the pixel locations along the cross-dispersion axis
            p3d_wavecal_calculate_mcurv, both, key[0], Key[1], linepos_low, linewidth = 10

            ;if the fit peaks / maxima are far from a measured peak then set them to !nan
            for j = 0, n_elements(linepos_low) - 1 do begin
              check = round(linepos_low[j]) - expected_pixel[i]
              if min(abs(check)) gt 10. then linepos_low[j]= !values.f_nan ;reject the calibration line if there's no measured line within 10 pixels       
            endfor 
            
            ;Get an accurate starting position for tracing the line from the upper section
            Fit = mpfitpeak(findgen(2.*key_search_range+1.),both[expected_pixel[i] - key_search_range:expected_pixel[i] + key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row], a)
            location = a[1] ;use the fit peak instead of the max peak
            key = [location + expected_pixel[i] - key_search_range, mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row]   
            
            ;For each line, get the pixel locations along the cross-dispersion axis
            p3d_wavecal_calculate_mcurv, both, key[0], Key[1], linepos_high, linewidth = 10
            
            ;if the fit peaks / maxima are far from a measured peak then set them to !nan
            for j = 0, n_elements(linepos_high) - 1 do begin
              check = round(linepos_high[j]) - expected_pixel[i]
              if min(abs(check)) gt 10. then linepos_high[j]= !values.f_nan ;reject the calibration line if there's no measured line within 20 pixels       
            endfor 

            line_pixel_positions[i, *] = [linepos_low[0:center],linepos_high[center+1:*]]
        endfor  

        ;correct these line positions using Gaussian fitting at every row.
          refrow = round(mean([lower_section_lower_row,lower_section_upper_row]) - lower_section_lower_row)
          i = 1
          REPEAT BEGIN ;grow the reference row until it's finite and gets all the lines, that might be a lot though!
             lowerpos = line_pixel_positions[*, refrow-i:refrow+i]
             i++
          ENDREP UNTIL finite(total(median(lowerpos, dimension=2)))
          p3d_wavecal_correct_maskpos, both, median(lowerpos, dimension=2),  10, lower_out, inlines, FWHM = 2.7, method = 'Gaussian', refrow = refrow  
          refrow = round(mean([upper_section_lower_row,upper_section_upper_row]) - lower_section_lower_row)
          i = 1
          REPEAT BEGIN ;grow the reference row until it's finite and gets all the lines, that might be a lot though!
             upperpos = line_pixel_positions[*, refrow-i:refrow+i]
             i++
          ENDREP UNTIL finite(total(median(upperpos, dimension=2)))
          p3d_wavecal_correct_maskpos, both, median(upperpos, dimension=2), 10, Upper_out, inlines, FWHM = 2.7, method = 'Gaussian', refrow = refrow  
          line_pixel_positions = [[lower_out[*, 0:center]], [upper_out[*, center+1:*]]]

        ;for each arc lamp emission line, fit the line positions in y with a simple polynomial, replace with the fit line pixel position.
        for i = 0, N_elements(lines_in_range) - 1 do begin
          x = findgen(n_elements(line_pixel_positions[i, *]))
          y = line_pixel_positions[i, *]
          reject = where((x gt MASK_COORDS_r[0] - lower_section_lower_row) and (x lt MASK_COORDS_r[1] - lower_section_lower_row), /null, COMPLEMENT = keep) ;reject points under the mask
          x[reject] = !values.F_nan & y[reject] = !values.F_nan
          reject = where(abs(y - expected_pixel[i]) gt 10., /null, COMPLEMENT = keep) ;reject points more than 10 pixels from home
          x[reject] = !values.F_nan & y[reject] = !values.F_nan 
          coeff = ROBUST_POLY_FIT(x[where(finite(y) eq 1, /null)], y[where(finite(y) eq 1, /null)], 2, yfit, SIG)   
          x = findgen(n_elements(line_pixel_positions[i, *]))
          yfit = POLY(X, Coeff)
          if sig lt 1 then line_pixel_positions[i, *] = yfit else line_pixel_positions[i, *] = !values.F_NaN  ;reject fits with high residuals
          blank_dummy[y, indgen(upper_section_upper_row - lower_section_lower_row + 1)] = 10000. ;mark the territory
          dummy[round(line_pixel_positions[i, *]), indgen(upper_section_upper_row - lower_section_lower_row + 1)] = 10000. ;mark the territory         
        endfor
        if keyword_set(debug) then begin
          window, 4, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
          tv, bytscl(dummy, 0, 100)
          window, 5, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
          tv, bytscl(blank_dummy, 0, 100)    
        endif
        
        ;truncate to only arc lines with a decent polynomial fit
        keep = where(finite(total(line_pixel_positions, 2)), /NULL)
        lines_in_range = lines_in_range[keep]
        line_pixel_positions = line_pixel_positions[keep,*]         
        
        Print, 'Rectifying red channel using', N_elements(lines_in_range), ' lines... '
        p3d_wavecal_fit_maskwd, both, N_elements(lines_in_range), line_pixel_positions, lines_in_range, 1, 3, -1, warp_coefficients, debug = debug, chisq = chisq, deadfibers = deadfibers, sdeadfibers = sdeadfibers    
        template_both_r = dblarr(2098, N_elements(warp_coefficients.m[0,*]))
        for k = 0L, N_elements(warp_coefficients.m[0,*]) - 1L do begin
          xinv = (2d0*dindgen(2098) - total(warp_coefficients.m[*,k])) / (warp_coefficients.m[1L,k] - warp_coefficients.m[0L,k])
          Template_both_r[*,k] = poly(xinv, warp_coefficients.c[k,*])
        endfor

        ;Rectify & make a dispersion mask
        P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_r, crval, cdelt, both_out, dout 
        dispersion_mask_r = both
        for k = 0L, n_elements(both_out[0,*])-1 do begin
          dispersion_mask_r[0:N_elements(both_out[*,k])-1, k] = findgen(N_elements(both_out[*,k]))*cdelt + crval
        endfor
        full_frame_dispersion_mask = make_array(2098, 1078, /float, value = !values.f_nan)
        full_frame_dispersion_mask[*, lower_section_lower_row:upper_section_upper_row] = dispersion_mask_r
        dispersion_mask_r = temporary(full_frame_dispersion_mask)   
        offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges          
        red_arcs[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = both_out
        
        if keyword_set(debug) then begin
          window, 0, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
          tv, bytscl(both, 0, 1000)
          window, 1, xs = 2098, ys = upper_section_upper_row - lower_section_lower_row + 1
          tv, bytscl(red_arcs[*,lower_section_lower_row:upper_section_upper_row], 0, 1000) 
        endif

        ;Find the sigma line width, and it's first order slope with wavelength 
        s = size(both_out, /dimensions)
        tv, bytscl(both_out)
        wid = []
        for i = 0, n_elements(lines_in_range)-1 do begin
          junk = min(abs(dispersion_mask_r[*,500] - lines_in_range[i]), location)
          save_widths = []
          for j = 0, s[1]-1 do begin
            if (j gt lower_section_upper_row-lower_section_lower_row) and (j gt upper_section_lower_row-lower_section_lower_row) then continue ;skip the mask
            Fit = mpfitpeak(findgen(2.*key_search_range+1.), both_out[location - key_search_range:location + key_search_range, j], a)  
            save_widths = [save_widths, a[2]]
          endfor 
          wid = [wid, median( save_widths )]
        endfor
        Lwidth_coeffs_r = ROBUST_POLY_FIT( lines_in_range, wid, 1, yfit)
        if keyword_set(debug) then begin
          window, 2, xs = 600, ys = 400
          plot, lines_in_range, wid, /ynozero, psym=4, xtitle = 'Angstroms', ytitle = 'Linewidth (sigma / pixels)' 
          oplot, lines_in_range, yfit
        endif

        ;unwarp and wavelength calibrate the science frames:  
        Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*r.BSCRJR.fits'), count = n_images) 
        for i = 0, n_elements(Images_r)-1 do begin
          orig_torus_image_r = MRDFITS(Images_r[i], 0, torus_header, /Dscale, /silent)  
          err_orig_torus_image_r = MRDFITS(Images_r[i], 1, /Dscale, /silent)  
          orig_torus_image_r = rotate(orig_torus_image_r, 5) ;flip everything so that wavelength increase with pixel number
          err_orig_torus_image_r = rotate(err_orig_torus_image_r, 5) 
          both = orig_torus_image_r[*, lower_section_lower_row:upper_section_upper_row]
          err_both = err_orig_torus_image_r[*, lower_section_lower_row:upper_section_upper_row]
          P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_r, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both    
          offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
          orig_torus_image_r[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = both_out 
          err_orig_torus_image_r[0: 2097 - offset[0], lower_section_lower_row:upper_section_upper_row] = err_both_out
          window, 1, xs = 2098, ys = 1098
          tv, bytscl(orig_torus_image_r, 0, .15)   
          ;rdpix, dispersion_mask_r
 
          ;write the dispersion parameters to the fits header. 
          sxaddpar, torus_header, 'DISPERS', string(cdelt) , 'Calculated Dispersion in A / pix'
          sxaddpar, torus_header, 'WAVELE1', string(crval) , 'Wavelength 1st pix disp mask (A)'
          sxaddpar, torus_header, 'W0_COEF', string(Lwidth_coeffs_r[0]), 'Sigma line width coeff 0'
          sxaddpar, torus_header, 'W1_COEF', string(Lwidth_coeffs_r[1]), 'Sigma line width coeff 1'
          Filename = STRMID(Images_r[i], 0, strpos(Images_r[i], '.fits' )) + 'Rec.fits' 
          MWRFITS, orig_torus_image_r, Filename, torus_header, /CREATE, /silent ;/create overwrites
          MWRFITS, err_orig_torus_image_r, Filename, /silent   ;Append the fits file -> extension 1
          MWRFITS, dispersion_mask_r, Filename, /silent ;Append the fits file -> extension 2
        endfor
        save, Template_both_r, filename = strcompress(Directory + 'Processed/Template_both_r.sav')    
endif ;part eq 3

if part eq 4 then begin ;Mask Neutral Density Calibration, Can skip if the same mask was used on the previous calibration
     Print, 'Data for computing mask throughput is NOT flat-fielded. Basically asssumes a physically perfect slit.
     Print, '. . . no way around this assumption unless flats are taken without the mask on.' 
     restore, strcompress(Directory + 'Processed\Template_both_b.sav')
     restore, strcompress(Directory + 'Processed\Template_both_r.sav')
   
     ;find the mask edges with a flat field
      mask_coords_b = intarr(2) & mask_coords_r = intarr(2)
      blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 0, flat_header_b, /Dscale, /silent )
      err_blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 1, /Dscale, /silent )
      red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 0, flat_header_r, /Dscale, /silent )
      err_red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 1, /Dscale, /silent )
      free_spectral_range_b = [float(sxpar(flat_header_b, 'DISPWC')) - float(sxpar(flat_header_b, 'DISPDW'))*2098./2., $
                               float(sxpar(flat_header_b, 'DISPWC')) + float(sxpar(flat_header_b, 'DISPDW'))*2098./2.] ;wavelength coverage in Angstroms
      free_spectral_range_r = [float(sxpar(flat_header_r, 'DISPWC')) - float(sxpar(flat_header_r, 'DISPDW'))*2098./2., $
                                float(sxpar(flat_header_r, 'DISPWC')) + float(sxpar(flat_header_r, 'DISPDW'))*2098./2.] ;wavelength coverage in Angstroms
      profile_r = total(red_flat, 1)
      profile_b = total(blue_flat, 1)
      mask_r = where(profile_r lt mean(profile_r))
      mask_r = mask_r[where((mask_r gt 350) and (mask_r lt 800))]
      mask_coords_r = [min(mask_r), max(mask_r)]
      mask_b = where(profile_b lt mean(profile_b))
      mask_b = mask_b[where((mask_b gt 200) and (mask_b lt 700))]
      mask_coords_b = [min(mask_b), max(mask_b)]
      aperture_b = where(blue_flat gt mean(blue_flat))
        dummy = blue_flat       
        y_aperture_ind = where( total(blue_flat,1) gt mean(total(blue_flat,1)), /NULL )
        ind = array_indices(blue_flat, aperture_b) 
        dummy[0:max(ind[0,*]),y_aperture_ind] = -100
        aperture_b = where(dummy eq -100, /NULL, complement = anti_aperture_blue) ;a much better aperture.
      aperture_r = where(red_flat gt mean(red_flat), /NULL, complement = anti_aperture_red)
      aperture_coords_b = ARRAY_INDICES(blue_flat, aperture_b)  
      lower_section_lower_row_b = min(aperture_coords_b[1,*]) 
      lower_section_upper_row_b = min(mask_b) - 1
      upper_section_lower_row_b = max(mask_b) + 1 
      upper_section_upper_row_b = max(aperture_coords_b[1,*]) 
      aperture_r = where(red_flat gt mean(red_flat))
      aperture_coords_r = ARRAY_INDICES(red_flat, aperture_r)  
      lower_section_lower_row_r = min(aperture_coords_r[1,*]) 
      lower_section_upper_row_r = min(mask_r) - 1
      upper_section_lower_row_r = max(mask_r) + 1 
      upper_section_upper_row_r = max(aperture_coords_r[1,*])   
      
      ;NORMALIZE THE FLAT FIELDS to a value of 1. (NO FLAT CORRECTION UNDER THE MASK)
        red_flat[ anti_aperture_red ] = !values.F_NaN
        mean_spec = mean(red_flat, dimension = 2, /NAN)
        red_flat[ anti_aperture_red ] = -1.
        p3d_tool_flatfield, red_flat, red_flat_norm, err_red_flat_norm, din = err_red_flat, /NORMALIZE, spec = mean_spec, /ndisp
        red_flat_norm[ anti_aperture_red ] = 1. & err_red_flat_norm[ anti_aperture_red ] = 0;err_red_flat[ anti_aperture_red ]
        red_flat = temporary(red_flat_norm) & err_red_flat = temporary(err_red_flat_norm)
        
        blue_flat[ anti_aperture_blue ] = !values.F_NaN
        mean_spec = mean(blue_flat, dimension = 2, /NAN)
        blue_flat[ anti_aperture_blue ] = -1.
        p3d_tool_flatfield, blue_flat, blue_flat_norm, err_blue_flat_norm, din = err_blue_flat, /NORMALIZE, spec = mean_spec, /ndisp
        blue_flat_norm[ anti_aperture_blue ] = 1. & err_blue_flat_norm[ anti_aperture_blue ] = 0;err_blue_flat[ anti_aperture_blue ]
        blue_flat = temporary(blue_flat_norm) & err_blue_flat = temporary(err_blue_flat_norm)

      ;LOAD A DISPERSION MASK IN EACH 
      Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*b.BSCRJRRec.fits'), count = n_images)           
      dispersion_mask_b = MRDFITS(Images_b[0], 2, /Dscale, /silent)
      Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*r.BSCRJRRec.fits'), count = n_images)           
      dispersion_mask_r = MRDFITS(Images_r[0], 2, /Dscale, /silent)

      ;Not all nights have an on-mask off mask calibration done, if calibration data doesn't exists, run this step for the closest night that does. It's same mask after all.
      if ((Jupiter_Off_Mask eq !NULL) and (Jupiter_On_Mask eq !NULL)) then begin
        current_directory = directory
        directory = Mask_Cal_Directory
        GOTO, load_directory
      endif

   ;==========================Blue Channel===========================================   

     ;find the center of the disk in an off mask exposure
     off_mask = MRDFITS(directory+'Processed\'+Jupiter_Off_Mask+'b.BSCR.fits', 0, off_mask_header, /Dscale, /silent)
     err_off_mask = MRDFITS(directory+'Processed\'+Jupiter_Off_Mask+'b.BSCR.fits', 1, /Dscale, /silent) 
     off_exptime = float(sxpar(off_mask_header, 'EXPTIME'))
     ;Flat field and rectify it to be consistent with the data that it will be used in correcting      
       off_mask_FF = off_mask / blue_flat 
       err_off_mask = abs(off_mask_FF) * sqrt( (err_off_mask/off_mask)^2. + (err_blue_flat/blue_flat)^2. )
       off_mask = temporary(off_mask_FF)
       both = off_mask[*, lower_section_lower_row_b : upper_section_upper_row_b]
       err_both = err_off_mask[*, lower_section_lower_row_b : upper_section_upper_row_b]
       P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_b, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both   
       offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
       off_mask[offset[0]:*, lower_section_lower_row_b : upper_section_upper_row_b] = both_out    
       err_off_mask[offset[0]:*, lower_section_lower_row_b : upper_section_upper_row_b] = err_both_out
     
     ;get Jupiter's instantaneous radius in pixels
      cspice_utc2ET, sxpar(off_mask_header, 'DATE-OBS'), ET
      torus_time = ET + float(sxpar(off_mask_header, 'EXPTIME'))/2. ;seconds past J2000       
      cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time, 'IAU_Jupiter', 'LT+S', 'Earth', Sub_Earth, trgepc, srfvec ;get sub observer point in catesian IAU Jupiter coords
      cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
      cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', ptarg, ltime
      R_j = 206264.806* atan(radii[0] / norm(ptarg))         
      R_J_per_pixel = float(sxpar(off_mask_header, 'PIXSCAL2')) / R_j 
       
  ;log the pixel center of jupiter's disc.
    dummy = off_mask
    under_mask = [lower_section_lower_row_b, upper_section_upper_row_b] ; be sure your only searching regions well inside the aperture    

    ;get a signal to noise estimate under the mask for every spectral bin
    y = total(off_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
    dy = total(err_off_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
    x = findgen(2098)

    ;Determine rebinning needed to have a uniform signal to noise
    new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
    Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
    for j = 0, N_elements(new_y) - 1 do begin 
      x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
      compress = where(BIN_map eq j, /null)
      y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
      junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
      fit = jjgaussfit(x, y, a)
      if (total(finite(a)) eq 5) then begin
        disk = where(y gt .1*a[0], /Null)               
        if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the ceter of a region 10% the peak differ by >20 pixels
        disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
        if (abs(median(x[disk]) - a[1]) gt 20.) then continue
        Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]
      endif
    endfor 
    keep = where(finite(Jupiter_centroid)) 
    coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
    off_center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's center
    if keyword_set(debug) then begin
      dummy = off_mask
      dummy[findgen(2098), round(off_center)] = 0
      window, 0, xs = 2098, ys = 1000
      tv, bytscl(dummy, 0., 500) ;inspect the center of Jupiter looks centered
      window, 1
      plot, new_x[keep], Jupiter_centroid[keep], psym=1;, yr = [590,610]
      oplot, off_center 
    endif

   ;find the center of the disk in an on mask exposure
   on_mask = MRDFITS(directory+'Processed\'+Jupiter_On_Mask+'b.BSCR.fits', 0, on_mask_header, /Dscale, /silent)       
   err_on_mask = MRDFITS(directory+'Processed\'+Jupiter_On_Mask+'b.BSCR.fits', 1, /Dscale, /silent) 

   on_exptime = float(sxpar(on_mask_header, 'EXPTIME'))
   ;Flat field and rectify it to be consistent with the data that it will be used in correcting
     on_mask_FF = on_mask / blue_flat 
     err_on_mask = abs(on_mask_FF) * sqrt( (err_on_mask/on_mask)^2. + (err_blue_flat/blue_flat)^2. )
     on_mask = temporary(on_mask_FF)
     both = on_mask[*, lower_section_lower_row_b : upper_section_upper_row_b]
     err_both = err_on_mask[*, lower_section_lower_row_b : upper_section_upper_row_b]
     P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_b, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both   
     onset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
     on_mask[onset[0]:*, lower_section_lower_row_b : upper_section_upper_row_b] = both_out    
     err_on_mask[onset[0]:*, lower_section_lower_row_b : upper_section_upper_row_b] = err_both_out   

   ;log the pixel center of jupiter's disc.
    dummy = on_mask
    under_mask = [mask_coords_b[0] + 20., mask_coords_b[1] - 20] ; be sure your only searching regions well inside the mask     
    
    ;get a signal to noise estimate under the mask for every spectral bin
    y = total(on_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
    dy = total(err_on_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
    x = findgen(2098)
    
    ;Determine rebinning needed to have a uniform signal to noise
    new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)

    Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
    for j = 0, N_elements(new_y) - 1 do begin 
      x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
      compress = where(BIN_map eq j, /null)
      y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
      junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
      fit = jjgaussfit(x, y, a, LIMITamp = [0.,1.e10], LIMITcen = [min(x)+10.,max(x)-10.], FIXsl = 0., guess = [max(y), mean(x), (2./R_J_per_pixel)/3., 0., 0.])
      if (total(finite(a)) eq 5) then begin
        disk = where(y gt .1*a[0], /Null)               
        if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the ceter of a region 10% the peak differ by >20 pixels
        disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
        if (abs(median(x[disk]) - a[1]) gt 20.) then continue
        Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]
      endif
    endfor 
    keep = where(finite(Jupiter_centroid)) 
    coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
    on_center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's center
    if keyword_set(debug) then begin
      dummy = on_mask
      dummy[findgen(2098), round(on_center)] = 0
      window, 0, xs = 2098, ys = 1000
      tv, bytscl(dummy, 0., 50) ;inspect the center of Jupiter looks centered
      window, 1
      plot, new_x[keep], Jupiter_centroid[keep], psym=1, yr = under_mask
      oplot, on_center 
    endif

     ;ratio the on mask to off mask to get the throughput using the total over jupiter's disk
     on_mask = on_mask * off_exptime / on_exptime & err_on_mask = err_on_mask * off_exptime / on_exptime
     cumulate_on = fltarr(2098) & cumulate_off = fltarr(2098)
     err_cumulate_on = fltarr(2098) & err_cumulate_off = fltarr(2098)
     for i = 0, 2097 do begin
        cumulate_on[i] = total(on_mask[i, round(on_center[i])-20 : round(on_center[i])+20], /nan)
        cumulate_off[i] = total(off_mask[i, round(off_center[i])-20 : round(off_center[i])+20], /nan)
        err_cumulate_on[i] =  sqrt( total(err_on_mask[i, round(on_center[i])-20 : round(on_center[i])+20]^2., /nan))
        err_cumulate_off[i] = sqrt( total(err_off_mask[i, round(off_center[i])-20 : round(off_center[i])+20]^2., /nan))
     endfor 
     Mask_throughput_b = cumulate_on / cumulate_off
     err_Mask_throughput_b = abs(Mask_throughput_b) * sqrt((err_cumulate_on/cumulate_on)^2. + (err_cumulate_off/cumulate_off)^2.)
     reject = where(((Mask_throughput_b gt .001) or (Mask_throughput_b lt .00001)), /null, COMPLEMENT = keep)
     Mask_throughput_b[reject] = !values.F_nan
     Mask_throughput_b[2048:*] = !values.F_nan
     reject = where(((dispersion_mask_b lt free_spectral_range_b[0]-100.) or (dispersion_mask_b gt free_spectral_range_b[1]+100.)), /null, COMPLEMENT = keep)
     dispersion_mask_b[reject] = !values.F_nan
     if keyword_set(debug) then begin
       window, 2, title = 'Mask Throughput (Blue Channel)'     
       cgplot, dispersion_mask_b[findgen(2098), round(on_center)], Mask_throughput_b, $, ERR_YLOW = err_Mask_throughput_b, ERR_YHIGH = err_Mask_throughput_b, $
         ystyle = 1, psym = 10, xrange = [free_spectral_range_b[0], free_spectral_range_b[1]], xstyle = 1, xtitle = 'Angstroms', ytitle = 'ND Mask Throughput in Blue Channel', color = 'Sky blue', err_thick = .1, err_width = 0, yrange = [.00006,.001]  
       cgplot, dispersion_mask_b[findgen(2098), round(on_center)], SMOOTH( Mask_throughput_b, 100, /EDGE_MIRROR, /NAN), /overplot
     endif 
     residual = abs(SMOOTH( Mask_throughput_b, 100, /EDGE_MIRROR, /NAN) - Mask_throughput_b)
     Mask_throughput_b = SMOOTH( Mask_throughput_b, 100, /EDGE_MIRROR, /NAN)
     err_Mask_throughput_b = rebin([robust_sigma(residual)], 2098)
     WL_Mask_throughput_b = dispersion_mask_b[findgen(2098), round(on_center)]    

   ;==========================Red Channel============================================    
     ;find the center of the disk in an off mask exposure
     off_mask = MRDFITS(directory+'Processed\'+Jupiter_Off_Mask+'r.BSCR.fits', 0, off_mask_header, /Dscale, /silent)
     ;FLIP and rectify it to be consistent with the data that it will be used in correcting  
       off_mask_FF = off_mask / red_flat  
       err_off_mask = abs(off_mask_FF) * sqrt( (err_off_mask/off_mask)^2. + (err_red_flat/red_flat)^2. )
       off_mask = temporary(off_mask_FF) 
       off_mask = rotate(off_mask, 5) & err_off_mask = rotate(err_off_mask, 5) 
       both = off_mask[*, lower_section_lower_row_r : upper_section_upper_row_r]
       err_both = err_off_mask[*, lower_section_lower_row_r : upper_section_upper_row_r]
       P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_r, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both   
       offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
       off_mask[offset[0]:*, lower_section_lower_row_r : upper_section_upper_row_r] = both_out    
       err_off_mask[offset[0]:*, lower_section_lower_row_r : upper_section_upper_row_r] = err_both_out     

     ;log the pixel center of jupiter's disc.
      dummy = off_mask
      under_mask = [lower_section_lower_row_r, upper_section_upper_row_r] ; be sure your only searching regions well inside the aperture    
  
      ;get a signal to noise estimate under the mask for every spectral bin
      y = total(off_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
      dy = total(err_off_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
      x = findgen(2098)
  
      ;Determine rebinning needed to have a uniform signal to noise
      new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
      Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
      for j = 0, N_elements(new_y) - 1 do begin 
        x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
        compress = where(BIN_map eq j, /null)
        y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
        junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
        fit = jjgaussfit(x, y, a)
        if (total(finite(a)) eq 5) then begin
          disk = where(y gt .1*a[0], /Null)               
          if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the ceter of a region 10% the peak differ by >20 pixels
          disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
          if (abs(median(x[disk]) - a[1]) gt 20.) then continue
          Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]
        endif
      endfor 
      keep = where(finite(Jupiter_centroid)) 
      coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
      off_center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's center
      if keyword_set(debug) then begin
        dummy = off_mask
        dummy[findgen(2098), round(off_center)] = 0
        window, 0, xs = 2098, ys = 1000
        tv, bytscl(dummy, 0., 500) ;inspect the center of Jupiter looks centered
        window, 1
        plot, new_x[keep], Jupiter_centroid[keep], psym=1;, yr = [590,610]
        oplot, off_center 
      endif
     
     dummy = off_mask
     dummy[findgen(2098), round(off_center)] = 0.
     s = size(dummy[*, lower_section_lower_row_r:upper_section_upper_row_r], /dimensions)
     if keyword_set(debug) then begin
       window, 0, xs = s[0], ys = s[1]
       tv, bytscl(dummy[*, lower_section_lower_row_r:upper_section_upper_row_r], 0, 1000)
       off_exptime = float(sxpar(off_mask_header, 'EXPTIME'))
     endif
     ;find the center of the disk in an on mask exposure
     on_mask = MRDFITS(directory+'Processed\'+Jupiter_On_Mask+'r.BSCR.fits', 0, on_mask_header, /Dscale, /silent)      
     on_exptime = float(sxpar(on_mask_header, 'EXPTIME'))
     ;FLIP and rectify it to be consistent with the data that it will be used in correcting  
       on_mask_FF = on_mask / red_flat  
       err_on_mask = abs(on_mask_FF) * sqrt( (err_on_mask/on_mask)^2. + (err_red_flat/red_flat)^2. )
       on_mask = temporary(on_mask_FF) 
       on_mask = rotate(on_mask, 5) & err_on_mask = rotate(err_on_mask, 5) 
       both = on_mask[*, lower_section_lower_row_r : upper_section_upper_row_r]
       err_both = err_on_mask[*, lower_section_lower_row_r : upper_section_upper_row_r]
       P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_r, crval, cdelt, both_out, err_both_out, /drizzle, dstack = err_both   
       onset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
       on_mask[onset[0]:*, lower_section_lower_row_r : upper_section_upper_row_r] = both_out    
       err_on_mask[onset[0]:*, lower_section_lower_row_r : upper_section_upper_row_r] = err_both_out   
       
     ;log the pixel center of jupiter's disc.
      dummy = on_mask
      under_mask = [mask_coords_r[0] + 20., mask_coords_r[1] - 20] ; be sure your only searching regions well inside the mask     
      
      ;get a signal to noise estimate under the mask for every spectral bin
      y = total(on_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
      dy = total(err_on_mask[*, under_mask[0]:under_mask[1]], 2, /NAN)
      x = findgen(2098)
      
      ;Determine rebinning needed to have a uniform signal to noise
      new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
  
      Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
      for j = 0, N_elements(new_y) - 1 do begin 
        x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
        compress = where(BIN_map eq j, /null)
        y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
        junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
        fit = jjgaussfit(x, y, a, LIMITamp = [0.,1.e10], LIMITcen = [min(x)+10.,max(x)-10.], FIXsl = 0., guess = [max(y), mean(x), (2./R_J_per_pixel)/3., 0., 0.])
        if (total(finite(a)) eq 5) then begin
          disk = where(y gt .1*a[0], /Null)               
          if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the ceter of a region 10% the peak differ by >20 pixels
          disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
          if (abs(median(x[disk]) - a[1]) gt 20.) then continue
          Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0]
        endif
      endfor 
      keep = where(finite(Jupiter_centroid)) 
      coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 2)
      on_center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's center
      if keyword_set(debug) then begin
        dummy = on_mask
        dummy[findgen(2098), round(on_center)] = 0
        window, 0, xs = 2098, ys = 1000
        tv, bytscl(dummy, 0., 5000) ;inspect the center of Jupiter looks centered
        window, 1
        plot, new_x[keep], Jupiter_centroid[keep], psym=1, yr = under_mask
        oplot, on_center 
      endif

     ;ratio of on mask to off mask to ge the throughput, total over jupiter's disk
     on_mask = on_mask * off_exptime / on_exptime & err_on_mask = err_on_mask * off_exptime / on_exptime
     cumulate_on = fltarr(2098)
     cumulate_off = fltarr(2098)
     for i = 0, 2097 do begin
        cumulate_on[i] = total(on_mask[i, round(on_center[i])-20 : round(on_center[i])+20], /nan)
        cumulate_off[i] = total(off_mask[i, round(off_center[i])-20 : round(off_center[i])+20], /nan)
        err_cumulate_on[i] =  sqrt( total(err_on_mask[i, round(on_center[i])-20 : round(on_center[i])+20]^2., /nan))
        err_cumulate_off[i] = sqrt( total(err_off_mask[i, round(off_center[i])-20 : round(off_center[i])+20]^2., /nan))
     endfor 
     Mask_throughput_r = cumulate_on / cumulate_off
     err_Mask_throughput_r = abs(Mask_throughput_r) * sqrt((err_cumulate_on/cumulate_on)^2. + (err_cumulate_off/cumulate_off)^2.)
     reject = where(((Mask_throughput_r gt .002) or (Mask_throughput_r lt .0005)), /null, COMPLEMENT = keep)
     Mask_throughput_r[reject] = !values.F_nan
     Mask_throughput_r[0:47] = !values.F_nan
     reject = where(((dispersion_mask_r lt free_spectral_range_r[0]-100.) or (dispersion_mask_r gt free_spectral_range_r[1]+100.)), /null, COMPLEMENT = keep)
     dispersion_mask_r[reject] = !values.F_nan
     if keyword_set(debug) then begin
       window, 2, title = 'Mask Throughput (Red Channel)'     
       cgplot, dispersion_mask_r[findgen(2098), round(on_center)], Mask_throughput_r, $, ERR_YLOW = err_Mask_throughput_r, ERR_YHIGH = err_Mask_throughput_r, $
         ystyle = 1, psym = 10, xrange = [free_spectral_range_r[0], free_spectral_range_r[1]], xstyle = 1, xtitle = 'Angstroms', ytitle = 'ND Mask Throughput in Red Channel', color = 'Red', err_thick = .1, err_width = 0, yrange = [.0006,.0018]  
       cgplot, dispersion_mask_r[findgen(2098), round(on_center)], SMOOTH( Mask_throughput_r, 100, /EDGE_MIRROR, /NAN), /overplot
     endif    
     residual = abs(SMOOTH( Mask_throughput_r, 100, /EDGE_MIRROR, /NAN) - Mask_throughput_r)
     Mask_throughput_r = SMOOTH( Mask_throughput_r, 100, /EDGE_MIRROR, /NAN)
     err_Mask_throughput_r = rebin([robust_sigma(residual)], 2098)
     WL_Mask_throughput_r = dispersion_mask_r[findgen(2098), round(on_center)]      
     
     if keyword_set(current_directory) then directory = current_directory ;the mask calibration directory is not necessarily the same as the current => switch back
     save, Mask_throughput_b, err_Mask_throughput_b, WL_Mask_throughput_b, Filename = directory + 'Processed\Mask_throughput_b.sav'
     save, Mask_throughput_r, err_Mask_throughput_r, WL_Mask_throughput_r, Filename = directory + 'Processed\Mask_throughput_r.sav'
endif     

if part eq 5 then begin ;absolute brightness calibration + spatial calibration           
  If ((directory eq 'D:\DATA\Apache Point Data\UT131107\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT131109\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT131217\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT131224\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT140107\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT140209\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT140214\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT140216\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT140219\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT150327\')) then begin
        restore, 'C:\IDL\Io\Apache_Point_Programs\Profiles\master_Mask_throughputs1.sav'  ;use the Master mask throughput, see: ND_Mask_consistency.pro 
        master_throughput_b = master_throughput_b1 & err_master_throughput_b = err_master_throughput_b1 
        master_throughput_r = master_throughput_r1 & err_master_throughput_r = err_master_throughput_r1
  endif
  If ((directory eq 'D:\DATA\Apache Point Data\UT141030\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT141101\')) then begin
        restore, 'C:\IDL\Io\Apache_Point_Programs\Profiles\master_Mask_throughputs3.sav'  ;use the Master mask throughput, see: ND_Mask_consistency.pro 
        master_throughput_b = master_throughput_b3 & err_master_throughput_b = err_master_throughput_b3 
        master_throughput_r = master_throughput_r3 & err_master_throughput_r = err_master_throughput_r3
  endif         
  If ((directory eq 'D:\DATA\Apache Point Data\UT160406\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT160414\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT160415\')) then begin
        restore, 'C:\IDL\Io\Apache_Point_Programs\Profiles\master_Mask_throughputs4.sav'  ;use the Master mask throughput, see: ND_Mask_consistency.pro 
        master_throughput_b = master_throughput_b4 & err_master_throughput_b = err_master_throughput_b4 
        master_throughput_r = master_throughput_r4 & err_master_throughput_r = err_master_throughput_r4
  endif                    
  If ((directory eq 'D:\DATA\Apache Point Data\UT170310\') or $ these def used the same mask, but still good to check for cal errors!
      (directory eq 'D:\DATA\Apache Point Data\UT170325\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170331\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170511\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170519\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170704\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170707\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT170714\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT180309\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT180320\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT180322\') or $
      (directory eq 'D:\DATA\Apache Point Data\UT180330\')) then begin
        restore, 'C:\IDL\Io\Apache_Point_Programs\Profiles\master_Mask_throughputs5.sav'  ;use the Master mask throughput, see: ND_Mask_consistency.pro 
        master_throughput_b = master_throughput_b5 & err_master_throughput_b = err_master_throughput_b5 
        master_throughput_r = master_throughput_r5 & err_master_throughput_r = err_master_throughput_r5
  endif  
  restore, strcompress(Directory + 'Processed/Template_both_b.sav')
  restore, strcompress(Directory + 'Processed/Template_both_r.sav')      
  blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 0, torus_header, /Dscale, /silent )
  red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 0, torus_header, /Dscale, /silent )
  profile_r = total(red_flat, 1)
  profile_b = total(blue_flat, 1)
  aperture_r = where(red_flat gt mean(red_flat))
  aperture_b = where(blue_flat gt mean(blue_flat))
    dummy = blue_flat       
    y_aperture_ind = where( total(blue_flat,1) gt mean(total(blue_flat,1)), /NULL )
    ind = array_indices(blue_flat, aperture_b) 
    dummy[0:max(ind[0,*]),y_aperture_ind] = -100
    aperture_b = where(dummy eq -100, /NULL, complement = anti_aperture_blue) ;a much better aperture.
  mask_r = where(profile_r lt mean(profile_r))
  mask_r = mask_r[where((mask_r gt 350) and (mask_r lt 800))]
  mask_coords_r = [min(mask_r), max(mask_r)]
  mask_b = where(profile_b lt mean(profile_b))
  mask_b = mask_b[where((mask_b gt 200) and (mask_b lt 700))]
  mask_coords_b = [min(mask_b), max(mask_b)]    
  aperture_coords_b = ARRAY_INDICES(blue_flat, aperture_b)  
  lower_section_lower_row_b = min(aperture_coords_b[1,*]) 
  lower_section_upper_row_b = min(mask_b) - 1
  upper_section_lower_row_b = max(mask_b) + 1 
  upper_section_upper_row_b = max(aperture_coords_b[1,*]) 
  aperture_coords_r = ARRAY_INDICES(red_flat, aperture_r)  
  lower_section_lower_row_r = min(aperture_coords_r[1,*]) 
  lower_section_upper_row_r = min(mask_r) - 1
  upper_section_lower_row_r = max(mask_r) + 1 
  upper_section_upper_row_r = max(aperture_coords_r[1,*])   
               
    ;absolute brightness: Data Thief'ed from Woodman et al. 1979.  
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot1.txt', F='A,A', WL, Albedo, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless     
      Woodman_WL = float([WL, WL_2]) ;STITCH THESE TOGETHER
      Woodman_Albedo = Float([albedo , albedo_2]);STITCH THESE TOGETHER
      
    ;absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
      READCOL,'C:\IDL\Io\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
      READCOL,'C:\IDL\Io\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless     
      Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2]) ;STITCH THESE TOGETHER
      Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ;STITCH THESE TOGETHER   
      Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
      Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]
      
    ;compare the two           
      cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', $
          ytitle = 'I/F Reflectivity'          
      cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
      cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
      cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'  
      
    ;make an informed choice via scaling  
      Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
      Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"
    
    ;Woodman only gives I/F to get the output light from Jupiter, multiply this by a solar spectrum
      READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
      start  = where(WL_nm eq '299.100')
      WL_nm = float(WL_nm[start:*]) 
      flux = float(flux[start:*])
      ; change flux units from W/m^2/nm to photons/(cm^2 s A)
      ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)    
      conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
      flux = flux * conversion 
      WL_A = temporary(WL_nm) * 10. ;wavelength from nm into angstroms 
      VACTOAIR, WL_A, WL_A_Air ;Vacuum to air wavelength conversion        
      WL_A = temporary(WL_A_Air)
      ;cross-checked this result against Huebner et al. 1992. 
   
     ;=========================================BLUE Channel===============================================
     Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*b.BSCRJRRec.fits'), count = n_images)           
     for i = 0, n_elements(Images_b)-1 do begin
        
        ;read in images and pull the timestamps mid-exposure
          image = MRDFITS(Images_b[i], 0, torus_header, /Dscale, /silent)
          err_image = MRDFITS(Images_b[i], 1, /Dscale, /silent)
          dispersion_mask = MRDFITS(Images_b[i], 2, /silent)
          cspice_utc2ET, sxpar(torus_header, 'DATE-OBS'), ET
          torus_time = ET + float(sxpar(torus_header, 'EXPTIME'))/2. ;seconds past J2000              

        ;-------------write a few things to the fits headers----------------- 
  
          ;The central meridian longitude: System III sub-Earth longitude accounting for light time
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time, 'IAU_Jupiter', 'LT+S', 'Earth', Sub_Earth, trgepc, srfvec ;get sub observer point in catesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Earth - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, SysIII_LONGITUDE, SysIII_LATITUDE, opgalt ;SysIII_LONGITUDE is the CML in radians      

          ;The Io Phase Angle, in degrees longitude from solar midnight 
              
            ;Get APO's J2000 position
            Observatory, 'APO', obs
            lon = obs.longitude 
            lat = obs.latitude 
            alt = obs.altitude
            cspice_bodvrd, 'EARTH', 'RADII', 3, abc
            equatr =  abc[0]
            polar  =  abc[2]
            f =  ( equatr - polar ) / equatr ;Earth Flatness
            cspice_georec, lon, lat, alt, equatr, f, epos
            cspice_pxform, 'IAU_EARTH', 'J2000', torus_time, rotate
            jpos = transpose(rotate) # epos ;APO's position in J2000
           
            ;Get the light time from Earth, and the Sun-Jupiter and Jupiter-Io vectors accounting for this light time
            cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Junk, ltime ;get earth to Jupiter Light time
            cspice_spkpos, 'Jupiter', torus_time-Ltime, 'J2000', 'None', 'Sun', Sun_Jupiter, junk ;get sun-jupiter vector
            cspice_spkpos, 'Io', torus_time-Ltime, 'J2000', 'None', 'Jupiter', Jupiter_Io, junk ;get sun-jupiter vector

            ;Find a vector from Jupiter's center to it's north pole in the J2000 frame 
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii 
            re = radii[0]   ; equatorial radius
            rp = radii[2]   ; polar radius
            f = (re-rp)/re  ; Jupiter flatness
            North_SYS3 = [0.,0.,1.]*rp ;location of the North Pole in system III
            cspice_pxform, 'IAU_Jupiter', 'J2000', torus_time - ltime, To_J2000 ; Find system III to J2000 rotation matrix for Jupiter
            cspice_mxv, To_J2000, North_SYS3, North_J2000 ; Rotate to J2000
    
            ;Project the Sun-Jupiter vector onto Jupiter's equatorial plane. 
            rot_angle = !dpi/2. - cspice_vsep( Sun_Jupiter, North_J2000 )
            rot_axis = crossp(Sun_Jupiter, North_J2000) ;take the cross-product of the vectors: Jovian North Pole and Sun_Jupiter >> this sets up the rotation axis
            cspice_vrotv, Sun_Jupiter, rot_axis, -rot_angle, Sun_Jupiter_equatorial ;Sun-Jupiter vector projected into Jupiter's equatorial plane
        
            ;Get Io's angular separation from a vector pointing towards equatorial midnight
            Separation = cspice_vsep( Sun_Jupiter_equatorial, Jupiter_Io ) * !radeg

            ;Separation only goes 0-180, fix it for clockwise rotations from midnight 0-360 degrees
            cspice_twovec, North_J2000, 3, Sun_Jupiter_equatorial, 2, xform ;this xform makes Z = Jovian North Pole, Y = Jovian Equator at Midnight, X = right hand completion 
            cspice_mxv, xform, Jupiter_Io, test ;x is dawnward, -x is duskward
            if test[0] lt 0. then IO_PA = Separation else IO_PA = 360. - Separation ;fixed, this result verified against http://www.igpp.ucla.edu/public/pds/wsmyth/ to good accuracy (<0.01 degree)
              
          ;Find the system III longitude of the subsolar point and of the Sun (CML is sub-Earth)
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time - ltime, 'IAU_Jupiter', 'None', 'Sun', Sub_sun, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Sun - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, SysIII_SS_LONGITUDE, SysIII_SS_LATITUDE, opgalt       
            
          ;Find the system III longitude and latitude of Io 
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Io - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt    
          
          ;get Jupiter's instantaneous radius in pixels
            cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', ptarg, ltime
            R_j = 206264.806* atan(re / norm(ptarg))         
            R_J_per_pixel = float(sxpar(torus_header, 'PIXSCAL2')) / R_j 

          ;write all these things to the Fits headers  

            sxaddpar, torus_header, 'EPHEM_TIME', string(torus_time, FORMAT='(F11.1)'), 'Mid-exposure seconds past J2000' ;write to header, accuracy is set to 1/10th of a second
            sxaddpar, torus_header, 'SSP_SYSIII', string(SysIII_SS_LONGITUDE*!radeg), 'System III sub Sun mid-exp' ;write to header 
            sxaddpar, torus_header, 'CML_SYSIII', string(SysIII_LONGITUDE*!radeg), 'System III sub Earth mid-exp' ;write to header   
            sxaddpar, torus_header, 'Io_SYSIII', string(Io_SysIII*!radeg), 'System III sub Io mid-exp' ;write to header  
            sxaddpar, torus_header, 'Io_PA', string(Io_PA), 'Io Phase angle, 0 = Jovian Midnight' ;write to header  
            sxaddpar, torus_header, 'RJ_PER_PIX', string(R_J_per_pixel), 'Jovian radii per pixel' ;write to header     

          if keyword_set(debug) then begin   
            cspice_ET2UTC, torus_time, 'C', 0, UTC 
            print, 'At Earth UTC time: ', UTC
            print, 'Sun''s sub-Jovian System III Longitude =        ', SysIII_SS_LONGITUDE * !radeg
            print, 'Earth''s sub-Jovian System III Longitude (CML) =', SysIII_LONGITUDE * !radeg
            print, 'Io''s sub-Jovian System III Longitude =         ', Io_SysIII * !radeg
            print, 'Io''s Solar Phase angle (0 = Midnight) =        ', IO_PA                        ;checked against http://www.igpp.ucla.edu/public/pds/wsmyth/
          endif   
      ;--------------End writing geometry to the fits headers------------------------------
       
      ;log the pixel center of jupiter's disc.
          dummy = image
          under_mask = [mask_coords_b[0] + 20., mask_coords_b[1] - 20] ; be sure your only searching regions well inside the mask     
          
          ;get a signal to noise estimate under the mask for every spectral bin
          y = total(image[*, under_mask[0]:under_mask[1]], 2, /NAN)
          dy = total(err_image[*, under_mask[0]:under_mask[1]], 2, /NAN)
          x = findgen(2098)

          ;Determine rebinning needed to have a uniform signal to noise
          new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
          Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
          deriv_image = fltarr(N_elements(new_x), under_mask[1] - under_mask[0] + 1)
          s = size(deriv_image, /dimensions)
          for j = 0, N_elements(new_y) - 1 do begin 
            x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
            compress = where(BIN_map eq j, /null)
            y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
            junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,10.], LIMITcen = [min(x)+10.,max(x)-10.], FIXsl = 0., guess = [max(y), mean(x), (2./R_J_per_pixel)/3., 0., 0.])
            if (total(finite(a)) eq 5) then begin
              disk = where(y gt .1*a[0], /Null)               
              if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the center of a region 10% the peak differ by >20 pixels
              disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
              if (abs(median(x[disk]) - a[1]) gt 20.) then continue
              Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0] ;ROUGH CENTROID
            endif
            deriv_image[j,*] = deriv(y,x)   
          endfor   
          deriv_image = 1./abs(deriv_image)
          top = fltarr(s[0]) & bottom =  fltarr(s[0])
          for j = 0, N_elements(new_y) - 1 do begin 
            x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
            y = reform(deriv_image[j,*]) 
            top_guess = Jupiter_centroid[j] - under_mask[0] + 1./R_j_per_pixel
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,.01], LIMITcen = [top_guess - 10., top_guess + 10.], FIXsl = 0., guess = [max(y), top_guess, 2., 0., 0.])
            top[j] = a[1]
            bottom_guess = Jupiter_centroid[j] - under_mask[0] - 1./R_j_per_pixel
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,.01], LIMITcen = [bottom_guess - 10., bottom_guess + 10.], FIXsl = 0., guess = [max(y), bottom_guess, 2., 0., 0.])
            bottom[j] = a[1]
          endfor
          better_centroid = (top+bottom)*.5
          if keyword_set(debug) then begin
            window, 0, xs = s[0], ys = s[1]
            deriv_image[findgen(s[0]), round(better_centroid)] = .1
            deriv_image[findgen(s[0]), round(better_centroid - 1./R_j_per_pixel)] = 0.
            deriv_image[findgen(s[0]), round(better_centroid + 1./R_j_per_pixel)] = 0.
            tv, bytscl(deriv_image, 0., 5.e-3)
          endif  
          Jupiter_centroid = better_centroid + under_mask[0] ;FINE TUNED centroid using the derivatives
          keep = where(finite(Jupiter_centroid) and (new_y gt 0.25), /NULL) 
          coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 1)
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 0)) then coeffs = [417., -0.0035] ;Hack
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 1)) then coeffs = [422., -0.0035] ;Hack
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 2)) then coeffs = [417., -0.0035] ;Hack  
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 3)) then coeffs = [415., -0.0035] ;Hack 
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 4)) then coeffs = [417., -0.0035] ;Hack  
          center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's Brightness center
          phase_angle = SysIII_LONGITUDE - SysIII_SS_LONGITUDE ; need a small correction since the disk isn't fully illuminated
          center = center + .5 * (phase_angle / abs(phase_angle) ) * (1. - cos(phase_angle)) / R_J_per_Pixel ;*Estimated* centroid WRT Jupiter's actual disk center because the disk isn't uniformly ulluminated due to phase angle

          if keyword_set(debug) then begin
            dummy = image
            dummy[findgen(2098), round(center)] = 0.
            dummy[findgen(2098), round(center+.9/R_J_per_pixel)] = 0.
            dummy[findgen(2098), round(center-.9/R_J_per_pixel)] = 0.
            dummy[*,under_mask[0]-20:under_mask[1]] = 50*dummy[*,under_mask[0]-20:under_mask[1]] 
            window, 0, xs = 2098, ys = 1000
            tv, bytscl(hist_equal(dummy, minv= 0., maxv = 1., /percent)) ;inspect the center of Jupiter looks centered
            window, 1
            plot, new_x[keep], Jupiter_centroid[keep], psym=1, yr = [mean(Jupiter_centroid) + [-10,10]]
            oplot, center 
          endif

        ;scale solar flux to Jupiter's instantaneous distance
          cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Sun', ptarg, ltime
          solar_distance = norm(ptarg) / 149597871.
          flux_at_jupiter = flux / solar_distance^2.   
          ;Albedo = INTERPOL(Woodman_Albedo, Woodman_WL, WL_A)          
          Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)   
          Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
          if keyword_set(debug) then begin
            window, 0, Title = 'Rayleighs per Angstrom: Center of Jupiter''s Disk' 
            plot, WL_A, Rayleighs_per_angstrom, xr = [5885,5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981
          endif

        ;Now take the disk center pixel counts and scale them to Jupiter's surface brightness
          ;convolve the "expected" jupiter spectrum by the instrumental line shape of DIS
          blue_arcs = MRDFITS(directory+'Processed/Master_blue_arcs.fits', 0, header, /Dscale, /silent )
          both = blue_arcs[*, lower_section_lower_row_b : upper_section_upper_row_b]
          P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_b, crval, cdelt, both_out, dout, /drizzle       
          offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
          blue_arcs[0: 2097 - offset[0], lower_section_lower_row_b : upper_section_upper_row_b] = both_out  
          convolve_kernel_line = 3888.648 ;A bright arc line that would make a nice convolution kernal
          kernel_indicies = where((dispersion_mask[*,300] gt convolve_kernel_line-5.) and (dispersion_mask[*,300] lt convolve_kernel_line+5.), /NULL)   
          kernel = mean([[blue_arcs[kernel_indicies, 298]],[blue_arcs[kernel_indicies, 299]], [blue_arcs[kernel_indicies, 300]], [blue_arcs[kernel_indicies, 301]],[blue_arcs[kernel_indicies, 302]]], dimension = 2, /NAN) ;Arbitarily pick pixel row 300 to sample the arc line profile.
          WL_Data = dispersion_mask[findgen(2098), center]
          Jupiter_brightness = INTERPOL(Rayleighs_per_angstrom, WL_A, WL_Data)  
          Jupiter_brightness = CONVOL(Jupiter_brightness, kernel, /EDGE_TRUNCATE, /NORMALIZE, /CENTER)
            
        ;Now the Doppler shift, note the following result matches Horizons. 
          cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_state, Jupiter_Sun_light_time ;Jupiter's J2000 position with respect to the Earth
          theta  = cspice_vsep(Jupiter_Sun_state[0:2], Jupiter_Sun_state[3:5])
          Jupiter_WRT_Sun = cos(theta) * sqrt(Jupiter_Sun_state[3]^2.+Jupiter_Sun_state[4]^2.+Jupiter_Sun_state[5]^2.) ;scalar projection of the relative velocity along the line of sight
          cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_state, Jupiter_Earth_light_time ;Jupiter's J2000 position with respect to the Earth
          theta  = cspice_vsep(Jupiter_Earth_state[0:2], Jupiter_Earth_state[3:5])
          Jupiter_WRT_Earth = cos(theta) * sqrt(Jupiter_Earth_state[3]^2.+Jupiter_Earth_state[4]^2.+Jupiter_Earth_state[5]^2.) ;scalar projection of the relative velocity along the line of sight
          shift_jupiter_reflectance = dispersion_mask[findgen(2098), round(center)] * ((Jupiter_WRT_Earth + Jupiter_WRT_Sun) / 299792.458)  
          Angstrom_per_pixel = float(SXPAR(torus_header, 'DISPERS'))
          Jupiter_brightness = smart_shift(Jupiter_brightness, mean(shift_jupiter_reflectance) / Angstrom_per_pixel, /interp, Missing = !values.f_nan) ;shift Jupiter's expected brightness to match it's Doppler Shift at the measurement    
       
        ;Create an image as if there were no mask, convert it from DN to rayleighs per angstrom units
          Mask_throughput_b = interpol(master_throughput_b, WL_master_throughput_b, dispersion_mask[indgen(2098),round(center)], /NaN)
          err_Mask_throughput_b = interpol(err_master_throughput_b, WL_master_throughput_b, dispersion_mask[indgen(2098),round(center)], /NaN)
          no_mask = rebin(Mask_throughput_b, 2098, (mask_coords_b[1] - mask_coords_b[0]) + 1)
          err_no_mask = rebin(err_Mask_throughput_b, 2098, (mask_coords_b[1] - mask_coords_b[0]) + 1)
          dummy = image & err_dummy = err_image  
          no_mask_DN = dummy[*, mask_coords_b[0]:mask_coords_b[1]] / no_mask
          err_no_mask_DN = abs(no_mask_DN) * sqrt( (err_dummy[*, mask_coords_b[0]:mask_coords_b[1]] / dummy[*, mask_coords_b[0]:mask_coords_b[1]] )^2. + (err_no_mask / no_mask)^2. ) 
          no_mask_DN = no_mask_DN*100. & err_no_mask_DN = err_no_mask_DN*100. ;The mask was already divided by 100 to make it visually appealling --> get that 100x back and correct for the mask density.  
          dummy[*, mask_coords_b[0]:mask_coords_b[1]] = temporary(no_mask_DN)
          err_dummy[*, mask_coords_b[0]:mask_coords_b[1]] =  temporary(err_no_mask_DN)
          ;Now the mask region has the acounts that would have been attained without a mask
                            
        ;Obtain the DN in a sample a region at Jupiter's center 
          sample_width = 5. ;Get the appropriate chunk          
          x_sample = rebin(findgen(2098), 2098, sample_width) 
          y_sample = rebin(round(center), 2098, sample_width) + rebin( transpose(indgen(sample_width) - fix(sample_width/2.)), 2098, sample_width) 
          counts = median(dummy[x_sample, y_sample], dimension = 2)
          err_counts = sqrt(!pi/2.) * sqrt(total(err_dummy[x_sample, y_sample]^2., 2, /NAN) / sample_width)

        ;Slight misalignments in wavelength will make noisey results ---> shift to co-align the measured and predicted doppler wells of calcium H and K
          near_min = where((dispersion_mask[indgen(2098),round(center)] gt 3925.) and (dispersion_mask[indgen(2098),round(center)] lt 3942.), /NULL) ;Calcium H region predicted
          junk = min(Jupiter_brightness[min(near_min):max(near_min)], well_pH)
          near_min = where((dispersion_mask[indgen(2098),round(center)] gt 3963.) and (dispersion_mask[indgen(2098),round(center)] lt 3976.), /NULL) ;Calcium K region predicted
          junk = min(Jupiter_brightness[min(near_min):max(near_min)], well_pK)
          near_min = where((dispersion_mask[indgen(2098),round(center)] gt 3925.) and (dispersion_mask[indgen(2098),round(center)] lt 3942.), /NULL) ;Calcium H region measured
          junk = min(counts[min(near_min):max(near_min)], well_mH)      
          near_min = where((dispersion_mask[indgen(2098),round(center)] gt 3963.) and (dispersion_mask[indgen(2098),round(center)] lt 3976.), /NULL) ;Calcium K region measured
          junk = min(counts[min(near_min):max(near_min)], well_mK)
          if keyword_set(debug) then print, 'To match Ca H & K Fraunhofer wells, shifted predicted Jupiter spectrum by', round(mean([well_mH - well_pH, well_mK - well_pK]))  
          Jupiter_brightness = shift(Jupiter_brightness, round(mean([well_mH - well_pH, well_mK - well_pK])))

          dispersion_mask  = smart_shift(dispersion_mask, round(mean([well_mH - well_pH, well_mK - well_pK])), 0, missing = !values.F_Nan) 

        ;Since there's no terrestrial absorption in the Kurucz, 2005 spectrum, scaling by Jupiter_brightness / Counts also corrects for telluric absorption.                         
          SCALE_TO_RAYLEIGHS = Jupiter_brightness * Angstrom_per_pixel / Counts ;expected brightness in (R / DN) = (Rayleighs / Angstrom) * (Angstrom / Pixel) / measured DN                
          err_SCALE_TO_RAYLEIGHS = SCALE_TO_RAYLEIGHS * (err_counts/counts)
            
;        ;Now some fiddly careful data rejection and interpolation  
          rel_err = abs(err_SCALE_TO_RAYLEIGHS / SCALE_TO_RAYLEIGHS)
          anomalies = (rel_err - medsmooth(rel_err, 100)) / medsmooth(rel_err, 100)
          sig = robust_sigma(anomalies)
          vector = SCALE_TO_RAYLEIGHS
          bad = where((abs(anomalies) gt 4.*sig), nbad, COMPLEMENT=good, NCOMPLEMENT=ngood, /NULL)
          IF nbad GT 0 && ngood GT 1 THEN vector[bad] = INTERPOL(vector[good], good, bad)
          anomalies = (vector - medsmooth(vector,100)) / medsmooth(vector,100)
          sig = robust_sigma(anomalies)
          bad = Where((abs(anomalies) gt 4.*sig), nbad, COMPLEMENT=good, NCOMPLEMENT=ngood, /NULL)
          IF nbad GT 0 && ngood GT 1 THEN vector[bad] = INTERPOL(vector[good], good, bad)
          unsmoothed = SCALE_TO_RAYLEIGHS
          SCALE_TO_RAYLEIGHS = medsmooth(SCALE_TO_RAYLEIGHS,15) ;this is the minimum to remove fraunhoffer mismatch lines
                  ;y = vector
                  ;x = findgen(2098)
                  ;x256 = interpol(x, 128, /NAN) ; interpolate to a lower resolution
                  ;y256 = interpol(y, 128, /NAN)
                  ;SPLINE_P, x256, y256, x2, y2 ; spline-interpolate back to original grid
                  ;fit = interpol(y2, x2, findgen(2098), /NAN)
                  ;window, 0, xs = 600, ys = 400
                  ;cgplot, x, y, psym = 3, yr = [300, 700], color = 'blue';, xr = [100,2000] 
                  ;cgplot, x256, y256, color = 'red', /overplot
                  ;cgplot, findgen(2098), fit, color = 'green', /overplot
                  ;cgplot, findgen(2098), medsmooth(vector, 25), color = 'orange', /overplot   
        if keyword_set(debug) then Print, 'Note that Chianti predicts a 3722A / 6312 SIII ratio is 0.3141 that is independent of electron density or temperature.          
        ;Put the images into Rayleigh units    
          Rayleighs = image * rebin(SCALE_TO_RAYLEIGHS, 2098, 1078)
          err_Rayleighs = abs(Rayleighs) * sqrt( (err_image / image )^2. + ( rebin(err_SCALE_TO_RAYLEIGHS / SCALE_TO_RAYLEIGHS, 2098, 1078))^2. ) 
          image = temporary(Rayleighs)
          err_image = temporary(err_Rayleighs)
          image[*,mask_coords_b[0]:mask_coords_b[1]] = 50.*image[*,mask_coords_b[0]:mask_coords_b[1]] ;increase mask region brightness by 50.
       
        ;Inspection    
          if keyword_set(debug) then begin
            window, 0, title = 'Estimated pre-smask counts and propagated error'
            cgplot, counts, err_yhigh = err_counts, err_ylow = err_counts
            dummy[findgen(2098), round(center)] = 0.
            window, 1, xs = 2098, ys = 1000, title = 'Jupiter Brightness Inspection: DN/s Estimate without Mask' 
            tv, bytscl(dummy, 0, 1.e4) ; display the image in DN / s
            window, 2 ;Compare the expected counts, with the measured counts
            cgPLOT, dispersion_mask[findgen(2098), round(center)], Jupiter_brightness, XRANGE = [3920., 3980.], color = 'blue', Title = 'Inspection: MR/A post calibration (Red), Theoretical Jupiter Brightness in MR/A (blue)'
            cgplot, dispersion_mask[findgen(2098), round(center)], counts * SCALE_TO_RAYLEIGHS / Angstrom_per_pixel, charsize = 1.5, color = 'red', psym = 10, /overplot                  
            cgplot, dispersion_mask[findgen(2098), round(center)], counts * 1500. / Angstrom_per_pixel, charsize = 1.5, color = 'orange', psym = 10, /overplot     
            window, 3, xs = 2098, ys = 1000, title = 'Torus Emission Inspection: Rayleighs Per pixel' 
            tv, bytscl(image, 0., 50.) 
            dummy = dummy * rebin(SCALE_TO_RAYLEIGHS, 2098, 1078) 
            window, 4, xs = 2098, ys = 1000, title = 'Jupiter Brightness Inspection: MegaRayleighs Per Angstrom' 
            tv, bytscl(dummy / (Angstrom_per_pixel *  1.e6), 0., 5.) ;display image in MegaRayleighs / Angstrom
            window, 5, title = 'Rayleighs / Pixel / count'     
            cgplot, dispersion_mask[findgen(2098), round(center)], SCALE_TO_RAYLEIGHS, yrange = [30,4500], ystyle = 1, xrange = [dispersion_mask[0, 500], dispersion_mask[2047, 500]], xstyle = 1,$ 
              ytitle = 'Rayleighs x Seconds / DN', xtitle = cgsymbol('Angstrom'), thick = 2.       
            cgplot, dispersion_mask[findgen(2098), round(center)], unsmoothed, /overplot, color = 'blue', psym = 2, SYMSIZE = .1  
          ;stop
          endif   
       
       ;write the final image calibrated into rayleighs to a fits file
       ;1st Extension the error in the calibrated image
       ;2nd Extension is a structure with the dispersion mask and center pixel location of Jupiter's Disk
          Filename = STRMID(Images_b[i], 0, strpos(Images_b[i], 'BSCRJRRec.fits' )) + '_CAL.fits' 
          Metadata_b = ['center','dispersion_mask'] 
          metadata_structure_b = CREATE_STRUCT(Name = 'Metadata_b', metadata_b, center, dispersion_mask)
          MWRFITS, image, Filename, torus_header, /CREATE, /silent ;/create overwrites
          MWRFITS, err_image, Filename, /silent ;Append the fits file -> extension 1
          MWRFITS, metadata_structure_b, Filename, /silent ;Append the fits file -> extension 2
     endfor
         
     ;=========================================RED Channel===============================================
     Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*r.BSCRJRRec.fits'), count = n_images)           
     for i = 0, n_elements(Images_r)-1 do begin
        
        ;read in the images and mid-exposure times
          image = MRDFITS(Images_r[i], 0, torus_header, /Dscale, /silent)
          err_image = MRDFITS(Images_r[i], 1, /Dscale, /silent)
          dispersion_mask = MRDFITS(Images_r[i], 2, /silent)
          cspice_utc2ET, sxpar(torus_header, 'DATE-OBS'), ET
          torus_time = ET + float(sxpar(torus_header, 'EXPTIME'))/2. ;seconds past J2000       
          
        ;-------------write a few things to the fits headers----------------- 
  
          ;The central meridian longitude: System III sub-Earth longitude accounting for light time
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time, 'IAU_Jupiter', 'LT+S', 'Earth', Sub_Earth, trgepc, srfvec ;get sub observer point in catesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Earth - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, SysIII_LONGITUDE, SysIII_LATITUDE, opgalt ;SysIII_LONGITUDE is the CML in radians      

          ;The Io Phase Angle, in degrees longitude from solar midnight 
              
            ;Get APO's J2000 position
            Observatory, 'APO', obs
            lon = obs.longitude 
            lat = obs.latitude 
            alt = obs.altitude
            cspice_bodvrd, 'EARTH', 'RADII', 3, abc
            equatr =  abc[0]
            polar  =  abc[2]
            f =  ( equatr - polar ) / equatr ;Earth Flatness
            cspice_georec, lon, lat, alt, equatr, f, epos
            cspice_pxform, 'IAU_EARTH', 'J2000', torus_time, rotate
            jpos = transpose(rotate) # epos ;APO's position in J2000
           
            ;Get the light time from Earth, and the Sun-Jupiter and Jupiter-Io vectors accounting for this light time
            cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Junk, ltime ;get earth to Jupiter Light time
            cspice_spkpos, 'Jupiter', torus_time-Ltime, 'J2000', 'None', 'Sun', Sun_Jupiter, junk ;get sun-jupiter vector
            cspice_spkpos, 'Io', torus_time-Ltime, 'J2000', 'None', 'Jupiter', Jupiter_Io, junk ;get sun-jupiter vector

            ;Find a vector from Jupiter's center to it's north pole in the J2000 frame 
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii 
            re = radii[0]   ; equatorial radius
            rp = radii[2]   ; polar radius
            f = (re-rp)/re  ; Jupiter flatness
            North_SYS3 = [0.,0.,1.]*rp ;location of the North Pole in system III
            cspice_pxform, 'IAU_Jupiter', 'J2000', torus_time - ltime, To_J2000 ; Find system III to J2000 rotation matrix for Jupiter
            cspice_mxv, To_J2000, North_SYS3, North_J2000 ; Rotate to J2000
    
            ;Project the Sun-Jupiter vector onto Jupiter's equatorial plane. 
            rot_angle = !dpi/2. - cspice_vsep( Sun_Jupiter, North_J2000 )
            rot_axis = crossp(Sun_Jupiter, North_J2000) ;take the cross-product of the vectors: Jovian North Pole and Sun_Jupiter >> this sets up the rotation axis
            cspice_vrotv, Sun_Jupiter, rot_axis, -rot_angle, Sun_Jupiter_equatorial ;Sun-Jupiter vector projected into Jupiter's equatorial plane
        
            ;Get Io's angular separation from a vector pointing towards equatorial midnight
            Separation = cspice_vsep( Sun_Jupiter_equatorial, Jupiter_Io ) * !radeg

            ;Separation only goes 0-180, fix it for clockwise rotations from midnight 0-360 degrees
            cspice_twovec, North_J2000, 3, Sun_Jupiter_equatorial, 2, xform ;this xform makes Z = Jovian North Pole, Y = Jovian Equator at Midnight, X = right hand completion 
            cspice_mxv, xform, Jupiter_Io, test ;x is dawnward, -x is duskward
            if test[0] lt 0. then IO_PA = Separation else IO_PA = 360. - Separation ;fixed, this result verified against http://www.igpp.ucla.edu/public/pds/wsmyth/ to good accuracy (<0.01 degree)
              
          ;Find the system III longitude of the subsolar point and of Io (CML is sub-Earth)
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time - ltime, 'IAU_Jupiter', 'None', 'Sun', Sub_sun, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Sun - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, SysIII_SS_LONGITUDE, SysIII_SS_LATITUDE, opgalt       
          
          ;Find the system III longitude and latitude of Io 
            cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', torus_time - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
            cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
            re = radii[0]
            rp = radii[2]
            f = (re-rp)/re
            obspos = Sub_Io - srfvec
            cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt  
                 
          ;get Jupiter's instantaneous radius in pixels
            cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', ptarg, ltime
            R_j = 206264.806* atan(re / norm(ptarg))         
            R_J_per_pixel = float(sxpar(torus_header, 'PIXSCAL2')) / R_j 

          ;write all these things to the Fits headers   
            sxaddpar, torus_header, 'EPHEM_TIME', string(torus_time, FORMAT='(F11.1)'), 'Mid-exposure seconds past J2000' ;write to header
            sxaddpar, torus_header, 'SSP_SYSIII', string(SysIII_SS_LONGITUDE*!radeg), 'System III sub Sun mid-exp' ;write to header 
            sxaddpar, torus_header, 'CML_SYSIII', string(SysIII_LONGITUDE*!radeg), 'System III sub Earth mid-exp' ;write to header   
            sxaddpar, torus_header, 'Io_SYSIII', string(Io_SysIII*!radeg), 'System III sub Io mid-exp' ;write to header  
            sxaddpar, torus_header, 'Io_PA', string(Io_PA), 'Io Phase angle, 0 = Jovian Midnight' ;write to header  
            sxaddpar, torus_header, 'RJ_PER_PIX', string(R_J_per_pixel), 'Jovian radii per pixel' ;write to header     

          if keyword_set(debug) then begin   
            cspice_ET2UTC, torus_time, 'C', 0, UTC 
            print, 'At Earth UTC time: ', UTC
            print, 'Sun''s sub-Jovian System III Longitude =        ', SysIII_SS_LONGITUDE * !radeg
            print, 'Earth''s sub-Jovian System III Longitude (CML) =', SysIII_LONGITUDE * !radeg
            print, 'Io''s sub-Jovian System III Longitude =         ', Io_SysIII * !radeg
            print, 'Io''s Solar Phase angle (0 = Midnight) =        ', IO_PA                        ;checked against http://www.igpp.ucla.edu/public/pds/wsmyth/
          endif   
      ;--------------End writing geometry to the fits headers------------------------------         

        ;log the pixel center of jupiter's disc.
              dummy = image
              under_mask = [mask_coords_r[0] + 20., mask_coords_r[1] - 20] ; be sure your only searching regions well inside the mask     
              
              ;get a signal to noise estimate under the mask for every spectral bin
              y = total(image[*, under_mask[0]:under_mask[1]], 2, /NAN)
              dy = total(err_image[*, under_mask[0]:under_mask[1]], 2, /NAN)
              x = findgen(2098)

          ;Determine rebinning needed to have a uniform signal to noise
          new_y = REGROUP(y, x, dy = dy, SNR = max(y/dy), bin_x = new_x, BIN_dy = new_dy, BIN_map = BIN_map)
          Jupiter_centroid = make_array(N_elements(new_y), /float, value = !values.F_nan)
          deriv_image = fltarr(N_elements(new_x), under_mask[1] - under_mask[0] + 1)
          s = size(deriv_image, /dimensions)
          for j = 0, N_elements(new_y) - 1 do begin 
            x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
            compress = where(BIN_map eq j, /null)
            y = total(dummy[compress, under_mask[0]:under_mask[1]], 1, /Nan)
            junk = max(medsmooth(y, 25), peak) ;the peak of the median smoothed profile
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,10.], LIMITcen = [min(x)+10.,max(x)-10.], FIXsl = 0., guess = [max(y), mean(x), (2./R_J_per_pixel)/3., 0., 0.])
            if (total(finite(a)) eq 5) then begin
              disk = where(y gt .1*a[0], /Null)               
              if (disk eq !NULL) then continue ;reject of the Gaussian centroid and the pixels at the center of a region 10% the peak differ by >20 pixels
              disk = disk[where(((disk gt a[1] - 1.2/R_J_per_pixel) and (disk lt 1.2/R_J_per_pixel + a[1])), /NULL)]
              if (abs(median(x[disk]) - a[1]) gt 20.) then continue
              Jupiter_centroid[j] = median(x[disk], /even) + under_mask[0] ;ROUGH CENTROID
            endif
            deriv_image[j,*] = deriv(y,x)   
          endfor   
          deriv_image = 1./abs(deriv_image)
          top = fltarr(s[0]) & bottom =  fltarr(s[0])
          for j = 0, N_elements(new_y) - 1 do begin 
            x = findgen(n_elements(dummy[0, under_mask[0]:under_mask[1]]))
            y = reform(deriv_image[j,*]) 
            top_guess = Jupiter_centroid[j] - under_mask[0] + 1./R_j_per_pixel
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,.02], LIMITcen = [top_guess - 10., top_guess + 10.], FIXsl = 0., guess = [max(y), top_guess, 2., 0., 0.])
            top[j] = a[1]
            bottom_guess = Jupiter_centroid[j] - under_mask[0] - 1./R_j_per_pixel
            fit = jjgaussfit(x, y, a, LIMITamp = [0.,.02], LIMITcen = [bottom_guess - 10., bottom_guess + 10.], FIXsl = 0., guess = [max(y), bottom_guess, 2., 0., 0.])
            bottom[j] = a[1]
          endfor
          better_centroid = (top+bottom)*.5
          if keyword_set(debug) then begin
            window, 0, xs = s[0], ys = s[1]
            deriv_image[findgen(s[0]), round(better_centroid)] = .1
            deriv_image[findgen(s[0]), round(better_centroid - 1./R_j_per_pixel)] = 0.
            deriv_image[findgen(s[0]), round(better_centroid + 1./R_j_per_pixel)] = 0.
            tv, bytscl(deriv_image, 0., 1.e-2)            
          endif  
          Jupiter_centroid = better_centroid + under_mask[0] ;FINE TUNED centroid using the derivatives
          keep = where(finite(Jupiter_centroid) and (new_y gt 0.25), /NULL) 

          coeffs = robust_POLY_FIT( new_x[keep], Jupiter_centroid[keep], 1)
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 0)) then coeffs = [553, -0.005] ;Hack
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 1)) then coeffs = [555, -0.005] ;Hack
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 2)) then coeffs = [552., -0.005] ;Hack  
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 3)) then coeffs = [551., -0.0065] ;Hack 
          if ((directory eq 'D:\DATA\Apache Point Data\UT131109\') and (i eq 4)) then coeffs = [552., -0.005] ;Hack  
          center = poly(findgen(2098), coeffs) ;center is now the pixel location of jupiter's center
          phase_angle = SysIII_LONGITUDE - SysIII_SS_LONGITUDE ; need a small correction since the disk isn't fully illuminated
          center = center + .5 * (phase_angle / abs(phase_angle) ) * (1. - cos(phase_angle)) / R_J_per_Pixel

          if keyword_set(debug) then begin
            dummy = image
            dummy[findgen(2098), round(center)] = 0.
            dummy[findgen(2098), round(center)+.9/R_J_per_pixel] = 0.
            dummy[findgen(2098), round(center)-.9/R_J_per_pixel] = 0.
            ;dummy[*,under_mask[0]:under_mask[1]] = 5*dummy[*,under_mask[0]:under_mask[1]] 
            window, 0, xs = 2098, ys = 1000
            tv, bytscl(dummy, 0., 1.e-1) ;inspect the center of Jupiter looks centered
            window, 1
            plot, new_x[keep], Jupiter_centroid[keep], psym=1, yr = [mean(Jupiter_centroid) + [-20, 20]]
            oplot, center 
            print, coeffs
          endif
        
        ;scale solar flux to Jupiter's instantaneous distance
          cspice_spkpos, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Sun', ptarg, ltime
          solar_distance = norm(ptarg) / 149597871.
          flux_at_jupiter = flux / solar_distance^2.   
          ;Albedo = INTERPOL(Woodman_Albedo, Woodman_WL, WL_A)          
          Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)    
          Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
          if keyword_set(debug) then begin
            window, 0, Title = 'Rayleighs per Angstrom: Center of Jupiter''s Disk' 
            plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981
            ;plot, WL_A, Rayleighs_per_angstrom, xr = [6700, 6750], charsize = 2 ;compare to ???
            ;stop
          endif

        ;Now take the disk center pixel counts and scale them to Jupiter's surface brightness  
          ;convolve the "expected" jupiter spectrum by the instrumental line shape of DIS
            red_arcs = MRDFITS(directory+'Processed/Master_red_arcs.fits', 0, header, /Dscale, /silent )
            red_arcs = rotate(red_arcs, 5) ;FLIP and rectify it to be consistent with the data 
            both = red_arcs[*, lower_section_lower_row_r : upper_section_upper_row_r]
            P3D_WAVECAL_DISPERSION_CORRECTION, both, Template_both_r, crval, cdelt, both_out, dout, /drizzle       
            offset = size(both, /dimensions) - size(both_out, /dimensions) ;rectified image looses some pixels at the edges
            red_arcs[0:2097-offset[0], lower_section_lower_row_r : upper_section_upper_row_r] = both_out  
            convolve_kernel_line = 6598.9528 ;A bright arc line that would make a nice convolution kernal
            kernel_indicies = where((dispersion_mask[*,300] gt convolve_kernel_line-5.) and (dispersion_mask[*,300] lt convolve_kernel_line+5.), /NULL) ;Arbitarily pick pixel row 300 to sample the arc line profile.          
            kernel = mean([[red_arcs[kernel_indicies, 298]],[red_arcs[kernel_indicies, 299]], [red_arcs[kernel_indicies, 300]], [red_arcs[kernel_indicies, 301]],[red_arcs[kernel_indicies, 302]]], dimension = 2, /NAN)         
            WL_Data = dispersion_mask[findgen(2098), center]
            Jupiter_brightness = INTERPOL(Rayleighs_per_angstrom, WL_A, WL_Data)  
            Jupiter_brightness = CONVOL(Jupiter_brightness, kernel, /EDGE_TRUNCATE, /NORMALIZE, /CENTER)

          ;Now add the round trip Doppler shift of Jupiter, note the following result matches Horizons. 
            cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_state, Jupiter_Sun_light_time ;Jupiter's J2000 position with respect to the Earth
            theta  = cspice_vsep(Jupiter_Sun_state[0:2], Jupiter_Sun_state[3:5])
            Jupiter_WRT_Sun = cos(theta) * sqrt(Jupiter_Sun_state[3]^2.+Jupiter_Sun_state[4]^2.+Jupiter_Sun_state[5]^2.) ;scalar projection of the relative velocity along the line of sight
            cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_state, Jupiter_Earth_light_time ;Jupiter's J2000 position with respect to the Earth
            theta  = cspice_vsep(Jupiter_Earth_state[0:2], Jupiter_Earth_state[3:5])
            Jupiter_WRT_Earth = cos(theta) * sqrt(Jupiter_Earth_state[3]^2.+Jupiter_Earth_state[4]^2.+Jupiter_Earth_state[5]^2.) ;scalar projection of the relative velocity along the line of sight
            shift_jupiter_reflectance = dispersion_mask[findgen(2098), round(center)] * ((Jupiter_WRT_Earth + Jupiter_WRT_Sun) / 299792.458)      
            Angstrom_per_pixel = float(SXPAR(torus_header, 'DISPERS'))
            Jupiter_brightness = smart_shift(Jupiter_brightness, mean(shift_jupiter_reflectance) / Angstrom_per_pixel, /interp, Missing = !values.f_nan) ;shift Jupiter's expected brightness to match it's Doppler Shift at the measurement

          ;Create an image as if there were no mask, convert it from DN to rayleighs per angstrom units
            Mask_throughput_r = interpol(master_throughput_r, WL_master_throughput_r, dispersion_mask[indgen(2098),round(center)], /NaN)
            err_Mask_throughput_r = interpol(err_master_throughput_r, WL_master_throughput_r, dispersion_mask[indgen(2098),round(center)], /NaN)
            no_mask = rebin(Mask_throughput_r, 2098, (mask_coords_r[1] - mask_coords_r[0]) + 1)
            err_no_mask = rebin(err_Mask_throughput_r, 2098, (mask_coords_r[1] - mask_coords_r[0]) + 1)
            dummy = image & err_dummy = err_image  
            no_mask_DN = dummy[*, mask_coords_r[0]:mask_coords_r[1]] / no_mask
            err_no_mask_DN = abs(no_mask_DN) * sqrt( (err_dummy[*, mask_coords_r[0]:mask_coords_r[1]] / dummy[*, mask_coords_r[0]:mask_coords_r[1]] )^2. + (err_no_mask / no_mask)^2. ) 
            no_mask_DN = no_mask_DN*100. & err_no_mask_DN = err_no_mask_DN*100. ;The mask was already divided by 100 to make it visually appealling --> get that 100x back and correct for the mask density.  
            dummy[*, mask_coords_r[0]:mask_coords_r[1]] = temporary(no_mask_DN)
            err_dummy[*, mask_coords_r[0]:mask_coords_r[1]] =  temporary(err_no_mask_DN)
            ;Now the mask region of the dummy array has the counts that would have been attained without a mask
                               
          ;Obtain the DN in a sample a region at Jupiter's center 
            sample_width = 5. ;Get the appropriate chunk          
            x_sample = rebin(findgen(2098), 2098, sample_width) 
            y_sample = rebin(round(center), 2098, sample_width) + rebin( transpose(indgen(sample_width) - fix(sample_width/2.)), 2098, sample_width) 
            counts = median(dummy[x_sample, y_sample], dimension = 2)
            err_counts = sqrt(!pi/2.) * sqrt(total(err_dummy[x_sample, y_sample]^2., 2, /NAN) / sample_width)
      
          ;Slight misalignments in wavelength will make noisey results ---> shift to co-align the measured and predicted doppler well of H alpha
            near_min = where((dispersion_mask[indgen(2098),round(center)] gt 6558.) and (dispersion_mask[indgen(2098),round(center)] lt 6568.), /NULL) ;H alpha region predicted
            junk = min(Jupiter_brightness[min(near_min):max(near_min)], well_pH)
            near_min = where((dispersion_mask[indgen(2098),round(center)] gt 6558.) and (dispersion_mask[indgen(2098),round(center)] lt 6568.), /NULL) ;Calcium H region measured
            junk = min(counts[min(near_min):max(near_min)], well_mH)      
            if keyword_set(debug) then print, 'To match Halpha Fraunhofer wells, shifted predicted Jupiter spectrum by', well_mH - well_pH
            Jupiter_brightness = shift(Jupiter_brightness, well_mH - well_pH)
      
          ;BETA Shift the disperion mask by this same increment, it was never accurate to better than a pixel anyhow! 
            dispersion_mask  = smart_shift(dispersion_mask, well_mH - well_pH, 0, missing = !values.F_Nan) 
      
          ;Since there's no terrestrial absorption in the Kurucz, 2005 spectrum, scaling by Jupiter_brightness / Counts also corrects for telluric absorption.                         
            SCALE_TO_RAYLEIGHS = Jupiter_brightness * Angstrom_per_pixel / Counts ;expected brightness in (R / DN) = (Rayleighs / Angstrom) * (Angstrom / Pixel) / measured DN                
            err_SCALE_TO_RAYLEIGHS = SCALE_TO_RAYLEIGHS * (err_counts/counts)
            
          ;Now some fiddly careful data rejection and interpolation  
            rel_err = abs(err_SCALE_TO_RAYLEIGHS / SCALE_TO_RAYLEIGHS)
            anomalies = (rel_err - medsmooth(rel_err, 100)) / medsmooth(rel_err, 100)
            sig = robust_sigma(anomalies)
            vector = SCALE_TO_RAYLEIGHS
            bad = where((abs(anomalies) gt 4.*sig), nbad, COMPLEMENT=good, NCOMPLEMENT=ngood, /NULL)
            IF nbad GT 0 && ngood GT 1 THEN vector[bad] = INTERPOL(vector[good], good, bad)
            anomalies = (vector - medsmooth(vector, 100)) / medsmooth(vector, 100)
            sig = robust_sigma(anomalies)
            bad = Where((abs(anomalies) gt 4.*sig), nbad, COMPLEMENT=good, NCOMPLEMENT=ngood, /NULL)
            IF nbad GT 0 && ngood GT 1 THEN vector[bad] = INTERPOL(vector[good], good, bad)
            unsmoothed = SCALE_TO_RAYLEIGHS
            SCALE_TO_RAYLEIGHS = medsmooth(SCALE_TO_RAYLEIGHS, 15) ;tried a few things, this seems the best result
   
          ;Put the images into Rayleigh units    
            Rayleighs = image * rebin(SCALE_TO_RAYLEIGHS, 2098, 1078)
            err_Rayleighs = abs(Rayleighs) * sqrt( (err_image / image )^2. + ( rebin(err_SCALE_TO_RAYLEIGHS / SCALE_TO_RAYLEIGHS, 2098, 1078))^2. ) 
            image = temporary(Rayleighs)
            err_image = temporary(err_Rayleighs)
            image[*,mask_coords_r[0]:mask_coords_r[1]] = 5.*image[*,mask_coords_r[0]:mask_coords_r[1]] ;increase mask region brightness by 5.

          ;Inspection    
          if keyword_set(debug) then begin
            window, 0, title = 'Estimated pre-smask counts and propagated error'
            cgplot, counts, err_yhigh = err_counts, err_ylow = err_counts
            dummy[findgen(2098), round(center)] = 0.
            window, 1, xs = 2098, ys = 1000, title = 'Jupiter Brightness Inspection: DN/s Estimate without Mask' 
            tv, bytscl(dummy, 0, 1.e4) ; display the image in DN / s
            window, 2, title = 'Jovian Reflected H Alpha Wavelength (A) = ' + string(6562.801 + shift_jupiter_reflectance[690])
            ;cgPLOT, dispersion_mask[findgen(2098), round(center)], Jupiter_brightness, XRANGE = [6555, 6570.], color = 'blue', Title = 'Inspection: MR/A post calibration (Red), Theoretical Jupiter Brightness in MR/A (blue)'
            cgPLOT, dispersion_mask[findgen(2098), round(center)], Jupiter_brightness, XRANGE = [5850, 6870.], color = 'blue', yr = [2.e6,7.e6],$
               Title = 'Inspection: R/A post calibration (Red), Theoretical Jupiter Brightness in R/A (blue)', xtitle = cgsymbol('Angstrom')          
            cgplot, dispersion_mask[findgen(2098), round(center)], counts*SCALE_TO_RAYLEIGHS/Angstrom_per_pixel, charsize = 1.5, color = 'red', psym = 10, /overplot                  
            window, 3, xs = 2098, ys = 1000, title = 'Torus Emission Inspection: Rayleighs Per pixel' 
            tv, bytscl(image, 0., 100.) 
            dummy = dummy * rebin(SCALE_TO_RAYLEIGHS, 2098, 1078) 
            window, 4, xs = 2098, ys = 1000, title = 'Jupiter Brightness Inspection: MegaRayleighs Per Angstrom' 
            tv, bytscl(dummy / (Angstrom_per_pixel *  1.e6), 0., 5.) 
            window, 5, title = 'Rayleighs / Pixel / count'     
            cgplot, dispersion_mask[findgen(2098), round(center)], SCALE_TO_RAYLEIGHS, yrange = [350,620], ystyle = 1, xrange = [dispersion_mask[0, 500], dispersion_mask[2047, 500]], xstyle = 1,$ 
              ytitle = 'Rayleighs x Seconds / DN', xtitle = cgsymbol('Angstrom'), thick = 2.   
            cgplot, dispersion_mask[findgen(2098), round(center)], unsmoothed, /overplot, color = 'red', psym = 2, SYMSIZE = .1    
          
            ;to get MR/A back from the calibrated image for jupiter brightness test
            cross_check = image
            cross_check[*,mask_coords_r[0]:mask_coords_r[1]] = cross_check[*,mask_coords_r[0]:mask_coords_r[1]] * 100./(5.*No_mask) * (1./(Angstrom_per_pixel*1.e6))
            window, 6, xs = 2098, ys = 1000, title = 'Double Check' 
            tv, bytscl(cross_check, 0, 6)
          endif     

         ;write the final image calibrated into rayleighs to a fits file
         ;1st Extension is the error
         ;2nd Extension is a structure witht the dispersion mask and center pixel location of Jupiter's Disk
            Filename = STRMID(Images_r[i], 0, strpos(Images_r[i], 'BSCRJRRec.fits' )) + '_CAL.fits' 
            Metadata_r = ['center','dispersion_mask'] 
            metadata_structure_r = CREATE_STRUCT(Name = 'Metadata_r', metadata_r, center, dispersion_mask)
            MWRFITS, image, Filename, torus_header, /CREATE, /silent ;/create overwrites
            MWRFITS, err_image, Filename, /silent  ;Append the fits file -> extension 1
            MWRFITS, metadata_structure_r, Filename, /silent  ;Append the fits file with this structure -> extension 2
     endfor 
endif

if part eq 6 then begin ;get all data products from the calibrated images!  
   Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*b._Cal.fits'), count = n_images) 
   Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*r._Cal.fits'), count = n_images) 
   
   ;Get Jovian and other contants
     cspice_bodvrd, 'Jupiter', "RADII", 5, radii  
     R_J = radii[0]                               ;Jovian radii, in KM
     cspice_bodvrd, 'Jupiter', "PM", 5, PM        ;Jovian rotation rate in deg per day
     period_J = 24.*3600. / (PM[1] / 360.)        ;Jovian period, in s                 
     c = cspice_clight()             ;speed of light (km / s)
   
   ;setup MPFIT Gaussian parameters
    sample_width = 7. 
    guess = [0.,sample_width,1.,0.,0.]
    A = guess
    parinfo = replicate( {value: 0., fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 5)
    parinfo[0].limited = 1b                                 ;limit amplitude 
    parinfo[0].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS. 
    parinfo[1].limited = 1b                                 ;limit centroid
    parinfo[1].limits = [sample_width-3., sample_width+3.]  ;centroid limits in pixels
    parinfo[2].fixed = 1b                                   ;limit linewidth sigma, it's wavelength dependent so revalue this for each emission line    
    parinfo[3].fixed = 1b                                   ;Zero dc
    parinfo[3].limited = 1b                                 ;limit DC Background     
    parinfo[4].fixed = 1b                                   ;limit slope to zero
 
   ;setup MPFIT Gaussian parameters for triplet line fitting:
    sample_width_Blue_triplet = 10 
    parinfo_blue_triplet = replicate( { value: 0., fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 5)
    parinfo_blue_triplet[0].limited = 1b                                 ;limit 3722 amplitude 
    parinfo_blue_triplet[0].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS. 
    parinfo_blue_triplet[1].limited = 1b                                 ;limit 3726 amplitude 
    parinfo_blue_triplet[1].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS. 
    parinfo_blue_triplet[2].limited = 1b                                 ;limit 3729 amplitude 
    parinfo_blue_triplet[2].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS.       
    parinfo_blue_triplet[3].limited = 1b                                 ;limit centroid
    parinfo_blue_triplet[3].limits = [sample_width_Blue_triplet-2., sample_width_Blue_triplet+2.]  ;4069 centroid limits in pixels
    parinfo_blue_triplet[4].fixed = 1b                                   ;Fix sigma   
    parinfo_blue_triplet[4].value = 0.                                   ;Set sigma based on the width of lines measured in arclamp frames, this is image dependent so it's set in the for loops    
    ;parinfo_blue_triplet[5].fixed = 1b                                   ;zero DC background
                   
   ;setup MPFIT Gaussian parameters for doublet line fitting:
    sample_width_Blue_doublet = 20 
    parinfo_Blue_doublet = replicate( {value: 0., fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 4)
    parinfo_Blue_doublet[0].limited = 1b                                 ;limit 4069 amplitude 
    parinfo_Blue_doublet[0].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS. 
    parinfo_Blue_doublet[1].limited = 1b                                 ;limit 4076 amplitude 
    parinfo_Blue_doublet[1].limits = [-100., 1.e4]                        ;positive AND A LITTLE NEGATIVE, THIS DOESN'T AFFECT WHETHER THERE'S ZERO EMISSION FAR FROM THE TORUS.      
    parinfo_Blue_doublet[2].limited = 1b                                 ;limit centroid
    parinfo_Blue_doublet[2].limits = [sample_width_Blue_doublet-3., sample_width_Blue_doublet+3.]  ;3726 centroid limits in pixels
    parinfo_Blue_doublet[3].fixed = 1b                                   ;Fix sigma   
    parinfo_Blue_doublet[3].value = 0.                                   ;Set sigma based on the width of lines measured in arclamp frames, this is image dependent so it's set in the for loops  
    ;parinfo_Blue_doublet[4].fixed = 1b                                   ;zero DC background     
    
    ;-------------------Red Channel----------------------------
    max_dawn_r = 0 & max_dusk_r = 0 ;the furthest pixels with any data of all the images considered

    for i = 0, n_elements(Images_r)-1 do begin
        
        ;Calculate the maximum distance in RJ where any data is available
          image = MRDFITS(Images_r[i], 0, torus_header_r, /silent)
          err_image = MRDFITS(Images_r[i], 1, /silent)
          dummy = image & err_dummy = err_image
          metadata_r = MRDFITS(Images_r[i], 2, /silent)
          Dispersion_mask_r = metadata_r.dispersion_mask 
          center_r = metadata_r.center
          R_J_per_pixel_r = float(sxpar(torus_header_r, 'RJ_PER_P')) 
          dispers_r = float(sxpar(torus_header_r, 'DISPERS'))  
          dispers = float(sxpar(torus_header_r, 'DISPERS')) 
          aperture_coords = ARRAY_INDICES(image, aperture_r)  
          distance_extent_dawn = abs(max(aperture_coords[1,*]) - center_r[1024]) * R_J_per_pixel_r - 0.2 ;Distance in RJ out to which torus data is present 
          distance_extent_dusk = abs(min(aperture_coords[1,*]) - center_r[1024]) * R_J_per_pixel_r - 0.2 ;Distance in RJ out to which torus data is present 
          linewidths_r = poly(lines_r, float([sxpar(torus_header_r, 'W0_COEF'), sxpar(torus_header_r, 'W1_COEF')])) 
 
        ;Calculate Jupiter's radial velocity WRT Earth (Earth center not observatory)
          torus_time = double(sxpar(torus_header_r, 'EPHEM_TI')) 
          cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_state, Jupiter_Earth_light_time ;Jupiter's J2000 position with respect to the Earth
          theta  = cspice_vsep(Jupiter_Earth_state[0:2], Jupiter_Earth_state[3:5])
          Jupiter_WRT_Earth = cos(theta) * sqrt(Jupiter_Earth_state[3]^2.+Jupiter_Earth_state[4]^2.+Jupiter_Earth_state[5]^2.) ;scalar projection of the relative velocity along the line of sight

        ;Calculate where the measured line centers are. Use the fit centroid position wavelength where the ribbon should be. For all lines
          CML = float(sxpar(torus_header_r, 'CML_SYSI'))
          CML_Dawn = (CML + 90. + 360.0) MOD 360.0 
          CML_Dusk = (CML - 90. + 360.0) MOD 360.0 
          ribbon_range_dawn = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          ribbon_range_dusk = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
          expected_ribbon_dawn = ribbon_range_dawn / R_J_per_pixel_r
          expected_ribbon_dusk = ribbon_range_dusk / R_J_per_pixel_r
          Line_centers = fltarr(N_elements(LINES_R))
          for n = 0, N_elements(LINES_R) - 1 do begin
            dawn_shift = (-(2. * !pi * R_J * ribbon_range_dawn * lines_r[n]/ period_J) / (c * dispers_r)) + (lines_r[n] * Jupiter_WRT_Earth) / (c * dispers_r) 
            dusk_shift = ( (2. * !pi * R_J * ribbon_range_dusk * lines_r[n]/ period_J) / (c * dispers_r)) + (lines_r[n] * Jupiter_WRT_Earth) / (c * dispers_r) 
            junk = min(abs(Dispersion_mask_r[indgen(2098), center_r] - Lines_r[n]), /NAN, origin) ;compute the line center pixel from the dispersion mask
            dawn_fit = total(Image[origin - sample_width + dawn_shift : origin + sample_width + dawn_shift, $ ;sum 10 pixels over the ribbon
                             center_r[[origin]] + expected_ribbon_dawn - 5: center_r[[origin]] + expected_ribbon_dawn + 5], 2) 
            dusk_fit = total(Image[origin - sample_width + dusk_shift : origin + sample_width + dusk_shift, $ ;sum 10 pixels over the ribbon
                             center_r[[origin]] - expected_ribbon_dusk - 5: center_r[[origin]] - expected_ribbon_dusk + 5], 2) 
            err_dawn_fit = sqrt(total(err_Image[origin - sample_width + dawn_shift : origin + sample_width + dawn_shift, $ ;sum 10 pixels over the ribbon
                                 center_r[[origin]] + expected_ribbon_dawn - 5: center_r[[origin]] + expected_ribbon_dawn + 5]^2, 2)) 
            err_dusk_fit = sqrt(total(err_Image[origin - sample_width + dusk_shift : origin + sample_width + dusk_shift, $ ;sum 10 pixels over the ribbon
                                 center_r[[origin]] - expected_ribbon_dusk - 5: center_r[[origin]] - expected_ribbon_dusk + 5]^2, 2)) 
            parinfo[2].value = linewidths_r[n] & guess[2] = linewidths_r[n]
            fit_dawn = mpfitpeak(findgen(n_elements(dawn_fit)), dawn_fit, A_dawn, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dawn_fit, ESTIMATES = guess, status = status, /NAN)
            fit_dusk = mpfitpeak(findgen(n_elements(dusk_fit)), dusk_fit, A_dusk, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dusk_fit, ESTIMATES = guess, status = status, /NAN)
                   
            ;REFINE THE LINE CENTER'S BASED ON THE FITS
            Line_centers[n] = origin + (mean([a_dawn[1],a_dusk[1]], /Nan) - sample_width)  
            if keyword_set(debug) then begin              
              if n eq 0 then window, 0, xs = 600, ys = 400
              wset, 0
              if n eq 0 then cgplot, dawn_fit, psym = 4, color = red_colors[n], yr = [-100, 1700] else cgplot, dawn_fit, psym = 4, color = red_colors[n], /overplot     
              ;print, 'shifted Dawn centroid of ', lines_r[n], ' line by ', (a_dawn[1]- sample_width) 
              cgplot, fit_dawn, color = red_colors[n], /overplot
              
              if n eq 0 then window, 1, xs = 600, ys = 400
              wset, 1
              if n eq 0 then cgplot, dusk_fit, psym = 4, color = red_colors[n], yr = [-100, 1700] else cgplot, dusk_fit, psym = 4, color = red_colors[n], /overplot     
              ;print, 'shifted Dusk centroid of ', lines_r[n], ' line by ', (a_dusk[1]- sample_width) 
              cgplot, fit_dusk, color = red_colors[n], /overplot
            endif
          endfor
        if keyword_set(debug) then begin    
          range = indgen(2098)
          y = total(Image[ min(range):max(range), mean(center_r[range]) + expected_ribbon_dawn - 15: mean(center_r[range]) + expected_ribbon_dawn + 15 ], 2, /NAN) 
          x = Dispersion_mask_r[range, mean(center_r[range])]
          wset, 0
          cgplot, x, y, xr = [6710, 6740]
          cgplot, Line_centers, make_array(N_elements(Line_centers), /float, value = 1000.), /overplot, psym =4
          cgplot, Dispersion_mask_r[round(Line_centers), center_r[round(Line_centers)]] + dawn_shift*dispers, make_array(N_elements(Line_centers), /float, value = 1000.), psym =4, /overplot
        endif
                 
        pixel_fraction = center_r[round(Line_centers)] - fix(center_r[round(Line_centers)]) ; the fractional pixel distance of the centroid

        ;--------------------dawn---------------------------        
        dawn_distance_array_r = findgen(distance_extent_dawn / R_J_per_pixel_r)
        dawn_distance_array_r_km = R_J * R_J_per_pixel_r * dawn_distance_array_r
        dawn_shift = -((2. * !pi * dawn_distance_array_r_km / period_J) / c)  + (Jupiter_WRT_Earth / c); just v over c, multiply by lambda for the wavelength shift
        dawn_line_profiles = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r)) & err_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r))
        dawn_line_profiles_offband = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r)) & err_dawn_line_profiles_offband = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r))
        fit_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r)) & err_fit_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_r), N_elements(lines_r))
        for j = 0, N_elements(lines_r)-1 do begin ;On band loop
          
          ;Straight up summing the counts to extract emissions, this is generally less accurate! 
          dawn_x_sample = rebin(round(Line_centers[j] + dawn_shift * lines_r[j] / dispers_r), n_elements(dawn_distance_array_r), sample_width) + $
                          rebin(transpose(indgen(sample_width) - fix(sample_width/2)), n_elements(dawn_distance_array_r), sample_width)
          dawn_y_sample = rebin(center_r[round(Line_centers[j])] + dawn_distance_array_r, n_elements(dawn_distance_array_r), sample_width)
          dawn_line_profiles[*,j] = total(image[dawn_x_sample, dawn_y_sample], 2) 
          err_dawn_line_profiles[*,j] = sqrt(total((err_image[dawn_x_sample, dawn_y_sample])^2., 2, /NAN))
          dummy[dawn_x_sample, dawn_y_sample] = !values.F_NaN & err_dummy[dawn_x_sample, dawn_y_sample] = !values.F_NaN ;never accidentally count these pixels again         
                
          ;fit the lines in individual frames, more accurate than summing if done right.
          parinfo[2].value = linewidths_r[j] & guess[2] = linewidths_r[j] ;use fixed linewidths interpolated from the arc frames 
          origin = round(Line_centers[j])
          for k = 0, n_elements(dawn_distance_array_r) - 1 do begin
            fitme = image[origin - sample_width + (Lines_r[j]*dawn_shift[k]/dispers_r) : origin + sample_width + (Lines_r[j]*dawn_shift[k]/dispers_r), fix(center_r[origin] + k)]
            err_fitme = err_image[origin - sample_width + (Lines_r[j]*dawn_shift[k]/dispers_r) : origin + sample_width + (Lines_r[j]*dawn_shift[k]/dispers_r), fix(center_r[origin] + k)]
            fit = mpfitpeak(findgen(n_elements(fitme)), fitme, A, PERROR = err_A, /POSITIVE, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_fitme, ESTIMATES = guess, $
              bestNORM = BESTNORM, DOF = DOF, COVAR = COVAR, /NAN)
            Fit_dawn_line_profiles[k,j] = A[0]*A[2]*SQRT(2*!DPI)  
            ;err_Fit_dawn_line_profiles[k,j] = Fit_dawn_line_profiles[k,j] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2. + abs(2.*A[2]*A[0]*covar[0,2]) ) ; uncertainties in width and height add in quadrature.
            err_Fit_dawn_line_profiles[k,j] = Fit_dawn_line_profiles[k,j] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2. ) ; uncertainties in width and height add in quadrature.
          endfor           
        endfor     
        for j = 0, N_elements(lines_r)-1 do begin ;Off band loop
          dawn_x_sample_offband = rebin(round(Line_centers[j] + dawn_shift * lines_r[j] / dispers_r), n_elements(dawn_distance_array_r), 4.*sample_width) + $
                                  rebin(transpose(indgen(4.*sample_width) - 2.*sample_width), n_elements(dawn_distance_array_r), 4.*sample_width)
          dawn_y_sample_offband = rebin(center_r[round(Line_centers[j])] + dawn_distance_array_r, n_elements(dawn_distance_array_r), 4.*sample_width)
          dawn_line_profiles_offband[*,j] = median(dummy[dawn_x_sample_offband, dawn_y_sample_offband], DIMENSION = 2, /EVEN) ;NaN's are non included to this is the median of adhacent spectral columns
          finite_elements = round(4.*sample_width * total(finite(err_dummy[dawn_x_sample_offband, dawn_y_sample_offband])) / N_elements(err_dummy[dawn_x_sample_offband, dawn_y_sample_offband]))
          err_dawn_line_profiles_offband[*,j] = sqrt(total((err_dummy[dawn_x_sample_offband, dawn_y_sample_offband])^2., 2, /NAN)) / float(finite_elements) ;average error per pixel along the dispersion axis
        endfor
        dawn_line_profiles = dawn_line_profiles - (dawn_line_profiles_offband * sample_width) 
        err_dawn_line_profiles = sqrt(err_dawn_line_profiles^2. + (err_dawn_line_profiles_offband * sample_width)^2.) 
        
        ;account for that partial pixel, consider the pixel center as the location where the brightness values apply  
        dawn_distance_array_r = rebin(dawn_distance_array_r, N_elements(dawn_distance_array_r), N_elements(lines_r))
        for j = 0, N_elements(lines_r)-1 do begin
          dawn_distance_array_r[*,j] = dawn_distance_array_r[*,j] + (.5 - pixel_fraction[j])
        endfor  
        
        dawn_distance_array_r = dawn_distance_array_r * R_J_per_pixel_r ;convert pixel distance array to an R_J array 
        max_dawn_r = max_dawn_r > N_elements(dawn_distance_array_r)

        ;--------------------dusk---------------------------        
        dusk_distance_array_r = findgen(distance_extent_dusk / R_J_per_pixel_r)
        dusk_distance_array_r_km = R_J * R_J_per_pixel_r * dusk_distance_array_r
        dusk_shift =  ((2. * !pi * dusk_distance_array_r_km / period_J) / c)  + (Jupiter_WRT_Earth / c)   ; just v over c, multiply by lambda for the wavelength shift
        dusk_line_profiles = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) & err_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) 
        dusk_line_profiles_offband = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) & err_dusk_line_profiles_offband = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) 
        fit_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) & err_fit_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_r), N_elements(lines_r)) 
        for j = 0, N_elements(lines_r)-1 do begin
          dusk_x_sample = rebin(round(Line_centers[j] + dusk_shift * lines_r[j] / dispers_r), n_elements(dusk_distance_array_r), sample_width) + $
                          rebin(transpose(indgen(sample_width) - fix(sample_width/2)), n_elements(dusk_distance_array_r), sample_width)
          dusk_y_sample = rebin(center_r[(Line_centers[j])] - dusk_distance_array_r, n_elements(dusk_distance_array_r), sample_width)
          dusk_line_profiles[*,j] = total(image[dusk_x_sample, dusk_y_sample], 2) 
          err_dusk_line_profiles[*,j] = sqrt(total((err_image[dusk_x_sample, dusk_y_sample])^2., 2, /NaN))
          dummy[dusk_x_sample, dusk_y_sample] = !values.F_NaN 
          parinfo[2].value = linewidths_r[j] & guess[2] = linewidths_r[j]
        
          ;fit the lines in individual frames
          origin = round(Line_centers[j])          
          for k = 0, n_elements(dusk_distance_array_r) - 1 do begin
            fitme = image[origin - sample_width + (Lines_r[j]*dusk_shift[k]/dispers_r) : origin + sample_width + (Lines_r[j]*dusk_shift[k]/dispers_r), fix(center_r[origin] - k)]
            err_fitme = err_image[origin - sample_width + (Lines_r[j]*dusk_shift[k]/dispers_r) : origin + sample_width + (Lines_r[j]*dusk_shift[k]/dispers_r), fix(center_r[origin] - k)]
            fit = mpfitpeak(findgen(n_elements(fitme)), fitme, A, PERROR = err_A, /POSITIVE, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_fitme, ESTIMATES = guess, $
              bestNORM = BESTNORM, DOF = DOF, COVAR = COVAR, /NAN)
            Fit_dusk_line_profiles[k,j] = A[0]*A[2]*SQRT(2*!DPI)  
            ;err_Fit_dusk_line_profiles[k,j] = Fit_dusk_line_profiles[k,j] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2. + abs(2.*A[2]*A[0]*covar[0,2]) ) ; uncertainties in width and height add in quadrature.                
            err_Fit_dusk_line_profiles[k,j] = Fit_dusk_line_profiles[k,j] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2. ) ; uncertainties in width and height add in quadrature. 
          endfor    
        endfor
        for j = 0, N_elements(lines_r)-1 do begin ;Off band loop for the straight sum portion, not used in line fitting
          dusk_x_sample_offband = rebin(round(Line_centers[j] + dusk_shift * lines_r[j] / dispers_r), n_elements(dusk_distance_array_r), 4.*sample_width) + $
                                  rebin(transpose(indgen(4.*sample_width) - 2.*sample_width), n_elements(dusk_distance_array_r), 4.*sample_width)
          dusk_y_sample_offband = rebin(center_r[round(Line_centers[j])] + dusk_distance_array_r, n_elements(dusk_distance_array_r), 4.*sample_width)
          dusk_line_profiles_offband[*,j] = median(dummy[dusk_x_sample_offband, dusk_y_sample_offband], DIMENSION = 2 , /EVEN)
          finite_elements = round(4.*sample_width * total(finite(err_dummy[dusk_x_sample_offband, dusk_y_sample_offband])) / N_elements(err_dummy[dusk_x_sample_offband, dusk_y_sample_offband]))
          err_dusk_line_profiles_offband[*,j] = sqrt(total((err_dummy[dusk_x_sample_offband, dusk_y_sample_offband])^2., 2, /NAN)) / float(finite_elements) ;average error per pixel along the dispersion axis
        endfor
        dusk_line_profiles = dusk_line_profiles - (dusk_line_profiles_offband * sample_width)
        err_dusk_line_profiles = sqrt(err_dusk_line_profiles^2. + (err_dusk_line_profiles_offband * sample_width)^2.) 

        ;account for that partial pixel, consider the pixel center as the location where the brightness values apply  
        dusk_distance_array_r = rebin(dusk_distance_array_r, N_elements(dusk_distance_array_r), N_elements(lines_r))
        for j = 0, N_elements(lines_r)-1 do begin
          dusk_distance_array_r[*,j] = dusk_distance_array_r[*,j] + (pixel_fraction[j] - .5)
        endfor  

        dusk_distance_array_r = dusk_distance_array_r * R_J_per_pixel_r ;convert pixel distance array to an R_J array    
        max_dusk_r = max_dusk_r > N_elements(dusk_distance_array_r)
        if keyword_set(debug) then begin   
          Window, 2, xs = 2098, ys = 1078
          tv, bytscl(dummy, 0, 50)        
        endif

        distance_profiles_structure_r = {dusk_distance_array_r:dusk_distance_array_r, dusk_line_profiles:dusk_line_profiles, err_dusk_line_profiles:err_dusk_line_profiles, $
                                         fit_dusk_line_profiles:fit_dusk_line_profiles, err_fit_dusk_line_profiles:err_fit_dusk_line_profiles, $
                                         dawn_distance_array_r:dawn_distance_array_r, dawn_line_profiles:dawn_line_profiles, err_dawn_line_profiles:err_dawn_line_profiles, $
                                         fit_dawn_line_profiles:fit_dawn_line_profiles, err_fit_dawn_line_profiles:err_fit_dawn_line_profiles}
        MWRFITS, image, Images_r[i], torus_header_r, /CREATE, /silent ;/create overwrites
        MWRFITS, err_image, Images_r[i], /silent ;Append the fits file -> extension 1  
        MWRFITS, metadata_r,Images_r[i], /silent ;Append the fits file -> extension 2      
        MWRFITS, distance_profiles_structure_r, Images_r[i], /silent ;Append the fits file -> extension 3
    endfor

    ;-------------------Blue Channel----------------------------  
    max_dawn_b = 0 & max_dusk_b = 0 ;the furthest pixels with any data of all the images considered   
    for i = 0, n_elements(Images_b)-1 do begin
   
      ;Calculate the maximum distance in RJ where any data is available
        image = MRDFITS(Images_b[i], 0, torus_header_b, /silent)
        err_image = MRDFITS(Images_b[i], 1, /silent)
        dummy = image & err_dummy = err_image
        metadata_b = MRDFITS(Images_b[i], 2, /silent)
        center_b = metadata_b.center
        Dispersion_mask_b = metadata_b.dispersion_mask 
        dispers_b = float(sxpar(torus_header_b, 'DISPERS')) 
        R_J_per_pixel_b = float(sxpar(torus_header_b, 'RJ_PER_P')) 
        dispers = float(sxpar(torus_header_b, 'DISPERS')) 
        date = strmid(sxpar(torus_header_b, 'DATE-OBS'),0,10)
        aperture_coords = ARRAY_INDICES(image , aperture_b)  
        distance_extent_dawn = abs(max(aperture_coords[1,*]) - center_b[1024]) * R_J_per_pixel_b - 0.2 ;Distance in RJ out to which torus data is present 
        distance_extent_dusk = abs(min(aperture_coords[1,*]) - center_b[1024]) * R_J_per_pixel_b - 0.2 ;Distance in RJ out to which torus data is present 
        linewidths_b = poly(lines_b, float([sxpar(torus_header_b, 'W0_COEF'), sxpar(torus_header_b, 'W1_COEF')]))  ; the sigma linewidth of emission lines in the blue channel
        parinfo_Blue_doublet[3].value = mean(linewidths_b[3:4])              ;Set sigma based on the width of lines measured in arclamp frames, ignore negligible line broadening
        parinfo_blue_triplet[4].value = mean(linewidths_b[0:2])              ;Set sigma based on the width of lines measured in arclamp frames, ignore negligible line broadening 
        initial_guess_Blue_triplet = [0., 0., 0., sample_width_Blue_triplet, mean(linewidths_b[0:2])]
        initial_guess_Blue_doublet = [0., 0.,     sample_width_Blue_doublet, mean(linewidths_b[3:4])]

      ;Calculate Jupiter's radial velocity WRT Earth (Earth center not observatory)
        torus_time = double(sxpar(torus_header_b, 'EPHEM_TI')) 
        cspice_spkezr, 'Jupiter', torus_time, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_state, Jupiter_Earth_light_time ;Jupiter's J2000 position with respect to the Earth
        theta  = cspice_vsep(Jupiter_Earth_state[0:2], Jupiter_Earth_state[3:5])
        Jupiter_WRT_Earth = cos(theta) * sqrt(Jupiter_Earth_state[3]^2.+Jupiter_Earth_state[4]^2.+Jupiter_Earth_state[5]^2.) ;scalar projection of the relative velocity along the line of sight

      ;Calculate where the measured line centers are. Use the fit centroid position wavelength where the ribbon should be. For all lines         
        CML = float(sxpar(torus_header_b, 'CML_SYSI'))
        CML_Dawn = (CML + 90. + 360.0) MOD 360.0 
        CML_Dusk = (CML - 90. + 360.0) MOD 360.0 
        ribbon_range_dawn = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
        ribbon_range_dusk = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
        expected_ribbon_dawn = ribbon_range_dawn / R_J_per_pixel_b
        expected_ribbon_dusk = ribbon_range_dusk / R_J_per_pixel_b
        Line_centers = fltarr(N_elements(Lines_b)) 
        if keyword_set(debug) then window, 0, xs = 600, ys = 400
        dawn_shift = (-(2. * !pi * R_J * ribbon_range_dawn * Lines_b/ period_J) / (c * dispers_b)) +  (Lines_b * Jupiter_WRT_Earth) / (c * dispers_b) ;doppler shift at the dawn ribbon for each emission line
        dusk_shift = ( (2. * !pi * R_J * ribbon_range_dusk * lines_b/ period_J) / (c * dispers_b)) +  (lines_b * Jupiter_WRT_Earth) / (c * dispers_b) ;doppler shift at the dusk ribbon for each emission line
        for n = 0, N_elements(Lines_b) - 1 do begin                  
            junk = min(abs(Dispersion_mask_b[indgen(2098), center_b] - Lines_b[n]), /NAN, origin) ;compute the line center pixel from the dispersion mask
            dawn_fit = total(Image[origin - sample_width + dawn_shift[n] : origin + sample_width + dawn_shift[n], $ ;sum 30 pixels over the ribbon
                             center_b[[origin]] + expected_ribbon_dawn - 15: center_b[[origin]] + expected_ribbon_dawn + 15], 2) 
            dusk_fit = total(Image[origin - sample_width + dusk_shift[n] : origin + sample_width + dusk_shift[n], $ ;sum 30 pixels over the ribbon
                             center_b[[origin]] - expected_ribbon_dusk - 15: center_b[[origin]] - expected_ribbon_dusk + 15], 2) 
            err_dawn_fit = sqrt(total(err_Image[origin - sample_width + dawn_shift[n] : origin + sample_width + dawn_shift[n], $ ;sum 30 pixels over the ribbon
                                 center_b[[origin]] + expected_ribbon_dawn - 15: center_b[[origin]] + expected_ribbon_dawn + 15]^2, 2)) 
            err_dusk_fit = sqrt(total(err_Image[origin - sample_width + dusk_shift[n] : origin + sample_width + dusk_shift[n], $ ;sum 30 pixels over the ribbon
                                 center_b[[origin]] - expected_ribbon_dusk - 15: center_b[[origin]] - expected_ribbon_dusk + 15]^2, 2)) 
            parinfo[2].value = linewidths_b[n] & guess[2] = linewidths_b[n]
            fit_dawn = mpfitpeak(findgen(n_elements(dawn_fit)), dawn_fit, A_dawn, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dawn_fit, ESTIMATES = guess, /NAN)
            fit_dusk = mpfitpeak(findgen(n_elements(dusk_fit)), dusk_fit, A_dusk, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dusk_fit, ESTIMATES = guess, /NAN)
      
          ;Refine line centers based on the measured ribbon pixel positions, assuming corotation 
          Line_centers[n] = origin; + (mean([a_dawn[1],a_dusk[1]], /Nan) - sample_width) 
                 
          if keyword_set(debug) then begin              
            if n eq 0 then window, 0, xs = 600, ys = 400, Title = 'Dawn Blue Line Positions'
            wset, 0
            if n eq 0 then cgplot, dawn_fit, psym = 4, color = blue_colors[n], yr = [-100, 2000] else cgplot, dawn_fit, psym = 4, color = blue_colors[n], /overplot     
            cgplot, fit_dawn, color = blue_colors[n], /overplot
            
            if n eq 0 then window, 1, xs = 600, ys = 400, Title = 'Dusk Blue Line Positions'
            wset, 1
            if n eq 0 then cgplot, dusk_fit, psym = 4, color = blue_colors[n], yr = [-100, 2000] else cgplot, dusk_fit, psym = 4, color = blue_colors[n], /overplot     
            cgplot, fit_dusk, color = blue_colors[n], /overplot
          endif
        endfor
        
        if keyword_set(debug) then begin ;check line centers!       
          range = indgen(2098)
          y = total(Image[ min(range):max(range), mean(center_b[range]) + expected_ribbon_dawn - 15: mean(center_b[range]) + expected_ribbon_dawn + 15 ], 2, /NAN) 
          x = Dispersion_mask_b[range, mean(center_b[range])]
          wset, 0
          cgplot, x, y, xr = [3705, 3750], psym = 10, yr = [-100, 1100], ystyle = 1
          cgplot, Dispersion_mask_b[round(Line_centers + dawn_shift), center_b[round(Line_centers + dawn_shift)]], make_array(N_elements(Line_centers), /float, value = 1000.), psym =4, /overplot
          wset, 1
          cgplot, x, y, xr = [4060, 4080], psym = 10, yr = [-100, 1100], ystyle = 1
          cgplot, Dispersion_mask_b[round(Line_centers + dawn_shift), center_b[round(Line_centers + dawn_shift)]], make_array(N_elements(Line_centers), /float, value = 1000.), psym =4, /overplot
        endif 
        pixel_fraction = center_b[round(Line_centers)] - fix(center_b[round(Line_centers)]) ; the fractional pixel distance of the centroid

        ;--------------------dawn---------------------------        
        dawn_distance_array_b = findgen(distance_extent_dawn / R_J_per_pixel_b)
        dawn_distance_array_b_km = R_J * R_J_per_pixel_b * dawn_distance_array_b
        dawn_shift = -((2. * !pi * dawn_distance_array_b_km / period_J) / c)  + (Jupiter_WRT_Earth / c); just v over c, multiply by lambda for the wavelength shift
        dawn_line_profiles = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b)) & err_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b)) 
        dawn_line_profiles_offband = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b)) & err_dawn_line_profiles_offband = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b)) 
        fit_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b)) & err_fit_dawn_line_profiles = fltarr(N_elements(dawn_distance_array_b), N_elements(lines_b))
        for j = 0, N_elements(lines_b)-1 do begin
          dawn_x_sample = rebin(round(Line_centers[j] + dawn_shift * lines_b[j] / dispers_b), n_elements(dawn_distance_array_b), sample_width) + rebin(transpose(indgen(sample_width) - fix(sample_width/2)), n_elements(dawn_distance_array_b), sample_width)
          dawn_y_sample = rebin(center_b[round(Line_centers[j])] + dawn_distance_array_b, n_elements(dawn_distance_array_b), sample_width)
          dawn_line_profiles[*,j] = total(image[dawn_x_sample, dawn_y_sample], 2) 
          err_dawn_line_profiles[*,j] = sqrt(total((err_image[dawn_x_sample, dawn_y_sample])^2., 2, /NaN)) 
          parinfo[2].value = linewidths_b[j] & guess[2] = linewidths_b[j]
          dummy[dawn_x_sample, dawn_y_sample] = !values.F_NaN 
        endfor
        for j = 0, N_elements(lines_b)-1 do begin ;Off band loop
          dawn_x_sample_offband = rebin(round(Line_centers[j] + dawn_shift * lines_b[j] / dispers_b), n_elements(dawn_distance_array_b), 4.*sample_width) + rebin(transpose(indgen(4.*sample_width) - 2.*sample_width), n_elements(dawn_distance_array_b), 4.*sample_width)
          dawn_y_sample_offband = rebin(center_b[round(Line_centers[j])] + dawn_distance_array_b, n_elements(dawn_distance_array_b), 4.*sample_width)            
          dawn_line_profiles_offband[*,j] = median(dummy[dawn_x_sample_offband, dawn_y_sample_offband], DIMENSION = 2 , /EVEN)
          finite_elements = round(4.*sample_width * total(finite(err_dummy[dawn_x_sample_offband, dawn_y_sample_offband])) / N_elements(err_dummy[dawn_x_sample_offband, dawn_y_sample_offband]))
          err_dawn_line_profiles_offband[*,j] = sqrt(total((err_dummy[dawn_x_sample_offband, dawn_y_sample_offband])^2., 2, /NAN)) / float(finite_elements) ;average error per pixel along the dispersion axis
        endfor
        dawn_line_profiles = dawn_line_profiles - (dawn_line_profiles_offband * sample_width)    
        err_dawn_line_profiles = sqrt(err_dawn_line_profiles^2. + (err_dawn_line_profiles_offband * sample_width)^2.)     
        
        for k = 0, n_elements(dawn_distance_array_b) - 1 do begin
          junk = min(abs(Dispersion_mask_b[indgen(2098), center_b] - Lines_b[1]), /NAN, origin) ;compute the line center pixel from the dispersion mask
          y = image[origin - sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[k]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[k]/dispers_b), center_b[origin] + k]
          err_y = err_image[origin - sample_width_Blue_triplet  + (Lines_b[1]*dawn_shift[k]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[k]/dispers_b), center_b[origin] + k]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_triplet', initial_guess_Blue_triplet, PERROR = err_p, funct=fa, STATUS=STATUS, parinfo = parinfo_Blue_triplet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan & err_P = make_array(5, value = !values.f_Nan)
          Fit_dawn_line_profiles[k,0] = P[0]*P[4]*SQRT(2*!DPI) 
          Fit_dawn_line_profiles[k,1] = P[1]*P[4]*SQRT(2*!DPI) 
          Fit_dawn_line_profiles[k,2] = P[2]*P[4]*SQRT(2*!DPI) 
          err_Fit_dawn_line_profiles[k,0] = P[0]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dawn_line_profiles[k,1] = P[1]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dawn_line_profiles[k,2] = P[2]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[2]/P[2])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          
          junk = min(abs(Dispersion_mask_b[indgen(2098), center_b] - Lines_b[3]), /NAN, origin) ;compute the line center pixel from the dispersion mask  
          y = image[origin - sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[k]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[k]/dispers_b), center_b[origin] + k]
          err_y = err_image[origin - sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[k]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[k]/dispers_b), center_b[origin] + k]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_doublet', initial_guess_Blue_doublet, PERROR = err_p, funct=fa, STATUS=STATUS, parinfo = parinfo_Blue_doublet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan & err_P = make_array(4, value = !values.f_Nan)
          Fit_dawn_line_profiles[k,3] = P[0]*P[3]*SQRT(2*!DPI) 
          Fit_dawn_line_profiles[k,4] = P[1]*P[3]*SQRT(2*!DPI) 
          err_Fit_dawn_line_profiles[k,3] = P[0]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[3]/P[3])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dawn_line_profiles[k,4] = P[1]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[3]/P[3])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
        endfor   

        ;account for that partial pixel, consider the pixel center as the location where the brightness values apply  
        ;convert pixel distance array to an R_J array
        dawn_distance_array_b = rebin(dawn_distance_array_b, N_elements(dawn_distance_array_b), N_elements(lines_b))
        for j = 0, N_elements(lines_b)-1 do begin
          dawn_distance_array_b[*,j] = (dawn_distance_array_b[*,j] + (.5 - pixel_fraction[j])) * R_J_per_pixel_b 
        endfor 
        max_dawn_b = max_dawn_b > N_elements(dawn_distance_array_b)
 
        ;--------------------dusk---------------------------
        dusk_distance_array_b = findgen(distance_extent_dusk / R_J_per_pixel_b)
        dusk_distance_array_b_km = R_J * R_J_per_pixel_b * dusk_distance_array_b
        dusk_shift = ((2. * !pi * dusk_distance_array_b_km / period_J) / c)  + (Jupiter_WRT_Earth / c)   ; just v over c, multiply by lambda for the wavelength shift
        dusk_line_profiles = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b)) & err_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b))
        dusk_line_profiles_offband = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b)) & err_dusk_line_profiles_offband = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b)) 
        fit_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b)) & err_fit_dusk_line_profiles = fltarr(N_elements(dusk_distance_array_b), N_elements(lines_b)) 
        for j = 0, N_elements(lines_b)-1 do begin
          dusk_x_sample = rebin(round(Line_centers[j] + dusk_shift * lines_b[j] / dispers_b), n_elements(dusk_distance_array_b), sample_width) + $
                          rebin(transpose(indgen(sample_width) - fix(sample_width/2)), n_elements(dusk_distance_array_b), sample_width)
          dusk_y_sample = rebin(center_b[round(Line_centers[j])] - dusk_distance_array_b, n_elements(dusk_distance_array_b), sample_width)
          dusk_line_profiles[*,j] = total(image[dusk_x_sample, dusk_y_sample], 2) 
          err_dusk_line_profiles[*,j] = sqrt(total((err_image[dusk_x_sample, dusk_y_sample])^2., 2, /NaN))
          parinfo[2].value = linewidths_b[j] & guess[2] = linewidths_b[j]
          dummy[dusk_x_sample, dusk_y_sample] = !values.F_NaN 
        endfor
        for j = 0, N_elements(lines_b)-1 do begin ;Off band loop
          dusk_x_sample_offband = rebin(round(Line_centers[j] + dusk_shift * lines_b[j] / dispers_b), n_elements(dusk_distance_array_b), 4.*sample_width) + $
                                  rebin(transpose(indgen(4.*sample_width) - 2.*sample_width), n_elements(dusk_distance_array_b), 4.*sample_width)
          dusk_y_sample_offband = rebin(center_b[round(Line_centers[j])] + dusk_distance_array_b, n_elements(dusk_distance_array_b), 4.*sample_width)
          dusk_line_profiles_offband[*,j] = median(dummy[dusk_x_sample_offband, dusk_y_sample_offband], DIMENSION = 2 , /EVEN) ;this is the average offband per pixel
          finite_elements = round(4.*sample_width * total(finite(err_dummy[dusk_x_sample_offband, dusk_y_sample_offband])) / N_elements(err_dummy[dusk_x_sample_offband, dusk_y_sample_offband]))
          err_dusk_line_profiles_offband[*,j] = sqrt(total((err_dummy[dusk_x_sample_offband, dusk_y_sample_offband])^2., 2, /NAN)) / float(finite_elements) ;average error per pixel along the dispersion axis
        endfor
        dusk_line_profiles = dusk_line_profiles - (dusk_line_profiles_offband * sample_width)  
        err_dusk_line_profiles = sqrt(err_dusk_line_profiles^2. + (err_dusk_line_profiles_offband * sample_width)^2.) 
        
        for k = 0, n_elements(dusk_distance_array_b) - 1 do begin
          junk = min(abs(Dispersion_mask_b[indgen(2098), center_b] - Lines_b[1]), /NAN, origin) ;compute the line center pixel from the dispersion mask
          y = image[origin - sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[k]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[k]/dispers_b), center_b[origin] - k]
          err_y = err_image[origin - sample_width_Blue_triplet  + (Lines_b[1]*dusk_shift[k]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[k]/dispers_b), center_b[origin] - k]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_triplet', initial_guess_Blue_triplet, PERROR = err_p, funct=fa, STATUS=STATUS, parinfo = parinfo_Blue_triplet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan & err_P = make_array(5, value = !values.f_Nan)
          Fit_dusk_line_profiles[k,0] = P[0]*P[4]*SQRT(2*!DPI) 
          Fit_dusk_line_profiles[k,1] = P[1]*P[4]*SQRT(2*!DPI) 
          Fit_dusk_line_profiles[k,2] = P[2]*P[4]*SQRT(2*!DPI) 
          err_Fit_dusk_line_profiles[k,0] = P[0]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dusk_line_profiles[k,1] = P[1]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dusk_line_profiles[k,2] = P[2]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[2]/P[2])^2. + (err_P[4]/P[4])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          
          junk = min(abs(Dispersion_mask_b[indgen(2098), center_b] - Lines_b[3]), /NAN, origin) ;compute the line center pixel from the dispersion mask  
          y = image[origin - sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[k]/dispers) : origin + sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[k]/dispers), center_b[origin] - k]
          err_y = err_image[origin - sample_width_Blue_doublet  + (Lines_b[3]*dusk_shift[k]/dispers) : origin + sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[k]/dispers), center_b[origin] - k]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_doublet', initial_guess_Blue_doublet, PERROR = err_p, funct=fa, STATUS=STATUS, parinfo = parinfo_Blue_doublet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan & err_P = make_array(4, value = !values.f_Nan)
          Fit_dusk_line_profiles[k,3] = P[0]*P[3]*SQRT(2*!DPI) 
          Fit_dusk_line_profiles[k,4] = P[1]*P[3]*SQRT(2*!DPI) 
          err_Fit_dusk_line_profiles[k,3] = P[0]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[3]/P[3])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
          err_Fit_dusk_line_profiles[k,4] = P[1]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[3]/P[3])^2.) ; uncertainties in width and height add in quadrature, neglect covariance.  
        endfor
        
        ;account for that partial pixel, consider the pixel center as the location where the brightness values apply  
        ;convert pixel distance array to an R_J array    
        dusk_distance_array_b = rebin(dusk_distance_array_b, N_elements(dusk_distance_array_b), N_elements(lines_b))
        for j = 0, N_elements(lines_b)-1 do begin
          dusk_distance_array_b[*,j] = (dusk_distance_array_b[*,j] + (pixel_fraction[j] - .5)) * R_J_per_pixel_b   
        endfor        
        max_dusk_b = max_dusk_b > N_elements(dusk_distance_array_b)        
        
        if keyword_set(debug) then begin     
          Window, 2, xs = 2098, ys = 1078, title = 'frame' + string(i)
          tv, bytscl(dummy, 0, 50)
          ;wait, 5.
        endif

        distance_profiles_structure_b = {dusk_distance_array_b:dusk_distance_array_b, dusk_line_profiles:dusk_line_profiles, err_dusk_line_profiles:err_dusk_line_profiles, $
                                         fit_dusk_line_profiles:fit_dusk_line_profiles, err_fit_dusk_line_profiles:err_fit_dusk_line_profiles, $
                                         dawn_distance_array_b:dawn_distance_array_b, dawn_line_profiles:dawn_line_profiles, err_dawn_line_profiles:err_dawn_line_profiles, $
                                         fit_dawn_line_profiles:fit_dawn_line_profiles, err_fit_dawn_line_profiles:err_fit_dawn_line_profiles}
        MWRFITS, image, Images_b[i], torus_header_b, /CREATE, /silent ;/create overwrites
        MWRFITS, err_image, Images_b[i], /silent  ;Append the fits file -> extension 1  
        MWRFITS, metadata_b, Images_b[i], /silent ;Append the fits file -> extension 2       
        MWRFITS, distance_profiles_structure_b, Images_b[i], /silent ;Append the fits file -> extension 3
    endfor

    ;-----------------------------------------------co_adding-------------------------------------------------------------
    ;
    ;     Shift, stack and average calibrated torus frames, based on Jupiter's Centroid and the measured line positions.
    ;
    ;------------------------------------------------------------------------------------------------------------------------
    
   img_b = MRDFITS(Images_b[key], 0, torus_header_b, /silent)
   metadata_b = MRDFITS(Images_b[key], 2, /silent)
   distance_profiles_b = MRDFITS(Images_b[key], 3, /silent)
   Dispersion_mask_b_key = metadata_b.Dispersion_mask
   center_b = metadata_b.center
   key_b = mean(metadata_b.center)
   dispers_b = float(sxpar(torus_header_b, 'DISPERS')) 
   
   img_r = MRDFITS(Images_r[key], 0, torus_header_r, /silent)
   metadata_r = MRDFITS(Images_r[key], 2, /silent)
   distance_profiles_r = MRDFITS(Images_r[key], 3, /silent)
   Dispersion_mask_r_Key = metadata_r.Dispersion_mask
   center_r = metadata_r.center
   key_r = mean(metadata_r.center)
   dispers_r = float(sxpar(torus_header_r, 'DISPERS')) 

   for j = 0, N_elements(CoAdd_frames[0,*]) - 1 do begin
    
      ;----------------------Tabulate Metadata---------------------
      CML = fltarr( n_elements(CoAdd_frames[*,j]) ) & CML_Dusk = fltarr( n_elements(CoAdd_frames[*,j]) ) & CML_Dawn = fltarr( n_elements(CoAdd_frames[*,j]) )
      observing_times = strarr( n_elements(CoAdd_frames[*,j]) ) 
      R_J_per_pixel_r = fltarr( n_elements(CoAdd_frames[*,j]) ) 
      R_J_per_pixel_b = fltarr( n_elements(CoAdd_frames[*,j]) )
      for i = 0, N_elements(CoAdd_frames[*,j]) - 1 do begin 
        torus_header_r = headfits(Images_r[CoAdd_frames[i,j]])
        torus_header_b = headfits(Images_b[CoAdd_frames[i,j]])
        R_J_per_pixel_r[i] = float(sxpar(torus_header_r, 'RJ_PER_P')) 
        R_J_per_pixel_b[i] = float(sxpar(torus_header_b, 'RJ_PER_P')) 
        CML[i] = float(sxpar(torus_header_r, 'CML_SYSI'))
        observing_times[i] = sxpar(torus_header_r, 'DATE-OBS')
        CML_Dawn[i] = (CML[i] + 90. + 360.0) MOD 360.0 
        CML_Dusk[i] = (CML[i] - 90. + 360.0) MOD 360.0   
      endfor
      R_J_per_pixel_r = mean(R_J_per_pixel_r) ;average platescale over the co-added images
      R_J_per_pixel_b = mean(R_J_per_pixel_b) ;average platescale over the co-added images
   
      ;get the apertures as to not include off aperture data in the combined image, under the mask is included in the final image
      blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 0, flat_header_b, /Dscale, /silent )
      red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 0, flat_header_r, /Dscale, /silent )
      err_blue_flat = MRDFITS(directory+'Processed\Blue_Master_flat.fits', 1, err_flat_header_b, /Dscale, /silent )
      err_red_flat = MRDFITS(directory+'Processed\Red_Master_flat.fits', 1, err_flat_header_r, /Dscale, /silent )   
      aperture_b = where(blue_flat gt mean(blue_flat)) 
        dummy = blue_flat       
        y_aperture_ind = where( total(blue_flat,1) gt mean(total(blue_flat,1)), /NULL )
        y_aperture_ind = indgen(max(y_aperture_ind) - min(y_aperture_ind ) + 1) + min(y_aperture_ind ) ;fill in the mask so that Jupiter doesn't disappear
        ind = array_indices(blue_flat, aperture_b) 
        dummy[0:max(ind[0,*]),y_aperture_ind] = -100
        aperture_b = where(dummy eq -100, /NULL, complement = anti_aperture_blue) ;a much better aperture.
      aperture_r = where(red_flat gt mean(red_flat)) 
        dummy = red_flat       
        y_aperture_ind = where( total(red_flat,1) gt mean(total(red_flat,1)), /NULL )
        y_aperture_ind = indgen(max(y_aperture_ind) - min(y_aperture_ind ) + 1) + min(y_aperture_ind ) ;fill in the mask so that Jupiter doesn't disappear
        ind = array_indices(red_flat, aperture_r) 
        dummy[0:max(ind[0,*]),y_aperture_ind] = -100
        aperture_r = where(dummy eq -100, /NULL, complement = anti_aperture_red) ;a much better aperture.
 
      ;-------------Combine Red Images for Equivalent Widths----------------------------
      bigarray = fltarr(2098, 1078, N_elements(CoAdd_frames[*,j]))
      err_bigarray = fltarr(2098, 1078, N_elements(CoAdd_frames[*,j]))
      for i = 0, N_elements(CoAdd_frames[*,j]) - 1 do begin
        img = MRDFITS(Images_r[CoAdd_frames[i,j]], 0, /silent)
        err_img = MRDFITS(Images_r[CoAdd_frames[i,j]], 1, /silent)
        metadata_r = MRDFITS(Images_r[CoAdd_frames[i,j]], 2, /silent) 
        Dispersion_mask_r = metadata_r.Dispersion_mask ;HACK, DOUBLE CHECK NEEDED that the DISPERSION MASK IS ALWAYS THE SAME over consecutive frames, BLUE ALSO!!!
        center = metadata_r.center
        img[anti_aperture_red] = !values.f_nan ;don't include points not in the slit aperture
        err_img[anti_aperture_red] = !values.f_nan ;don't include points not in the slit aperture
        bigarray[*,*,i] = smart_shift(img, 0, key_r - mean(metadata_r.center), /interp, Missing = !values.f_nan)
        err_bigarray[*,*,i] = smart_shift(err_img, 0, key_r - mean(metadata_r.center), /interp, Missing = !values.f_nan)
      endfor  
     
      ;replace outliers in the stack, careful not to replace the torus instrinsic variability
        sig = STDDEV(bigarray, DIMENSION=3, /DOUBLE, /NAN)
        stack_mean = MEAN(bigarray, DIMENSION=3, /double, /NAN) 
        deviation = bigarray - rebin(stack_mean, 2098, 1078, N_elements(CoAdd_frames[*,j]))
        bigarray[where(abs(deviation) gt 2.75*sig, /NULL)] = !values.F_NAN
        err_bigarray[where(abs(deviation) gt 2.75*sig, /NULL)] = !values.F_NAN
        pretty_r = MEAN(bigarray, DIMENSION=3, /double, /NAN) 
        err_pretty_r = sqrt(total(err_bigarray^2, 3, /NAN)) / N_elements(CoAdd_frames[*,j])
        if keyword_set(debug) then begin
          window, 0, xs = 2098, ys =1000
          tv, bytscl(stack_mean - pretty_r, -2, 2)
        endif

      dawn_dist_r_km = distance_profiles_r.dawn_distance_array_r * R_J
      dusk_dist_r_km = distance_profiles_r.dusk_distance_array_r * R_J
      dawn_shift = (-2.*!pi * dawn_dist_r_km / period_J) / c  ;V/C
      dusk_shift = ( 2.*!pi * dusk_dist_r_km / period_J) / c  ;V/C  
      dawn_avg_r = fltarr(N_elements(dawn_dist_r_km[*,0]), N_elements(Lines_r)) & err_dawn_avg_r = fltarr(N_elements(dawn_dist_r_km[*,0]), N_elements(Lines_r))
      dusk_avg_r = fltarr(N_elements(dusk_dist_r_km[*,0]), N_elements(Lines_r)) & err_dusk_avg_r = fltarr(N_elements(dusk_dist_r_km[*,0]), N_elements(Lines_r))
          
      ;----------Fine tune line center's-----------------
          expected_ribbon_dawn = 5.85 / R_J_per_pixel_r
          expected_ribbon_dusk = 5.57 / R_J_per_pixel_r
          Line_centers = Fltarr(N_elements(Lines_r))
          for n = 0, N_elements(Lines_r) - 1 do begin
            junk = min(abs(Dispersion_mask_r_key[indgen(2098), center_r] - Lines_r[n]), /NAN, origin) ;compute the line center pixel from the dispersion mask
            parinfo[2].value = linewidths_r[n] & guess[2] = linewidths_r[n]
            dawn_fit = total(Pretty_r[origin - sample_width + (Lines_r[n]*dawn_shift[expected_ribbon_dawn,n]/dispers_r) : origin + sample_width + (Lines_r[n]*dawn_shift[expected_ribbon_dawn,n]/dispers_r), $ ;sum 10 pixels over the ribbon
                             center_r[[origin]] + expected_ribbon_dawn - 5: center_r[[origin]] + expected_ribbon_dawn + 5], 2) 
            dusk_fit = total(Pretty_r[origin - sample_width + (Lines_r[n]*dusk_shift[expected_ribbon_dusk,n]/dispers_r) : origin + sample_width + (Lines_r[n]*dusk_shift[expected_ribbon_dusk,n]/dispers_r), $ ;sum 10 pixels over the ribbon
                             center_r[[origin]] - expected_ribbon_dusk - 5: center_r[[origin]] - expected_ribbon_dusk + 5], 2) 
            err_dawn_fit = sqrt(total(err_Pretty_r[origin - sample_width + (Lines_r[n]*dawn_shift[expected_ribbon_dawn,n]/dispers_r) : origin + sample_width + (Lines_r[n]*dawn_shift[expected_ribbon_dawn,n]/dispers_r), $ ;sum 10 pixels over the ribbon
                                 center_r[[origin]] + expected_ribbon_dawn - 5: center_r[[origin]] + expected_ribbon_dawn + 5]^2, 2)) 
            err_dusk_fit = sqrt(total(err_Pretty_r[origin - sample_width + (Lines_r[n]*dusk_shift[expected_ribbon_dusk,n]/dispers_r) : origin + sample_width + (Lines_r[n]*dusk_shift[expected_ribbon_dusk,n]/dispers_r), $ ;sum 10 pixels over the ribbon
                                 center_r[[origin]] - expected_ribbon_dusk - 5: center_r[[origin]] - expected_ribbon_dusk + 5]^2, 2)) 
            fit_dawn = mpfitpeak(findgen(n_elements(dawn_fit)), dawn_fit, A_dawn, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dawn_fit, ESTIMATES = guess, /NAN)
            fit_dusk = mpfitpeak(findgen(n_elements(dusk_fit)), dusk_fit, A_dusk, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_dusk_fit, ESTIMATES = guess, /NAN)
                              
            Line_centers[n] = origin ;+ (mean([a_dawn[1],a_dusk[1]], /Nan) - sample_width)  ;REFINE THE LINE CENTERS BASED ON THE FITS
            if keyword_set(debug) then begin              
              x = findgen(2*sample_width) - sample_width
              if n eq 0 then window, 2, xs = 600, ys = 400, title = 'Dawn-side Line Profiles at the Ribbon' 
              wset, 2
              if n eq 0 then cgplot, x, dawn_fit, psym = 10, color = red_colors[n], yr = [-100, 1700], xstyle=1 else cgplot, x, dawn_fit, psym = 10, color = red_colors[n], /overplot     
              ;print, 'shifted Dawn centroid of ', lines_r[n], ' line by ', (a_dawn[1]- sample_width) 
              cgplot, x, fit_dawn, color = red_colors[n], /overplot
              
              if n eq 0 then window, 3, xs = 600, ys = 400, title = 'Dusk-side Line Profiles at the Ribbon' 
              wset, 3
              if n eq 0 then cgplot, x, dusk_fit, psym = 10, color = red_colors[n], yr = [-100, 1700], xstyle=1 else cgplot, x, dusk_fit, psym = 10, color = red_colors[n], /overplot     
              ;print, 'shifted Dusk centroid of ', lines_r[n], ' line by ', (a_dusk[1]- sample_width) 
              cgplot, x, fit_dusk, color = red_colors[n], /overplot
            endif            
          endfor
        if keyword_set(debug) then begin ;check line centers!   
          range = indgen(2098)
          y = total(Pretty_r[ min(range):max(range), mean(center_r[range]) + expected_ribbon_dawn - 15: mean(center_r[range]) + expected_ribbon_dawn + 15 ], 2, /NAN) 
          x = Dispersion_mask_r[range, mean(center_r[range])]
          window, 0, xs = 600, ys = 400
          cgplot, x, y, xr = [6300, 6324], psym = 10, yr = [-100, 1100], ystyle = 1
          cgplot, Dispersion_mask_r[round(Line_centers + (Lines_r*dawn_shift[expected_ribbon_dawn,*]/dispers_r)), center_r[round(Line_centers + (Lines_r*dawn_shift[expected_ribbon_dawn,*]/dispers_r))]], $
                  make_array(N_elements(Line_centers), /float, value = 1000.), psym =4, /overplot
          window, 1, xs = 600, ys = 400
          cgplot, x, y, xr = [6710, 6740], psym = 10, yr = [-100, 2000], ystyle = 1
          cgplot, Dispersion_mask_r[round(Line_centers + (Lines_r*dawn_shift[expected_ribbon_dawn,*]/dispers_r)), center_r[round(Line_centers + (Lines_r*dawn_shift[expected_ribbon_dawn,*]/dispers_r))]], $
                  make_array(N_elements(Line_centers), /float, value = 1700.), psym =4, /overplot
        endif 
      ;--------------------------------------------------    
          
      for n = 0, N_elements(Lines_r) - 1 do begin
        parinfo[2].value = linewidths_r[n] & guess[2] = linewidths_r[n]
        junk = min(abs(Dispersion_mask_r_key[indgen(2098), key_r] - Lines_r[n]), /NAN, origin)
        for i = 0, n_elements(dawn_dist_r_km[*,n]) - 1 do begin
          fitme = pretty_r[origin - sample_width + (Lines_r[n]*dawn_shift[i]/dispers_r) : origin + sample_width + (Lines_r[n]*dawn_shift[i]/dispers_r), center_r[origin] + i]
          err_fitme = err_pretty_r[origin - sample_width + (Lines_r[n]*dawn_shift[i]/dispers_r) : origin + sample_width + (Lines_r[n]*dawn_shift[i]/dispers_r), center_r[origin] + i]
          fit = mpfitpeak(findgen(n_elements(fitme)), fitme, A, PERROR = err_A, /POSITIVE, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_fitme, ESTIMATES = guess, $
            bestNORM = BESTNORM, DOF = DOF, COVAR = COVAR, /NAN)
          dawn_avg_r[i,n] = A[0]*A[2]*SQRT(2*!DPI)  
          err_dawn_avg_r[i,n] = dawn_avg_r[i,n] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. ) 
        endfor 
        for i = 0, n_elements(dusk_dist_r_km[*,n]) - 1 do begin
          fitme = pretty_r[origin - sample_width + (Lines_r[n]*dusk_shift[i]/dispers_r) : origin + sample_width + (Lines_r[n]*dusk_shift[i]/dispers_r), center_r[origin] - i]
          err_fitme = err_pretty_r[origin - sample_width + (Lines_r[n]*dusk_shift[i]/dispers_r) : origin + sample_width + (Lines_r[n]*dusk_shift[i]/dispers_r), center_r[origin] - i]
          fit = mpfitpeak(findgen(n_elements(fitme)), fitme, A, PERROR = err_A, /POSITIVE, NTERMS = 5, parinfo=parinfo, MEASURE_ERRORS = err_fitme, ESTIMATES = guess, $
            bestNORM = BESTNORM, DOF = DOF, COVAR = COVAR, /NAN)
          dusk_avg_r[i,n] = A[0]*A[2]*SQRT(2*!DPI)
          err_dusk_avg_r[i,n] = dusk_avg_r[i,n] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. ) 
        endfor      
        if keyword_set(debug) then begin  
          if n eq 0 then window, 0, xs = 600, ys = 400
          if n eq 0 then cgplot, fitme, psym = 4, color = red_colors[n], yr = [-100, 1700] else cgplot, fitme, psym = 4, color = red_colors[n], /overplot     
          cgplot, fit, color = red_colors[n], /overplot
        endif
      endfor

      ;-------------Combine Blue Images for Equivalent Widths----------------------------
      bigarray = fltarr(2098, 1078, N_elements(CoAdd_frames[*,j]))
      err_bigarray = fltarr(2098, 1078, N_elements(CoAdd_frames[*,j]))
      for i = 0, N_elements(CoAdd_frames[*,j]) - 1 do begin
        img = MRDFITS(Images_b[CoAdd_frames[i,j]], 0, /silent)
        err_img = MRDFITS(Images_b[CoAdd_frames[i,j]], 1, /silent)
        metadata_b = MRDFITS(Images_b[CoAdd_frames[i,j]], 2, /silent)  
        center = metadata_b.center
        dispersion_mask_b = metadata_b.dispersion_mask
        img[anti_aperture_blue] = !values.f_nan ;don't include points not in the slit aperture
        err_img[anti_aperture_blue] = !values.f_nan ;don't include points not in the slit aperture
        bigarray[*,*,i] = smart_shift(img, 0, key_b - mean(metadata_b.center) ,/interp, Missing = !values.f_nan)
        err_bigarray[*,*,i] = smart_shift(err_img, 0, key_b - mean(metadata_b.center) ,/interp, Missing = !values.f_nan)
      endfor  
      
      ;replace outliers in the stack, careful not to replace the torus instrinsic variability
        sig = STDDEV(bigarray, DIMENSION=3, /DOUBLE, /NAN)
        stack_mean = MEAN(bigarray, DIMENSION=3, /double, /NAN) 
        deviation = bigarray - rebin(stack_mean, 2098, 1078, N_elements(CoAdd_frames[*,j]))
        bigarray[where(abs(deviation) gt 2.75*sig, /NULL)] = !values.F_NAN
        err_bigarray[where(abs(deviation) gt 2.75*sig, /NULL)] = !values.F_NAN
        pretty_b = MEAN(bigarray, DIMENSION=3, /NAN) 
        err_pretty_b = sqrt(total(err_bigarray^2, 3, /NAN)) / N_elements(CoAdd_frames[*,j])
        if keyword_set(debug) then begin
          window, 0, xs = 2098, ys =1000
          tv, bytscl(stack_mean - pretty_b, -2, 2)
        endif
        dummy = pretty_b

      dawn_dist_b_km = distance_profiles_b.dawn_distance_array_b * R_J
      dusk_dist_b_km = distance_profiles_b.dusk_distance_array_b * R_J
      dawn_shift = (-2. * !pi * dawn_dist_b_km / period_J) / c  
      dusk_shift = ( 2. * !pi * dusk_dist_b_km / period_J) / c    
      dawn_avg_b = fltarr(N_elements(dawn_dist_b_km[*,0]), N_elements(Lines_b)) & err_dawn_avg_b = fltarr(N_elements(dawn_dist_b_km[*,0]), N_elements(Lines_b))
      dusk_avg_b = fltarr(N_elements(dusk_dist_b_km[*,0]), N_elements(Lines_b)) & err_dusk_avg_b = fltarr(N_elements(dusk_dist_b_km[*,0]), N_elements(Lines_b))
      expected_ribbon_dawn = 5.85 / R_J_per_pixel_b
      expected_ribbon_dusk = 5.57 / R_J_per_pixel_b
              
      ;----------multiline fitting----Inspectrion only------------          
      if keyword_set(debug) then begin
              junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[1]), /NAN, origin) ;compute the line center pixel from the dispersion mask
              y = total(Pretty_b[origin - sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[expected_ribbon_dawn,1]/dispers_b) : $
                                 origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[expected_ribbon_dawn,1]/dispers_b), $ ;sum 10 pixels over the ribbon
                                 center_b[[origin]] + expected_ribbon_dawn - 5: center_b[[origin]] + expected_ribbon_dawn + 5], 2) 
              err_y = sqrt(total(err_Pretty_b[origin - sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[expected_ribbon_dawn,1]/dispers_b) : $
                                              origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[expected_ribbon_dawn,1]/dispers_b), $ ;sum 10 pixels over the ribbon
                                              center_b[[origin]] + expected_ribbon_dawn - 5: center_b[[origin]] + expected_ribbon_dawn + 5]^2, 2)) 
              fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
              p = mpfit( 'Fit_blue_triplet', initial_guess_Blue_triplet, PERROR = err_p, funct=fa, maxiter=100, $
                  STATUS = STATUS, parinfo = parinfo_blue_triplet, Npegged=Npegged, /quiet)       
              if status le 0 then P[*] = !values.f_Nan & err_P = make_array(5, value = !values.f_Nan)
              window, 0
              cgplot, fa.x, fa.y, psym = 10, ERR_YLOW = fa.err, ERR_YHigh = fa.err
              cgplot, fa.x, blue_triplet(p, fa.x), /overplot                                       

              junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[3]), /NAN, origin) ;compute the line center pixel from the dispersion mask
              y = total(Pretty_b[origin - sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[expected_ribbon_dawn,3]/dispers_b) : $
                                 origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[expected_ribbon_dawn,3]/dispers_b), $ ;sum 10 pixels over the ribbon
                                 center_b[[origin]] + expected_ribbon_dawn - 5: center_b[[origin]] + expected_ribbon_dawn + 5], 2) 
              err_y = sqrt(total(err_Pretty_b[origin - sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[expected_ribbon_dawn,3]/dispers_b) : $ 
                                 origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[expected_ribbon_dawn,3]/dispers_b), $ ;sum 10 pixels over the ribbon
                                 center_b[[origin]] + expected_ribbon_dawn - 5: center_b[[origin]] + expected_ribbon_dawn + 5]^2, 2)) 
              fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
              p = mpfit( 'Fit_blue_doublet', initial_guess_blue_doublet, PERROR = err_p, funct=fa, maxiter = 100, $
                  STATUS = STATUS, parinfo = parinfo_Blue_doublet, Npegged=Npegged, /quiet)       
              if status le 0 then P[*] = !values.f_Nan & err_P = make_array(4, value = !values.f_Nan)
              window, 1
              cgplot, fa.x, fa.y, psym = 10, ERR_YLOW = fa.err, ERR_YHigh = fa.err
              cgplot, fa.x, blue_doublet(p, fa.x), /overplot
              ;stop
        endif ;end debug inspection

      ;----------------multiline integrated flux----------------
        for i = 0, N_elements(dawn_avg_b[*,0]) - 1 do begin
          junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[1]), /NAN, origin) ;compute the line center pixel from the dispersion mask
          y = pretty_b[origin - sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[i,1]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[i,1]/dispers_b), center_b[origin] + i]
          err_y = err_pretty_b[origin - sample_width_Blue_triplet  + (Lines_b[1]*dawn_shift[i,1]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dawn_shift[i,1]/dispers_b), center_b[origin] + i]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_triplet', initial_guess_Blue_triplet, PERROR = err_p, funct=fa, STATUS = Status, parinfo = parinfo_Blue_triplet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan 
          if status le 0 then err_P = make_array(5, value = !values.f_Nan)
          dawn_avg_b[i,0] =     P[0]*P[4]*SQRT(2*!DPI) 
          dawn_avg_b[i,1] =     P[1]*P[4]*SQRT(2*!DPI) 
          dawn_avg_b[i,2] =     P[2]*P[4]*SQRT(2*!DPI) 
          err_dawn_avg_b[i,0] = P[0]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.    
          err_dawn_avg_b[i,1] = P[1]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance. 
          err_dawn_avg_b[i,2] = P[2]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[2]/P[2])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.      
          dummy[origin + (Lines_b[1]*dawn_shift[i,1]/dispers_b), center_b[origin] + i] = 1000. ;mark the locations
          
          junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[3]), /NAN, origin) ;compute the line center pixel from the dispersion mask  
          y = pretty_b[origin - sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[i,3]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[i,3]/dispers_b), center_b[origin] + i]
          err_y = err_pretty_b[origin - sample_width_Blue_doublet  + (Lines_b[3]*dawn_shift[i,3]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dawn_shift[i,3]/dispers_b), center_b[origin] + i]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_doublet', initial_guess_Blue_doublet, PERROR = err_p, funct=fa, STATUS = Status, parinfo = parinfo_Blue_doublet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan  
          if status le 0 then err_P = make_array(4, value = !values.f_Nan)
          dawn_avg_b[i,3] =     P[0]*P[3]*SQRT(2*!DPI) 
          dawn_avg_b[i,4] =     P[1]*P[3]*SQRT(2*!DPI) 
          err_dawn_avg_b[i,3] = P[0]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[3]/P[3])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.    
          err_dawn_avg_b[i,4] = P[1]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[3]/P[3])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance. 

          dummy[origin + (Lines_b[1]*dawn_shift[i,3]/dispers_b), center_b[origin] + i] = 1000. ;mark the locations  
        endfor 
        for i = 0, N_elements(dusk_avg_b[*,0]) - 1 do begin 
          junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[1]), /NAN, origin) ;compute the line center pixel from the dispersion mask
          y = pretty_b[origin - sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[i,1]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[i,1]/dispers_b), center_b[origin] - i]
          err_y = err_pretty_b[origin - sample_width_Blue_triplet  + (Lines_b[1]*dusk_shift[i,1]/dispers_b) : origin + sample_width_Blue_triplet + (Lines_b[1]*dusk_shift[i,1]/dispers_b), center_b[origin] - i]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_triplet', initial_guess_Blue_triplet, PERROR = err_p, funct=fa, STATUS = status, parinfo = parinfo_Blue_triplet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan 
          if status le 0 then err_P = make_array(5, value = !values.f_Nan)
          dusk_avg_b[i,0] =     P[0]*P[4]*SQRT(2*!DPI) 
          dusk_avg_b[i,1] =     P[1]*P[4]*SQRT(2*!DPI) 
          dusk_avg_b[i,2] =     P[2]*P[4]*SQRT(2*!DPI) 
          err_dusk_avg_b[i,0] = P[0]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.    
          err_dusk_avg_b[i,1] = P[1]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance. 
          err_dusk_avg_b[i,2] = P[2]*P[4]*SQRT(2*!DPI) * sqrt( (err_P[2]/P[2])^2. + (err_P[4]/P[4])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.      
          dummy[origin + (Lines_b[1]*dusk_shift[i,1]/dispers_b), center_b[origin] + i] = 1000. ;mark the locations
          
          junk = min(abs(Dispersion_mask_b_key[indgen(2098), center_b] - Lines_b[3]), /NAN, origin) ;compute the line center pixel from the dispersion mask
          y = pretty_b[origin - sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[i,3]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[i,3]/dispers_b), center_b[origin] - i]
          err_y = err_pretty_b[origin - sample_width_Blue_doublet  + (Lines_b[3]*dusk_shift[i,3]/dispers_b) : origin + sample_width_Blue_doublet + (Lines_b[3]*dusk_shift[i,3]/dispers_b), center_b[origin] - i]
          fa = {x:findgen(N_elements(y)), y:y, err:err_y} 
          p = mpfit( 'Fit_blue_doublet', initial_guess_Blue_doublet, PERROR = err_p, funct=fa, STATUS=Status, parinfo = parinfo_Blue_doublet, maxiter = 1000, /quiet) 
          if status le 0 then P[*] = !values.f_Nan 
          if status le 0 then err_P = make_array(4, value = !values.f_Nan) 
          dusk_avg_b[i,3] =     P[0]*P[3]*SQRT(2*!DPI) 
          dusk_avg_b[i,4] =     P[1]*P[3]*SQRT(2*!DPI) 
          err_dusk_avg_b[i,3] = P[0]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[0]/P[0])^2. + (err_P[3]/P[3])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance.    
          err_dusk_avg_b[i,4] = P[1]*P[3]*SQRT(2*!DPI) * sqrt( (err_P[1]/P[1])^2. + (err_P[3]/P[3])^2.)  ; uncertainties in width and height add in quadrature, neglect covariance. 

          dummy[origin + (Lines_b[1]*dusk_shift[i,3]/dispers_b), center_b[origin] + i] = 1000. ;mark the locations  
        endfor 
      if keyword_set(debug) then begin  
        window, 1, xs = 2000, ys = 1000
        tv, bytscl(dummy, 0, 50)
      endif  
         
  date = strmid(sxpar(torus_header_b, 'DATE-OBS'),0,10)
  R_J_per_pixel_b = float(sxpar(torus_header_b, 'RJ_PER_P'))
  R_J_per_pixel_r = float(sxpar(torus_header_r, 'RJ_PER_P'))
  
  err_dusk_avg_r = abs(err_dusk_avg_r)
  err_dusk_avg_b = abs(err_dusk_avg_b)
  err_dawn_avg_r = abs(err_dawn_avg_r)
  err_dawn_avg_b = abs(err_dawn_avg_b)
  
  ;SMOOTH in the spatial dimension. These go in the save files for the NIGHTLY co-added data, UNSMOOTHED already written to fits themselves 
    print, 'Smoothing 4069, 6716,6731 by (Rj)', smooth_bright_by
    print, 'Smoothing The Other Lines by (Rj)', smooth_by 
    dawn_avg_b[*,[0,1,2,4]] = Smooth(dawn_avg_b[*,[0,1,2,4]], [round(smooth_by / R_J_per_pixel_b), 0], /edge_mirror)
    dawn_avg_b[*,3] = Smooth(dawn_avg_b[*,3], round(smooth_bright_by / R_J_per_pixel_b), /edge_mirror)
    err_dawn_avg_b[*,[0,1,2,4]] = err_dawn_avg_b[*,[0,1,2,4]] / sqrt(round(smooth_by / R_J_per_pixel_b))
    err_dawn_avg_b[*,3] = err_dawn_avg_b[*,3] / sqrt(round(smooth_bright_by / R_J_per_pixel_b))
    dawn_avg_r[*,0] = Smooth(dawn_avg_r[*,0], round(smooth_by / R_J_per_pixel_r), /edge_mirror)
    dawn_avg_r[*,1:2] = Smooth(dawn_avg_r[*,1:2], [round(smooth_bright_by / R_J_per_pixel_r), 0], /edge_mirror) 
    err_dawn_avg_r[*,0] = err_dawn_avg_r[*,0] / sqrt(round(smooth_by / R_J_per_pixel_r))
    err_dawn_avg_r[*,1:2] = err_dawn_avg_r[*,1:2] / sqrt(round(smooth_bright_by / R_J_per_pixel_r))
    dusk_avg_b[*,[0,1,2,4]] = Smooth(dusk_avg_b[*,[0,1,2,4]], [round(smooth_by / R_J_per_pixel_b), 0], /edge_mirror)
    dusk_avg_b[*,3] = Smooth(dusk_avg_b[*,3], round(smooth_bright_by / R_J_per_pixel_b), /edge_mirror)
    err_dusk_avg_b[*,[0,1,2,4]] = err_dusk_avg_b[*,[0,1,2,4]] / sqrt(round(smooth_by / R_J_per_pixel_b))
    err_dusk_avg_b[*,3] = err_dusk_avg_b[*,3] / sqrt(round(smooth_bright_by / R_J_per_pixel_b))
    dusk_avg_r[*,0] = Smooth(dusk_avg_r[*,0], round(smooth_by / R_J_per_pixel_r), /edge_mirror)
    dusk_avg_r[*,1:2] = Smooth(dusk_avg_r[*,1:2], [round(smooth_bright_by / R_J_per_pixel_r), 0], /edge_mirror) 
    err_dusk_avg_r[*,0] = err_dusk_avg_r[*,0] / sqrt(round(smooth_by / R_J_per_pixel_r))
    err_dusk_avg_r[*,1:2] = err_dusk_avg_r[*,1:2] / sqrt(round(smooth_bright_by / R_J_per_pixel_r))
    
    dawn_search_peak = where((distance_profiles_r.dawn_distance_array_r gt 5.) and (distance_profiles_r.dawn_distance_array_r lt 7.))
    dusk_search_peak = where((distance_profiles_r.dusk_distance_array_r gt 5.) and (distance_profiles_r.dusk_distance_array_r lt 7.))
    plot_ymax = max([max(dawn_avg_r[dawn_search_peak,2]), max(dusk_avg_r[dusk_search_peak,2])])

  cgPS_Open, filename = Directory + strcompress('Processed\Figures\Profiles_'+strmid(sxpar(torus_header_r, 'DATE-OBS'), 0,10) + '_' + strcompress(j) + '.eps', /remove_all), /ENCAPSULATED, xsize = 16., ysize = 6.5
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.

        ;Setup 2 panels
        positions = cglayout([2,1], XGAP = 1., OYMargin=[7,5], OXMargin=[8,8] )

        ;Setup Axis
        cgplot, distance_profiles_r.dawn_distance_array_r, dawn_avg_r[*,0], xrange = [8.1,4.5], xstyle = 1, yrange = [0, 1.05*plot_ymax], ytitle = 'Rayleighs', $
                xtitle = 'Dawn-side Distance (R!DJupiter!N)', psym = 10, thick = 6., Title = Date + ' - Ansa ' + cgSymbol('lambda') + '!DIII!N = ' + string(CML_Dawn[0], $
                FORMAT = '(F5.1)') + cgSymbol('deg') + ' to ' + string(CML_Dawn[-1], FORMAT = '(F5.1)') + cgSymbol('deg'), /NODATA, position = positions[*,0], YStyle=1, /NoErase
        y = !y.crange
        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, 1.1*plot_ymax], /overplot, linestyle = 2.
          ribbon_range = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [1., 1.048*plot_ymax, 1.048*plot_ymax, 1.], /data, COLOR='yellow'

        cgplot, distance_profiles_r.dawn_distance_array_r[*,0], dawn_avg_r[*,0], ERR_YLOW = err_dawn_avg_r[*,0], ERR_YHIGH = err_dawn_avg_r[*,0], /overplot, psym = 10, thick = 6., color = Red_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dawn_distance_array_r[*,1], dawn_avg_r[*,1], ERR_YLOW = err_dawn_avg_r[*,1], ERR_YHIGH = err_dawn_avg_r[*,1], /overplot, psym = 10, thick = 6., color = Red_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dawn_distance_array_r[*,2], dawn_avg_r[*,2], ERR_YLOW = err_dawn_avg_r[*,2], ERR_YHIGH = err_dawn_avg_r[*,2], /overplot, psym = 10, thick = 6., color = Red_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,0], dawn_avg_b[*,0], ERR_YLOW = err_dawn_avg_b[*,0], ERR_YHIGH = err_dawn_avg_b[*,0], /overplot, psym = 10, thick = 6., color = Blue_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,1], dawn_avg_b[*,1], ERR_YLOW = err_dawn_avg_b[*,1], ERR_YHIGH = err_dawn_avg_b[*,1], /overplot, psym = 10, thick = 6., color = Blue_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,2], dawn_avg_b[*,2], ERR_YLOW = err_dawn_avg_b[*,2], ERR_YHIGH = err_dawn_avg_b[*,2], /overplot, psym = 10, thick = 6., color = Blue_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,3], dawn_avg_b[*,3], ERR_YLOW = err_dawn_avg_b[*,3], ERR_YHIGH = err_dawn_avg_b[*,3], /overplot, psym = 10, thick = 6., color = Blue_colors[3], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,4], dawn_avg_b[*,4], ERR_YLOW = err_dawn_avg_b[*,4], ERR_YHIGH = err_dawn_avg_b[*,4], /overplot, psym = 10, thick = 6., color = Blue_colors[4], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgText, 7.95, .34*(1.1*plot_ymax), 'SIII 3722'+cgSymbol('Angstrom'), color = cgcolor(Blue_colors[0])    
        cgText, 7.95, .42*(1.1*plot_ymax), 'OII 3726' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[1])
        cgText, 7.95, .5*(1.1*plot_ymax),  'OII 3729' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[2])  
        cgText, 7.95, .58*(1.1*plot_ymax), 'SII 4069' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[3])
        cgText, 7.95, .66*(1.1*plot_ymax), 'SII 4076' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[4])
        cgText, 7.95, .74*(1.1*plot_ymax), 'SIII 6312'+cgSymbol('Angstrom'), color = cgcolor(Red_colors[0])        
        cgText, 7.95, .82*(1.1*plot_ymax), 'SII 6716' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[1]) 
        cgText, 7.95, .9*(1.1*plot_ymax),  'SII 6731' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[2]) 
      
        ;Setup Axis
        cgplot, distance_profiles_r.dusk_distance_array_r[*,0], dusk_avg_r[*,0], xrange = [4.5,8.1], xstyle = 1, yrange = [0, 1.05*plot_ymax], $
                xtitle = 'Dusk-side Distance (R!DJupiter!N)', psym = 10, thick = 6., Title = Date + ' - Ansa ' + cgSymbol('lambda') + '!DIII!N = ' + string(CML_dusk[0], $
                FORMAT = '(F5.1)') + cgSymbol('deg') + ' to ' + string(CML_dusk[-1], FORMAT = '(F5.1)') + cgSymbol('deg'), /NODATA, YTickformat='(A1)', YStyle=1, Position = positions[*,1], /NoErase, /window;, Position=[0.075, 0.125, 0.875, 0.9]
        ;cgaxis, /Yaxis, yrange = [0, 1.05*plot_ymax], ytitle = 'Rayleighs', /window, /save

        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, 1.1*plot_ymax], /overplot, linestyle = 2.
          ribbon_range = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg);dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [1., 1.048*plot_ymax, 1.048*plot_ymax, 1.], /data, COLOR='yellow'
          
        cgplot, distance_profiles_r.dusk_distance_array_r[*,0], dusk_avg_r[*,0], ERR_YLOW = err_dusk_avg_r[*,0], ERR_YHIGH = err_dusk_avg_r[*,0], /overplot, psym = 10, thick = 6., color = Red_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dusk_distance_array_r[*,1], dusk_avg_r[*,1], ERR_YLOW = err_dusk_avg_r[*,1], ERR_YHIGH = err_dusk_avg_r[*,1], /overplot, psym = 10, thick = 6., color = Red_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dusk_distance_array_r[*,2], dusk_avg_r[*,2], ERR_YLOW = err_dusk_avg_r[*,2], ERR_YHIGH = err_dusk_avg_r[*,2], /overplot, psym = 10, thick = 6., color = Red_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,0], dusk_avg_b[*,0], ERR_YLOW = err_dusk_avg_b[*,0], ERR_YHIGH = err_dusk_avg_b[*,0], /overplot, psym = 10, thick = 6., color = Blue_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,1], dusk_avg_b[*,1], ERR_YLOW = err_dusk_avg_b[*,1], ERR_YHIGH = err_dusk_avg_b[*,1], /overplot, psym = 10, thick = 6., color = Blue_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,2], dusk_avg_b[*,2], ERR_YLOW = err_dusk_avg_b[*,2], ERR_YHIGH = err_dusk_avg_b[*,2], /overplot, psym = 10, thick = 6., color = Blue_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,3], dusk_avg_b[*,3], ERR_YLOW = err_dusk_avg_b[*,3], ERR_YHIGH = err_dusk_avg_b[*,3], /overplot, psym = 10, thick = 6., color = Blue_colors[3], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,4], dusk_avg_b[*,4], ERR_YLOW = err_dusk_avg_b[*,4], ERR_YHIGH = err_dusk_avg_b[*,4], /overplot, psym = 10, thick = 6., color = Blue_colors[4], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgText, 7.3, .34*(1.1*plot_ymax), 'SIII 3722'+cgSymbol('Angstrom'), color = cgcolor(Blue_colors[0])    
        cgText, 7.3, .42*(1.1*plot_ymax), 'OII 3726' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[1])
        cgText, 7.3, .5*(1.1*plot_ymax),  'OII 3729' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[2])  
        cgText, 7.3, .58*(1.1*plot_ymax), 'SII 4069' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[3])
        cgText, 7.3, .66*(1.1*plot_ymax), 'SII 4076' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[4])
        cgText, 7.3, .74*(1.1*plot_ymax), 'SIII 6312'+cgSymbol('Angstrom'), color = cgcolor(Red_colors[0])        
        cgText, 7.3, .82*(1.1*plot_ymax), 'SII 6716' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[1]) 
        cgText, 7.3, .9*(1.1*plot_ymax),  'SII 6731' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[2]) 
    cgPS_Close 
    set_plot,'WIN'

  ;*******************************************************************Dawn & Dusk Individual Plots**************************************************************************************
  cgPS_Open, filename = Directory + strcompress('Processed\Figures\Dawn_Profiles_'+strmid(sxpar(torus_header_r, 'DATE-OBS'), 0,10) + '_' + strcompress(j) + '_New.eps', /remove_all), /ENCAPSULATED, xsize = 8.5, ysize = 7.5
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.

        ;Setup Axis
        cgplot, distance_profiles_r.dawn_distance_array_r[*,0], dawn_avg_r[*,0], xrange = [4.5,8.1], xstyle = 1, yrange = [0, 1.05*max(dawn_avg_r[dawn_search_peak,2])], ytitle = 'Rayleighs', $
                xtitle = 'Dawn-side Distance (R!DJupiter!N)', psym = 10, thick = 6., Title = Date + ' - Ansa ' + cgSymbol('lambda') + '!DIII!N = ' + string(CML_Dawn[0], $
                FORMAT = '(F5.1)') + cgSymbol('deg') + ' to ' + string(CML_Dawn[-1], FORMAT = '(F5.1)') + cgSymbol('deg'), /NODATA
        
        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, 1.1*max(dawn_avg_r[dawn_search_peak,2])], /overplot, linestyle = 2.
          ribbon_range = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [1., 1.048*max(dawn_avg_r[dawn_search_peak,2]), 1.048*max(dawn_avg_r[dawn_search_peak,2]), 1.], /data, COLOR='yellow'
          
        cgplot, distance_profiles_r.dawn_distance_array_r[*,0], dawn_avg_r[*,0], ERR_YLOW = err_dawn_avg_r[*,0], ERR_YHIGH = err_dawn_avg_r[*,0], /overplot, psym = 10, thick = 6., color = Red_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dawn_distance_array_r[*,1], dawn_avg_r[*,1], ERR_YLOW = err_dawn_avg_r[*,1], ERR_YHIGH = err_dawn_avg_r[*,1], /overplot, psym = 10, thick = 6., color = Red_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dawn_distance_array_r[*,2], dawn_avg_r[*,2], ERR_YLOW = err_dawn_avg_r[*,2], ERR_YHIGH = err_dawn_avg_r[*,2], /overplot, psym = 10, thick = 6., color = Red_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,0], dawn_avg_b[*,0], ERR_YLOW = err_dawn_avg_b[*,0], ERR_YHIGH = err_dawn_avg_b[*,0], /overplot, psym = 10, thick = 6., color = Blue_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,1], dawn_avg_b[*,1], ERR_YLOW = err_dawn_avg_b[*,1], ERR_YHIGH = err_dawn_avg_b[*,1], /overplot, psym = 10, thick = 6., color = Blue_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,2], dawn_avg_b[*,2], ERR_YLOW = err_dawn_avg_b[*,2], ERR_YHIGH = err_dawn_avg_b[*,2], /overplot, psym = 10, thick = 6., color = Blue_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,3], dawn_avg_b[*,3], ERR_YLOW = err_dawn_avg_b[*,3], ERR_YHIGH = err_dawn_avg_b[*,3], /overplot, psym = 10, thick = 6., color = Blue_colors[3], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dawn_distance_array_b[*,4], dawn_avg_b[*,4], ERR_YLOW = err_dawn_avg_b[*,4], ERR_YHIGH = err_dawn_avg_b[*,4], /overplot, psym = 10, thick = 6., color = Blue_colors[4], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgText, 7.3, .34*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'SIII 3722'+cgSymbol('Angstrom'), color = cgcolor(Blue_colors[0])    
        cgText, 7.3, .42*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'OII 3726' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[1])
        cgText, 7.3, .5*(1.1*max(dawn_avg_r[dawn_search_peak,2])),  'OII 3729' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[2])  
        cgText, 7.3, .58*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'SII 4069' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[3])
        cgText, 7.3, .66*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'SII 4076' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[4])
        cgText, 7.3, .74*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'SIII 6312'+cgSymbol('Angstrom'), color = cgcolor(Red_colors[0])        
        cgText, 7.3, .82*(1.1*max(dawn_avg_r[dawn_search_peak,2])), 'SII 6716' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[1]) 
        cgText, 7.3, .9*(1.1*max(dawn_avg_r[dawn_search_peak,2])),  'SII 6731' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[2]) 
      cgPS_Close 
      set_plot,'WIN'      

  cgPS_Open, filename=Directory + strcompress('Processed\Figures\Dusk_Profiles_'+strmid(sxpar(torus_header_r, 'DATE-OBS'), 0,10) + '_' + strcompress(j) + '_New.eps', /remove_all), /ENCAPSULATED, xsize = 8.5, ysize = 7.5
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.
        
        ;Setup Axis
        cgplot, distance_profiles_r.dusk_distance_array_r, dusk_avg_r[*,0], xrange = [4.5,8.1], xstyle = 1, yrange = [0, 1.05*max(dusk_avg_r[dusk_search_peak,2])], $
                xtitle = 'Dusk-side Distance (R!DJupiter!N)', psym = 10, thick = 6., Title = Date + ' - Ansa ' + cgSymbol('lambda') + '!DIII!N = ' + string(CML_dusk[0], $
                FORMAT = '(F5.1)') + cgSymbol('deg') + ' to ' + string(CML_dusk[-1], FORMAT = '(F5.1)') + cgSymbol('deg'), /NODATA, YStyle=9, YTickformat='(A1)', Position=[0.075, 0.125, 0.875, 0.9]
        cgaxis, Yaxis = 1, ytitle = 'Rayleighs'
        
        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, 1.1*max(dusk_avg_r[dusk_search_peak,2])], /overplot, linestyle = 2.
          ribbon_range = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg);dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [1., 1.048*max(dusk_avg_r[dusk_search_peak,2]), 1.048*max(dusk_avg_r[dusk_search_peak,2]), 1.], /data, COLOR='yellow'
          
        cgplot, distance_profiles_r.dusk_distance_array_r[*,0], dusk_avg_r[*,0], ERR_YLOW = err_dusk_avg_r[*,0], ERR_YHIGH = err_dusk_avg_r[*,0], /overplot, psym = 10, thick = 6., color = Red_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dusk_distance_array_r[*,1], dusk_avg_r[*,1], ERR_YLOW = err_dusk_avg_r[*,1], ERR_YHIGH = err_dusk_avg_r[*,1], /overplot, psym = 10, thick = 6., color = Red_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_r.dusk_distance_array_r[*,2], dusk_avg_r[*,2], ERR_YLOW = err_dusk_avg_r[*,2], ERR_YHIGH = err_dusk_avg_r[*,2], /overplot, psym = 10, thick = 6., color = Red_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,0], dusk_avg_b[*,0], ERR_YLOW = err_dusk_avg_b[*,0], ERR_YHIGH = err_dusk_avg_b[*,0], /overplot, psym = 10, thick = 6., color = Blue_colors[0], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,1], dusk_avg_b[*,1], ERR_YLOW = err_dusk_avg_b[*,1], ERR_YHIGH = err_dusk_avg_b[*,1], /overplot, psym = 10, thick = 6., color = Blue_colors[1], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,2], dusk_avg_b[*,2], ERR_YLOW = err_dusk_avg_b[*,2], ERR_YHIGH = err_dusk_avg_b[*,2], /overplot, psym = 10, thick = 6., color = Blue_colors[2], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,3], dusk_avg_b[*,3], ERR_YLOW = err_dusk_avg_b[*,3], ERR_YHIGH = err_dusk_avg_b[*,3], /overplot, psym = 10, thick = 6., color = Blue_colors[3], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgplot, distance_profiles_b.dusk_distance_array_b[*,4], dusk_avg_b[*,4], ERR_YLOW = err_dusk_avg_b[*,4], ERR_YHIGH = err_dusk_avg_b[*,4], /overplot, psym = 10, thick = 6., color = Blue_colors[4], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgText, 7.3, .34*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'SIII 3722'+cgSymbol('Angstrom'), color = cgcolor(Blue_colors[0])    
        cgText, 7.3, .42*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'OII 3726' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[1])
        cgText, 7.3, .5*(1.1*max(dusk_avg_r[dusk_search_peak,2])),  'OII 3729' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[2])  
        cgText, 7.3, .58*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'SII 4069' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[3])
        cgText, 7.3, .66*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'SII 4076' +cgSymbol('Angstrom'), color = cgcolor(Blue_colors[4])
        cgText, 7.3, .74*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'SIII 6312'+cgSymbol('Angstrom'), color = cgcolor(Red_colors[0])        
        cgText, 7.3, .82*(1.1*max(dusk_avg_r[dusk_search_peak,2])), 'SII 6716' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[1]) 
        cgText, 7.3, .9*(1.1*max(dusk_avg_r[dusk_search_peak,2])),  'SII 6731' +cgSymbol('Angstrom'), color = cgcolor(Red_colors[2]) 
      cgPS_Close 
      set_plot,'WIN'
      loadct, 0

    ;***********************************************************************Image Plots**************************************************************************************
     cgPS_Open, filename = Directory + strcompress('Processed\Figures\Pretty_'+strmid(sxpar(torus_header_r, 'DATE-OBS'), 0,10) + '_' + strcompress(j) + '_r.eps', /remove_all), /ENCAPSULATED, xsize = 8.5, ysize = 4. 
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.      
        minvalue = 0
        maxvalue = 50
        Ycrop = [ceil(mean(CENTER_R - max(distance_profiles_r.dusk_distance_array_r[*,0]) / R_J_per_pixel_R)), $
                floor(mean(CENTER_R + max(distance_profiles_r.dawn_distance_array_r[*,0]) / R_J_per_pixel_R))]   
        if float(sxpar(torus_header_r, 'DISPWC')) gt 6500. then begin
          xstart = 6300.
          xfin = 7350.
        endif else begin
          xstart = 5800.
          xfin = 6810.
        endelse 
        junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - xstart), start_location)
        junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - xfin), fin_location)
        Xcrop  = [start_location, fin_location]
        
        ;Mark emission lines
          junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - 5889.9), Mark_Na_D2)
          junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - 5895.9), Mark_Na_D1)
          junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - lines_r[0]), Mark_6312)
          junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - lines_r[1]), Mark_6716)
          junk = min(abs(DISPERSION_MASK_r_KEY[indgen(2098), center] - lines_r[2]), Mark_6731)
        
        ;Mark centroid guide lines
          pretty_r[indgen(2098), round(CENTER_r)] = 0. 
          pretty_r[indgen(2098), round(CENTER_r + .9 / R_J_per_pixel_r)] = 0.
          pretty_r[indgen(2098), round(CENTER_r - .9 / R_J_per_pixel_r)] = 0. 
        
        cgImage, pretty_r[Xcrop[0]:Xcrop[1],Ycrop[0]:Ycrop[1]], Title='DIS Red Channel' + strmid(sxpar(torus_header_r, 'DATE-OBS'), 0,10), /KEEP_ASPECT_RATIO, stretch = 1, /Scale, $
            minvalue = minvalue, maxvalue = maxvalue, Position = [0.0375, 0.175, 0.9, 0.85], /Window, /save 
        if float(sxpar(torus_header_r, 'DISPWC')) gt 6500. then $
          cgColorbar, title = 'Rayleighs', /Vertical, Range=[MinValue, MaxValue], Position = [!x.window[1]+.10, !y.window[0], !x.window[1]+.13, !y.window[1]] $    
        else cgColorbar, title = 'Rayleighs', /Vertical, Range=[MinValue, MaxValue], Position = [!x.window[1]+.11, !y.window[0], !x.window[1]+.14, !y.window[1]]    
        
        ;format x axis     
          step = 100.
          xtickv = step * indgen(13) / dispers_r  + (start_location - Xcrop[0])         
          xticks = n_elements(xtickv)-1
          xtickname =  step * indgen(13) + DISPERSION_MASK_r_KEY[start_location, center[start_location]]
          xtickname = strcompress(string(round(xtickname), format = '(I)'),/remove_all)                   
         AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.5, Xthick=3., /data, Xtitle = cgSymbol('Angstrom'),$
            xticklen = -.01, xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname)
         AXIS, XAXIS = 1, XSTYLE = 1, color=0., charsize=1.5, Xthick=3., /data, xticklen = -.01, xticks = xticks, xtickv = xtickv, XTICKFORMAT='(A1)'   
          
       ;format y axis  
          step = 1.
          ytickv = (step * indgen(19) -9) * (1./R_J_PER_PIXEL_R) + CENTER_R[Xcrop[0]] - Ycrop[0]
          yticks = n_elements(ytickv)-1
          ytickname = (step * indgen(19) -9)                 
          ytickname = strcompress(string(abs(ytickname), format = '(I)'),/remove_all) 
         AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.5, ythick=3., Ytitle = STRCOMPRESS('Jovian Radii'), /data, yticklen = -.005, $
            yticks = yticks, ytickv = ytickv, ytickname = ytickname
         AXIS, YAXIS = 1, color=0., charsize=1.5, ythick=3., /data, yticklen = -.005, yticks = yticks, ytickv = ytickv, YTICKFORMAT='(A1)'   
        
        ;Mark lines 
         cgarrow, mean([Mark_Na_D1, Mark_Na_D2]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_Na_D2 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.      
         cgarrow, mean([Mark_Na_D1, Mark_Na_D2]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_Na_D1 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.    
         cgtext, mean([Mark_Na_D1, Mark_Na_D2]) - start_location - 20, Ycrop[1] - Ycrop[0]+30., 'Na', /data, charsize = 1.2
         cgarrow, Mark_6312 - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_6312 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.       
         cgtext, Mark_6312 - start_location - 15, Ycrop[1] - Ycrop[0]+30., 'S!U++!N', /data, charsize = 1.2
         cgarrow, mean([Mark_6731, Mark_6716]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_6716 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.      
         cgarrow, mean([Mark_6731, Mark_6716]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_6731 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.    
         cgtext, mean([Mark_6731, Mark_6716]) - start_location - 15, Ycrop[1] - Ycrop[0]+30., 'S!U+!N', /data, charsize = 1.2
        
        ;Mark Dawn side and Dusk side        
         cgtext, .815, .22, 'Dusk', /Normal, orientation = 90., charsize = 1.2    
         cgtext, .815, .7, 'Dawn', /Normal, orientation = 90., charsize = 1.2              
     cgPS_Close 

     cgPS_Open, filename = Directory + strcompress('Processed\Figures\Pretty_'+strmid(sxpar(torus_header_b, 'DATE-OBS'), 0,10) + '_' + strcompress(j) + '_b.eps', /remove_all), /ENCAPSULATED, xsize = 8.5, ysize = 4. 
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.        
        xstart = 3600.
        xfin = 4220.
        minvalue = 0
        maxvalue = 50
        Ycrop = [ceil(mean(CENTER_B - max(distance_profiles_b.dusk_distance_array_b) / R_J_per_pixel_b)), $
                floor(mean(CENTER_B + max(distance_profiles_b.dawn_distance_array_b) / R_J_per_pixel_b))]
        junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - xstart), start_location, /NAN)
        junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - xfin), fin_location, /NAN)
        Xcrop  = [start_location, fin_location]
        
        ;Mark emission lines
          junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - lines_b[0]), Mark_3722, /NAN)
          junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - lines_b[1]), Mark_3726, /NAN)
          junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - lines_b[2]), Mark_3729, /NAN)
          junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - lines_b[3]), Mark_4069, /NAN)
          junk = min(abs(DISPERSION_MASK_b_KEY[indgen(2098), center] - lines_b[4]), Mark_4076, /NAN)
        
        ;mark centroid
          pretty_b[indgen(2098), round(CENTER_b)] = 0. 
          pretty_b[indgen(2098), round(CENTER_b + .9 / R_J_per_pixel_b)] = 0.
          pretty_b[indgen(2098), round(CENTER_b - .9 / R_J_per_pixel_b)] = 0. 

        ;That bright mask region is distracting, scale it back
          show_pretty_b = pretty_b
          show_pretty_b[*, 312:530] = show_pretty_b[*, 312:530] / 4.
          show_pretty_b[*, 312:329] = show_pretty_b[*, 312:329] / 4.
        
        cgImage, show_pretty_b[Xcrop[0]:Xcrop[1],Ycrop[0]:Ycrop[1]], Title='DIS Red Channel' + strmid(sxpar(torus_header_b, 'DATE-OBS'), 0,10), /KEEP_ASPECT_RATIO, stretch = 1, /Scale, $
            minvalue = minvalue, maxvalue = maxvalue, Position = [0.0375, 0.175, 0.9, 0.925], /Window, /save 
            
        cgColorbar, title = 'Rayleighs', /Vertical, Range=[MinValue, MaxValue], Position = [!x.window[1]+.13, !y.window[0], !x.window[1]+.16, !y.window[1]]    
        
        ;format x axis            
          step = 100.            
          xtickv = step * indgen(13) / dispers_b  + (start_location - Xcrop[0])         
          xticks = n_elements(xtickv)-1
          xtickname =  step * indgen(13) + DISPERSION_MASK_b_KEY[start_location, center[start_location]]
          xtickname = strcompress(string(round(xtickname), format = '(I)'),/remove_all)                   
         AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.5, Xthick=3., /data, Xtitle = cgSymbol('Angstrom'),$
            xticklen = -.01, xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname)
         AXIS, XAXIS = 1, XSTYLE = 1, color=0., charsize=1.5, Xthick=3., /data, xticklen = -.01, xticks = xticks, xtickv = xtickv, XTICKFORMAT='(A1)'   
          
       ;format y axis  
          step = 1.
          ytickv = (step * indgen(19) -9) * (1./R_J_PER_PIXEL_b) + CENTER_b[Xcrop[0]] - Ycrop[0]
          yticks = n_elements(ytickv)-1
          ytickname = (step * indgen(19) -9)                
          ytickname = strcompress(string(abs(ytickname), format = '(I)'),/remove_all) 
         AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.5, ythick=3., Ytitle = STRCOMPRESS('Jovian Radii'), /data, yticklen = -.005, $
            yticks = yticks, ytickv = ytickv, ytickname = ytickname
         AXIS, YAXIS = 1, color=0., charsize=1.5, ythick=3., /data, yticklen = -.005, yticks = yticks, ytickv = ytickv, YTICKFORMAT='(A1)'             
     
       ;Mark lines 
         cgarrow, Mark_3722 - start_location - 30, Ycrop[1] - Ycrop[0] + 25., Mark_3722 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.      
         cgtext, Mark_3722 - start_location - 50, Ycrop[1] - Ycrop[0]+30., 'S!U++!N', /data, charsize = 1.2
         
         cgarrow, mean([Mark_3726, Mark_3729]) - start_location + 30, Ycrop[1] - Ycrop[0] + 25., Mark_3726 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.    
         cgarrow, mean([Mark_3726, Mark_3729]) - start_location + 30, Ycrop[1] - Ycrop[0] + 25., Mark_3729 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.       
         cgtext, mean([Mark_3726, Mark_3729]) - start_location + 20, Ycrop[1] - Ycrop[0]+30., 'O!U+!N', /data, charsize = 1.2

         cgarrow, mean([Mark_4069, Mark_4076]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_4069 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.      
         cgarrow, mean([Mark_4069, Mark_4076]) - start_location, Ycrop[1] - Ycrop[0] + 25., Mark_4076 - start_location, Ycrop[1] - Ycrop[0], /solid, /data, hsize = 0.    
         cgtext, mean([Mark_4069, Mark_4076]) - start_location - 13, Ycrop[1] - Ycrop[0]+30., 'S!U+!N', /data, charsize = 1.2
        
       ;Mark Dawn side and Dusk side        
         cgtext, .71, .27, 'Dusk', /Normal, orientation = 90., charsize = 1.2    
         cgtext, .71, .75, 'Dawn', /Normal, orientation = 90., charsize = 1.2        
     cgPS_Close 
     set_plot,'WIN'
    SYSIII = CML[CoAdd_frames]
    Times = observing_times[CoAdd_frames]
    Dawn_distance_array_r = distance_profiles_r.Dawn_DISTANCE_ARRAY_R 
    Dusk_distance_array_r = distance_profiles_r.DUSK_DISTANCE_ARRAY_R 
    Dawn_distance_array_b = distance_profiles_b.Dawn_DISTANCE_ARRAY_b 
    Dusk_distance_array_b = distance_profiles_b.DUSK_DISTANCE_ARRAY_b 
    save, Dawn_distance_array_b, Dawn_avg_b, err_dawn_avg_b, Dawn_distance_array_r, Dawn_avg_r, err_dawn_avg_r, $
          Dusk_distance_array_b, Dusk_avg_b, err_dusk_avg_b, Dusk_distance_array_r, Dusk_avg_r, err_dusk_avg_r, SYSIII, Times, filename = 'C:\IDL\Io\Apache_Point_Programs\Profiles\'+ strmid(directory, 26, 8) + '_' + strcompress(j, /remove_all) + '.sav'
    if keyword_set(debug) then begin
     window, 0, xs = 2098, ys = 1000
     tv, bytscl(pretty_b, 0, 40)
     window, 1, xs = 2098, ys = 1000
     tv, bytscl(pretty_r, 0, 40)
    endif 
   endfor
   MWRFITS, float(pretty_r), strcompress(Directory + 'Processed\Pretty_r.fits'), torus_header_r, /CREATE, /silent ;/create overwrites
   MWRFITS, float(pretty_b), strcompress(Directory + 'Processed\Pretty_b.fits'), torus_header_b, /CREATE, /silent ;/create overwrites
   stop
endif          

if part eq 6.5 then begin        
    ;-------------------------------------------Stack Plots----(Frame by Frame SII Profiles)------------------------------------------------------    
          
   ;Get Jovian and other contants
     cspice_bodvrd, 'Jupiter', "RADII", 5, radii  
     R_J = radii[0]                  ;Jovian radii, in KM
   
   Images_r = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*r._Cal.fits'), count = n_images) 
   Images_b = FILE_SEARCH(strcompress(Directory + 'Processed\' + '*b._Cal.fits'), count = n_images)       
   SMOOTH_BY = 0.1 ;rj 
    
   ;get the headers and dates
   UT = strarr(N_elements(Images_r))
   for i = 0, N_elements(Images_r) - 1 do begin
          img = MRDFITS(Images_r[i], 0, torus_header, /silent)
          ET = double(sxpar(torus_header, 'EPHEM_TI'))
          cspice_et2utc, et, 'c', 0, utcstr
          UT[i] = utcstr
   endfor

    ;Find the height of the topmost profile in the stack plot
        img = MRDFITS(Images_r[-1], 0, torus_header, /silent)
        profiles = MRDFITS(Images_r[-1], 3, /silent) 
        y = [profiles.FIT_DAWN_LINE_PROFILES, profiles.FIT_Dusk_LINE_PROFILES]
        x = [profiles.DAWN_DISTANCE_ARRAY_R[*,0], profiles.DUSK_DISTANCE_ARRAY_R[*,0]]
        range = where((abs(x) gt 5.) and (abs(x) lt 7.), /Null)
        height = max(y[range,2])

   cgPS_Open, filename = Directory + strcompress('Processed\Figures\Stack_plot_'+strmid(sxpar(torus_header, 'DATE-OBS'), 0,10) + '.eps', /remove_all), /ENCAPSULATED, xsize = 16., ysize = 16.  
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.
        ;step = 200.
        step = 400.

        ;Setup 2 panels
        positions = cglayout([2,1], XGAP = 1., OYMargin=[7,5], OXMargin=[12,3] )
        
        ;Setup Axis
        ytickv = findgen(N_elements(Images_r)+1)*step
        ytickname = make_array(N_elements(Images_r)+1, value = '0, '+ strcompress(fix(step)))
        ytickname[0] = '0'
        ytickname[-1] = strcompress(fix(step)) 
        ytickname = [ ytickname, strcompress(fix(2*step)) ]      
        cgplot, [5.91, 5.91], [0, 1.], xrange = [7.7,4.5001], xstyle = 1, yrange = [0, (N_elements(Images_r)-1)*step + height], ytitle = 'Rayleighs', $
                xtitle = 'Dawn-side Distance (R!DJupiter!N)', psym = 10, thick = 6., /NODATA, position = positions[*,0], YStyle=1, ytickv=ytickv, ytickname = ytickname, YTICKINTERVAL = step
        y_range = !y.crange
        cgText, 7.45, y_range[1]*1.01, strmid(UT[0], 0,11) + ' ' + strmid(UT[0], 12,9) + ' - ' + strmid(UT[i-1], 12,9) + '   S!U+!N 6731' + cgsymbol('Angstrom') + ' (Black) & 6716' + cgsymbol('Angstrom') + ' (Gray)', charsize = 3.5 
        
        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, (N_elements(Images_r)-1)*step + height], /overplot, linestyle = 2.

        for i = 0, N_elements(Images_r) - 1 do begin
          img = MRDFITS(Images_r[i], 0, torus_header, /silent)
          err_img = MRDFITS(Images_r[i], 1, /silent)
          profiles = MRDFITS(Images_r[i], 3, /silent)  
          R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
          CML = float(sxpar(torus_header, 'CML_SYSI'))
          CML_Dawn = (CML + 90. + 360.0) MOD 360.0 
          x = profiles.DAWN_DISTANCE_ARRAY_R
          y = profiles.FIT_DAWN_LINE_PROFILES
          ERR_y = profiles.err_FIT_DAWN_LINE_PROFILES
          Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
          err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))

          cgplot, x[*,1], y[*,1] + i*step, /overplot, psym = 10, thick = 6., color = 'Charcoal', ERR_YLOW = err_y[*,1], ERR_YHIGH = err_y[*,1], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
          cgplot, x[*,2], y[*,2] + i*step, /overplot, psym = 10, thick = 6., color = 'Black', ERR_YLOW = err_y[*,2], ERR_YHIGH = err_y[*,2], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0. 
          cgText, 7.65, i*step + 100.,  cgSymbol('lambda') + '!DIII!N = ' + string(CML_dawn, FORMAT = '(F5.1)') + cgSymbol('deg')
          cgplot, [0.,10.], [i*step,i*step], /overplot, psym = 10, thick = 1., linestyle = 3
        endfor   
        ;        for i = 0, N_elements(Images_b) - 1 do begin
        ;          img = MRDFITS(Images_b[i], 0, torus_header, /silent)
        ;          err_img = MRDFITS(Images_b[i], 1, /silent)
        ;          profiles = MRDFITS(Images_b[i], 3, /silent)  
        ;          R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
        ;          CML = float(sxpar(torus_header, 'CML_SYSI'))
        ;          CML_Dawn = (CML + 90. + 360.0) MOD 360.0 
        ;          x = profiles.DAWN_DISTANCE_ARRAY_b
        ;          y = profiles.FIT_DAWN_LINE_PROFILES
        ;          ERR_y = profiles.err_FIT_DAWN_LINE_PROFILES
        ;          Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
        ;          err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
        ;          
        ;          cgplot, x, y[*,1] + i*step, /overplot, psym = 10, thick = 6., color = blue_colors[1], ERR_YLOW = err_y[*,1], ERR_YHIGH = err_y[*,1], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
        ;          cgplot, x, y[*,2] + i*step, /overplot, psym = 10, thick = 6., color = blue_colors[2], ERR_YLOW = err_y[*,2], ERR_YHIGH = err_y[*,2], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0. 
        ;        endfor   

        ;Setup Axis
        cgplot, [5.91, 5.91], [0, 1.e3], xrange = [4.5001,7.7], xstyle = 1, yrange = [0, (N_elements(Images_r)-1)*step + height], $
                xtitle = 'Dusk-side Distance (R!DJupiter!N)', psym = 10, thick = 6., /NODATA, YTickformat='(A1)', YStyle=1, Position = positions[*,1], /NoErase, /window

        ;label Io and label S&T 95 range
          cgplot, [5.91, 5.91], [0, (N_elements(Images_r)-1)*step + height], /overplot, linestyle = 2
          
        for i = 0, N_elements(Images_r) - 1 do begin
          img = MRDFITS(Images_r[i], 0, torus_header, /silent)
          err_img = MRDFITS(Images_r[i], 1, /silent)
          profiles = MRDFITS(Images_r[i], 3, /silent)  
          R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
          CML = float(sxpar(torus_header, 'CML_SYSI'))
          CML_Dusk = (CML - 90. + 360.0) MOD 360.0 
          x = profiles.dusk_DISTANCE_ARRAY_R
          y = profiles.FIT_dusk_LINE_PROFILES
          ERR_y = profiles.err_FIT_dusk_LINE_PROFILES
          Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
          err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
        
          cgplot, x[*,1], y[*,1] + i*step, /overplot, psym = 10, thick = 6., color = 'Charcoal', ERR_YLOW = err_y[*,1], ERR_YHIGH = err_y[*,1], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
          cgplot, x[*,2], y[*,2] + i*step, /overplot, psym = 10, thick = 6., color = 'Black', ERR_YLOW = err_y[*,2], ERR_YHIGH = err_y[*,2], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
          cgText, 7.0, i*step + 100.,  cgSymbol('lambda') + '!DIII!N = ' + string(CML_dusk, FORMAT = '(F5.1)') + cgSymbol('deg') 
          cgplot, [0.,10.], [i*step,i*step], /overplot, psym = 10, thick = 1., linestyle = 3    
          print, i, '  ', Images_r[i]
        endfor 

;        for i = 0, N_elements(Images_b) - 1 do begin
;          img = MRDFITS(Images_b[i], 0, torus_header, /silent)
;          err_img = MRDFITS(Images_b[i], 1, /silent)
;          profiles = MRDFITS(Images_b[i], 3, /silent)  
;          R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
;          CML = float(sxpar(torus_header, 'CML_SYSI'))
;          CML_Dusk = (CML - 90. + 360.0) MOD 360.0 
;          x = profiles.dusk_DISTANCE_ARRAY_b
;          y = profiles.FIT_dusk_LINE_PROFILES
;          ERR_y = profiles.err_FIT_dusk_LINE_PROFILES
;          Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
;          err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
;        
;          cgplot, x, y[*,1] + i*step, /overplot, psym = 10, thick = 6., color = blue_colors[1], ERR_YLOW = err_y[*,1], ERR_YHIGH = err_y[*,1], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
;          cgplot, x, y[*,2] + i*step, /overplot, psym = 10, thick = 6., color = blue_colors[2], ERR_YLOW = err_y[*,2], ERR_YHIGH = err_y[*,2], /err_clip, ERR_THICK = 0., ERR_WIDTH = 0.
;        endfor         

    cgPS_Close 
    set_plot,'WIN'

;
;    ;ribbon fitting parameters
;    parinfo = replicate( { fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 5)
;    parinfo[0].limited = 1b                                 ;limit amplitude 
;    parinfo[0].limits = [10., 1.e3]                         ;positive only between 10 and 1000 Rayleighs
;    parinfo[1].limited = 1b                                 ;limit centroid
;    parinfo[1].limits = [5.3, 6.2]                          ;centroid limits is system III dependent, allocated later within the loop
;    parinfo[2].limited = 1b                                 ;limit sigma   
;    parinfo[2].limits = [.06, 0.3]                         ;limit sigma width in R_J, 0.05 to 0.5
;    parinfo[3].fixed = 1b                                   ;limit DC Background to zero
;    parinfo[4].fixed = 1b                                   ;limit slope to zero
;    
;    ribbon_centroids = fltarr(N_elements(Images_r), 2)
;    err_ribbon_centroids = fltarr(N_elements(Images_r), 2)
;    sysIII = fltarr(N_elements(Images_r), 2)
;    for i = 0, N_elements(Images_r) - 1 do begin
;      img = MRDFITS(Images_r[i], 0, torus_header, /silent)
;      err_img = MRDFITS(Images_r[i], 1, /silent)
;      profiles = MRDFITS(Images_r[i], 3, /silent)  
;      R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
;      CML = float(sxpar(torus_header, 'CML_SYSI'))
;      sysIII[i,*] = [(CML + 90. + 360.0) MOD 360.0, (CML - 90. + 360.0) MOD 360.0] ;[dawn,dusk] ribbon location in R_J from Schneider & Trauger's (1995) fit
;      avg_ribbon_dawn = 5.85 + 0.049*cos((sysIII[i,0] - 167.)/!Radeg)
;      avg_ribbon_dusk = 5.57 + 0.073*cos((sysIII[i,1] - 130.)/!Radeg) 
;      
;      x = profiles.DAWN_DISTANCE_ARRAY_R
;      y = profiles.FIT_DAWN_LINE_PROFILES
;      ERR_y = profiles.err_FIT_DAWN_LINE_PROFILES
;      Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
;      err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
;      Y = total(Y[*,1:2], 2)
;      err_y = sqrt(err_y[*,1]^2. + err_y[*,2]^2.)
;
;      ;trim x to the valid region of interest   
;      y = y[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      ERR_y = ERR_y[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      x = x[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      parinfo[1].limits = [avg_ribbon_dawn - .4, avg_ribbon_dawn + .4]  ;centroid limits in RJ
;      parinfo[0].limits = [.2*max(y), 1000.]  ;centroid limits in RJ
;      cold_peak = where(y eq max(y[where(x lt avg_ribbon_dawn - 0.4)]), /null)
;      ribbon_range = where(abs(x - avg_ribbon_dawn) lt 0.4, /NULL, complement = quadratic_range)
;      COEFF = ROBUST_POLY_FIT(X[cold_peak:*], Y[cold_peak:*], 2, YFIT)
;      residual = y - poly(x, coeff)
;      err_y[where(err_y eq 0., /NULL)] = mean(err_y) ;Hack, why are their zeros in the error arrays?
;      initial_guess = [max(residual), avg_ribbon_dawn, .12, 0., 0.]
;      fa = {x:double(x), y:double(residual), err:double(err_y)} 
;      p = mpfit('Gaussian_Ribbon', initial_guess, PERROR = err_p, funct=fa, maxiter=20, STATUS = Did_it_work, parinfo = parinfo, /Quiet, NPEGGED = NPEGGED) 
;      if ((Did_it_work gt 0) and (NPEGGED lt 2)) then begin
;        ribbon_centroids[i,0] = p[1]  
;        err_ribbon_centroids[i,0] = err_p[1] 
;        print, sysIII[i,0], p
;        cgplot, x, y, color = 'red', psym = 10
;        cgplot, X[cold_peak:*], poly(x[cold_peak:*], coeff), /overplot, psym = 5
;        cgplot, X, poly(x, coeff) + gaussian(X, P[0:2]), /overplot
;      endif else print, 'No Ribbon Found'     
;
;      x = profiles.Dusk_DISTANCE_ARRAY_R
;      y = profiles.FIT_Dusk_LINE_PROFILES
;      ERR_y = profiles.err_FIT_Dusk_LINE_PROFILES
;      Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
;      err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
;      Y = total(Y[*,1:2], 2)
;      err_y = sqrt(err_y[*,1]^2. + err_y[*,2]^2.)
;
;      ;trim x to the valid region of interest   
;      y = y[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      ERR_y = ERR_y[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      x = x[where((abs(x) gt 4.5) and (abs(x) lt 8.1), /Null)]
;      parinfo[1].limits = [avg_ribbon_dusk - .4, avg_ribbon_dusk + .4]  ;centroid limits in RJ
;      parinfo[0].limits = [.2*max(y), 1000.]  ;centroid limits in RJ
;      cold_peak = where(y eq max(y[where(x lt avg_ribbon_dusk - 0.4)]), /null)
;      ribbon_range = where(abs(x - avg_ribbon_dusk) lt 0.4, /NULL, complement = quadratic_range)
;      COEFF = ROBUST_POLY_FIT(X[cold_peak[0]:*], Y[cold_peak[0]:*], 2, YFIT)
;      residual = y - poly(x, coeff)
;      err_y[where(err_y eq 0., /NULL)] = mean(err_y) ;Hack, why are their zeros in the error arrays?
;      initial_guess = [max(residual), avg_ribbon_dusk, .12, 0., 0.]
;      fa = {x:double(x), y:double(residual), err:double(err_y)}
;      p = mpfit('Gaussian_Ribbon', initial_guess, PERROR = err_p, funct=fa, maxiter=20, STATUS = Did_it_work, parinfo = parinfo, /Quiet, NPEGGED = NPEGGED) 
;      if ((Did_it_work gt 0) and (NPEGGED lt 2)) then begin
;        ribbon_centroids[i,1] = p[1]  
;        err_ribbon_centroids[i,1] = err_p[1] 
;        print, sysIII[i,1], p
;        cgplot, x, y, color = 'red', psym = 10
;        cgplot, X[cold_peak[0]:*], poly(x[cold_peak[0]:*], coeff), /overplot, psym = 5
;        cgplot, X, poly(x, coeff) + gaussian(X, P[0:2]), /overplot
;      endif else print, 'No Ribbon Found' 
;    endfor  
;    
;   cgPS_Open, Directory + strcompress('Processed\Figures\Ribbon_Pos_'+strmid(sxpar(torus_header, 'DATE-OBS'), 0,10) + '.eps', /remove_all), /ENCAPSULATED, xsize = 8, ysize = 6  
;        !P.font=1
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;        !p.charsize = 2.
;        
;      cgplot, [sysIII[*,0], sysIII[*,0]+360.], [ribbon_centroids[*,0], ribbon_centroids[*,0]], $
;              ERR_YLOW = [err_ribbon_centroids[*,0], err_ribbon_centroids[*,0]] + 1.*R_J_per_pixel, ERR_YHIGH = [err_ribbon_centroids[*,0], err_ribbon_centroids[*,0]] + 1.*R_J_per_pixel, $
;              ERR_WIDTH = 0., ERR_thick = .1, /ERR_CLIP, xtitle = 'Ansae (Earth LOS)' + cgsymbol('lambda')+'!DIII!N', ytitle = 'Ribbon Position (R!DJupiter!N)', psym = 4, color = 'red', $
;              yr = [5.4, 6.1], xr = [0, 720], Aspect = .66
;      cgplot, [sysIII[*,1], sysIII[*,1]+360.], [ribbon_centroids[*,1], ribbon_centroids[*,1]], $
;              ERR_YLOW = [err_ribbon_centroids[*,1], err_ribbon_centroids[*,1]]+ 1.*R_J_per_pixel, ERR_YHIGH = [err_ribbon_centroids[*,1], err_ribbon_centroids[*,1]] + 1.*R_J_per_pixel, $
;              ERR_WIDTH = 0., ERR_thick = .1,/ERR_CLIP, psym = 4, color = 'blue', /overplot
;      cgplot, findgen(721), 5.85 + 0.049*cos((findgen(721) - 167.)/!Radeg), color = 'red', /overplot
;      cgplot, findgen(721), 5.57 + 0.073*cos((findgen(721) - 130.)/!Radeg), color = 'blue', /overplot
;  cgPS_Close  

endif        

If part eq 6.75 then begin
  
  ;Setup CHIANTI for use with the line ratios
    use_chianti, 'C:\ssw\packages\chianti', abund='C:\ssw\packages\chianti\abundance\sun_photospheric_2011_caffau.abund'
  
;  ;  ;Test chianti vs my cloudy type program:
;  EMISS = EMISS_CALC(16, 2, TEMP = 4., DENS=alog10(1600.), /NO_DE)
;  WL = emiss.lambda
;  result = emiss.em
;  vactoair, WL[[810,1198,1200]], test_wavelengths
;  print, test_wavelengths
;  print, result[1200] / result[1198], result[1200] / result[810], result[1198] / result[810]
;  ;density_ratios,'s_2', 4060., 6740., alog10(1600.), alog10(1601.), dens, ratio, desc, temp=1.e+4, /PHOTONS
;  ;temperature_ratios,'s_2', 4000., 6740., 4., 4.01, temps, ratio, desc, density=1.6e+3, /PHOTONS
;  Print, 'Chianti: 6731/6716 = 1.2816, 6731/4069 = 10.5199, 6716/4069 = 8.2084'
;  result = line_ratio('S_II', dens = 1600., temp = 10000.)
;  emissivity = result.emissivity / result.linewave ;change from units of ergs to (relative) photon units 
;  print, emissivity[1,0] / emissivity[2,0] ;6731 / 6716 ratio
;  print, emissivity[1,0] / emissivity[4,0] ;6731 / 4069 ratio
;  print, emissivity[2,0] / emissivity[4,0] ;6716 / 4069 ratio
;  print, Temperature sensitive ratio in cloudy
;stop

  k_B = 1.38064852e-23 ;m2 kg s-2 K-1
  J_per_eV = 1.60218e-19 ;joules per electronvolt
  
  for j = 0, N_elements(CoAdd_frames[0,*]) - 1 do begin
      restore, 'C:\IDL\Io\Apache_Point_Programs\Profiles\'+ strmid(directory, 26, 8) + '_' + strcompress(j, /remove_all) + '.sav'  
      
;     ;trim to the region of the torus
;      furthest = min([max(Dawn_distance_array_r), max(Dusk_distance_array_r), max(Dawn_distance_array_b), max(Dusk_distance_array_b)])
;      ;range = where(((Dawn_distance_array_r gt 4.5) and (Dawn_distance_array_r lt furthest - .1)), /NULL)
;      range = where(((Dawn_distance_array_r gt 4.5) and (Dawn_distance_array_r lt 8.1)), /NULL)
;      Dawn_avg_r = Dawn_avg_r[min(range):max(range), *]
;      err_Dawn_avg_r = err_Dawn_avg_r[min(range):max(range), *]
;      dawn_distance = Dawn_distance_array_r[min(range):max(range)]
;      Dawn_avg_b_new = fltarr(size( Dawn_avg_b[0:N_elements(range), *], /dimensions )) 
;      err_Dawn_avg_b_new = fltarr(size( Dawn_avg_b[0:N_elements(range), *], /dimensions )) 
;      
;    ;interpolate the blue channel to the same spatial scale as the red. 
;      for i = 0, n_elements(Dawn_avg_b[0,*])-1 do begin
;          Dawn_avg_b_new[*,i] = interpol( Dawn_avg_b[*,i], Dawn_distance_array_b, dawn_distance )
;          err_Dawn_avg_b_new[*,i] = interpol( err_Dawn_avg_b[*,i], Dawn_distance_array_b, dawn_distance )
;      endfor
;      Dawn_avg_b = Dawn_avg_b_new
;      err_Dawn_avg_b = err_Dawn_avg_b_new
;   
;     ;trim to the region of the torus
;      range = where(((dusk_distance_array_r gt 4.5) and (dusk_distance_array_r lt 8.1)), /NULL)
;      dusk_avg_r = dusk_avg_r[min(range):max(range), *]
;      err_dusk_avg_r = err_dusk_avg_r[min(range):max(range), *]
;      dusk_distance = dusk_distance_array_r[min(range):max(range)]
;      dusk_avg_b_new = fltarr(size( dusk_avg_b[0:N_elements(range), *], /dimensions )) 
;      err_dusk_avg_b_new = fltarr(size( dusk_avg_b[0:N_elements(range), *], /dimensions )) 
;
;    ;interpolate the blue channel to the same spatial scale as the red. 
;      for i = 0, n_elements(dusk_avg_b[0,*])-1 do begin
;          dusk_avg_b_new[*,i] = interpol( dusk_avg_b[*,i], dusk_distance_array_b, dusk_distance )
;          err_dusk_avg_b_new[*,i] = interpol( err_dusk_avg_b[*,i], dusk_distance_array_b, dusk_distance )
;      endfor
;      dusk_avg_b = dusk_avg_b_new
;      err_dusk_avg_b = err_dusk_avg_b_new
;       
;       
;       ;A flawed plan
;       
;        ;      ;rebin the data, using variable bin size to only include line strengths with a SNR > 3     
;        ;        SNR = 3
;        ;        for i = 0, 2 do begin
;        ;          if i eq 0 then continue ;skip the 6312A line
;        ;          y = REGROUP(dusk_avg_r[*,i], Dusk_Distance, dy = err_dusk_avg_r[*,i], SNR = SNR, bin_x = x, BIN_dy = dy, BIN_map = huh)
;        ;          print, max( huh )
;        ;        endfor
;        ;        
;        ;        for i = 0, 4 do begin
;        ;          if i lt 3 then continue ;skip the other lines line
;        ;          y = REGROUP(dusk_avg_b[*,i], Dusk_Distance, dy = err_dusk_avg_b[*,i], SNR = SNR, bin_x = x, BIN_dy = dy, BIN_map = huh)
;        ;          print, max( huh )
;        ;        endfor
;        ;      ;To actually use these ratios, they all must have the same x coefficients, search and make this so.  
;        
;          
;      Dawn_A6716_over_A6731 = Dawn_avg_r[*,1] / Dawn_avg_r[*,2] 
;      Dawn_A4069_over_A6731 = Dawn_avg_b[*,3] / Dawn_avg_r[*,2] 
;      Dawn_A4069_over_A6716 = Dawn_avg_b[*,3] / Dawn_avg_r[*,1] 
;      dusk_A6716_over_A6731 = dusk_avg_r[*,1] / dusk_avg_r[*,2] 
;      dusk_A4069_over_A6731 = dusk_avg_b[*,3] / dusk_avg_r[*,2]
;      dusk_A4069_over_A6716 = dusk_avg_b[*,3] / dusk_avg_r[*,1]
;
;      err_Dawn_A4069_over_A6716 = abs(Dawn_A4069_over_A6716) * sqrt((err_Dawn_avg_r[*,1] / Dawn_avg_r[*,1])^2. + (err_Dawn_avg_b[*,3] / Dawn_avg_b[*,3])^2.) 
;      err_dusk_A4069_over_A6716 = abs(dusk_A4069_over_A6716) * sqrt((err_dusk_avg_r[*,1] / dusk_avg_r[*,1])^2. + (err_dusk_avg_b[*,3] / dusk_avg_b[*,3])^2.) 
;      err_Dawn_A4069_over_A6731 = abs(Dawn_A4069_over_A6731) * sqrt((err_Dawn_avg_r[*,2] / Dawn_avg_r[*,2])^2. + (err_Dawn_avg_b[*,3] / Dawn_avg_b[*,3])^2.) 
;      err_Dawn_A6716_over_A6731 = abs(Dawn_A6716_over_A6731) * sqrt((err_Dawn_avg_r[*,2] / Dawn_avg_r[*,2])^2. + (err_Dawn_avg_r[*,1] / Dawn_avg_r[*,1])^2.) 
;      err_dusk_A4069_over_A6731 = abs(dusk_A4069_over_A6731) * sqrt((err_dusk_avg_r[*,2] / dusk_avg_r[*,2])^2. + (err_dusk_avg_b[*,3] / dusk_avg_b[*,3])^2.)  
;      err_dusk_A6716_over_A6731 = abs(dusk_A6716_over_A6731) * sqrt((err_dusk_avg_r[*,2] / dusk_avg_r[*,2])^2. + (err_dusk_avg_r[*,1] / dusk_avg_r[*,1])^2.)
;
;      ;rebin the data using the smallest bin which allows the a signal to noise of >SNR
;        SNR = 3 
;        
;        ;Hack doesn't yet account for different dawn and dusk distance scales     
;        dens_ratio  = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        dens_dist   = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        temp_ratio1 = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        temp_dist1  = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        temp_ratio2 = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        temp_dist2  = make_array(6, n_elements(Dawn_Distance), /Float, value = !values.F_nan)
;        y = REGROUP(Dawn_A6716_over_A6731, Dawn_Distance, dy = err_Dawn_A6716_over_A6731, SNR = SNR, bin_x = x, BIN_dy = dy)
;        dens_ratio [0,  0:N_elements(y)-1] = y
;        dens_ratio [1:2,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        dens_dist  [0:2,0:N_elements(y)-1] = transpose([[x], [x], [x]])       
;        y = REGROUP(Dusk_A6716_over_A6731, Dusk_Distance, dy = err_Dusk_A6716_over_A6731, SNR = SNR, bin_x = x, BIN_dy = dy)
;        dens_ratio [3,  0:N_elements(y)-1] = y
;        dens_ratio [4:5,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        dens_dist  [3:5,0:N_elements(y)-1] = transpose([[x], [x], [x]])  
;        y = REGROUP(Dawn_A4069_over_A6731, Dawn_Distance, dy = err_Dawn_A4069_over_A6731, SNR = SNR, bin_x = x, BIN_dy = dy)
;        temp_ratio1[0,  0:N_elements(y)-1] = y
;        temp_ratio1[1:2,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        temp_dist1 [0:2,0:N_elements(y)-1] = transpose([[x], [x], [x]])        
;        y = REGROUP(Dusk_A4069_over_A6731, Dusk_Distance, dy = err_Dusk_A4069_over_A6731, SNR = SNR, bin_x = x, BIN_dy = dy)
;        temp_ratio1[3,  0:N_elements(y)-1] = y
;        temp_ratio1[4:5,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        temp_dist1 [3:5,0:N_elements(y)-1] = transpose([[x], [x], [x]])
;        y = REGROUP(Dawn_A4069_over_A6716, Dawn_Distance, dy = err_Dawn_A4069_over_A6716, SNR = SNR, bin_x = x, BIN_dy = dy)
;        temp_ratio2[0,  0:N_elements(y)-1] = y
;        temp_ratio2[1:2,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        temp_dist2 [0:2,0:N_elements(y)-1] = transpose([[x], [x], [x]])        
;        y = REGROUP(Dusk_A4069_over_A6716, Dusk_Distance, dy = err_Dusk_A4069_over_A6716, SNR = SNR, bin_x = x, BIN_dy = dy)
;        temp_ratio2[3,  0:N_elements(y)-1] = y
;        temp_ratio2[4:5,0:N_elements(y)-1] = transpose([[dy], [dy]])
;        temp_dist2 [3:5,0:N_elements(y)-1] = transpose([[x], [x], [x]])
;
;        if keyword_set(debug) then begin
;          cgplot, Dawn_Distance, Dawn_A4069_over_A6731, ERR_YLOW = err_dawn_A4069_over_A6731, ERR_YHIGH = err_dawn_A4069_over_A6731, ERR_THICK = 2., ERR_WIDTH = 0., thick =4
;          cgplot, temp_dist1[0,*], temp_ratios[0,*], psym =10, ERR_YLOW = temp_ratios[2,*], ERR_YHIGH = temp_ratios[1,*], ERR_THICK = .5, ERR_WIDTH = 0., /overplot, color = 'red'
;        endif
;
;    ;Plot the line ratios
;      cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Profiles\Ratios--'+strmid(directory, 26, 8)+'.eps', /remove_all), /ENCAPSULATED, xsize = 16., ysize = 6.5
;        !P.font=1
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;        !p.charsize = 2.
;        
;        yrange = [0,1.6]
;
;        ;Setup 2 panels
;        positions = cglayout([2,1], XGAP = 1., OYMargin=[7,5], OXMargin=[8,8] )
;    
;        ;Setup Axis
;        cgplot, dawn_distance, Dawn_A4069_over_A6731, xrange = [8.1, 4.5], xstyle = 1, yrange = yrange, psym = 4, ystyle =1., ytitle = 'Measured Line Ratio', $
;                xtitle = 'Dawn-side Distance (R!DJupiter!N)', /NODATA, position = positions[*,0], /NoErase
;                
;        ;label Io and label S&T 95 range
;          CML_Dawn = (SYSIII + 90. + 360.0) MOD 360.0 
;          cgplot, [5.91, 5.91], yrange, /overplot, linestyle = 2.
;          ribbon_range = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
;          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [yrange[0], yrange[1], yrange[1], yrange[0]], /data, COLOR='yellow'
;        
;        cgText, 5.7, 1.4, '6716' + cgSymbol('Angstrom') + '/6731' +  cgSymbol('Angstrom'), color = 'Blue'
;        cgText, 5.7, 1.25,  '4069' + cgSymbol('Angstrom') + '/6716'+ cgSymbol('Angstrom'), color = 'red'
;        cgText, 5.7, 1.1,'4069' + cgSymbol('Angstrom') + '/6731' + cgSymbol('Angstrom'), color = 'dark green'
;                 
;;        cgplot, dawn_distance, dawn_A4069_over_A6731, ERR_YLOW = err_dawn_A4069_over_A6731, ERR_YHIGH = err_dawn_A4069_over_A6731, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='dark green'     
;;        cgplot, dawn_distance, dawn_A6716_over_A6731, ERR_YLOW = err_dawn_A6716_over_A6731, ERR_YHIGH = err_dawn_A6716_over_A6731, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='blue'
;;        cgplot, dawn_distance, dawn_A4069_over_A6716, ERR_YLOW = err_dawn_A4069_over_A6716, ERR_YHIGH = err_dawn_A4069_over_A6716, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='red' 
;    
;         cgplot, dens_dist [0,*], dens_ratio[0,*], ERR_YLOW = dens_ratio[1,*], ERR_YHIGH = dens_ratio[2,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='blue'     
;         cgplot, temp_dist1[0,*],temp_ratio1[0,*], ERR_YLOW =temp_ratio1[1,*], ERR_YHIGH =temp_ratio1[2,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='dark green'    
;         cgplot, temp_dist2[0,*],temp_ratio2[0,*], ERR_YLOW =temp_ratio2[1,*], ERR_YHIGH =temp_ratio1[2,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='red'    
;        
;        ;Setup Axis
;        cgplot, dusk_distance, dusk_A4069_over_A6731, xrange = [4.5, 8.1], xstyle = 1, yrange = yrange, psym = 4, ystyle =1., YTickformat='(A1)', $
;                xtitle = 'Dusk-side Distance (R!DJupiter!N)', /NODATA, position = positions[*,1], /NoErase        
;         
;        ;label Io and label S&T 95 range
;          CML_Dusk = (SYSIII - 90. + 360.0) MOD 360.0 
;          cgplot, [5.91, 5.91], yrange, /overplot, linestyle = 2.
;          ribbon_range = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
;          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [yrange[0], yrange[1], yrange[1], yrange[0]], /data, COLOR='yellow'
;
;;        cgplot, dusk_distance, dusk_A4069_over_A6731, ERR_YLOW = err_dusk_A4069_over_A6731, ERR_YHIGH = err_dusk_A4069_over_A6731, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='dark green'     
;;        cgplot, dusk_distance, dusk_A6716_over_A6731, ERR_YLOW = err_dusk_A6716_over_A6731, ERR_YHIGH = err_dusk_A6716_over_A6731, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='blue'
;;        cgplot, dusk_distance, dusk_A4069_over_A6716, ERR_YLOW = err_dusk_A4069_over_A6716, ERR_YHIGH = err_dusk_A4069_over_A6716, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 10, thick=3, COLOR='red' 
;        
;        cgplot, dens_dist [3,*], dens_ratio[3,*], ERR_YLOW = dens_ratio[4,*], ERR_YHIGH = dens_ratio[5,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='blue'     
;        cgplot, temp_dist1[3,*],temp_ratio1[3,*], ERR_YLOW =temp_ratio1[4,*], ERR_YHIGH =temp_ratio1[5,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='dark green'    
;        cgplot, temp_dist2[3,*],temp_ratio2[3,*], ERR_YLOW =temp_ratio2[4,*], ERR_YHIGH =temp_ratio1[5,*], /err_clip, ERR_THICK = 1., ERR_WIDTH = 0., /overplot, psym = 16, thick=3, COLOR='red'  
;      cgPS_Close     
;
;
;  temp_ratios = [[Dawn_A4069_over_A6731], $ ;arrays all the temp ratios 
;                 [Dawn_A4069_over_A6731 - err_Dawn_A4069_over_A6731], $
;                 [Dawn_A4069_over_A6731 + err_Dawn_A4069_over_A6731], $
;                 [Dusk_A4069_over_A6731], $
;                 [Dusk_A4069_over_A6731 - err_dusk_A4069_over_A6731], $ 
;                 [Dusk_A4069_over_A6731 + err_dusk_A4069_over_A6731]] 
;                 
;  temp_ratios2 = [[Dawn_A4069_over_A6716], $ ;arrays all the temp ratios 
;                 [Dawn_A4069_over_A6716 - err_Dawn_A4069_over_A6716], $
;                 [Dawn_A4069_over_A6716 + err_Dawn_A4069_over_A6716], $
;                 [Dusk_A4069_over_A6716], $
;                 [Dusk_A4069_over_A6716 - err_dusk_A4069_over_A6716], $ 
;                 [Dusk_A4069_over_A6716 + err_dusk_A4069_over_A6716]] 
;                 
;  den_ratios =  [[Dawn_A6716_over_A6731], $  ;arrays all the density ratios
;                 [Dawn_A6716_over_A6731 - err_Dawn_A6716_over_A6731], $
;                 [Dawn_A6716_over_A6731 + err_Dawn_A6716_over_A6731], $
;                 [Dusk_A6716_over_A6731], $
;                 [Dusk_A6716_over_A6731 - err_dusk_A6716_over_A6731], $ 
;                 [Dusk_A6716_over_A6731 + err_dusk_A6716_over_A6731]] 
;;        
;
;      ;setup parameter space distribution in density and temperature
;        space = 40. ; size of either axis of the density temp parameter space
;        ;core = asin(1.0-2.0*findgen(space)/space)  ;gives random numbers evenly distributed between -90 and 90 degrees (in radians) 
;        ;core = (core + (!Pi/2.)) / (!pi/2.)
;        ;test_dens = core*2500. + 10.
;        ;test_temps = core*40000. + 2000.
;        ;test_temps = test_temps[sort(test_temps)]
;        ;test_dens = test_dens[sort(test_dens)]        
;        test_dens = 10.^ ((findgen(space+1)/space)*3 + 1.5) ;100 - 10,000 per cc 
;        test_temps= 10.^ ((findgen(space+1)/space)*2 + 3.5) ;~2,000 - 200,000 Kelvin
;        Chianti_dens_ratio = fltarr(n_elements(test_dens)) 
;        Chianti_temp_ratio = fltarr(n_elements(test_temps))
;        window, 0, xs = 500, ys = 400
;        cghistoplot, test_dens
;        
;        window, 1, xs = 500, ys = 400
;        cghistoplot, test_temps * K_b / J_per_Ev 
;
;  ;the arrays that results will be written to:
;    N_e_arr = den_ratios & T_e_arr = den_ratios
;
;    ;CLOUDY type method
;    print, 'Building CLOUDY N_e vs T_e parameter Space. This will take a while...', Systime()  
;    dens_table = fltarr(N_elements(test_dens),N_elements(test_temps)) 
;    temp_table = fltarr(N_elements(test_dens),N_elements(test_temps)) 
;    temp_table2 = fltarr(N_elements(test_dens),N_elements(test_temps))
;    ;make a 2d Lookup table for each ratio
;    for ind_i = 0, N_elements(test_dens) - 1 do begin
;      for ind_j = 0, N_elements(test_temps) - 1 do begin
;        result = line_ratio('S_II', dens = test_dens[ind_i], temp = test_temps[ind_j])
;        emissivity = result.emissivity / result.linewave ;change from units of ergs to (relative) photon units 
;        dens_table[ind_i, ind_j] = emissivity[2,0] / emissivity[1,0]  ;6716 / 6731 ratio
;        temp_table[ind_i, ind_j] = emissivity[4,0] / emissivity[1,0] ;4069/6731 ratio
;        temp_table2[ind_i, ind_j] =  emissivity[4,0] / emissivity[2,0]  ;4069/6716 ratio
;      endfor
;    endfor
;    window, 0
;    tv, bytscl(dens_table, 0, 2)
;    window, 1
;    tv, bytscl(temp_table, 2, 5)
;    window, 2
;    tv, bytscl(temp_table2, 2, 5)
;    save, temp_table, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_temp_table.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6731A/4069A SII lines'
;    save, dens_table, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_dens_table.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6731A/6716A SII lines'   
;    save, temp_table2, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_temp_table2.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6716A/4069A SII lines'   
;    print, Systime() 
;
;
;    ;Chianti type method
;    print, 'Building CHIANTI N_e vs T_e parameter Space. This will take an even longer while...', Systime()  
;    dens_table = fltarr(N_elements(test_dens),N_elements(test_temps)) 
;    temp_table = fltarr(N_elements(test_dens),N_elements(test_temps)) 
;    temp_table2 = fltarr(N_elements(test_dens),N_elements(test_temps))
;    ;make a 2d Lookup table for each ratio
;    for ind_i = 0, N_elements(test_dens) - 1 do begin
;      for ind_j = 0, N_elements(test_temps) - 1 do begin
;        result = EMISS_CALC(16, 2, TEMP = alog10(test_temps[ind_j]), DENS=alog10(test_dens[ind_i]), /NO_DE, /QUIET)
;        emissivity = result.em
;        dens_table[ind_i, ind_j] =  emissivity[1198] / emissivity[1200] ;6731 / 6716 ratio
;        temp_table[ind_i, ind_j] =  emissivity[810] / emissivity[1200]  ;6731 / 4069 ratio
;        temp_table2[ind_i, ind_j] = emissivity[810] / emissivity[1198]  ;6716 / 4069 ratio
;      endfor    
;      print, Systime()  
;    endfor
;       
;    window, 0
;    tv, bytscl(dens_table, 0, 2)
;    window, 1
;    tv, bytscl(temp_table, 2, 5)
;    window, 2 
;    tv, bytscl(temp_table2, 2, 5)
;    save, temp_table, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_temp_table.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6731A/4069A SII lines'
;    save, dens_table, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_dens_table.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6731A/6716A SII lines'   
;    save, temp_table2, test_dens, test_temps, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_temp_table2.sav', $
;        Description = 'X-axis = test_dens, Y-Axis = test_temps, Value = Photon emissivity ratio of 6716A/4069A SII lines'   


    ;Test chianti vs my cloudy type program:
    CHIANTI = 1
    if keyword_set(chianti) then begin
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_dens_table.sav'
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_temp_table.sav'
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_temp_table2.sav'
    endif else begin
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_dens_table.sav'
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_temp_table.sav'
      restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CLOUDY_temp_table2.sav'
    endelse
        
    cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Profiles\Figure_F.eps', /remove_all), /ENCAPSULATED, xsize = 6, ysize = 4.5
      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 2.
      
       image = dens_table
       overlay = temp_table
       overlay2 = dens_table
     
       minValue = .2 ;Floor(Min(image))
       maxValue = 1.6 ;Ceil(Max(image))
       minOverlayValue = Floor(Min(overlay))
       maxOverlayValue = Ceil(Max(overlay))
       nLevels = 10
       xtitle = 'Log(N!De-!N) cm!U-3!N'
       ytitle = 'Log(T!De-!N) K'
       position =   [0.05, 0.225, 0.90, 0.975]
       cbposition = [0.85, 0.225, 0.9, 0.975]
       ;cbTitle = '6731' + cgsymbol('Angstrom') + ' / 6716' + cgsymbol('Angstrom') + ' Ratio'
       cbTitle = '6716' + cgsymbol('Angstrom') + ' / 6731' + cgsymbol('Angstrom') + ' Ratio'
       
       ; Set up colors for contour plot.
       cgLoadCT, 33, CLIP=[30,255]
       
       ; Display the image on the display. Keep its aspect ratio.
       cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, XTitle=xtitle, YTitle=ytitle, /AXES, /fit_inside, position = position, $
           XRange=[min(alog10(test_dens)),Max(alog10(test_dens))], YRange=[min(alog10(test_temps)),Max(alog10(test_temps))], /Keep_Aspect

       ;contourLevels = [.0001,.305d] ;ribbon ratio 
       contourLevels = [.0001, 0.313561d] ;ribbon ratio 
       cgContour, overlay, Levels=contourLevels, /OnImage, Color='charcoal', C_LABELS = [1,0], C_LINESTYLE = 1 
       ;contourLevels2 = [.0001,0.7561d] ;ribbon ratio 
       contourLevels2 = [.0001, 0.760868d] ;ribbon ratio 
       cgContour, overlay2, Levels=contourLevels2, /OnImage, Color='charcoal', C_LABELS = [1,0], C_LINESTYLE = 1 
  
       ;contourLevels = [.0001,.001,.01,.05,.15,.3,.5,.75,1., 1.5]
       contourLevels = [.0001,.001,.01,.05,.15,.25,.35,.5,.75, 1., 1.5]
       contourLevels = string(contourLevels, Format = '(F10.3)')
       cgContour, overlay, Levels=contourLevels, /OnImage, Color='charcoal'
;       cgText, 1.6, 4.8, '4069' + cgsymbol('Angstrom') + ' / 6731' + cgsymbol('Angstrom'), Color='charcoal', charsize =1.9 
;       cgText, 1.6, 4.65, 'contours', Color='charcoal', charsize =1.9 
       cgText, 1.6, 5.3, '4069' + cgsymbol('Angstrom') + ' / 6731' + cgsymbol('Angstrom'), Color='charcoal', charsize =1.9 
       cgText, 1.6, 5.15, 'contours', Color='charcoal', charsize =1.9  
       
       ; Draw the color bar. 
       cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
          Title=cbTitle, TLocation='Right', /Vertical  
    cgPS_Close
stop

    ;bracket the values coming in to the line fits so that they only lie in the allowable range of line ratios
    ;den_ratios = min(dens_table) > den_ratios < max(dens_table)
    ;temp_ratios = min(temp_table) > temp_ratios < max(temp_table)

    ;fix the x,y guess index:
;    junk = min(abs(test_temps - 10000.), temp_guess)
;    junk = min(abs(test_dens  - 1600.), dens_guess)
    junk = min(abs(test_temps - 40000.), temp_guess)
    junk = min(abs(test_dens  - 2600.), dens_guess)

    ;new method
    window, 0, xs = 600, ys = 400
    window, 1, xs = 600, ys = 400
    window, 2, xs = 600, ys = 400
    window, 3, xs = 600, ys = 400
    N_loc = make_array(n_elements(N_e_arr[*,0]), 6, value = dens_guess)
    T_loc = make_array(n_elements(N_e_arr[*,0]), 6, value = temp_guess)
    T2_loc= make_array(n_elements(N_e_arr[*,0]), 6, value = temp_guess)
    T_Mean= make_array(n_elements(N_e_arr[*,0]), 6, value = temp_guess)
    T_Mean_Mean = make_array(n_elements(N_e_arr[*,0]), 6, value = 0)
    N_Mean_Mean = make_array(n_elements(N_e_arr[*,0]), 6, value = 0)
    
   READCOL,'C:\IDL\Io\Bagenal_94_Fig1b.txt', F='A,A', Radius_bagenal, T_e_bagenal, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
      T_e_bagenal = float(T_e_bagenal) / (K_b / J_per_Ev)  ;convert from eV to Kelvin 
      Radius_bagenal = float(Radius_bagenal)  
    
    n_iterations = 15.
    for iter = 0, n_iterations - 1 do begin
      old = N_loc[*,0]
        for ind_i = 0, 5 do begin  
          for ind_j = 0, n_elements(N_e_arr[*,0]) - 1 do begin

           ;find nearest dens 6731/6716
            junk = min(abs(dens_table[*, T_Mean[ind_j,ind_i]] - den_ratios[ind_j, ind_i]), N_loc_j)
            n_loc[ind_j, ind_i] = N_loc_j
            
           ;find nearest temp 6731/6716
            junk = min(abs(temp_table[N_loc[ind_j,ind_i], *] - temp_ratios[ind_j, ind_i]), T_loc_j)
            T_loc[ind_j, ind_i] = T_loc_j
            
           ;find nearest temp 6716/4069
            junk = min(abs(temp_table2[N_loc[ind_j,ind_i], *] - temp_ratios2[ind_j, ind_i]), T_loc_j)
            T2_loc[ind_j, ind_i] = T_loc_j
    
          endfor
          wset, 0
          cgplot, Dawn_distance, test_dens[N_loc[*,0]], xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_dens), max(test_dens)], $
              ytitle = 'Electron Density (cm!U-3!N) (Line of sight average)'
          wset, 1
;          cgplot, Dawn_distance, test_temps[T_loc[*,0]], xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_temps), max(test_temps)], $
;              ytitle = 'Electron Temp (K) (Line of sight average)'     
;          cgplot, Dawn_distance, test_temps[T2_loc[*,0]], linestyle = 2., /overplot      
           cgplot, Dawn_distance, test_temps[T_Mean[*,0]], xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_temps), max(test_temps)], $
             ytitle = 'Electron Temp (K) (Line of sight average)'   
        endfor       
        ;T_Mean = make_array(n_elements(N_e_arr[*,0]), 6, value = temp_guess) 
        T_Mean = round((T_loc + T2_loc) / 2.)
        T_Mean_Mean = T_Mean_Mean + T_Mean
        N_Mean_Mean = N_Mean_Mean + N_Loc
        ;wait, 1.5
        print, total(abs(N_loc[*,0] - old)) ;test for convergence
    endfor
    T_Mean_Mean = T_Mean_Mean / n_iterations
    N_Mean_Mean = N_Mean_Mean / n_iterations
    wset, 2
    cgplot, Dawn_distance, test_dens[N_Mean_Mean[*,0]], xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_dens), max(test_dens)], $
       ytitle = 'Electron Density (cm!U-3!N) (Line of sight average)'
    wset, 3
    cgplot, Dawn_distance, test_temps[T_Mean_mean[*,0]], xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_temps), max(test_temps)], $
       ytitle = 'Electron Temp (K) (Line of sight average)'   
    
    ;allocate the final answers
    for ind_i = 0, 5 do begin  
      T_e_arr[*,ind_i] = test_temps[T_Mean_Mean[*,ind_i]]
      N_e_arr[*,ind_i] = test_dens[N_Mean_Mean[*,ind_i]]
    endfor

      ;Above code calculated them for high and low line ratios + err. Get the error back without this offset:

      N_e_arr[*,1] = abs(N_e_arr[*,1] - N_e_arr[*,0]) 
      N_e_arr[*,2] = abs(N_e_arr[*,2] - N_e_arr[*,0])
      T_e_arr[*,1] = abs(T_e_arr[*,1] - T_e_arr[*,0])
      T_e_arr[*,2] = abs(T_e_arr[*,2] - T_e_arr[*,0])
      N_e_arr[*,4] = abs(N_e_arr[*,4] - N_e_arr[*,3]) 
      N_e_arr[*,5] = abs(N_e_arr[*,5] - N_e_arr[*,3])
      T_e_arr[*,4] = abs(T_e_arr[*,4] - T_e_arr[*,3])
      T_e_arr[*,5] = abs(T_e_arr[*,5] - T_e_arr[*,3])
      T_e_arr = T_e_Arr * K_b / J_per_Ev 
      
      ;window, 0, xs = 600, ys = 400
      ;plot, Dawn_distance, N_e_arr[*,0], xrange = [4.5, 8.1]
      ;oplot,  Dawn_distance, N_e_arr[*,0] + N_e_arr[*,1], linestyle = 2.
      ;oplot,  Dawn_distance, N_e_arr[*,0] + N_e_arr[*,2], linestyle = 3

      N_range = [100,4500]
      cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Profiles\Ne_and_Te_CHIANTI--'+strmid(directory, 26, 8)+'.eps', /remove_all), /ENCAPSULATED, xsize = 16., ysize = 6.5
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 2.

        ;Setup 2 panels
        positions = cglayout([2,1], XGAP = 1., OYMargin=[7,5], OXMargin=[10,6] )
    
        ;Setup Axis
        cgplot, Dawn_distance, N_e_arr[*,0], xrange = [8.1, 4.5], xstyle = 1, yrange = N_range, ytitle = 'N!De-!N  (cm!U-3!N) (Line of sight average)', ystyle =9, $
                xtitle = 'Dawn-side Distance (R!DJupiter!N)', /NODATA, position = positions[*,0], /NoErase     
        ;label Io and label S&T 95 range
          CML_Dawn = (SYSIII + 90. + 360.0) MOD 360.0 
          cgplot, [5.91, 5.91], N_range, /overplot, linestyle = 2.
          ribbon_range = 5.85 + 0.049*cos((CML_Dawn - 167.)/!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [N_range[0],N_range[1],N_range[1],N_range[0]], /data, COLOR='yellow'
        cgplot, Dawn_distance, N_e_arr[*,0], err_yhigh = N_e_arr[*,2], err_ylow = N_e_arr[*,1], /overplot, psym = 10, thick = 4., /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.   
        cgAxis, YAxis=1.0, /Save, Color='Red', YTickformat = '(A1)', YRange = [0, 7], /ystyle
        cgOPlot, Dawn_distance, T_e_arr[*,0], err_yhigh = T_e_arr[*,1], err_ylow = T_e_arr[*,2], Color='Red', psym = 10, thick = 4., /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
    
        ;Setup Axis
        cgplot, dusk_distance, N_e_arr[*,3], xrange = [4.5, 8.1], xstyle = 1, yrange = N_range, ystyle =9, YTickformat='(A1)', $
                xtitle = 'Dusk-side Distance (R!DJupiter!N)', /NODATA, position = positions[*,1], /NoErase
        ;label Io and label S&T 95 range
          CML_Dusk = (SYSIII - 90. + 360.0) MOD 360.0 
          cgplot, [5.91, 5.91], N_range, /overplot, linestyle = 2.
          ribbon_range = 5.57 + 0.073*cos((CML_Dusk - 130.)/!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
          cgColorFill, [min(ribbon_range), min(ribbon_range), max(ribbon_range), max(ribbon_range)], [N_range[0],N_range[1],N_range[1],N_range[0]], /data, COLOR='yellow'
        cgplot, Dusk_distance, N_e_arr[*,3], err_yhigh = N_e_arr[*,5], err_ylow = N_e_arr[*,4], /overplot, psym = 10, thick = 4., /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
        cgAxis, YAxis=1.0, /Save, Color='red', YTitle='T!De-!N (eV) (Line of sight average)', YRange=[0, 7], /ystyle
        cgOPlot, Dusk_distance, T_e_arr[*,3], err_yhigh = T_e_arr[*,4], err_ylow = T_e_arr[*,5], Color='red', psym = 10, thick = 4, /err_clip, ERR_THICK = 1., ERR_WIDTH = 0.
      cgPS_Close     
  endfor
 stop
endif  


















;   WIDTH = 2.5                ;search / integration width in angstroms
;     
;   
;   ;define the arrays to hold the data products [value, error value]
;   Ribbon_brightness_Dawn_6731 = fltarr(n_elements(Images_r),2)
;   Ribbon_brightness_Dusk_6731 = fltarr(n_elements(Images_r),2)
;   Ribbon_brightness_Dawn_6716 = fltarr(n_elements(Images_r),2)
;   Ribbon_brightness_Dusk_6716 = fltarr(n_elements(Images_r),2)
;   Ribbon_brightness_Dawn_6312 = fltarr(n_elements(Images_r),2)
;   Ribbon_brightness_Dusk_6312 = fltarr(n_elements(Images_r),2)
;   warmT_brightness_Dawn_6731 = fltarr(n_elements(Images_r),2)
;   warmT_brightness_Dusk_6731 = fltarr(n_elements(Images_r),2)
;   warmT_brightness_Dawn_6716 = fltarr(n_elements(Images_r),2)
;   warmT_brightness_Dusk_6716 = fltarr(n_elements(Images_r),2)
;   
;   Ribbon_brightness_Dawn_4069 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dusk_4069 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dawn_4076 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dusk_4076 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dawn_3726 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dusk_3726 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dawn_3729 = fltarr(n_elements(Images_b),2)
;   Ribbon_brightness_Dusk_3729 = fltarr(n_elements(Images_b),2)
;
;   dawn_SYSIII= fltarr(n_elements(Images_r))
;   dusk_SYSIII= fltarr(n_elements(Images_r))
;   dawn_ratio = fltarr(n_elements(Images_r),2) 
;   dusk_ratio = fltarr(n_elements(Images_r),2)  
;   warmT_dawn_ratio = fltarr(n_elements(Images_r),2)
;   warmT_dusk_ratio = fltarr(n_elements(Images_r),2)  
;   dawn_Temp_ratio = fltarr(min([n_elements(Images_r), n_elements(Images_b)]),2) 
;   dusk_Temp_ratio = fltarr(min([n_elements(Images_r), n_elements(Images_b)]),2)  
   
   ;---------------------------------------------Ribbon Coordinates & Line Distance Profiles----------------------------------------------
;    ;Compare to expected distance 
;    phi = 4.*!pi*findgen(100) / 100.  
;    dawn_ribbon = 5.85 + 0.049*cos(phi - 167./!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
;    dusk_ribbon = 5.57 + 0.073*cos(phi - 130./!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
;    window, 1
;    plot, phi*!radeg, dawn_ribbon, yr = [5.4,6.], xr = [0.,720.], xstyle = 1
;    oplot, phi*!radeg, dusk_ribbon
;
;   ;------------------------Blue Data Products---------------------------------
;   for i = 0, n_elements(Images_b)-1 do begin
;      image_b = MRDFITS(Images_b[i], 0, torus_header_b, /Dscale, /silent)
;      metadata_b = MRDFITS(Images_b[i], 1)
;      Dispersion_mask_b = metadata_b.dispersion_mask     
;       
;      CML = float(sxpar(torus_header_b, 'CML_SYSI'))
;      RJ_PER_Pixel = float(sxpar(torus_header_b, 'RJ_PER_P'))
;      dawn_ribbon_index_b = metadata_b.center + ((5.85 + 0.049*cos(((CML - 90.) - 167.)/!Radeg)) / RJ_PER_Pixel)
;      dusk_ribbon_index_b = metadata_b.center - ((5.57 + 0.073*cos(((CML + 90.) - 130.)/!Radeg)) / RJ_PER_Pixel)     
;      dawn_shift = -(2. * !pi * (5.85 + 0.049*cos(((CML - 90.) - 167.)/!Radeg)) * R_J / period_J) / c ; =  delta lambda / lambda at the ribbon position in Schneider & Trauger, 1995
;      dusk_shift = (2. * !pi * (5.57 + 0.073*cos(((CML + 90.) - 130.)/!Radeg)) * R_J / period_J) / c  ; =  delta lambda / lambda at the ribbon position in Schneider & Trauger, 1995
;
;      SII_4069_dawn = where((Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] gt 4068.60 + 4068.60*dawn_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] lt 4068.60 + 4068.60*dawn_shift + width))
;      SII_4069_dusk = where((Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] gt 4068.60 + 4068.60*dusk_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] lt 4068.60 + 4068.60*dusk_shift + width))
;      SII_4076_dawn = where((Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] gt 4076.35 + 4076.35*dawn_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] lt 4076.35 + 4076.35*dawn_shift + width))
;      SII_4076_dusk = where((Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] gt 4076.35 + 4076.35*dusk_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] lt 4076.35 + 4076.35*dusk_shift + width))
;      ;width = 1.5 ;for the very close together lines
;      OII_3726_dawn = where((Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] gt 3726.04 + 3726.04*dawn_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] lt 3726.04 + 3726.04*dawn_shift + width))
;      OII_3726_dusk = where((Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] gt 3726.04 + 3726.04*dusk_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] lt 3726.04 + 3726.04*dusk_shift + width)) 
;      OII_3729_dawn = where((Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] gt 3728.80 + 3728.80*dawn_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dawn_ribbon_index_b] lt 3728.80 + 3728.80*dawn_shift + width))
;      OII_3729_dusk = where((Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] gt 3728.80 + 3728.80*dusk_shift - width) and $
;                            (Dispersion_mask_b[indgen(2098),dusk_ribbon_index_b] lt 3728.80 + 3728.80*dusk_shift + width))                                  
;      ;width = 2.5                      
;      Ribbon_brightness_dawn_4069[i,0] = (total(image_b[SII_4069_dawn, round(dawn_ribbon_index_b[SII_4069_dawn[5]]) - 1.: round(dawn_ribbon_index_b[SII_4069_dawn[5]]) + 1.]) / 3.); - offset_dawn
;      Ribbon_brightness_dusk_4069[i,0] = (total(image_b[SII_4069_dusk, round(dusk_ribbon_index_b[SII_4069_dusk[5]]) - 1.: round(dusk_ribbon_index_b[SII_4069_dusk[5]]) + 1.]) / 3.); - offset_dusk
;      Ribbon_brightness_dawn_4069[i,1] = robust_sigma(image_b[SII_4069_dawn - 20., round(dawn_ribbon_index_b[SII_4069_dawn[5]]) - 1.: round(dawn_ribbon_index_b[SII_4069_dawn[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dusk_4069[i,1] = robust_sigma(image_b[SII_4069_dusk - 20., round(dusk_ribbon_index_b[SII_4069_dusk[5]]) - 1.: round(dusk_ribbon_index_b[SII_4069_dusk[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dawn_3726[i,0] = (total(image_b[OII_3726_dawn, round(dawn_ribbon_index_b[OII_3726_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3726_dawn[5]]) + 1.]) / 3.); - offset_dawn
;      Ribbon_brightness_dusk_3726[i,0] = (total(image_b[OII_3726_dusk, round(dusk_ribbon_index_b[OII_3726_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3726_dusk[5]]) + 1.]) / 3.); - offset_dusk
;      Ribbon_brightness_dawn_3726[i,1] = robust_sigma(image_b[OII_3726_dawn - 20., round(dawn_ribbon_index_b[OII_3726_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3726_dawn[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dusk_3726[i,1] = robust_sigma(image_b[OII_3726_dusk - 20., round(dusk_ribbon_index_b[OII_3726_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3726_dusk[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dawn_3729[i,0] = (total(image_b[OII_3729_dawn, round(dawn_ribbon_index_b[OII_3729_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3729_dawn[5]]) + 1.]) / 3.); - offset_dawn
;      Ribbon_brightness_dusk_3729[i,0] = (total(image_b[OII_3729_dusk, round(dusk_ribbon_index_b[OII_3729_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3729_dusk[5]]) + 1.]) / 3.); - offset_dusk
;      Ribbon_brightness_dawn_3729[i,1] = robust_sigma(image_b[OII_3729_dawn - 20., round(dawn_ribbon_index_b[OII_3729_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3729_dawn[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dusk_3729[i,1] = robust_sigma(image_b[OII_3729_dusk - 20., round(dusk_ribbon_index_b[OII_3729_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3729_dusk[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      
;      ;inspection: white is line sampling regions
;      dummy = image_b
;      dummy[SII_4069_dawn, round(dawn_ribbon_index_b[SII_4069_dawn[5]]) - 1.: round(dawn_ribbon_index_b[SII_4069_dawn[5]]) + 1.] = 1.e4
;      dummy[SII_4069_dusk, round(dusk_ribbon_index_b[SII_4069_dusk[5]]) - 1.: round(dusk_ribbon_index_b[SII_4069_dusk[5]]) + 1.] = 1.e4             
;      dummy[SII_4076_dawn, round(dawn_ribbon_index_b[SII_4076_dawn[5]]) - 1.: round(dawn_ribbon_index_b[SII_4076_dawn[5]]) + 1.] = 1.e4
;      dummy[SII_4076_dusk, round(dusk_ribbon_index_b[SII_4076_dusk[5]]) - 1.: round(dusk_ribbon_index_b[SII_4076_dusk[5]]) + 1.] = 1.e4
;      dummy[OII_3726_dawn, round(dawn_ribbon_index_b[OII_3726_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3726_dawn[5]]) + 1.] = 1.e4
;      dummy[OII_3726_dusk, round(dusk_ribbon_index_b[OII_3726_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3726_dusk[5]]) + 1.] = 1.e4             
;      dummy[OII_3729_dawn, round(dawn_ribbon_index_b[OII_3729_dawn[5]]) - 1.: round(dawn_ribbon_index_b[OII_3729_dawn[5]]) + 1.] = 1.e4
;      dummy[OII_3729_dusk, round(dusk_ribbon_index_b[OII_3729_dusk[5]]) - 1.: round(dusk_ribbon_index_b[OII_3729_dusk[5]]) + 1.] = 1.e4                 
;      window, 0, xs = 2098, ys = 1078
;      tv, bytscl(dummy, 0, 50)
;      stop
;   endfor
;
;   ;-----------------------------------Red Data Products--------------------------------------------------------------
;   for i = 0, n_elements(Images_r)-1 do begin
;      image = MRDFITS(Images_r[i], 0, torus_header, /Dscale, /silent)
;      metadata = MRDFITS(Images_r[i], 1)
;
;    ;retrieve the all line brightnesses in rayleighs
;      CML = float(sxpar(torus_header, 'CML_SYSI'))
;      RJ_PER_Pix = float(sxpar(torus_header, 'RJ_PER_P'))
;      dawn_ribbon_index = metadata.center + ((5.85 + 0.049*cos(((CML - 90.) - 167.)/!Radeg)) / RJ_PER_Pix)
;      dusk_ribbon_index = metadata.center - ((5.57 + 0.073*cos(((CML + 90.) - 130.)/!Radeg)) / RJ_PER_Pix)      
;      dawn_warmT_index = dawn_ribbon_index + (.5 / RJ_PER_Pix) 
;      dusk_warmT_index = dusk_ribbon_index + (.5 / RJ_PER_Pix) 
;      Dispersion_mask = metadata.dispersion_mask 
;      dawn_shift = -(2. * !pi * (5.85 + 0.049*cos(((CML - 90.) - 167.)/!Radeg)) * R_J / period_J) / c ; =  delta lambda / lambda at the ribbon position in Schneider & Trauger, 1995
;      dusk_shift = (2. * !pi * (5.57 + 0.073*cos(((CML + 90.) - 130.)/!Radeg)) * R_J / period_J) / c  ; =  delta lambda / lambda at the ribbon position in Schneider & Trauger, 1995
;  
;      SII_6731_dawn = where((Dispersion_mask[indgen(2098),dawn_ribbon_index] gt 6730.815 + 6730.815*dawn_shift - width) and $
;                            (Dispersion_mask[indgen(2098),dawn_ribbon_index] lt 6730.815 + 6730.815*dawn_shift + width))
;      SII_6731_dusk = where((Dispersion_mask[indgen(2098),dusk_ribbon_index] gt 6730.815 + 6730.815*dusk_shift - width) and  $
;                            (Dispersion_mask[indgen(2098),dusk_ribbon_index] lt 6730.815 + 6730.815*dusk_shift + width))
;      SII_6716_dusk = where((Dispersion_mask[indgen(2098),dusk_ribbon_index] gt 6716.440 + 6716.440*dusk_shift - width) and $
;                            (Dispersion_mask[indgen(2098),dusk_ribbon_index] lt 6716.440 + 6716.440*dusk_shift + width))
;      SII_6716_dawn = where((Dispersion_mask[indgen(2098),dawn_ribbon_index] gt 6716.440 + 6716.440*dawn_shift - width) and $
;                            (Dispersion_mask[indgen(2098),dawn_ribbon_index] lt 6716.440 + 6716.440*dawn_shift + width))                                   
;      between_dawn = where((Dispersion_mask[indgen(2098),dawn_ribbon_index] gt 6723.63 + 6723.63*dawn_shift - width) and $
;                           (Dispersion_mask[indgen(2098),dawn_ribbon_index] lt 6723.63 + 6723.63*dawn_shift + width))
;      between_dusk = where((Dispersion_mask[indgen(2098),dusk_ribbon_index] gt 6723.63 + 6723.63*dusk_shift - width) and $
;                           (Dispersion_mask[indgen(2098),dusk_ribbon_index] lt 6723.63 + 6723.63*dusk_shift + width))
;      SIII_6312_dusk = where((Dispersion_mask[indgen(2098),dusk_ribbon_index] gt 6312.060 + 6312.060*dusk_shift - width) and $
;                             (Dispersion_mask[indgen(2098),dusk_ribbon_index] lt 6312.060 + 6312.060*dusk_shift + width))
;      SIII_6312_dawn = where((Dispersion_mask[indgen(2098),dawn_ribbon_index] gt 6312.060 + 6312.060*dawn_shift - width) and $
;                             (Dispersion_mask[indgen(2098),dawn_ribbon_index] lt 6312.060 + 6312.060*dawn_shift + width))    
;                             
;      ;Error Estimates
;      offset_dawn = total(image[between_dawn, round(dawn_ribbon_index[between_dawn[5]]) - 1.: round(dawn_ribbon_index[between_dawn[5]]) + 1.]) / 3. ;AVERAGE BRIGHTNESS BETWEEN THE LINES, subtract this from the emission as a DC offset
;      offset_dusk = total(image[between_dusk, round(dusk_ribbon_index[between_dusk[5]]) - 1.: round(dusk_ribbon_index[between_dusk[5]]) + 1.]) / 3. ;AVERAGE BRIGHTNESS BETWEEN THE LINES, subtract this from the emission as a DC offset
;      err_dawn = sqrt(2)*robust_sigma(image[between_dawn, round(dawn_ribbon_index[between_dawn[5]]) - 1.: round(dawn_ribbon_index[between_dawn[5]]) + 1.])
;      err_dusk = sqrt(2)*robust_sigma(image[between_dusk, round(dusk_ribbon_index[between_dusk[5]]) - 1.: round(dusk_ribbon_index[between_dusk[5]]) + 1.])
;      ;rough error is the standard deviation of the pixels integrated added in quadrature with that of the offset
;
;      Ribbon_brightness_Dawn_6731[i,0] = (total(image[SII_6731_dawn, round(dawn_ribbon_index[SII_6731_dawn[5]]) - 1.: round(dawn_ribbon_index[SII_6731_dawn[5]]) + 1.]) / 3.) - offset_dawn
;      Ribbon_brightness_dusk_6731[i,0] = (total(image[SII_6731_dusk, round(dusk_ribbon_index[SII_6731_dusk[5]]) - 1.: round(dusk_ribbon_index[SII_6731_dusk[5]]) + 1.]) / 3.) - offset_dusk
;      Ribbon_brightness_dawn_6716[i,0] = (total(image[SII_6716_dawn, round(dawn_ribbon_index[SII_6716_dawn[5]]) - 1.: round(dawn_ribbon_index[SII_6716_dawn[5]]) + 1.]) / 3.) - offset_dawn
;      Ribbon_brightness_dusk_6716[i,0] = (total(image[SII_6716_dusk, round(dusk_ribbon_index[SII_6716_dusk[5]]) - 1.: round(dusk_ribbon_index[SII_6716_dusk[5]]) + 1.]) / 3.) - offset_dusk
;      Ribbon_brightness_dawn_6312[i,0] = (total(image[SIII_6312_dawn, round(dawn_ribbon_index[SIII_6312_dawn[5]]) - 1.: round(dawn_ribbon_index[SIII_6312_dawn[5]]) + 1.]) / 3.)
;      Ribbon_brightness_dusk_6312[i,0] = (total(image[SIII_6312_dusk, round(dusk_ribbon_index[SIII_6312_dusk[5]]) - 1.: round(dusk_ribbon_index[SIII_6312_dusk[5]]) + 1.]) / 3.)
;      
;      Ribbon_brightness_dawn_6731[i,1] = err_dawn
;      Ribbon_brightness_dusk_6731[i,1] = err_dusk
;      Ribbon_brightness_dawn_6716[i,1] = err_dawn
;      Ribbon_brightness_dusk_6716[i,1] = err_dusk
;      Ribbon_brightness_dawn_6312[i,1] = err_dawn
;      Ribbon_brightness_dusk_6312[i,1] = err_dusk
;      Ribbon_brightness_dawn_6312[i,1] = robust_sigma(image_b[SIII_6312_dawn - 20., round(dawn_ribbon_index_b[SIII_6312_dawn[5]]) - 1.: round(dawn_ribbon_index_b[SIII_6312_dawn[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      Ribbon_brightness_dusk_6312[i,1] = robust_sigma(image_b[SIII_6312_dusk - 20., round(dusk_ribbon_index_b[SIII_6312_dusk[5]]) - 1.: round(dusk_ribbon_index_b[SIII_6312_dusk[5]]) + 1.]) ;for now call the error the sigma of this box 20 pixel's blueward
;      
;      warmT_brightness_Dawn_6731[i,0] = (total(image[SII_6731_dawn, round(dawn_warmT_index[SII_6731_dawn[5]]) - 1.: round(dawn_warmT_index[SII_6731_dawn[5]]) + 1.]) / 3.) - offset_dawn
;      warmT_brightness_dusk_6731[i,0] = (total(image[SII_6731_dusk, round(dusk_warmT_index[SII_6731_dusk[5]]) - 1.: round(dusk_warmT_index[SII_6731_dusk[5]]) + 1.]) / 3.) - offset_dusk
;      warmT_brightness_dawn_6716[i,0] = (total(image[SII_6716_dawn, round(dawn_warmT_index[SII_6716_dawn[5]]) - 1.: round(dawn_warmT_index[SII_6716_dawn[5]]) + 1.]) / 3.) - offset_dawn
;      warmT_brightness_dusk_6716[i,0] = (total(image[SII_6716_dusk, round(dusk_warmT_index[SII_6716_dusk[5]]) - 1.: round(dusk_warmT_index[SII_6716_dusk[5]]) + 1.]) / 3.) - offset_dusk
;
;      dawn_SYSIII[i] = (360. + CML - 90.) mod 360.
;      dusk_SYSIII[i] = (CML + 90.) mod 360.
;      dusk_ratio[i,0] = Ribbon_brightness_dusk_6716[i,0] / Ribbon_brightness_dusk_6731[i,0]
;      dawn_ratio[i,0] = Ribbon_brightness_dawn_6716[i,0] / Ribbon_brightness_dawn_6731[i,0]
;      dusk_ratio[i,1] = dusk_ratio[i,0] * sqrt( (Ribbon_brightness_dusk_6716[i,1]/Ribbon_brightness_dusk_6716[i,0])^2. + (Ribbon_brightness_dusk_6731[i,1]/Ribbon_brightness_dusk_6731[i,0])^2.)
;      dawn_ratio[i,1] = dawn_ratio[i,0] * sqrt( (Ribbon_brightness_dawn_6716[i,1]/Ribbon_brightness_dawn_6716[i,0])^2. + (Ribbon_brightness_dawn_6731[i,1]/Ribbon_brightness_dawn_6731[i,0])^2.)
;      warmT_dusk_ratio[i,0] = warmT_brightness_dusk_6716[i,0] / warmT_brightness_dusk_6731[i,0]
;      warmT_dawn_ratio[i,0] = warmT_brightness_dawn_6716[i,0] / warmT_brightness_dawn_6731[i,0]
;      warmT_dusk_ratio[i,1] = dusk_ratio[i,0] * sqrt( (warmT_brightness_dusk_6716[i,1]/warmT_brightness_dusk_6716[i,0])^2. + (warmT_brightness_dusk_6731[i,1]/warmT_brightness_dusk_6731[i,0])^2.)
;      warmT_dawn_ratio[i,1] = dawn_ratio[i,0] * sqrt( (warmT_brightness_dawn_6716[i,1]/warmT_brightness_dawn_6716[i,0])^2. + (warmT_brightness_dawn_6731[i,1]/warmT_brightness_dawn_6731[i,0])^2.)
;   
;      ;inspection: white is line sampling region.
;      dummy = image
;      dummy[SII_6731_dawn, round(dawn_ribbon_index[SII_6731_dawn[5]]) - 1.: round(dawn_ribbon_index[SII_6731_dawn[5]]) + 1.] = 1.e4
;      dummy[SII_6731_dusk, round(dusk_ribbon_index[SII_6731_dusk[5]]) - 1.: round(dusk_ribbon_index[SII_6731_dusk[5]]) + 1.] = 1.e4
;      dummy[SII_6716_dawn, round(dawn_ribbon_index[SII_6716_dawn[5]]) - 1.: round(dawn_ribbon_index[SII_6716_dawn[5]]) + 1.] = 1.e4
;      dummy[SII_6716_dusk, round(dusk_ribbon_index[SII_6716_dusk[5]]) - 1.: round(dusk_ribbon_index[SII_6716_dusk[5]]) + 1.] = 1.e4 
;      dummy[SIII_6312_dawn, round(dawn_ribbon_index[SIII_6312_dawn[5]]) - 1.: round(dawn_ribbon_index[SIII_6312_dawn[5]]) + 1.] = 1.e4
;      dummy[SIII_6312_dusk, round(dusk_ribbon_index[SIII_6312_dusk[5]]) - 1.: round(dusk_ribbon_index[SIII_6312_dusk[5]]) + 1.] = 1.e4          
;      dummy[between_dawn, round(dawn_ribbon_index[between_dawn[5]]) - 1.: round(dawn_ribbon_index[between_dawn[5]]) + 1.] = 250.  
;      dummy[between_dusk, round(dusk_ribbon_index[between_dusk[5]]) - 1.: round(dusk_ribbon_index[between_dusk[5]]) + 1.] = 250.                 
;      window, 0, xs = 2098, ys = 1078
;      tv, bytscl(dummy, 0, 100)
;   endfor    

      ;test_densities = 10^(findgen(500.)/100.)
      ;line_ratio_10000K = fltarr(N_elements(test_densities))
      ;line_ratio_25000K = fltarr(N_elements(test_densities))
      ;line_ratio_50000K = fltarr(N_elements(test_densities))
      ;line_ratio_75000K = fltarr(N_elements(test_densities))
      ;line_ratio_100000K = fltarr(N_elements(test_densities))
      ;for i = 0, N_elements(test_densities)-1 do begin
      ;    result = line_ratio('S_II', dens = test_densities[i], temp = 1d4, matrix = 100.)
      ;    Line_ratio_10000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
      ;    result = line_ratio('S_II', dens = test_densities[i], temp = 2.5d4, matrix = 100.)
      ;    Line_ratio_25000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
      ;    result = line_ratio('S_II', dens = test_densities[i], temp = 5d4, matrix = 100.)
      ;    Line_ratio_50000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
      ;    result = line_ratio('S_II', dens = test_densities[i], temp = 7.5d4, matrix = 100.)
      ;    Line_ratio_75000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
      ;    result = line_ratio('S_II', dens = test_densities[i], temp = 1d5, matrix = 100.)
      ;    Line_ratio_100000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
      ;endfor  
      ;Dawn_densities = INTERPOL(test_densities, Line_ratio_50000K, dawn_ratio[*,0])
      ;Dusk_densities = INTERPOL(test_densities, Line_ratio_50000K, dusk_ratio[*,0])  
      ;warmT_Dawn_densities = INTERPOL(test_densities, Line_ratio_100000K, warmT_dawn_ratio[*,0])
      ;warmT_Dusk_densities = INTERPOL(test_densities, Line_ratio_100000K, warmT_dusk_ratio[*,0])   
      ;
;      ;==================================================Brightness=====================================================================
;          ;Compare to expected distance 
;
;      dawn_ribbon = 5.85 + 0.049*cos(Dawn_SYSIII - 167./!Radeg) ;dawn ribbon location in R_J from Schneider & Trauger's (1995) fit
;      dusk_ribbon = 5.57 + 0.073*cos(Dusk_SYSIII - 130./!Radeg) ;dusk ribbon location in R_J from Schneider & Trauger's (1995) fit
;
;      Ribbon_brightness_Dawn_r = fltarr(n_elements(images_r), n_elements(lines_r))
;      err_Ribbon_brightness_Dawn_r = fltarr(n_elements(images_r), n_elements(lines_r))
;      for i = 0, n_elements(images_r)-1 do begin
;          near = Min(Abs((R_J_per_pixel * dawn_distance_array_r) - dawn_ribbon[i]), index)
;          Ribbon_brightness_Dawn_r[i,*] = mean(dawn_r[ i, index + [-2,-1,0,1,2], * ], dimension = 2) ;get the average emission over a 5 pixel (.1 RJ range)
;          err_Ribbon_brightness_Dawn_r[i,*] = stddev(dawn_r[ i, index + [-2,-1,0,1,2], * ], dimension = 2)
;      endfor
;      Ribbon_brightness_dusk_r = fltarr(n_elements(images_r), n_elements(lines_r))
;      err_Ribbon_brightness_dusk_r = fltarr(n_elements(images_r), n_elements(lines_r))
;      for i = 0, n_elements(images_r)-1 do begin
;          near = Min(Abs((R_J_per_pixel * dusk_distance_array_r) - dusk_ribbon[i]), index)
;          Ribbon_brightness_dusk_r[i,*] = median(dusk_r[ i, index + [-2,-1,0,1,2], * ], dimension = 2) ;get the average emission over a 5 pixel (.1 RJ range)
;          err_Ribbon_brightness_dusk_r[i,*] = stddev(dusk_r[ i, index + [-2,-1,0,1,2], * ], dimension = 2)
;      endfor
;      Ribbon_brightness_Dawn_b = fltarr(n_elements(images_b), n_elements(lines_b))
;      err_ribbon_brightness_Dawn_b = fltarr(n_elements(images_b), n_elements(lines_b))
;      for i = 0, n_elements(images_b)-1 do begin
;          near = Min(Abs((R_J_per_pixel * dawn_distance_array_b) - dawn_ribbon[i]), index)
;          Ribbon_brightness_Dawn_b[i,*] = mean(dawn_b[ i, index + [-2,-1,0,1,2], * ], dimension = 2) ;get the average emission over a 5 pixel (.1 RJ range)
;          err_ribbon_brightness_Dawn_b[i,*] = stddev(dawn_b[ i, index + [-2,-1,0,1,2], * ], dimension = 2)
;      endfor
;      Ribbon_brightness_dusk_b = fltarr(n_elements(images_b), n_elements(lines_b))
;      err_ribbon_brightness_dusk_b = fltarr(n_elements(images_b), n_elements(lines_b))
;      for i = 0, n_elements(images_b)-1 do begin
;          near = Min(Abs((R_J_per_pixel * dusk_distance_array_b) - dusk_ribbon[i]), index)
;          Ribbon_brightness_dusk_b[i,*] = median(dusk_b[ i, index + [-2,-1,0,1,2], * ], dimension = 2) ;get the average emission over a 5 pixel (.1 RJ range)
;          err_ribbon_brightness_dusk_b[i,*] = stddev(dusk_b[ i, index + [-2,-1,0,1,2], * ], dimension = 2)
;      endfor
;      
;      Print, 'the Errors here are not propogated, but estimated from the data, bit of a hack!'
;      cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Brightness'+strmid(sxpar(torus_header, 'DATE-OBS'), 0,10) +'.eps'), /ENCAPSULATED, xsize = 8.5, ysize = 6.9
;        !P.font=1
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;        !p.charsize = 2.
;;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,0], xtitle = cgSymbol('lambda') + '!DIII!N Longitude', ytitle = 'Rayleighs', xrange = [0, 360], yrange = [000, 220], SYMSIZE = 2, ERR_THICK = 5, $
;;          PSYM = CGSYMCAT(14), Title = 'Ribbon Surface Brightness - 7 Nov. 2013', ERR_YLOW = err_Ribbon_brightness_Dawn_r[*,0], ERR_YHIGH = err_Ribbon_brightness_Dawn_r[*,0], color = cgcolor('Orange Red')          
;;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,1], ERR_YLOW = err_Ribbon_brightness_Dawn_r[*,1],  ERR_YHIGH = err_Ribbon_brightness_Dawn_r[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Red')  
;;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,2], ERR_YLOW = err_Ribbon_brightness_Dawn_r[*,2],  ERR_YHIGH = err_Ribbon_brightness_Dawn_r[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Dark Red')  
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,0], ERR_YLOW = err_Ribbon_brightness_dusk_r[*,0] + [0,0,0,-26],  ERR_YHIGH = err_Ribbon_brightness_dusk_r[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Orange Red')  
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,1], ERR_YLOW = err_Ribbon_brightness_dusk_r[*,1],  ERR_YHIGH = err_Ribbon_brightness_dusk_r[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Red')  
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,2], ERR_YLOW = err_Ribbon_brightness_dusk_r[*,2],  ERR_YHIGH = err_Ribbon_brightness_dusk_r[*,2] + [-3,0,0,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Dark Red') 
;;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,0], ERR_YLOW = err_ribbon_brightness_dawn_b[*,0],  ERR_YHIGH = err_ribbon_brightness_dawn_b[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Violet')  
;;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,1], ERR_YLOW = err_ribbon_brightness_dawn_b[*,1],  ERR_YHIGH = err_ribbon_brightness_dawn_b[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Blue Violet')  
;;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,2], ERR_YLOW = err_ribbon_brightness_dawn_b[*,2],  ERR_YHIGH = err_ribbon_brightness_dawn_b[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Dodger blue') 
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,0], ERR_YLOW = err_ribbon_brightness_dusk_b[*,0],  ERR_YHIGH = err_ribbon_brightness_dusk_b[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Violet')  
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,1], ERR_YLOW = err_ribbon_brightness_dusk_b[*,1],  ERR_YHIGH = err_ribbon_brightness_dusk_b[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Blue Violet')   
;;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,2], ERR_YLOW = err_ribbon_brightness_dusk_b[*,2],  ERR_YHIGH = err_ribbon_brightness_dusk_b[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Dodger blue') 
;;        
;;        
;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,0], xtitle = cgSymbol('lambda') + '!DIII!N Longitude', ytitle = 'Rayleighs', xrange = [0, 360], yrange = [000, 220], SYMSIZE = 2, ERR_THICK = 5, $
;          PSYM = CGSYMCAT(14), Title = 'Ribbon Surface Brightness - 7 Nov. 2013', color = cgcolor('Orange Red')          
;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Red')  
;        cgPlot, Dawn_SYSIII, Ribbon_brightness_Dawn_r[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Dark Red')  
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Orange Red')  
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Red')  
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_r[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Dark Red') 
;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Violet')  
;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Blue Violet')  
;        cgPlot, dawn_SYSIII, Ribbon_brightness_dawn_b[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(14), color = cgcolor('Dodger blue') 
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,0], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Violet')  
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,1], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Blue Violet')   
;        cgPlot, dusk_SYSIII, Ribbon_brightness_dusk_b[*,2], SYMSIZE = 2, ERR_THICK = 5, /OVERPLOT, PSYM = CGSYMCAT(16), color = cgcolor('Dodger blue') 
;        
;        ;legend
;        cgText, 280, 180, 'SII 6731'+cgSymbol('Angstrom'), color = cgcolor('Dark Red')  
;        cgText, 280, 165, 'SII 6716'+cgSymbol('Angstrom'), color = cgcolor('Red') 
;        cgText, 280, 150, 'SIII 6312'+cgSymbol('Angstrom'), color = cgcolor('Orange Red') 
;        cgText, 280, 135, 'SII 4069'+cgSymbol('Angstrom'), color = cgcolor('Dodger blue')
;        cgText, 280, 120, 'OII 3729'+cgSymbol('Angstrom'), color = cgcolor('Blue Violet') 
;        cgText, 280, 105, 'OII 3726'+cgSymbol('Angstrom'), color = cgcolor('Violet')
;      
;        !p.charsize =1.  
;      cgPS_Close 
;      set_plot,'WIN'
;
;      ;==================================================Electron Temperture & Density =====================================================================
;      
;      ;search density and temperature
;      temp_ratios_dawn = Ribbon_brightness_dawn_r[*,2] / Ribbon_brightness_dawn_b[*,2] ;SII 6731 / SII 4069
;      dens_ratios_dawn = Ribbon_brightness_dawn_r[*,1] / Ribbon_brightness_dawn_r[*,2] ;SII 6716 / SII 6731
;      temp_ratios_dusk = Ribbon_brightness_dusk_r[*,2] / Ribbon_brightness_dusk_b[*,2] ;SII 6731 / SII 4069
;      dens_ratios_dusk = Ribbon_brightness_dusk_r[*,1] / Ribbon_brightness_dusk_r[*,2] ;SII 6716 / SII 6731
;      
;      Print, 'This step is manual convergence until I get amoeba up and running'
;      ;initial guesses: density first:
;      dens_line_ratio_50000K = fltarr(n_elements(test_densities))
;      temp = 1.82d4
;      for i = 0, N_elements(test_densities) - 1 do begin
;        result = line_ratio('S_II', dens = test_densities[i], temp = temp, matrix = 100.)
;        dens_line_ratio_50000K[i] = result.emissivity[2,0] / result.emissivity[1,0] 
;      endfor 
;      print, result.linewave[2,0], '   /', result.linewave[1,0] 
;      for j = 0, n_elements(dens_ratios_dusk) -1 do begin
;         near = Min(Abs(dens_line_ratio_50000K - dens_ratios_dusk[j]), dusk_index)
;         print, 'density at', temp, 'K =', test_densities[dusk_index] 
;      endfor
;      
;      test_temps = Findgen(500)*200. + 200.
;      temp_line_ratio_3000cc = fltarr(n_elements(test_temps)) 
;      dens = 870.     
;      for i = 0, N_elements(test_temps) - 1 do begin
;        result = line_ratio('S_II', dens = dens, temp = test_temps[i])
;        emissivity = result.emissivity
;        linewave = result.linewave
;        temp_line_ratio_3000cc[i] = emissivity[1,0] / emissivity[4,0] 
;      endfor 
;      print, linewave[1,0], '   /', linewave[4,0] 
;      for j = 0, n_elements(temp_ratios_dusk) -1 do begin
;         near = Min(Abs(temp_line_ratio_3000cc - temp_ratios_dusk[j]), dusk_index)
;         print, 'temp at ',dens, 'cc =', test_temps[dusk_index] 
;      endfor
;
;      ;results of manual fits
;      dawn_dens = [1380., 1621., 1000., 1348]
;      dawn_temp = [14600., 14600., 13200., 10400.]
;      
;      dusk_dens = [1548., 1948., 2041., 871.]
;      dusk_temp = [9400., 10800., 11200., 18200.]
;
;      err_dens_ratios_dusk = dens_ratios_dusk * sqrt((err_ribbon_brightness_dusk_r[*,1] / Ribbon_brightness_dusk_r[*,1])^2. + (err_ribbon_brightness_dusk_r[*,2] / Ribbon_brightness_dusk_r[*,2])^2.) 
;      err_dens_ratios_dawn = dens_ratios_dawn * sqrt((err_ribbon_brightness_dawn_r[*,1] / Ribbon_brightness_dawn_r[*,1])^2. + (err_ribbon_brightness_dawn_r[*,2] / Ribbon_brightness_dawn_r[*,2])^2.) 
;      ;err_dens_ratios_dusk = dens_ratios_dusk * sqrt((err_ribbon_brightness_dusk_r[i,1] / Ribbon_brightness_dusk_r[i,1])^2. + (err_ribbon_brightness_dusk_r[i,2] / Ribbon_brightness_dusk_r[i,2])^2.) 
;
;      ;==================================================Electron Density=====================================================================
;      cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Electron_density'+strmid(sxpar(torus_header, 'DATE-OBS'), 0,10)+'.eps'), /ENCAPSULATED, xsize = 8.5, ysize = 7.5
;        !P.font=1
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;        !p.charsize = 2.
;        cgPlot, test_densities, Line_ratio_10000K, /xlog, xtitle = 'n!De-!N (cm!U-3!N)', ytitle = '6716' + cgSymbol('Angstrom') +' / 6731'+ cgSymbol('Angstrom') + ' Line Ratio', xrange = [1.e1, 1.e5], yrange = [0.3, 1.5]
;        cgPlot, test_densities, Line_ratio_10000K, COLOR = cgcolor('Blue'), /overplot  
;        cgPlot, test_densities, Line_ratio_25000K, COLOR = cgcolor('Green'), /overplot
;        cgPlot, test_densities, Line_ratio_50000K, COLOR = cgcolor('Yellow'), /overplot 
;        cgPlot, test_densities, Line_ratio_75000K, COLOR = cgcolor('Orange'), /overplot   
;        cgPlot, test_densities, Line_ratio_100000K, COLOR = cgcolor('Red'), /OVERPLOT
;              
;;        cgPlot, dawn_dens, dens_ratios_dawn, ERR_YLOW = err_dens_ratios_dawn, ERR_Yhigh = err_dens_ratios_dawn, /OVERPLOT, SYMSIZE = 1.5, Psym = CGSYMCAT(14) 
;;        cgPlot, dusk_dens, dens_ratios_dusk, ERR_YLOW = err_dens_ratios_dusk, ERR_Yhigh = err_dens_ratios_dusk, /OVERPLOT, SYMSIZE = 1.5, Psym = CGSYMCAT(16)        
;
;        cgPlot, dawn_dens, dens_ratios_dawn, /OVERPLOT, SYMSIZE = 1.5, Psym = CGSYMCAT(14) 
;        cgPlot, dusk_dens, dens_ratios_dusk, /OVERPLOT, SYMSIZE = 1.5, Psym = CGSYMCAT(16)
;        ;cgPlot, Dusk_densities, dusk_ratio[*,0], ERR_YLOW = dusk_ratio[*,1], ERR_Yhigh = dusk_ratio[*,1], psym = 6, /OVERPLOT, SYMSIZE = 1.5
;;        cgPlot, warmT_Dawn_densities, warmT_dawn_ratio[*,0], ERR_YLOW = warmT_dawn_ratio[*,1], ERR_Yhigh = warmT_dawn_ratio[*,1], psym = 3, /OVERPLOT, SYMSIZE = 1.5 
;;        cgPlot, warmT_Dusk_densities, warmT_dusk_ratio[*,0], ERR_YLOW = warmT_dusk_ratio[*,1], ERR_Yhigh = warmT_dusk_ratio[*,1], psym = 7, /OVERPLOT, SYMSIZE = 1.5    
;        
;        ;legend
;        XYouts, 22, .5, 'Dawn Ribbon (5.8-5.9 R!DJ!N)'
;        XYouts, 20, .4, 'Dusk Ribbon (5.5-5.65 R!DJ!N)'
;        cgPlot, 1300, .515, psym = CGSYMCAT(14) , /OVERPLOT, SYMSIZE = 2. 
;        cgPlot, 1300, .415, psym = CGSYMCAT(16), /OVERPLOT, SYMSIZE = 2.        
;        XYouts, 1.e4, 1.2, '10,000K T!De-!N'
;        XYouts, 1.e4, 1.25, '25,000K T!De-!N'
;        XYouts, 1.e4, 1.3, '50,000K T!De-!N'
;        XYouts, 1.e4, 1.35, '75,000K T!De-!N'
;        XYouts, 1.e4, 1.4, '100,000K T!De-!N'
;        cgPlot, [7000,9000], [1.22,1.22], COLOR = cgcolor('Blue'), /OVERPLOT
;        cgPlot, [7000,9000], [1.27,1.27], COLOR = cgcolor('Green'), /OVERPLOT
;        cgPlot, [7000,9000], [1.32,1.32], COLOR = cgcolor('Yellow'), /OVERPLOT
;        cgPlot, [7000,9000], [1.37,1.37], COLOR = cgcolor('Orange'), /OVERPLOT
;        cgPlot, [7000,9000], [1.42,1.42], COLOR = cgcolor('Red'), /OVERPLOT         
;        !p.charsize = 1.  
;      cgPS_Close 
;      set_plot,'WIN'   

;      cgPS_Open, filename=strcompress('C:\IDL\Io\Apache_Point_Programs\Distance_Profile.eps'), /ENCAPSULATED, xsize = 8.5, ysize = 6.9
;        !P.font=1
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;        !p.charsize = 2.
;
;        !p.charsize = 1.  
;      cgPS_Close 
;      set_plot,'WIN'

            ;One day I will learn how to properly correct for the Earth's rotational velocity at a surface point, for now ignore it.
            ;Horizons says it's .02 km/s difference on 2013-11-07T11:32:40.412, likely this so small because this is when Jupiter is near zenith
                  ;;Now get the Apache Point location at these times. . .    
                  ;cspice_bodvrd, 'EARTH', 'RADII', 3, abc
                  ;equatr =  abc[0]
                  ;polar  =  abc[2]
                  ;f =  ( equatr - polar ) / equatr
                  ;cspice_georec, lon, lat, observatory_altitude, equatr, f, epos ;McDonald WRT Earth:
                  ;cspice_pxform, 'IAU_EARTH', 'J2000', torus_time, rotate
                  ;jpos = transpose(rotate) # epos ;APO position WRT Geocenter in J2000  
                  ;
                  ;;Get the velocity of the APO surface point WRT geocenter in the J2000 Frame.            
                  ;target = 'Jupiter'
                  ;obsctr = 'EARTH'
                  ;obsref = 'ITRF93'
                  ;obsepc    =  0.0D
                  ;obspos =  jpos
                  ;outref = 'DSS-14_TOPO'
                  ;abcorr = 'CN+S'     
                  ;refloc = 'OBSERVER'
                  ;cspice_spkcpo, target, et, torus_time, refloc, abcorr, obspos, obsctr, obsref, state0, lt0




If part eq 7 then begin

  ;Setup CHIANTI for use with the line ratios
  ;use_chianti, 'C:\ssw\packages\chianti', abund='C:\ssw\packages\chianti\abundance\sun_photospheric_2011_caffau.abund'
  E_field_shift = 0.25 ;offset between the dawn and dusk ansae in R_j.
  cooler_dusk = 0. ;1.1 ;The Yoshioka et al. 2014 data shows it a bit cooler on the dawn-side, call it 1.1 eV cooler, and convert to Kelvin
  smooth_by = .1 ;R_J
  
  k_B = 1.38064852e-23 ;m2 kg s-2 K-1
  J_per_eV = 1.60218e-19 ;joules per electronvolt
  
  restore, FILENAME = 'C:\IDL\Io\Apache_Point_Programs\CHIANTI_dens_table.sav'

  ;Interpolate everything to the same grid
  grid_step = .02 ;spacing in R_J
  x_grid = findgen((8.1-4.5)/grid_step)*grid_step + 4.5  ;an even grid from 4.5 to 8.1 R_J with .02 R_J spacing
  
  READCOL,'C:\IDL\Io\Voyager_1_data.txt', F='X,X,X,X,A,X,X,A,X,X,X,X,X,X,X,X,X,X,X,X,X,X,A', V1_radius, V1_N_e, V1_T_e, STRINGSKIP = '#', /Silent 
  k_B = 1.38064852e-23 ;m2 kg s-2 K-1
  J_per_eV = 1.60218e-19 ;joules per electronvolt
  V1_radius = float(V1_radius[2:*])
  V1_T_e = float(V1_T_e[2:*]) / (K_b / J_per_Ev) 
  V1_N_e = float(V1_N_e[2:*])

  Images_r = FILE_SEARCH(strcompress(directory + 'Processed\' + '*r._Cal.fits'), count = n_images) 

    for i = 0, N_elements(Images_r) - 1 do begin
      img = MRDFITS(Images_r[i], 0, torus_header, /silent)
      err_img = MRDFITS(Images_r[i], 1, /silent)
      metadata_r = MRDFITS(Images_r[i], 2, /silent)
      profiles = MRDFITS(Images_r[i], 3, /silent)  
  
      R_J_per_pixel = float(sxpar(torus_header, 'RJ_PER_P'))
      CML = float(sxpar(torus_header, 'CML_SYSI'))
      IO_PA = float(sxpar(torus_header, 'IO_PA'))
      IO_SYSIII = float(sxpar(torus_header, 'IO_SYSII'))
      SSP_SYSIII = float(sxpar(torus_header, 'SSP_SYSI '))
      ET = double(sxpar(torus_header, 'EPHEM_TI '))

      ;------------------Dusk N_e Fits----------------------------------------

      x = profiles.Dusk_DISTANCE_ARRAY_R
      y = profiles.FIT_Dusk_LINE_PROFILES
      ERR_y = profiles.err_FIT_Dusk_LINE_PROFILES
      Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
      err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
      
      ;Interpolate everything to the same grid
      I_Y = fltarr(N_elements(x_grid),N_elements(lines_r))
      err_I_Y = fltarr(N_elements(x_grid),N_elements(lines_r))
      For j = 0, N_elements(lines_r)-1 do begin
        I_Y[*,j] = interpol(Y[*,j], x[*,j], x_grid, /NaN, /spline)
        err_I_Y[*,j] = interpol(err_Y[*,j], x[*,j], x_grid, /NaN, /spline)
      endfor
         
      SII_ratio = I_Y[*,1] / I_Y[*,2]
      ;don't forget to include covariance term in the error propagation
      err_SII_ratio = abs(SII_ratio) * sqrt( (err_I_y[*,1]/I_y[*,1])^2. + (err_I_y[*,2]/I_y[*,2])^2. - 2. * CORRELATE(I_Y[*,1], I_Y[*,2]) * err_I_y[*,1] * err_I_y[*,2] / (I_y[*,1]*I_y[*,2]) )

      ;Get T_e from Bagenal (1994)
      T_e = interpol(V1_T_e, V1_radius, x_grid, /NAN)
      junk = where(finite(T_e), /Null, complement = extrapolate) ;indices to extrapolate in where there's not Bagenal measurement
      T_e[extrapolate] = min(T_e, /NAN) ;set them to the min T_e here
      N_e_dusk = fltarr(N_elements(T_e))
      N_e_dusk_hi = fltarr(N_elements(T_e))
      N_e_dusk_lo = fltarr(N_elements(T_e))
      
      ;Get N_e from the chianti lookup table
      for j = 0, n_elements(T_e)-1 do begin
        junk = min(abs( test_temps - T_e[j] ), T_e_index ) 
        if ( (SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                   N_e_dusk[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
        endif else N_e_dusk[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j], /spline) 
      endfor
        ;Find the error
        for j = 0, n_elements(T_e)-1 do begin
          junk = min(abs( test_temps - T_e[j] ), T_e_index )    
          if ( (SII_ratio[j]+err_SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j]+err_SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                     N_e_dusk_lo[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
          endif else N_e_dusk_lo[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j]+err_SII_ratio[j], /spline) 
        endfor
        for j = 0, n_elements(T_e)-1 do begin
          junk = min(abs( test_temps - T_e[j] ), T_e_index ) 
          if ( (SII_ratio[j]-err_SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j]-err_SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                     N_e_dusk_hi[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
          endif else N_e_dusk_hi[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j]-err_SII_ratio[j], /spline) 
        endfor
        
;        ;if N_e is out of bounds in the parameter space then reject it. If an error bar is out of bounds, then set to the bound 
;        N_e_dusk [ where( (N_e_dusk ge max(test_dens)) or (N_e_dusk le min(test_dens)), /NULL ) ] = !values.F_Nan
;        N_e_dusk_lo [ where( N_e_dusk_lo le min(test_dens), /NULL ) ] = min(test_dens)
;        N_e_dusk_hi [ where( N_e_dusk_hi ge max(test_dens), /NULL ) ] = max(test_dens)
        
        ;Now get the error from the actual density value
        N_e_dusk_lo = N_e_dusk - N_e_dusk_lo
        N_e_dusk_hi = N_e_dusk_hi - N_e_dusk 

      if keyword_set(debug) then begin
        window, 0
        cgplot, x_grid, N_e_dusk, xrange = [4.5, 8.1], xstyle = 1, yrange = [0., 8.e3], thick = 2, $
         ytitle = 'Electron Density (cm!U-3!N) (Line of sight average)', ERR_YLOW = N_e_dusk_lo, ERR_YHIGH = N_e_dusk_hi, /ERR_CLIP, $
         err_thick =.1,  ERR_WIDTH = 0.
        window, 1
        cgplot, x_grid, T_e, xrange = [4.5, 8.1], xstyle = 1, yrange = [min(test_temps), max(test_temps)], $
         ytitle = 'Electron Temp (K) (Line of sight average)'   
      endif  

      ;------------------Dawn N_e Fits----------------------------------------
      x = profiles.DAWN_DISTANCE_ARRAY_R
      y = profiles.FIT_DAWN_LINE_PROFILES
      ERR_y = profiles.err_FIT_DAWN_LINE_PROFILES
      Y[*,1:2] = Smooth(Y[*,1:2], [round( smooth_by / R_J_per_pixel), 0], /edge_mirror) 
      err_Y[*,1:2] = err_Y[*,1:2] / sqrt(round(smooth_by / R_J_per_pixel))
      
      ;Interpolate everything to the same grid
      I_Y = fltarr(N_elements(x_grid),N_elements(lines_r))
      err_I_Y = fltarr(N_elements(x_grid),N_elements(lines_r))
      For j = 0, N_elements(lines_r)-1 do begin
        I_Y[*,j] = interpol(Y[*,j], x[*,j], x_grid, /NaN, /spline)
        err_I_Y[*,j] = interpol(err_Y[*,j], x[*,j], x_grid, /NaN, /spline)
      endfor
      
      SII_ratio = I_Y[*,1] / I_Y[*,2]
      ;don't forget to include covariance term in the error propagation
      err_SII_ratio = abs(SII_ratio) * sqrt( (err_I_y[*,1]/I_y[*,1])^2. + (err_I_y[*,2]/I_y[*,2])^2. - 2. * CORRELATE(I_Y[*,1], I_Y[*,2]) * err_I_y[*,1] * err_I_y[*,2] / (I_y[*,1]*I_y[*,2]) )

      ;Get T_e from Bagenal (1994)
      T_e = interpol(V1_T_e, V1_radius + E_field_shift, x_grid, /NAN)
      T_e = T_e - cooler_dusk*1.16045221e4 
      junk = where(finite(T_e), /Null, complement = extrapolate) ;indices to extrapolate in where there's not Bagenal measurement
      T_e[extrapolate] = min(T_e, /NAN) ;set them to the min T_e here
      N_e_dawn = fltarr(N_elements(T_e))
      N_e_dawn_hi = fltarr(N_elements(T_e))
      N_e_dawn_lo = fltarr(N_elements(T_e))
      
      ;Get N_e from the chianti lookup table
      for j = 0, n_elements(T_e)-1 do begin
        junk = min(abs( test_temps - T_e[j] ), T_e_index ) 
        if ( (SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                   N_e_dawn[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
        endif else N_e_dawn[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j], /spline) 
      endfor
        ;Find the error
        for j = 0, n_elements(T_e)-1 do begin
          junk = min(abs( test_temps - T_e[j] ), T_e_index )    
          if ( (SII_ratio[j]+err_SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j]+err_SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                     N_e_dawn_lo[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
          endif else N_e_dawn_lo[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j]+err_SII_ratio[j], /spline) 
        endfor
        for j = 0, n_elements(T_e)-1 do begin
          junk = min(abs( test_temps - T_e[j] ), T_e_index ) 
          if ( (SII_ratio[j]-err_SII_ratio[j] gt max(dens_table[*, T_e_index])) or (SII_ratio[j]-err_SII_ratio[j] lt min(dens_table[*, T_e_index])) ) then begin
                     N_e_dawn_hi[j] = !values.F_NaN ;bound it to the parameter space, set to NaN when out of bounds  
          endif else N_e_dawn_hi[j] = interpol(test_dens, dens_table[*, T_e_index], SII_ratio[j]-err_SII_ratio[j], /spline) 
        endfor
        
;        ;If N_e is out of bounds in the parameter space then reject it. If an error bar is out of bounds, then set to the bound 
;        N_e_dawn [ where( (N_e_dawn ge max(test_dens)) or (N_e_dawn le min(test_dens)), /NULL ) ] = !values.F_Nan
;        N_e_dawn_lo [ where( N_e_dawn_lo le min(test_dens), /NULL ) ] = min(test_dens)
;        N_e_dawn_hi [ where( N_e_dawn_hi ge max(test_dens), /NULL ) ] = max(test_dens)
        
        ;Now get the error from the actual density value
        N_e_dawn_lo = N_e_dawn - N_e_dawn_lo
        N_e_dawn_hi = N_e_dawn_hi - N_e_dawn 

      if keyword_set(debug) then begin
        window, 0
        cgplot, x_grid, N_e_dawn, xrange = [8.1, 4.5], xstyle = 1, yrange = [0., 8.e3], $
         ytitle = 'Electron Density (cm!U-3!N) (Line of sight average)', ERR_YLOW = N_e_dawn_lo, ERR_YHIGH = N_e_dawn_hi, /ERR_CLIP, thick =4
        window, 1
        cgplot, x_grid, T_e, xrange = [8.1, 4.5], xstyle = 1, yrange = [min(test_temps), max(test_temps)], $
         ytitle = 'Electron Temp (K) (Line of sight average)'   
      endif

      ;Can only update these tags if they exist already, otherwise create them
      update = {N_e_dawn:N_e_dawn, N_e_dawn_lo:N_e_dawn_lo, N_e_dawn_hi:N_e_dawn_hi, N_e_dusk:N_e_dusk, N_e_dusk_lo:N_e_dusk_lo, N_e_dusk_hi:N_e_dusk_hi, X_grid:X_grid}
      if TAG_EXIST(profiles, 'N_e_dawn') then begin
        struct_assign, update, profiles, /verbose, /NOZERO
        updated_structure = profiles
      endif else begin
        Updated_Structure = CREATE_STRUCT(profiles, update)
      endelse
      ;Write it into the fits files
      MWRFITS, img, Images_r[i], torus_header, /CREATE, /silent ;/create overwrites
      MWRFITS, err_img, Images_r[i], /silent  ;Append the fits file -> extension 1  
      MWRFITS, metadata_r, Images_r[i], /silent ;Append the fits file -> extension 2  
      MWRFITS, Updated_Structure, Images_r[i], /silent ;Append the fits file -> extension 3
    endfor ;loop of images in this directory 
endif ;part eq 7     
  print, 'DONE with everything!!!!'
end
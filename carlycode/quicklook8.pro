
skipgetdata=0
if skipgetdata eq 0 then begin

   offline=0
   
;Red
snr=215.0
saturation=0.47
lambda=617.7863e-9 ;m

;Blue
snr=123.0
saturation=0.16
lambda=490.78823e-9 ;m

rd=31.0
gain=58.6
h=6.626e-22

c=299792458.0 ;m/s
r=2.37e-9
r=r/gain
r=r*h*c/lambda

a=(-1.0)
b=snr^2.0 ;SNR^2.0
c=(snr^2.0)*(rd^2.0)

val1=(-1.0)*b
val2=sqrt((b^2.0)-(4*a*c))
val3=(2*a)

ans=(val1-val2)/val3
print, ans/gain/saturation

etc_dn=ans/gain ;dn for I/F=1.0

SPICEdir='/Users/carlyhowett1/Desktop/isis3/data/science_spice/'
;spicesetupfile=SPICEdir+'15sci.tm'
spicesetupfile=SPICEdir+'15sci_rhr.tm'
cspice_furnsh,spicesetupfile
;cspice_furnsh,'/Users/carlyhowett1/Desktop/isis3/data/science_spice/kernels/pck/pck00010.tpc'

;INPUTS
;temp_filename='ralph_temps.txt'
;readcol, temp_filename, temp_met, cdhtempwarmrad, cdhtempcoldrad, cdhtempbulk, cdhtempsc, cdhtempelecbox, format='d, i, f, i, i,i', skipline=1
;temp_julian=met2julian(temp_met)

;date_label = LABEL_DATE(DATE_FORMAT = ['%D/%N/%Y'])
;setupplot, 'ralph_temps.ps'
;plot, temp_julian, cdhtempwarmrad, xtitle='Date', ytitle='CDH Warm',XTICKFORMAT = 'LABEL_DATE', psym=1
;plot, temp_julian, cdhtempcoldrad, xtitle='Date', ytitle='CDH Cold',XTICKFORMAT = 'LABEL_DATE'
;plot, temp_julian, cdhtempbulk, xtitle='Date', ytitle='CDH Bulk',XTICKFORMAT = 'LABEL_DATE'
;plot, temp_julian, cdhtempsc, xtitle='Date', ytitle='CDH Sc/ft',XTICKFORMAT = 'LABEL_DATE'
;plot, temp_julian, cdhtempelecbox, xtitle='Date', ytitle='CDH Electronic Box',XTICKFORMAT = 'LABEL_DATE'
;device,/close

SPICEdir='/Users/carlyhowett1/Desktop/isis3/data/science_spice/'
spicesetupfile=SPICEdir+'15sci_rhr.tm'
cspice_furnsh,spicesetupfile

runbasephote=1
first=1

bsfilename_pluto_red=strcompress('basphot_Pluto_red.txt', /remove_all)
bsfilename_pluto_blue=strcompress('basphot_Pluto_blue.txt', /remove_all)
bsfilename_pluto_nir=strcompress('basphot_Pluto_nir.txt', /remove_all)
bsfilename_pluto_ch4=strcompress('basphot_Pluto_ch4.txt', /remove_all)

bsfilename_charon_red=strcompress('basphot_Charon_red.txt', /remove_all)
bsfilename_charon_blue=strcompress('basphot_Charon_blue.txt', /remove_all)
bsfilename_charon_nir=strcompress('basphot_Charon_nir.txt', /remove_all)
bsfilename_charon_ch4=strcompress('basphot_Charon_ch4.txt', /remove_all)

if offline eq 0 then begin
   openw, 1, bsfilename_pluto_red, width=250
   openw, 2, bsfilename_pluto_blue, width=250
   openw, 3, bsfilename_pluto_nir, width=250
   openw, 4, bsfilename_pluto_ch4, width=250
   
   for n=1, 4 do printf, n,  'name', ' ',' reqid', '                  ', 'shfilename', '                       ', 'midtime', ' ', 'Sub-Obs Lon (deg RHR)', ' ', 'Sub-Obs Lat (deg)', ' ', 'flux (photons/s)',' ', 'flux_err (ph/s)', ' ', 'fwhm', '      ', 'mag', '         ', 'max', '           ', 'skymean', '          ', 'skyerr', '       ', 'skysig', '       ', 'xcen', '         ', 'ycen'

   openw, 11, bsfilename_charon_red, width=250
   openw, 12, bsfilename_charon_blue, width=250
   openw, 13, bsfilename_charon_nir, width=250
   openw, 14, bsfilename_charon_ch4, width=250
   for n=11, 14 do printf, n, 'name', ' ',' reqid', '                  ', 'shfilename', '                       ', 'midtime', ' ', 'Sub-Obs Lon (deg RHR)', ' ', 'Sub-Obs Lat (deg)', ' ', 'flux (photons/s)',' ', 'flux_err (ph/s)', ' ', 'fwhm', '      ', 'mag', '         ', 'max', '           ', 'skymean', '          ', 'skyerr', '       ', 'skysig', '       ', 'xcen', '         ', 'ycen'
endif



level1=0
dir='/Users/carlyhowett1/New_Horizons_Data/MVIC/Level2/pluto/'
if level1 eq 1 then dir='/Users/carlyhowett1/New_Horizons_Data/MVIC/Level1/pluto/'
hidir='/Users/carlyhowett1/New_Horizons_Data/MVIC/Level2/hicad/'

runfirst=0
maxobs=0
minobs=0

minobs=30
minobs=0
maxobs=56


setupplot, 'aperature_info.ps'
loadct, 13

if offline eq 0 then openw, 25, 'Charon_summary.txt', width=250

for obs=minobs, maxobs do begin
   if obs ne 51 then begin
   load_data_v2, obs, file0, file1, file2, file3, c_red, c_blue, c_nir, c_ch4, c_red_charon, c_blue_charon, c_nir_charon, c_ch4_charon, bad, onb
   

;#########


   if obs le 30 then begin
      pradius=4.0
      cradius=3.0
      
   endif
   
   if obs gt 30 then begin
      pradius=5.0
      cradius=4.0
   endif 
   
   if obs gt 50 then begin
      pradius=6.0
      cradius=5.0
   endif
   

;Basphote Inputs
   gain=58.6               ;e-/DN (Digital Number, or least significant bit ;nd
;if pradius eq -1 then pradius=5.0             
;if cradius eq -1 then cradius=3.0

   rdnoise=10.0
   
   pluto_charon_seperation = (((c_red(0)-c_red_charon(0))^2.0) + ((c_red(1)-c_red_charon(1))^2.0))^0.5
   print, pluto_charon_seperation
   
   wiggleroom=5.0
   pluto_charon_seperation=pluto_charon_seperation+wiggleroom
   
   
   maxn=3
   
   
;-------------------------------
;Main Loop through filters
;-------------------------------


for n=0, maxn do begin
   if n eq 0 then begin
      c_int=c_red
      c_charon=c_red_charon
      filename=file0
      pivot_lambda_um=617.79863e-3
      name='Red'
   endif
      
   if n eq 1 then begin
      c_int=c_blue
      c_charon=c_blue_charon
      filename=file1
      pivot_lambda_um=490.78823e-3
      name='Blue'
   endif
      
   if n eq 2 then begin
      c_int=c_nir
      c_charon=c_nir_charon
      filename=file2
      name='NIR'
   endif
   
   if n eq 3 then begin
      c_int=c_ch4
      c_charon=c_ch4_charon
      filename=file3
      name='CH4'
   endif

   ;stop
   shfilename=file_basename(filename, '.fit')
   img0=readfits(filename, header)
   get_subscft_latlon, header, plat, plon, clat, clon
   
   if plon lt 0 then plon=360+plon
   if clon lt 0 then clon=360+clon

   ;clon=clon-180 ;to make it from -180 to 180
   
   midtime=sxpar(header, 'SPCUTCAL')
   reqid=sxpar(header, 'REQID')
   exp_time=sxpar(header, 'EXPTIME')
   rdnoise=sxpar(header, 'READNOI')
   gain=sxpar(header, 'GAIN')
   range=sxpar(header, 'SPCTRANG')
   sidev=sxpar(header, 'SIDE')
   metv=sxpar(header, 'MET')
   metv_header=sxpar(header, 'MET_HK')

   spcft_sunv=sxpar(header, 'SPCSSCRN') ;To Sun Center
   spcft_earthv=sxpar(header, 'SPCESCRN') ;To Earth Center
   targ_sunv=sxpar(header, 'SPCTSORN') ; [km] Sun range to Target center     
   targ_earthv=sxpar(header, 'SPCTEORN') ;Earth observer to target center

   PSOLAR = sxpar(header, 'PSOLAR')
   PJUPITER = sxpar(header, 'PJUPITER')
   PPHOLUS = sxpar(header, 'PPHOLUS')
   PPLUTO = sxpar(header, 'PPLUTO')
   PCHARON = sxpar(header, 'PCHARON')
   
   RSOLAR = sxpar(header, 'RSOLAR')
   RJUPITER = sxpar(header, 'RJUPITER')
   RPHOLUS = sxpar(header, 'RPHOLUS')
   RPLUTO = sxpar(header, 'RPLUTO')
   RCHARON = sxpar(header, 'RCHARON')
   
  


   ;PCHARON_RED =  8.0606e+13
   ;PPHOLUS_RED =  8.3189e+13
   ;PJUPITER_RED = 8.5762e+13
   ;PSOLAR_RED =   8.0836e+13
   ;PPLUTO_RED =   8.0748e+13
   ;PCHARON_BLUE =   2.063e+13
   ;PPHOLUS_BLUE =   2.1424e+13
   ;PJUPITER_BLUE =  2.048e+13
   ;PSOLAR_BLUE =    2.0685e+13
   ;PPLUTO_BLUE =    2.0974e+13
   ;PCHARON_NIR =    1.0959e+14
   ;PPHOLUS_NIR =    1.0634e+14
   ;PJUPITER_NIR =   1.7801e+14
   ;PSOLAR_NIR =     1.096e+14
   ;PPLUTO_NIR =     1.1041e+14
   ;PCHARON_CH4 =   2.6702e+13
   ;PPHOLUS_CH4 =   2.6578e+13
   ;PJUPITER_CH4 =  6.3653e+13
   ;PSOLAR_CH4 =    2.6703e+13
   ;PPLUTO_CH4 =    2.6872e+13
   ;PCHARON_PAN1 =   2.5192e+14
   ;PPHOLUS_PAN1 =   2.5427e+14
   ;PJUPITER_PAN1 =  2.1761e+14
   ;PSOLAR_PAN1 =    2.5341e+14
   ;PPLUTO_PAN1 =    2.4377e+14
   ;PCHARON_PAN2 =   2.4398e+14
   ;PPHOLUS_PAN2 =   2.4626e+14
   ;PJUPITER_PAN2 =  2.1076e+14
   ;PSOLAR_PAN2 =    2.4543e+14
   ;PPLUTO_PAN2 =    2.3609e+14



   
   
   ;subsollon=sxpar(header, 'SPCTSOLO')
   ;subsclon=sxpar(header, 'SPCTSCLO')
   ;phasev=subsollon-subsclon


   abcorr = 'LT'                ;only correcting for light time
   frame='J2000'
   target='Pluto'
   
   obscode = '-98' ; New Horizons is the observer
   utctime=midtime
   CSPICE_STR2ET, utcTime, et
   CSPICE_ILLUM, target, et, abcorr, obscode, [0.,0.,0.], phasev, solar, emisn
   phasev=phasev*180./!pi

   target='Charon'
   CSPICE_ILLUM, target, et, abcorr, obscode, [0.,0.,0.], phasev_charon, solar, emisn
   phasev_charon=phasev_charon*180./!pi

   cspice_bodvrd,target, "RADII", 5, radii
   re = radii[0]
   rp = radii[2]
   f = ( re-rp)/re
   cspice_subslr, 'Near point: ellipsoid',target, $
                  et,'IAU_'+target, 'LT+S',obscode, $
                  ssolar, trgepc, srfvec
   cspice_recpgr,target, ssolar, re, f, sunloni, sunlati, sunalt
   sunloni=sunloni/!dtor
   sunloni=sunloni-180.0
   sublati=sunlati/!dtor
   
   if n eq 0 then printf, 25, midtime, metv_header, exp_time, range, sidev, phasev_charon, clon, clat, sunloni, sunlati
   print, sunloni, sunlati
   
   
   xloc=c_int(0)
   yloc=c_int(1)
   logfile='phot.txt'
   
   year=strmid(midtime, 0, 4)
   month=strmid(midtime, 5, 2)
   day=strmid(midtime, 8, 2)
   hour=strmid(midtime, 11,2)
   min=strmid(midtime, 14, 2)
   sec=strmid(midtime, 17, 2)
   
   jtime=julday(month, day, year, hour, min, sec)

   midtime_string=midtime
   midtime=jtime

   print, exp_time


   sky_values=get_sky(pluto_charon_seperation, metv_header, n)
   sky1=sky_values(0)
   sky2=sky_values(1)
   charon_sky1=sky_values(2)
   charon_sky2=sky_values(3)
  
   !P.Multi=[0,2,2]
   plot_sky_values, filename, img0, sky1, sky2, charon_sky1, charon_sky2, c_int, c_charon, pradius, cradius
   
   
   print, 'sky1:', sky1
   print, 'sky2:', sky2
   ;Destripe

   destripe=0
   if destripe eq 1 then begin
      nlines=n_elements(img0(0,*))
      newimg=img0
      newimg(*,*)=0.0
      medout=fltarr(nlines)
      for na=0, nlines-1 do begin
         vals=reform(img0(*,na))
         allvals=vals

         tmp=where(vals gt -23 and vals lt 2000)
         vals=vals(tmp)

         tmp=where(vals ne -1.0)
         if tmp(0) ne -1 then begin
            vals=vals(tmp)
            
            med=median(vals)
            
            h=where(abs(vals-med) lt 5)
            
            newmed=median(vals(h))
         endif

         if tmp(0) eq -1 then begin
            newmed=median(allvals)
         endif
         
         
         medout(na)=newmed
            
            
         newimg(*,na)=allvals-newmed
         
      endfor

      img0=newimg
   endif
   


   ;---------------------
   ;Running basphote
   ;---------------------

   basphote_v80, gain, img0, exp_time, xloc, yloc, pradius, sky1, sky2, logfile, $
                 rdnoise=rdnoise,$
                 err=err, flerr=flerr, flux=flux, $
                 fwhm=fwhm, onchip=onchip, outjd=outjd, mag=mag, max=max, $
                 skymean=skymean, skyerr=skyerr, skysig=skysig,$
                 xcen=xcen, ycen=ycen
   
   flux_counts=flux*exp_time/gain
   flux_counts_err=flerr*exp_time/gain
   print, flux_counts
   
   if n eq 0 then skylim=[-1.5, 2.0]
   ;if n eq 1 then skylim=[-1.5, 1.0]
   if n eq 1 then skylim=[-1.5, 1.5]
   if n eq 2 then skylim=[-6, 2]
   if n eq 3 then skylim=[-3, 3]
   

   if skymean lt skylim(0) then stop
   if skymean gt skylim(1) then stop
   print, flux, flerr           ;
   
   if flux lt 0 then begin
      print, 'Negative flux in filename:'
      print, filename
      stop
   endif
   
   
   
   apersky=sky1
   APER, img0, xloc, yloc, aper_flux_out, aper_flerr_out, apersky, skyerr, 1, pradius, [sky1,sky2],[-32767,8000],/exact,/flux
   print, aper_flux_out         ;counts

   
   if flux_counts lt 1e-13 then stop
   ;stop
   if n eq 0 then begin
      onbefore_sm=onb
      range_sm=range
      side_sm=sidev
      phase_pluto_sm=phasev
      phase_charon_sm=phasev_charon
      subsollon_sm=sunloni
      subsollat_sm=sunlati


      spcft_sunsm=spcft_sunv
      spcft_earthsm=spcft_earthv
      targ_earthsm=targ_earthv
      targ_sunsm=targ_sunv
      
      psolar_red=psolar
      pjupiter_red=pjupiter
      ppholus_red=ppholus
      ppluto_red=ppluto
      pcharon_red=pcharon
     
      rsolar_red=rsolar
      rjupiter_red=rjupiter
      rpholus_red=rpholus
      rpluto_red=rpluto
      rcharon_red=rcharon

      
      
      met_sm=metv
      met_header_sm=metv_header
      midjulian_sm=midtime
      reqid_sm=reqid
      pluto_radius_sm=pradius
      charon_radius_sm=cradius
      pluto_lon_sm=plon
      pluto_lat_sm=plat
      red_sky1_sm=sky1
      red_sky2_sm=sky2
      charon_lon_sm=clon
      charon_lat_sm=clat
      red_aper_flux_sm=aper_flux_out
      red_fwhm_sm=fwhm
      red_flux_sm=flux
      red_flux_err_sm=flerr
      red_counts_sm=flux_counts
      red_counts_err_sm=flux_counts_err
      red_err_sm=err
      red_mag_sm=mag
      red_max_sm=max
      red_skymean_sm=skymean
      red_skyerr_sm=skyerr
      red_skysig_sm=skysig
      red_xcen_sm=xcen
      red_ycen_sm=ycen
      red_exptime_sm=exp_time
      red_midtime_sm=midtime
   endif
         
   if n eq 1 then begin

      psolar_blue=psolar
      pjupiter_blue=pjupiter
      ppholus_blue=ppholus
      ppluto_blue=ppluto
      pcharon_blue=pcharon
     
      rsolar_blue=rsolar
      rjupiter_blue=rjupiter
      rpholus_blue=rpholus
      rpluto_blue=rpluto
      rcharon_blue=rcharon

      
      blue_fwhm_sm=fwhm
      blue_sky1_sm=sky1
      blue_sky2_sm=sky2
      blue_flux_sm=flux
      blue_flux_err_sm=flerr
      blue_aper_flux_sm=aper_flux_out
      blue_err_sm=err
      blue_counts_sm=flux_counts
      blue_counts_err_sm=flux_counts_err
      blue_mag_sm=mag
      blue_max_sm=max
      blue_skymean_sm=skymean
      blue_skyerr_sm=skyerr
      blue_skysig_sm=skysig
      blue_xcen_sm=xcen
      blue_ycen_sm=ycen
      blue_exptime_sm=exp_time
      blue_midtime_sm=midtime
   endif
   
   if n eq 2 then begin
      psolar_nir=psolar
      pjupiter_nir=pjupiter
      ppholus_nir=ppholus
      ppluto_nir=ppluto
      pcharon_nir=pcharon
     
      rsolar_nir=rsolar
      rjupiter_nir=rjupiter
      rpholus_nir=rpholus
      rpluto_nir=rpluto
      rcharon_nir=rcharon

      
      nir_fwhm_sm=fwhm
      nir_sky1_sm=sky1
      nir_sky2_sm=sky2
      nir_flux_sm=flux
      nir_flux_err_sm=flerr
      nir_counts_sm=flux_counts
      nir_counts_err_sm=flux_counts_err
      nir_aper_flux_sm=aper_flux_out
      nir_err_sm=err
      nir_mag_sm=mag
      nir_max_sm=max
      nir_skymean_sm=skymean
      nir_skyerr_sm=skyerr
      nir_skysig_sm=skysig
      nir_xcen_sm=xcen
      nir_ycen_sm=ycen
      nir_exptime_sm=exp_time
      nir_midtime_sm=midtime
   endif
   
   if n eq 3 then begin

      psolar_ch4=psolar
      pjupiter_ch4=pjupiter
      ppholus_ch4=ppholus
      ppluto_ch4=ppluto
      pcharon_ch4=pcharon
     
      rsolar_ch4=rsolar
      rjupiter_ch4=rjupiter
      rpholus_ch4=rpholus
      rpluto_ch4=rpluto
      rcharon_ch4=rcharon


      
      ch4_fwhm_sm=fwhm
      ch4_sky1_sm=sky1
      ch4_sky2_sm=sky2
      ch4_flux_sm=flux
      ch4_flux_err_sm=flerr
      ch4_counts_sm=flux_counts
      ch4_counts_err_sm=flux_counts_err
      ch4_aper_flux_sm=aper_flux_out
      ch4_err_sm=err
      ch4_mag_sm=mag
      ch4_max_sm=max
      ch4_skymean_sm=skymean
      ch4_skyerr_sm=skyerr
      ch4_skysig_sm=skysig
      ch4_xcen_sm=xcen
      ch4_ycen_sm=ycen
      ch4_exptime_sm=exp_time
      ch4_midtime_sm=midtime
   endif
   
   
   unitnum=1+n
   
   printf, unitnum, name, ' ', reqid, ' ', shfilename, ' ', midtime, ' ', plon, ' ', plat, ' ', flux, ' ', flerr, ' ', fwhm, ' ', mag, ' ', max, ' ', skymean(0), ' ', skyerr(0), ' ', skysig(0), ' ', xcen, ' ', ycen
   
   xloc=c_charon(0)
   yloc=c_charon(1)
   print, name, ' Charon'

   
   APER, img0, xloc, yloc, aper_flux_out, aper_flerr_out, apersky, skyerr, 1, cradius, [6,10],[-32767,8000],/exact,/flux
   print, aper_flux_out         ;counts

   
   
   basphote, gain, img0, exp_time, xloc, yloc, cradius, charon_sky1, charon_sky2, logfile, $
             rdnoise=rdnoise,$
             err=err, flerr=flerr, flux=flux, $
             fwhm=fwhm, onchip=onchip, outjd=outjd, mag=mag, max=max, $
             skymean=skymean, skyerr=skyerr, skysig=skysig,$
             xcen=xcen, ycen=ycen
   
   ;flux is in phot/sec
   flux_counts=flux*exp_time/gain
   flux_counts_err=flerr*exp_time/gain

   ;if obs eq 50 then stop

   
   cunitnum=11+n
   printf, cunitnum, ' ', reqid, ' ', shfilename, ' ', midtime, ' ', clon, ' ', clat,  ' ', flux, ' ', flerr, ' ', fwhm, ' ', mag, ' ', max, ' ', skymean, ' ', skyerr, ' ', skysig, ' ', xcen, ' ', ycen
   
   
   if n eq 0 then begin
      charon_red_fwhm_sm=fwhm
      charon_red_flux_sm=flux
      charon_red_flux_err_sm=flerr
      charon_red_counts_sm=flux_counts
      charon_red_counts_err_sm=flux_counts_err
      charon_red_err_sm=err
      charon_red_mag_sm=mag
      charon_red_max_sm=max
      charon_red_skymean_sm=skymean
      charon_red_skyerr_sm=skyerr
      charon_red_skysig_sm=skysig
      charon_red_xcen_sm=xcen
      charon_red_ycen_sm=ycen
   endif
   
   if n eq 1 then begin
      charon_blue_fwhm_sm=fwhm
      charon_blue_flux_sm=flux
      charon_blue_flux_err_sm=flerr
      charon_blue_counts_sm=flux_counts
      charon_blue_counts_err_sm=flux_counts_err
      charon_blue_err_sm=err
      charon_blue_mag_sm=mag
      charon_blue_max_sm=max
      charon_blue_skymean_sm=skymean
      charon_blue_skyerr_sm=skyerr
      charon_blue_skysig_sm=skysig
      charon_blue_xcen_sm=xcen
      charon_blue_ycen_sm=ycen
   endif
   
   if n eq 2 then begin
      charon_nir_fwhm_sm=fwhm
      charon_nir_flux_sm=flux
      charon_nir_flux_err_sm=flerr
      charon_nir_counts_sm=flux_counts
      charon_nir_counts_err_sm=flux_counts_err
      charon_nir_err_sm=err
      charon_nir_mag_sm=mag
      charon_nir_max_sm=max
      charon_nir_skymean_sm=skymean
      charon_nir_skyerr_sm=skyerr
      charon_nir_skysig_sm=skysig
      charon_nir_xcen_sm=xcen
      charon_nir_ycen_sm=ycen
   endif
   
   if n eq 3 then begin
      charon_ch4_fwhm_sm=fwhm
      charon_ch4_flux_sm=flux
      charon_ch4_flux_err_sm=flerr
      charon_ch4_counts_sm=flux_counts
      charon_ch4_counts_err_sm=flux_counts_err
      charon_ch4_err_sm=err
      charon_ch4_mag_sm=mag
      charon_ch4_max_sm=max
      charon_ch4_skymean_sm=skymean
      charon_ch4_skyerr_sm=skyerr
      charon_ch4_skysig_sm=skysig
      charon_ch4_xcen_sm=xcen
      charon_ch4_ycen_sm=ycen
   endif

endfor



if first eq 1 then begin
   onbefore=onbefore_sm
   sc_range=range_sm
   pluto_radius=pluto_radius_sm
   charon_radius=charon_radius_sm
   pluto_lat=pluto_lat_sm
   pluto_lon=pluto_lon_sm
   phase_pluto=phase_pluto_sm
   phase_charon=phase_charon_sm
   subsollon=subsollon_sm
   subsollat=subsollat_sm


   spcft_sun=spcft_sunsm
   spcft_earth=spcft_earthsm
   targ_earth=targ_earthsm
   targ_sun=targ_sunsm
   
   charon_lat=charon_lat_sm
   charon_lon=charon_lon_sm
   side=side_sm
   met=met_sm
   met_header=met_header_sm
   midjulian=midjulian_sm
   reqid_o=reqid_sm
   red_fwhm=red_fwhm_sm
   red_sky1=red_sky1_sm
   red_sky2=red_sky2_sm
   red_flux=red_flux_sm
   red_flux_err=red_flux_err_sm
   red_counts=red_counts_sm
   red_counts_err=red_counts_err_sm
   red_aper_flux=red_aper_flux_sm
   red_err=red_err_sm
   red_mag=red_mag_sm
   red_max=red_max_sm
   red_skymean=red_skymean_sm
   red_skyerr=red_skyerr_sm
   red_skysig=red_skysig_sm
   red_xcen=red_xcen_sm
   red_ycen=red_ycen_sm
   red_exptime=red_exptime_sm
   red_midtime=red_midtime_sm

   blue_fwhm=blue_fwhm_sm 
   blue_sky1=blue_sky1_sm
   blue_sky2=blue_sky2_sm
   blue_flux=blue_flux_sm
   blue_flux_err=blue_flux_err_sm
   blue_counts=blue_counts_sm
   blue_counts_err=blue_counts_err_sm
   blue_aper_flux=blue_aper_flux_sm
   blue_err=blue_err_sm
   blue_mag=blue_mag_sm
   blue_max=blue_max_sm
   blue_skymean=blue_skymean_sm
   blue_skyerr=blue_skyerr_sm
   blue_skysig=blue_skysig_sm
   blue_xcen=blue_xcen_sm   
   blue_ycen=blue_ycen_sm
   blue_exptime=blue_exptime_sm
   blue_midtime=blue_midtime_sm
   
   nir_fwhm=nir_fwhm_sm
   nir_sky1=nir_sky1_sm
   nir_sky2=nir_sky2_sm
   nir_flux=nir_flux_sm
   nir_flux_err=nir_flux_err_sm
   nir_counts=nir_counts_sm
   nir_counts_err=nir_counts_err_sm
   nir_aper_flux=nir_aper_flux_sm
   nir_err=nir_err_sm
   nir_mag= nir_mag_sm
   nir_max=nir_max_sm
   nir_skymean=nir_skymean_sm
   nir_skyerr=nir_skyerr_sm
   nir_skysig= nir_skysig_sm
   nir_xcen=nir_xcen_sm
   nir_ycen=nir_ycen_sm
   nir_exptime=nir_exptime_sm
   nir_midtime=nir_midtime_sm
   
   ch4_fwhm=ch4_fwhm_sm
   ch4_sky1=ch4_sky1_sm
   ch4_sky2=ch4_sky2_sm
   ch4_flux=ch4_flux_sm
   ch4_flux_err=ch4_flux_err_sm
   ch4_counts=ch4_counts_sm
   ch4_counts_err=ch4_counts_err_sm
   ch4_aper_flux=ch4_aper_flux_sm
   ch4_err=ch4_err_sm
   ch4_mag=ch4_mag_sm
   ch4_max=ch4_max_sm
   ch4_skymean=ch4_skymean_sm
   ch4_skyerr=ch4_skyerr_sm
   ch4_skysig=ch4_skysig_sm
   ch4_xcen=ch4_xcen_sm
   ch4_ycen=ch4_ycen_sm
   ch4_exptime=ch4_exptime_sm
   ch4_midtime=ch4_midtime_sm
      
   charon_red_fwhm=charon_red_fwhm_sm
   charon_red_flux=charon_red_flux_sm
   charon_red_flux_err=charon_red_flux_err_sm
   charon_red_counts=charon_red_counts_sm
   charon_red_counts_err=charon_red_counts_err_sm
   charon_red_err=charon_red_err_sm
   charon_red_mag=charon_red_mag_sm
   charon_red_max=charon_red_max_sm
   charon_red_skymean=charon_red_skymean_sm
   charon_red_skyerr=charon_red_skyerr_sm
   charon_red_skysig=charon_red_skysig_sm
   charon_red_xcen=charon_red_xcen_sm
   charon_red_ycen=charon_red_ycen_sm

   charon_blue_fwhm=charon_blue_fwhm_sm
   charon_blue_flux=charon_blue_flux_sm
   charon_blue_flux_err=charon_blue_flux_err_sm
   charon_blue_counts=charon_blue_counts_sm
   charon_blue_counts_err=charon_blue_counts_err_sm
   charon_blue_err=charon_blue_err_sm
   charon_blue_mag=charon_blue_mag_sm
   charon_blue_max=charon_blue_max_sm
   charon_blue_skymean=charon_blue_skymean_sm
   charon_blue_skyerr=charon_blue_skyerr_sm
   charon_blue_skysig=charon_blue_skysig_sm
   charon_blue_xcen=charon_blue_xcen_sm   
   charon_blue_ycen=charon_blue_ycen_sm
   
   charon_nir_fwhm=charon_nir_fwhm_sm
   charon_nir_flux=charon_nir_flux_sm
   charon_nir_flux_err=charon_nir_flux_err_sm
   charon_nir_counts=charon_nir_counts_sm
   charon_nir_counts_err=charon_nir_counts_err_sm
   charon_nir_err=charon_nir_err_sm
   charon_nir_mag= charon_nir_mag_sm
   charon_nir_max=charon_nir_max_sm
   charon_nir_skymean=charon_nir_skymean_sm
   charon_nir_skyerr=charon_nir_skyerr_sm
   charon_nir_skysig= charon_nir_skysig_sm
   charon_nir_xcen=charon_nir_xcen_sm
   charon_nir_ycen=charon_nir_ycen_sm
   
   charon_ch4_fwhm=charon_ch4_fwhm_sm
   charon_ch4_flux=charon_ch4_flux_sm
   charon_ch4_flux_err=charon_ch4_flux_err_sm
   charon_ch4_counts=charon_ch4_counts_sm
   charon_ch4_counts_err=charon_ch4_counts_err_sm
   charon_ch4_err=charon_ch4_err_sm
   charon_ch4_mag=charon_ch4_mag_sm
   charon_ch4_max=charon_ch4_max_sm
   charon_ch4_skymean=charon_ch4_skymean_sm
   charon_ch4_skyerr=charon_ch4_skyerr_sm
   charon_ch4_skysig=charon_ch4_skysig_sm
   charon_ch4_xcen=charon_ch4_xcen_sm
   charon_ch4_ycen=charon_ch4_ycen_sm
endif
 


if first eq 0 then begin
   onbefore=[onbefore, onbefore_sm]
   phase_pluto=[phase_pluto, phase_pluto_sm]
   phase_charon=[phase_charon, phase_charon_sm]
   subsollon=[subsollon,subsollon_sm]
   subsollat=[subsollat,subsollat_sm]
   sc_range=[sc_range, range_sm]
   
   spcft_sun=[spcft_sun, spcft_sunsm]
   spcft_earth=[spcft_earth, spcft_earthsm]
   targ_sun=[targ_sun, targ_sunsm]
   targ_earth=[targ_earth, targ_earthsm]

   
   pluto_radius=[pluto_radius, pluto_radius_sm]
   charon_radius=[charon_radius, charon_radius_sm]
   pluto_lat=[pluto_lat, pluto_lat_sm]
   pluto_lon=[pluto_lon, pluto_lon_sm]
   charon_lat=[charon_lat, charon_lat_sm]
   charon_lon=[charon_lon, charon_lon_sm]
   side=[side, side_sm]
   met=[met, met_sm]
   met_header=[met_header, met_header_sm]
   midjulian=[midjulian, midjulian_sm]
   reqid_o=[reqid_o, reqid_sm]
   red_sky1=[red_sky1, red_sky1_sm]
   red_sky2=[red_sky2, red_sky2_sm]
   red_fwhm=[red_fwhm,red_fwhm_sm]
  
   red_flux=[red_flux, red_flux_sm]
   red_flux_err=[red_flux_err, red_flux_err_sm]
   red_counts=[red_counts, red_counts_sm]
   red_counts_err=[red_counts_err,red_counts_err_sm]
   red_aper_flux=[red_aper_flux,red_aper_flux_sm]
   red_err=[red_err,red_err_sm]
   red_mag=[red_mag,red_mag_sm]
   red_max=[red_max,red_max_sm]
   red_skymean=[red_skymean,red_skymean_sm]
   red_skyerr=[red_skyerr,red_skyerr_sm]
   red_skysig=[red_skysig,red_skysig_sm]
   red_xcen=[red_xcen, red_xcen_sm]
   red_ycen=[red_ycen, red_ycen_sm]
   red_exptime=[red_exptime, red_exptime_sm]
   red_midtime=[red_midtime, red_midtime_sm]

   blue_sky1=[blue_sky1, blue_sky1_sm]
   blue_sky2=[blue_sky2, blue_sky2_sm]
   blue_fwhm=[blue_fwhm, blue_fwhm_sm]
   blue_flux=[blue_flux, blue_flux_sm]
   blue_flux_err=[blue_flux_err,blue_flux_err_sm]
   blue_counts=[blue_counts, blue_counts_sm]
   blue_counts_err=[blue_counts_err,blue_counts_err_sm]
   blue_aper_flux=[blue_aper_flux,blue_aper_flux_sm]
   blue_err=[blue_err,blue_err_sm]
   blue_mag=[blue_mag,blue_mag_sm]
   blue_max=[blue_max,blue_max_sm]
   blue_skymean=[blue_skymean,blue_skymean_sm]
   blue_skyerr=[blue_skyerr,blue_skyerr_sm]
   blue_skysig=[blue_skysig,blue_skysig_sm]
   blue_xcen=[ blue_xcen,blue_xcen_sm ]  
   blue_ycen=[blue_ycen,blue_ycen_sm]
   blue_exptime=[blue_exptime, blue_exptime_sm]
   blue_midtime=[blue_midtime, blue_midtime_sm]
   
   nir_sky1=[nir_sky1, nir_sky1_sm]
   nir_sky2=[nir_sky2, nir_sky2_sm]
   nir_fwhm=[nir_fwhm,nir_fwhm_sm]
   nir_flux=[nir_flux, nir_flux_sm]
   nir_flux_err=[nir_flux_err,nir_flux_err_sm]
   nir_counts=[nir_counts, nir_counts_sm]
   nir_counts_err=[nir_counts_err,nir_counts_err_sm]
   nir_aper_flux=[nir_aper_flux,nir_aper_flux_sm]
   nir_err=[nir_err,nir_err_sm]
   nir_mag=[nir_mag, nir_mag_sm]
   nir_max=[nir_max,nir_max_sm]
   nir_skymean=[nir_skymean,nir_skymean_sm]
   nir_skyerr=[nir_skyerr,nir_skyerr_sm]
   nir_skysig=[nir_skysig, nir_skysig_sm]
   nir_xcen=[nir_xcen,nir_xcen_sm]
   nir_ycen=[nir_ycen,nir_ycen_sm]
   nir_exptime=[nir_exptime, nir_exptime_sm]
   nir_midtime=[nir_midtime, nir_midtime_sm]
   
   ch4_sky1=[ch4_sky1, ch4_sky1_sm]
   ch4_sky2=[ch4_sky2, ch4_sky2_sm]
   ch4_fwhm=[ch4_fwhm,ch4_fwhm_sm]
   ch4_flux=[ch4_flux, ch4_flux_sm]
   ch4_flux_err=[ch4_flux_err,ch4_flux_err_sm]
   ch4_counts=[ch4_counts, ch4_counts_sm]
   ch4_counts_err=[ch4_counts_err,ch4_counts_err_sm]
   ch4_aper_flux=[ch4_aper_flux,ch4_aper_flux_sm]
   ch4_err=[ch4_err,ch4_err_sm]
   ch4_mag=[ch4_mag,ch4_mag_sm]
   ch4_max=[ch4_max,ch4_max_sm]
   ch4_skymean=[ch4_skymean,ch4_skymean_sm]
   ch4_skyerr=[ch4_skyerr,ch4_skyerr_sm]
   ch4_skysig=[ch4_skysig,ch4_skysig_sm]
   ch4_xcen=[ch4_xcen,ch4_xcen_sm]
   ch4_ycen=[ch4_ycen,ch4_ycen_sm]
   ch4_exptime=[ch4_exptime, ch4_exptime_sm]
   ch4_midtime=[ch4_midtime, ch4_midtime_sm]
   
   charon_red_fwhm=[charon_red_fwhm,charon_red_fwhm_sm] 
   charon_red_flux=[charon_red_flux, charon_red_flux_sm]
   charon_red_flux_err=[charon_red_flux_err,charon_red_flux_err_sm]
   charon_red_counts=[charon_red_counts, charon_red_counts_sm]
   charon_red_counts_err=[charon_red_counts_err,charon_red_counts_err_sm]
   charon_red_err=[charon_red_err,charon_red_err_sm]
   charon_red_mag=[charon_red_mag,charon_red_mag_sm]
   charon_red_max=[charon_red_max,charon_red_max_sm]
   charon_red_skymean=[charon_red_skymean,charon_red_skymean_sm]
   charon_red_skyerr=[charon_red_skyerr,charon_red_skyerr_sm]
   charon_red_skysig=[charon_red_skysig,charon_red_skysig_sm]
   charon_red_xcen=[charon_red_xcen,charon_red_xcen_sm]
   charon_red_ycen=[charon_red_ycen,charon_red_ycen_sm]

   charon_blue_fwhm=[charon_blue_fwhm,charon_blue_fwhm_sm] 
   charon_blue_flux=[charon_blue_flux, charon_blue_flux_sm]
   charon_blue_flux_err=[charon_blue_flux_err,charon_blue_flux_err_sm]
   charon_blue_counts=[charon_blue_counts, charon_blue_counts_sm]
   charon_blue_counts_err=[charon_blue_counts_err,charon_blue_counts_err_sm]
   charon_blue_err=[charon_blue_err,charon_blue_err_sm]
   charon_blue_mag=[charon_blue_mag,charon_blue_mag_sm]
   charon_blue_max=[charon_blue_max,charon_blue_max_sm]
   charon_blue_skymean=[charon_blue_skymean, charon_blue_skymean_sm]
   charon_blue_skyerr=[charon_blue_skyerr,charon_blue_skyerr_sm]
   charon_blue_skysig=[charon_blue_skysig,charon_blue_skysig_sm]
   charon_blue_xcen=[charon_blue_xcen,charon_blue_xcen_sm]   
   charon_blue_ycen=[charon_blue_ycen,charon_blue_ycen_sm]
   
   charon_nir_fwhm=[charon_nir_fwhm,charon_nir_fwhm_sm]
   charon_nir_flux=[charon_nir_flux, charon_nir_flux_sm]
   charon_nir_flux_err=[charon_nir_flux_err,charon_nir_flux_err_sm]
   charon_nir_counts=[charon_nir_counts, charon_nir_counts_sm]
   charon_nir_counts_err=[charon_nir_counts_err,charon_nir_counts_err_sm]
   charon_nir_err=[charon_nir_err,charon_nir_err_sm]
   charon_nir_mag= [charon_nir_mag, charon_nir_mag_sm]
   charon_nir_max=[charon_nir_max,charon_nir_max_sm]
   charon_nir_skymean=[charon_nir_skymean, charon_nir_skymean_sm]
   charon_nir_skyerr=[charon_nir_skyerr, charon_nir_skyerr_sm]
   charon_nir_skysig= [charon_nir_skysig,charon_nir_skysig_sm]
   charon_nir_xcen=[charon_nir_xcen,charon_nir_xcen_sm]
   charon_nir_ycen=[charon_nir_ycen,charon_nir_ycen_sm]
   
   charon_ch4_fwhm=[charon_ch4_fwhm, charon_ch4_fwhm_sm]
   charon_ch4_flux=[charon_ch4_flux, charon_ch4_flux_sm]
   charon_ch4_flux_err=[charon_ch4_flux_err,charon_ch4_flux_err_sm]
   charon_ch4_counts=[charon_ch4_counts, charon_ch4_counts_sm]
   charon_ch4_counts_err=[charon_ch4_counts_err,charon_ch4_counts_err_sm]
   charon_ch4_err=[charon_ch4_err,charon_ch4_err_sm]
   charon_ch4_mag=[charon_ch4_mag,charon_ch4_mag_sm]
   charon_ch4_max=[charon_ch4_max,charon_ch4_max_sm]
   charon_ch4_skymean=[charon_ch4_skymean,charon_ch4_skymean_sm]
   charon_ch4_skyerr=[charon_ch4_skyerr,charon_ch4_skyerr_sm]
   charon_ch4_skysig=[ charon_ch4_skysig,charon_ch4_skysig_sm]
   charon_ch4_xcen=[charon_ch4_xcen,charon_ch4_xcen_sm]
   charon_ch4_ycen=[charon_ch4_ycen,charon_ch4_ycen_sm]
endif


first=0
endif
   

   
endfor
   
close, /all
device,/close

if runbasephote eq 0 then savefilename='quicklook_results_aper.idlsave'
if runbasephote eq 1 then savefilename='quicklook_results_basphote.idlsave'

np=n_elements(side)

save, filename='quicklook8.idlsave', /all

endif

;##########################################################################



;From Alex's "Radiometric Calibration Keywords" document
;PCHARON_RED =  8.0606e+13
;PPHOLUS_RED =  8.3189e+13
;PJUPITER_RED = 8.5762e+13
;PSOLAR_RED =   8.0836e+13
;PPLUTO_RED =   8.0748e+13
;PCHARON_BLUE =   2.063e+13
;PPHOLUS_BLUE =   2.1424e+13
;PJUPITER_BLUE =  2.048e+13
;PSOLAR_BLUE =    2.0685e+13
;PPLUTO_BLUE =    2.0974e+13
;PCHARON_NIR =    1.0959e+14
;PPHOLUS_NIR =    1.0634e+14
;PJUPITER_NIR =   1.7801e+14
;PSOLAR_NIR =     1.096e+14
;PPLUTO_NIR =     1.1041e+14
;PCHARON_CH4 =   2.6702e+13
;PPHOLUS_CH4 =   2.6578e+13
;PJUPITER_CH4 =  6.3653e+13
;PSOLAR_CH4 =    2.6703e+13
;PPLUTO_CH4 =    2.6872e+13
;PCHARON_PAN1 =   2.5192e+14
;PPHOLUS_PAN1 =   2.5427e+14
;PJUPITER_PAN1 =  2.1761e+14
;PSOLAR_PAN1 =    2.5341e+14
;PPLUTO_PAN1 =    2.4377e+14
;PCHARON_PAN2 =   2.4398e+14
;PPHOLUS_PAN2 =   2.4626e+14
;PJUPITER_PAN2 =  2.1076e+14
;PSOLAR_PAN2 =    2.4543e+14
;PPLUTO_PAN2 =    2.3609e+14



calib_red_flux=red_counts/(red_exptime*PPLUTO_RED)
calib_blue_flux=blue_counts/(blue_exptime*PPLUTO_BLUE)
calib_nir_flux=nir_counts/(nir_exptime*PPLUTO_NIR)
calib_ch4_flux=ch4_counts/(ch4_exptime*PPLUTO_CH4)

calib_charon_red_flux=charon_red_counts/(red_exptime*PCHARON_RED)
calib_charon_blue_flux=charon_blue_counts/(blue_exptime*PCHARON_BLUE)
calib_charon_nir_flux=charon_nir_counts/(nir_exptime*PCHARON_NIR)
calib_charon_ch4_flux=charon_ch4_counts/(ch4_exptime*PCHARON_CH4)

calib_charon_red_flux_err=charon_red_counts_err/(red_exptime*PCHARON_RED)
calib_charon_blue_flux_err=charon_blue_counts_err/(blue_exptime*PCHARON_BLUE)
calib_charon_nir_flux_err=charon_nir_counts_err/(nir_exptime*PCHARON_NIR)
calib_charon_ch4_flux_err=charon_ch4_counts_err/(ch4_exptime*PCHARON_CH4)


save, filename='quicklook8_prenewscrange.idlsave', /all

openw, 1, 'quicklook8_prenewscrange.txt', width=1000
  printf, 1, 'Mid Time Julian, MET, Spacecraft-Sun (km), Spacecraft-Target (km), Side, Red Flux, Red Flux Error, Blue Flux, Blue Flux Error, NIR Flux, NIR Flux Error, CH4 Flux, CH4 Flux Error'
  
for n=0, np-1 do begin
   printf, 1, string(midjulian(n), format='(F20.6)'), ' ', met(n), ' ', $
           spcft_sun(n), ' ', sc_range(n), ' ', side(n), ' ', $
           red_flux(n), ' ', red_flux_err(n), ' ', $
           blue_flux(n), ' ', blue_flux_err(n), ' ', $
           nir_flux(n), ' ', nir_flux_err(n), ' ', $
           ch4_flux(n), ' ', ch4_flux_err(n)
endfor
close, 1


stop


;most distance spacecraft observation
sc_range0=max(sc_range)
new_sc_range=sc_range/sc_range0

corr_calib_red_flux=calib_red_flux*(new_sc_range^2.0)
corr_calib_blue_flux=calib_blue_flux*(new_sc_range^2.0)
corr_calib_nir_flux=calib_nir_flux*(new_sc_range^2.0)
corr_calib_ch4_flux=calib_ch4_flux*(new_sc_range^2.0)


corr_charon_calib_red_flux=calib_charon_red_flux*(new_sc_range^2.0)
corr_charon_calib_blue_flux=calib_charon_blue_flux*(new_sc_range^2.0)
corr_charon_calib_nir_flux=calib_charon_nir_flux*(new_sc_range^2.0)
corr_charon_calib_ch4_flux=calib_charon_ch4_flux*(new_sc_range^2.0)


corr_charon_calib_red_flux_err=calib_charon_red_flux_err*(new_sc_range^2.0)
corr_charon_calib_blue_flux_err=calib_charon_blue_flux_err*(new_sc_range^2.0)
corr_charon_calib_nir_flux_err=calib_charon_nir_flux_err*(new_sc_range^2.0)
corr_charon_calib_ch4_flux_err=calib_charon_ch4_flux_err*(new_sc_range^2.0)


set_plot, 'x'
!P.Multi=[0,1,1]
setdecomposedstate, 0
dummy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])

h=sort(midjulian)
bp=sort(pluto_lon)
bc=sort(charon_lon)


loadct, 0, ncolors=256, /silent
plot, midjulian(h), charon_nir_skymean, charsize=2, yrange=[-4, 4], $
      ytitle='Sky Mean', xtitle='Observation number', $
      XTICKFORMAT='LABEL_DATE'

for n=0, np-1 do if side(n) eq 1 then oplot, [midjulian(h(n)),midjulian(h(n))], !y.crange, linestyle=1
loadct, 13, ncolors=256, /silent
oplot, midjulian(h), charon_red_skymean, color=255, psym=-1
oplot, midjulian(h), charon_blue_skymean, color=100, psym=-1
oplot, midjulian(h), charon_nir_skymean, color=240, psym=-1
oplot, midjulian(h), charon_ch4_skymean, color=140, psym=-1


plot, corr_calib_blue_flux, charsize=2
oplot, corr_calib_red_flux, color=255
oplot, corr_calib_blue_flux, color=100
oplot, corr_calib_nir_flux, color=240
oplot, corr_calib_ch4_flux, color=140

loadct, 0, ncolors=256, /silent
plot, midjulian(h), corr_calib_blue_flux(h), charsize=2, $
      XTICKFORMAT='LABEL_DATE'
loadct, 13, ncolors=256, /silent
oplot, midjulian(h), corr_calib_red_flux(h), color=255
oplot, midjulian(h), corr_calib_blue_flux(h), color=100
oplot, midjulian(h), corr_calib_nir_flux(h), color=240
oplot, midjulian(h), corr_calib_ch4_flux(h), color=140



plot, corr_charon_calib_blue_flux, charsize=2
oplot, corr_charon_calib_red_flux, color=255
oplot, corr_charon_calib_blue_flux, color=100
oplot, corr_charon_calib_nir_flux, color=240
oplot, corr_charon_calib_ch4_flux, color=140

loadct, 0, ncolors=256, /silent
plot, midjulian(h), corr_charon_calib_blue_flux(h), charsize=2, $
      XTICKFORMAT='LABEL_DATE'
loadct, 13, ncolors=256, /silent
oplot, midjulian(h), corr_charon_calib_red_flux(h), color=255
oplot, midjulian(h), corr_charon_calib_blue_flux(h), color=100
oplot, midjulian(h), corr_charon_calib_nir_flux(h), color=240
oplot, midjulian(h), corr_charon_calib_ch4_flux(h), color=140


chlon=charon_lon(bc)
chside=side(bc)
rflux=corr_charon_calib_red_flux(bc)
bflux=corr_charon_calib_blue_flux(bc)
nflux=corr_charon_calib_nir_flux(bc)
cflux=corr_charon_calib_ch4_flux(bc)

rfluxerr=corr_charon_calib_red_flux_err(bc)
bfluxerr=corr_charon_calib_blue_flux_err(bc)
nfluxerr=corr_charon_calib_nir_flux_err(bc)
cfluxerr=corr_charon_calib_ch4_flux_err(bc)

good=where(chside eq 0) 
bad=where(chside eq 1)


fitted_lons=findgen(360)
nterms=2

fourfit, chlon/360.0, rflux, rfluxerr, nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,rfluxfit
fourfit, chlon/360.0, bflux, bfluxerr, nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,bfluxfit
fourfit, chlon/360.0, nflux, nfluxerr, nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,nfluxfit
fourfit, chlon/360.0, cflux, cfluxerr, nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,cfluxfit

fourfit, chlon(good)/360.0, rflux(good), rfluxerr(good), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,rgfluxfit
fourfit, chlon(good)/360.0, bflux(good), bfluxerr(good), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,bgfluxfit
fourfit, chlon(good)/360.0, nflux(good), nfluxerr(good), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,ngfluxfit
fourfit, chlon(good)/360.0, cflux(good), cfluxerr(good), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,cgfluxfit

fourfit, chlon(bad)/360.0, rflux(bad), rfluxerr(bad), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,rbfluxfit
fourfit, chlon(bad)/360.0, bflux(bad), bfluxerr(bad), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,bbfluxfit
fourfit, chlon(bad)/360.0, nflux(bad), nfluxerr(bad), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,nbfluxfit
fourfit, chlon(bad)/360.0, cflux(bad), cfluxerr(bad), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,cbfluxfit

stop



loadct, 0, ncolors=256, /silent
plot,charon_lon(bc), corr_charon_calib_blue_flux(bc), charsize=2, /nodata
loadct, 13, ncolors=256, /silent
oplot, charon_lon(bc), corr_charon_calib_red_flux(bc), color=255, psym=1
oplot, charon_lon(bc), corr_charon_calib_blue_flux(bc), color=100, psym=1
oplot, charon_lon(bc), corr_charon_calib_nir_flux(bc), color=240, psym=1
oplot, charon_lon(bc), corr_charon_calib_ch4_flux(bc), color=140, psym=1

for n=0, np-1 do oplot, [charon_lon(bc(n)), charon_lon(bc(n))], $
                        [corr_charon_calib_red_flux(bc(n))-corr_charon_calib_red_flux_err(bc(n)), corr_charon_calib_red_flux(bc(n))+corr_charon_calib_red_flux_err(bc(n))], color=255

for n=0, np-1 do oplot, [charon_lon(bc(n)), charon_lon(bc(n))], $
                        [corr_charon_calib_blue_flux(bc(n))-corr_charon_calib_blue_flux_err(bc(n)), corr_charon_calib_blue_flux(bc(n))+corr_charon_calib_blue_flux_err(bc(n))], color=100

for n=0, np-1 do oplot, [charon_lon(bc(n)), charon_lon(bc(n))], $
                        [corr_charon_calib_nir_flux(bc(n))-corr_charon_calib_nir_flux_err(bc(n)), corr_charon_calib_nir_flux(bc(n))+corr_charon_calib_nir_flux_err(bc(n))], color=240

for n=0, np-1 do oplot, [charon_lon(bc(n)), charon_lon(bc(n))], $
                        [corr_charon_calib_ch4_flux(bc(n))-corr_charon_calib_ch4_flux_err(bc(n)), corr_charon_calib_ch4_flux(bc(n))+corr_charon_calib_ch4_flux_err(bc(n))], color=140

oplot, fitted_lons, rfluxfit, color=255
oplot, fitted_lons, bfluxfit, color=100
oplot, fitted_lons, nfluxfit, color=240
oplot, fitted_lons, cfluxfit, color=140


loadct, 0, ncolors=256, /silent
plot,charon_lon(bc), corr_charon_calib_red_flux(bc), charsize=2, /nodata, $
     yrange=[8e-13, 13e-13], ystyle=1, xrange=[0, 360], xstyle=1, title='Red'
loadct, 13, ncolors=256, /silent
;oplot, charon_lon(bc), corr_charon_calib_red_flux(bc), color=255, psym=1

for n=0, np-1 do if side(bc(n)) eq 1 then oplot, [charon_lon(bc(n)),charon_lon(bc(n))], [corr_charon_calib_red_flux(bc(n)),corr_charon_calib_red_flux(bc(n))], psym=4, color=255
for n=0, np-1 do if side(bc(n)) eq 0 then oplot, [charon_lon(bc(n)),charon_lon(bc(n))], [corr_charon_calib_red_flux(bc(n)),corr_charon_calib_red_flux(bc(n))], psym=6

for n=0, np-1 do oplot, [charon_lon(bc(n)), charon_lon(bc(n))], $
                        [corr_charon_calib_red_flux(bc(n))-corr_charon_calib_red_flux_err(bc(n)), corr_charon_calib_red_flux(bc(n))+corr_charon_calib_red_flux_err(bc(n))], color=255

oplot, fitted_lons, rfluxfit, color=255
oplot, fitted_lons, rgfluxfit, color=255, linestyle=1

;for n=0, np-1 do if side(bc(n)) eq 1 then oplot, [charon_lon(bc(n)), charon_lon(bc(n))], !y.crange, linestyle=1


;if side(n) eq 0 then symv=4
;if side(n) eq 1 then symv=6

     


plot, fitted_lons, rfluxfit,/nodata
oplot, fitted_lons, rfluxfit, color=255
oplot, fitted_lons, bfluxfit, color=100
oplot, fitted_lons, nfluxfit, color=240
oplot, fitted_lons, cfluxfit, color=140

oplot, fitted_lons, rgfluxfit, color=255, linestyle=1
oplot, fitted_lons, bgfluxfit, color=100, linestyle=1
oplot, fitted_lons, ngfluxfit, color=240, linestyle=1
oplot, fitted_lons, cgfluxfit, color=140, linestyle=1






fourfit, charon_lon(h)/360.0, corr_charon_blue_flux(h), corr_charon_blue_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_blue_fit_flux









   
stop

corr_red_flux=red_flux*(new_sc_range^2.0)
corr_blue_flux=blue_flux*(new_sc_range^2.0)
corr_nir_flux=nir_flux*(new_sc_range^2.0)
corr_ch4_flux=ch4_flux*(new_sc_range^2.0)

corr_red_flux_err=red_flux_err*(new_sc_range^2.0)
corr_blue_flux_err=blue_flux_err*(new_sc_range^2.0)
corr_nir_flux_err=nir_flux_err*(new_sc_range^2.0)
corr_ch4_flux_err=ch4_flux_err*(new_sc_range^2.0)


corr_charon_red_flux=charon_red_flux*(new_sc_range^2.0)
corr_charon_blue_flux=charon_blue_flux*(new_sc_range^2.0)
corr_charon_nir_flux=charon_nir_flux*(new_sc_range^2.0)
corr_charon_ch4_flux=charon_ch4_flux*(new_sc_range^2.0)

corr_charon_red_flux_err=charon_red_flux_err*(new_sc_range^2.0)
corr_charon_blue_flux_err=charon_blue_flux_err*(new_sc_range^2.0)
corr_charon_nir_flux_err=charon_nir_flux_err*(new_sc_range^2.0)
corr_charon_ch4_flux_err=charon_ch4_flux_err*(new_sc_range^2.0)

stop

corr_red_mag=-2.5*alog10(corr_red_flux)
corr_blue_mag=-2.5*alog10(corr_blue_flux)
corr_nir_mag=-2.5*alog10(corr_nir_flux)
corr_ch4_mag=-2.5*alog10(corr_ch4_flux)

corr_red_mag_err_max=-2.5*alog10(corr_red_flux+corr_red_flux_err)
corr_blue_mag_err_max=-2.5*alog10(corr_blue_flux+corr_blue_flux_err)
corr_nir_mag_err_max=-2.5*alog10(corr_nir_flux+corr_nir_flux_err)
corr_ch4_mag_err_max=-2.5*alog10(corr_ch4_flux+corr_ch4_flux_err)

corr_red_mag_err_min=-2.5*alog10(corr_red_flux-corr_red_flux_err)
corr_blue_mag_err_min=-2.5*alog10(corr_blue_flux-corr_blue_flux_err)
corr_nir_mag_err_min=-2.5*alog10(corr_nir_flux-corr_nir_flux_err)
corr_ch4_mag_err_min=-2.5*alog10(corr_ch4_flux-corr_ch4_flux_err)



corr_charon_red_mag=-2.5*alog10(corr_charon_red_flux)
corr_charon_blue_mag=-2.5*alog10(corr_charon_blue_flux)
corr_charon_nir_mag=-2.5*alog10(corr_charon_nir_flux)
corr_charon_ch4_mag=-2.5*alog10(corr_charon_ch4_flux)

corr_charon_red_mag_err_max=-2.5*alog10(corr_charon_red_flux+corr_charon_red_flux_err)
corr_charon_blue_mag_err_max=-2.5*alog10(corr_charon_blue_flux+corr_charon_blue_flux_err)
corr_charon_nir_mag_err_max=-2.5*alog10(corr_charon_nir_flux+corr_charon_nir_flux_err)
corr_charon_ch4_mag_err_max=-2.5*alog10(corr_charon_ch4_flux+corr_charon_ch4_flux_err)

corr_charon_red_mag_err_min=-2.5*alog10(corr_charon_red_flux-corr_charon_red_flux_err)
corr_charon_blue_mag_err_min=-2.5*alog10(corr_charon_blue_flux-corr_charon_blue_flux_err)
corr_charon_nir_mag_err_min=-2.5*alog10(corr_charon_nir_flux-corr_charon_nir_flux_err)
corr_charon_ch4_mag_err_min=-2.5*alog10(corr_charon_ch4_flux-corr_charon_ch4_flux_err)


nterms=2
h=sort(pluto_lon)
!P.Multi=[0,1,1]


fitted_lons=findgen(361)
fourfit, pluto_lon(h)/360.0, corr_red_flux(h), corr_red_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,red_fit_flux

fourfit, pluto_lon(h)/360.0, corr_blue_flux(h), corr_blue_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,blue_fit_flux

fourfit, pluto_lon(h)/360.0, corr_nir_flux(h), corr_nir_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,nir_fit_flux

fourfit, pluto_lon(h)/360.0, corr_ch4_flux(h), corr_ch4_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,ch4_fit_flux


fourfit, charon_lon(h)/360.0, corr_charon_blue_flux(h), corr_charon_blue_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_blue_fit_flux

fourfit, charon_lon(h)/360.0, corr_charon_red_flux(h), corr_charon_red_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_red_fit_flux

fourfit, charon_lon(h)/360.0, corr_charon_nir_flux(h), corr_charon_nir_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_nir_fit_flux

fourfit, charon_lon(h)/360.0, corr_charon_ch4_flux(h), corr_charon_ch4_flux_err(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_ch4_fit_flux

set_plot, 'x'
plot, pluto_lon, corr_red_flux, psym=1
oplot,fitted_lons, red_fit_flux

plot, pluto_lon, corr_blue_flux, psym=1
oplot,fitted_lons, blue_fit_flux

plot, pluto_lon, corr_nir_flux, psym=1
oplot,fitted_lons, nir_fit_flux

plot, pluto_lon, corr_ch4_flux, psym=1
oplot,fitted_lons, ch4_fit_flux

plot, charon_lon, corr_charon_red_flux, psym=1
oplot, fitted_lons, charon_red_fit_flux

plot, charon_lon, corr_charon_blue_flux, psym=1
oplot, fitted_lons, charon_blue_fit_flux

fourfit, pluto_lon(h)/360.0, corr_red_mag(h), corr_red_mag_err_max(h)-corr_red_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,red_fit_mag
fourfit, pluto_lon(h)/360.0, corr_blue_mag(h), corr_blue_mag_err_max(h)-corr_blue_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,blue_fit_mag
fourfit, pluto_lon(h)/360.0, corr_nir_mag(h), corr_nir_mag_err_max(h)-corr_nir_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,nir_fit_mag
fourfit, pluto_lon(h)/360.0, corr_ch4_mag(h), corr_ch4_mag_err_max(h)-corr_ch4_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c, ch4_fit_mag

fourfit, charon_lon(h)/360.0, corr_charon_red_mag(h), corr_charon_red_mag_err_max(h)-corr_charon_red_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_red_fit_mag
fourfit, charon_lon(h)/360.0, corr_charon_blue_mag(h), corr_charon_blue_mag_err_max(h)-corr_charon_blue_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c, charon_blue_fit_mag
fourfit, charon_lon(h)/360.0, corr_charon_nir_mag(h), corr_charon_nir_mag_err_max(h)-corr_charon_nir_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c,charon_nir_fit_mag
fourfit, charon_lon(h)/360.0, corr_charon_ch4_mag(h), corr_charon_ch4_mag_err_max(h)-corr_charon_ch4_mag(h), nterms, c, csig, yfit=yfit
fourfunc, fitted_lons/360.0, c, charon_ch4_fit_mag




ha=where(side eq 0)
nha=n_elements(ha)
hb=where(side eq 1)
nhb=n_elements(hb)
date_label = LABEL_DATE(DATE_FORMAT = ['%D/%N/%Y'])


setupplot, 'sides_vs_time.ps'
plot, red_midtime, side, psym=1, $
      XTICKFORMAT = 'LABEL_DATE', $
      xtitle='Time', ytitle='Side', $
      xthick=3, ythick=3, thick=3, $
      yrange=[-0.5, 1.5], ystyle=1
device,/close

stop

setupplot, 'pluto_charon_quicklook_lightcurves.ps'
loadct, 13
!P.Multi=[0,2,2]

;************************************************
;Pluto
;************************************************
for na=0, maxn do begin

   if na eq 0 then begin
      x=pluto_lon
      y=corr_red_flux
      yerr=corr_red_flux_err
      yrange=[3,4.5]*1e4
      title='Red'
      fitv=red_fit_flux

      
      if panframe eq 1 then yrange=[9.5,12]*1e4
   endif
   
   if na eq 1 then begin
      x=pluto_lon
      y=corr_blue_flux
      yerr=corr_blue_flux_err
      yrange=[7.5,10.5]*1e3
      yrange=[6.5,11.5]*1e3
      title='Blue'
      fitv=blue_fit_flux
   endif

   if na eq 2 then begin
      x=pluto_lon
      y=corr_nir_flux
      yerr=corr_nir_flux_err
      yrange=[2,4]*1e4
      title='NIR'
      fitv=nir_fit_flux
   endif

   if na eq 3 then begin
      x=pluto_lon
      y=corr_ch4_flux
      yerr=corr_ch4_flux_err
      yrange=[5000,9000]
      title='CH4'
      fitv=ch4_fit_flux
   endif
   
   nlon=n_elements(y)
   
   plot, x, y, psym=1, $
         title=title, $
         ytitle='Flux', $
         yrange=yrange, ystyle=1, $
         xthick=3, ythick=3, charthick=3, charsize=1.0, $
         xtitle='Sub-Observer Longitude (RHR)',  $
         xrange=[0, 360], xticks=12, xstyle=1

   for n=0, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3
      oplot, [x(n), x(n)], [y(n)-yerr(n), y(n)+yerr(n)], thick=3
   endfor

   ;count=0
   ;colorv=[40,100,220,255]
   ;for n=nlon-1-3, nlon-1 do begin
   ;   if side(n) eq 0 then symv=4
   ;   if side(n) eq 1 then symv=6
   ;   
   ;   oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=colorv(count)
   ;   oplot, [x(n), x(n)], [y(n)-yerr(n), y(n)+yerr(n)], thick=3, color=colorv(count)
   ;   ;print, n, count;;

;      count=count+1
;   endfor

   ;cmin=min(tfheater)
   ;cmax=max(tfheater)
   ;dt=256.0/(cmax-cmin)
   ;for n=0, nlon-1 do begin
   ;   col=(tfheater(n)-cmin)*dt
   ;   if side(n) eq 0 then symv=4
   ;   if side(n) eq 1 then symv=6
   ;   oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=col
   ;   
   ;endfor
   oplot, fitted_lons, fitv
   ;XYOUTS, 30, 1.18*1e5, ' 6 Hours after heater on', color=250
   ; XYOUTS, 30, 1.17*1e5, '15 Hours after heater on', color=230
    XYOUTS, 30, 1.15*1e5, 'Side 0 - Diamonds'
    XYOUTS, 30, 1.14*1e5, 'Side 1 - Squares'
endfor

;************************************************
;Charon
;************************************************
for na=0, maxn do begin
   x=charon_lon
   nlon=n_elements(charon_lon)
   
   if na eq 0 then begin
      y=corr_charon_red_flux
      yerr=corr_charon_red_flux_err
      yrange=[4000,6000]
      title='Charon Red'
      fitv=charon_red_fit_flux
      if panframe eq 1 then yrange=[1.3, 1.5]*1e4
   endif
   
   if na eq 1 then begin
      y=corr_charon_blue_flux
      yerr=corr_charon_blue_flux_err
      yrange=[0,3000]
      title='Charon Blue'
      fitv=charon_blue_fit_flux
   endif

   if na eq 2 then begin
      y=corr_charon_nir_flux
      yerr=corr_charon_nir_flux_err
      yrange=[2000,5000]
      title='Charon NIR'
      fitv=charon_nir_fit_flux
   endif

   if na eq 3 then begin
      y=corr_charon_ch4_flux
      yerr=corr_charon_ch4_flux_err
      yrange=[500,1500]
      title='Charon CH4'
      fitv=charon_ch4_fit_flux
   endif
   
   nlon=n_elements(y)
   
   plot, x, y, psym=1, $
         title=title, $
         ytitle='Flux', $
         yrange=yrange, ystyle=1, $
         xthick=3, ythick=3, charthick=3, charsize=1.0, $
         xtitle='Sub-Observer Longitude (RHR)',  $
         xrange=[0, 360], xticks=12, xstyle=1

   for n=0, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3
      oplot, [x(n), x(n)], [y(n)-yerr(n), y(n)+yerr(n)], thick=3
   endfor

   ;count=0
   ;colorv=[40,100,220,255]
   ;for n=nlon-1-3, nlon-1 do begin
   ;   if side(n) eq 0 then symv=4
   ;   if side(n) eq 1 then symv=6
      
   ;   oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=colorv(count)
   ;   oplot, [x(n), x(n)], [y(n)-yerr(n), y(n)+yerr(n)], thick=3, color=colorv(count)
   ;   count=count+1
   ;endfor
   
   ;cmin=min(tfheater)
   ;cmax=max(tfheater)
   ;dt=256.0/(cmax-cmin)
   for n=0, nlon-1 do begin
                                ;   col=(tfheater(n)-cmin)*dt
      col=255
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=col
     
   endfor
  
   
   oplot, fitted_lons, fitv, thick=3


   

   
endfor


;CgColorbar, POSITION=[0.05,0.1,0.08,0.45],$
;            /right,$
;            /VERTICAL,$
;            charthick=2, $
;            charsize=1.5, $
;            title='Time since Heater Off (mins)', $
;            Range=[cmin, cmax], ncolors=256


device,/close

;************************************************
;Pluto - Residuals
;************************************************
for na=0, 3 do begin
   x=pluto_lon
   nx=n_elements(pluto_lon)
   smflux=fltarr(nx)


   if na eq 0 then begin
      for p=0, nx-1 do smflux(p)=interpol(red_fit_flux, fitted_lons, pluto_lon(p))
  
      y=corr_red_flux-smflux
      yerrmin=y-corr_red_flux_err
      yerrmax=y+corr_red_flux_err
      title='Red'
      yrange=[-3000,3000]
   endif
   
   if na eq 1 then begin
      for p=0, nx-1 do smflux(p)=interpol(blue_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_blue_flux-smflux
      yerrmin=y-corr_blue_flux_err
      yerrmax=y+corr_blue_flux_err
      title='Blue'
      yrange=[-2000,2000]
   endif

   if na eq 2 then begin
      for p=0, nx-1 do smflux(p)=interpol(nir_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_nir_flux-smflux
      yerrmin=y-corr_nir_flux_err
      yerrmax=y+corr_nir_flux_err
      title='NIR'
      yrange=[-0.5,0.5]*1e4
   endif

   if na eq 3 then begin
      for p=0, nx-1 do smflux(p)=interpol(ch4_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_ch4_flux-smflux
      yerrmin=y-corr_ch4_flux_err
      yerrmax=y+corr_ch4_flux_err
      title='CH4'
      yrange=[-1000,1000]
   endif
   
   nlon=n_elements(y)
   
   plot, x, y, psym=1, $
         title=title, $
         ytitle='Flux Difference From Trend Line', $
         yrange=yrange, ystyle=1, $
         xthick=3, ythick=3, charthick=3, charsize=1.0, $
         xtitle='Sub-Observer Longitude (RHR)',  $
         xrange=[0, 360], xticks=12, xstyle=1

   for n=0, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3
      oplot, [x(n), x(n)], [yerrmin(n), yerrmax(n)], thick=3
   endfor

   count=0
   colorv=[40,100,220,255]
   for n=nlon-1-3, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=colorv(count)
      oplot, [x(n), x(n)], [yerrmin(n), yerrmax(n)], thick=3, color=colorv(count)
      ;print, n, count

      count=count+1
   endfor
endfor


;************************************************
;Charon
;************************************************
for na=0, 3 do begin
   x=charon_lon
   nx=n_elements(charon_lon)
   smflux=fltarr(nx)
   
   if na eq 0 then begin
      for p=0, nx-1 do smflux(p)=interpol(charon_red_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_charon_red_flux-smflux
      yerr_min=y-corr_charon_red_flux_err
      yerr_max=y+corr_charon_red_flux_err
      yrange=[4000,6000]
      title='Charon Red'
   endif
   
   if na eq 1 then begin
      for p=0, nx-1 do smflux(p)=interpol(charon_blue_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_charon_blue_flux-smflux
      yerr_min=y-corr_charon_blue_flux_err
      yerr_max=y+corr_charon_blue_flux_err
      title='Charon Blue'
   endif

   if na eq 2 then begin
      for p=0, nx-1 do smflux(p)=interpol(charon_nir_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_charon_nir_flux-smflux
      yerr_min=y-corr_charon_nir_flux_err
      yerr_max=y+corr_charon_nir_flux_err
      title='Charon NIR'
   endif

   if na eq 3 then begin
      for p=0, nx-1 do smflux(p)=interpol(charon_ch4_fit_flux, fitted_lons, pluto_lon(p))
      y=corr_charon_ch4_flux-smflux
      yerr_min=y-corr_charon_ch4_flux_err
      yerr_max=y+corr_charon_ch4_flux_err
      title='Charon CH4'
   endif
   
   nlon=n_elements(y)
   
   plot, x, y, psym=1, $
         title=title, $
         ytitle='Flux Difference From Trend Line', $
         yrange=[-1000, 1000], ystyle=1, $
         xthick=3, ythick=3, charthick=3, charsize=1.0, $
         xtitle='Sub-Observer Longitude (RHR)',  $
         xrange=[0, 360], xticks=12, xstyle=1

   for n=0, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3
      oplot, [x(n), x(n)], [yerr_min(n), yerr_max(n)], thick=3
   endfor

   count=0
   colorv=[40,100,220,255]
   for n=nlon-1-3, nlon-1 do begin
      if side(n) eq 0 then symv=4
      if side(n) eq 1 then symv=6
      
      oplot, [x(n), x(n)], [y(n), y(n)], psym=symv, thick=3, color=colorv(count)
      oplot, [x(n), x(n)], [yerr_min(n), yerr_max(n)], thick=3, color=colorv(count)
      count=count+1
   endfor
endfor

!P.multi=[0,1,2]

Plot, pluto_lon, red_skymean, xtitle='Longitude', ytitle='Sky Mean', /nodata, yrange=[-5.0, 5]
oplot, pluto_lon, red_skymean, color=255, psym=1, thick=3
for n=0, nlon-1 do oplot, [pluto_lon(n), pluto_lon(n)], [red_skymean(n)-red_skysig(n), red_skymean(n)+red_skysig(n)], color=255
oplot, pluto_lon, blue_skymean, color=100, psym=1, thick=3
for n=0, nlon-1 do oplot, [pluto_lon(n), pluto_lon(n)], [blue_skymean(n)-blue_skysig(n), blue_skymean(n)+blue_skysig(n)], color=100
oplot, pluto_lon, nir_skymean, color=230, psym=1, thick=3
for n=0, nlon-1 do oplot, [pluto_lon(n), pluto_lon(n)], [nir_skymean(n)-nir_skysig(n), nir_skymean(n)+nir_skysig(n)], color=230
oplot, pluto_lon, ch4_skymean, color=140, psym=1, thick=3
for n=0, nlon-1 do oplot, [pluto_lon(n), pluto_lon(n)], [ch4_skymean(n)-ch4_skysig(n), ch4_skymean(n)+ch4_skysig(n)], color=140



date_label = LABEL_DATE(DATE_FORMAT = ['%D/%N/%Y'])
plot, red_midtime, red_skymean, xtitle='Date', ytitle='Sky Mean',XTICKFORMAT = 'LABEL_DATE', /nodata, yrange=[-6,6]
oplot, red_midtime, red_skymean, color=255, psym=1
for n=0, nlon-1 do oplot, [red_midtime(n), red_midtime(n)], [red_skymean(n)-red_skysig(n), red_skymean(n)+red_skysig(n)], color=255
oplot, red_midtime, blue_skymean, color=100, psym=1
for n=0, nlon-1 do oplot, [blue_midtime(n), blue_midtime(n)], [blue_skymean(n)-blue_skysig(n), blue_skymean(n)+blue_skysig(n)], color=100
oplot, red_midtime, nir_skymean, color=230, psym=1
for n=0, nlon-1 do oplot, [nir_midtime(n), nir_midtime(n)], [nir_skymean(n)-nir_skysig(n), nir_skymean(n)+nir_skysig(n)], color=230
oplot, red_midtime, ch4_skymean, color=140, psym=1
for n=0, nlon-1 do oplot, [ch4_midtime(n), ch4_midtime(n)], [ch4_skymean(n)-ch4_skysig(n), ch4_skymean(n)+ch4_skysig(n)], color=140


date_label = LABEL_DATE(DATE_FORMAT = ['%D/%N/%Y'])
plot, red_midtime, red_fwhm, xtitle='Date', ytitle='FWHM',XTICKFORMAT = 'LABEL_DATE', /nodata, yrange=[1,4.5], ystyle=1
oplot, red_midtime, red_fwhm, color=255, psym=-1
oplot, red_midtime, blue_fwhm, color=100, psym=-1
oplot, red_midtime, nir_fwhm, color=230, psym=-1
oplot, red_midtime, ch4_fwhm, color=140, psym=-1

device,/close


;blue_yfit_sm = GAUSSFIT(fitx,fity, coeff2, NTERMS=nterms2)



save, filename='phot_fit.idlsave', $
      fitted_lons, red_fit_mag, blue_fit_mag, nir_fit_mag, ch4_fit_mag, $
      charon_red_fit_mag, charon_blue_fit_mag, charon_nir_fit_mag, charon_ch4_fit_mag,$
      pluto_lon, pluto_lat, charon_lon, charon_lat, $
      corr_red_mag, corr_red_mag_err_min, corr_red_mag_err_max, $
      corr_blue_mag, corr_blue_mag_err_min, corr_blue_mag_err_max, $
      corr_charon_red_mag, corr_charon_red_mag_err_min, corr_charon_red_mag_err_max, $
      corr_charon_blue_mag, corr_charon_blue_mag_err_min, corr_charon_blue_mag_err_max, red_midtime


;Save output files
save, filename=savefilename, $
      onbefore, $
      temp_filename, temp_met, cdhtempwarmrad, cdhtempcoldrad, cdhtempbulk, cdhtempsc, cdhtempelecbox, temp_julian, $
      fitted_lons, red_fit_mag, blue_fit_mag, charon_red_fit_mag, charon_blue_fit_mag, $
      sc_range, pluto_lat, pluto_lon, charon_lat, charon_lon, $
      side, met, $
      nterms, red_fit_flux, blue_fit_flux, nir_fit_flux, ch4_fit_flux, $
      red_fwhm, red_flux, red_flux_err, red_err, $
      red_mag, red_max, red_skymean, red_skyerr, red_skysig, red_xcen, $
      red_ycen, red_exptime, red_midtime, $
      blue_fwhm, blue_flux, blue_flux_err, blue_err, $
      blue_mag, blue_max, blue_skymean, blue_skyerr, blue_skysig, blue_xcen, $
      blue_ycen, blue_exptime, blue_midtime, $
      nir_fwhm, nir_flux, nir_flux_err, nir_err, $
      nir_mag, nir_max, nir_skymean, nir_skyerr, nir_skysig, nir_xcen, $
      nir_ycen, nir_exptime, nir_midtime, $
      ch4_fwhm, ch4_flux, ch4_flux_err, ch4_err, $
      ch4_mag, ch4_max, ch4_skymean, ch4_skyerr, ch4_skysig, ch4_xcen, $
      ch4_ycen, ch4_exptime, ch4_midtime, $
      corr_red_flux, corr_blue_flux, corr_nir_flux, corr_ch4_flux, $
      corr_red_flux_err, corr_blue_flux_err, corr_nir_flux_err, corr_ch4_flux_err, $
      corr_red_mag, corr_blue_mag, corr_nir_mag, corr_ch4_mag, $
      corr_red_mag_err_min, corr_blue_mag_err_min, corr_nir_mag_err_min ,corr_ch4_mag_err_min, $
      corr_red_mag_err_max, corr_blue_mag_err_max, corr_nir_mag_err_max ,corr_ch4_mag_err_max, $
      charon_red_fwhm, charon_red_flux, charon_red_flux_err, charon_red_err, $
      charon_red_mag, charon_red_max, charon_red_skymean, charon_red_skyerr, charon_red_skysig, charon_red_xcen, charon_red_ycen, $
      charon_blue_fwhm, charon_blue_flux, charon_blue_flux_err, charon_blue_err, $
      charon_blue_mag, charon_blue_max, charon_blue_skymean, charon_blue_skyerr,charon_blue_skysig, charon_blue_xcen, charon_blue_ycen, $
      charon_nir_fwhm, charon_nir_flux, charon_nir_flux_err,charon_nir_err, $
      charon_nir_mag, charon_nir_max, charon_nir_skymean, charon_nir_skyerr, charon_nir_skysig, charon_nir_xcen, charon_nir_ycen, $
      charon_ch4_fwhm, charon_ch4_flux, charon_ch4_flux_err, charon_ch4_err, $
      charon_ch4_mag, charon_ch4_max, charon_ch4_skymean, charon_ch4_skyerr, charon_ch4_skysig, charon_ch4_xcen, charon_ch4_ycen, $
      corr_charon_red_flux, corr_charon_blue_flux, corr_charon_nir_flux, corr_charon_ch4_flux, $
      corr_charon_red_flux_err, corr_charon_blue_flux_err, corr_charon_nir_flux_err, corr_charon_ch4_flux_err, $
      corr_charon_red_mag, corr_charon_blue_mag, corr_charon_nir_mag, corr_charon_ch4_mag, $
      corr_charon_red_mag_err_min, corr_charon_blue_mag_err_min, corr_charon_nir_mag_err_min ,corr_charon_ch4_mag_err_min, $
      corr_charon_red_mag_err_max, corr_charon_blue_mag_err_max, corr_charon_nir_mag_err_max ,corr_charon_ch4_mag_err_max

end


;M6/M& star, choose an M6/7 star with good signal
;Solar star HD (Solar keyword) 


;-----------------
;Load MVIC throughputs
;-----------------
   throughputdir='/Users/carlyhowett1/Desktop/New_Horizons/Flight/Commissioning/newCHI/'
   red_qq=readfits(throughputdir + 'RED_NEW_CHI.FITS')
   lambda_A=reform(red_qq(0,*)*1e4)
   red_chi=reform(red_qq(1,*))
   
   blu_qq=readfits(throughputdir + 'BLUE_NEW_CHI.FITS')
   lambda_A=reform(blu_qq(0,*)*1e4)
   blu_chi=reform(blu_qq(1,*))
   
   nir_qq=readfits(throughputdir + 'NIR_NEW_CHI.FITS')
   lambda_A=reform(nir_qq(0,*)*1e4)
   nir_chi=reform(nir_qq(1,*))
   
   ch4_qq=readfits(throughputdir + 'CH4_NEW_CHI.FITS')
   lambda_A=reform(ch4_qq(0,*)*1e4)
   ch4_chi=reform(ch4_qq(1,*))
   
   pan_qq=readfits(throughputdir + 'PAN_1_NEW_CHI.FITS')
   lambda_A=reform(pan_qq(0,*)*1e4)
   pan_chi=reform(pan_qq(1,*))
   
;Plot the throughputs
   set_plot, 'x'
   plot, lambda_A, red_chi
   oplot, lambda_A, blu_chi
   oplot, lambda_A, nir_chi
   oplot, lambda_A, ch4_chi
   oplot, lambda_A, pan_chi




;readcol, '/Users/carlyhowett1/Desktop/Hubble/UVIS/F555W_UVIS_throughput.csv', num, wave1, chi1, $
;         wave2, chi2, format='f,f,f,f,f', skipline=1

;Extends beyond Blue filter coverage in MVIC
readcol, '/Users/carlyhowett1/Desktop/Hubble/ACS/wfc_F435W.dat', wave1, chi1, format='f,f' ;B

;Within the MVIC wavelength coverage
;readcol, '/Users/carlyhowett1/Desktop/Hubble/ACS/wfc_F555W.dat', wave1, chi1, format='f,f' ;V

oplot, wave1, chi1


filter_throughput=interpol(chi1, wave1, lambda_A)



target='charon'

if target eq 'pluto' then filename='/Users/carlyhowett1/Desktop/New_Horizons/Software/howett/spectrum/pluto_spectrum.dat'
if target eq 'charon' then  filename='/Users/carlyhowett1/Desktop/New_Horizons/Software/howett/spectrum/charon_spectrum.dat'
if target eq 'pholus' then filename='/Users/carlyhowett1/Desktop/New_Horizons/Software/howett/spectrum/pholus_spectrum.dat'
if target eq 'jupiter' then filename='/Users/carlyhowett1/Desktop/New_Horizons/Software/howett/spectrum/jupiter_spectrum.dat'
readcol,filename,wave_microns,flux ;in microns
target_spec0=flux
target_wave=wave_microns*1e4
target_spectrum=interpol(flux, target_wave, lambda_A)



;----------------------------------------
;Get Kurucz model spectrum for Vega
;----------------------------------------
;Assuming Vega's (α Lyrae) temperature is  9602 ± 180 (from Kinman and
;Castelli, 2002)

star_temp=9602D0
VALID_KURUCZ2, kurucz_t, logg, logZ=0       ;Get kurucz temperature array
neart_index = indx2(kurucz_t, star_temp)    ;find nearest index to star_temp[i]

;OUTPUT from KURUCZ
;	WAVE - Wavelength array (Angstroms)
;	FLAM - (F Lambda) Flux (ergs/cm2/s/A).  To convert to flux at Earth, 
;		multiply by !PI*(R/d)^2  where R is the stellar radius, and d 
;		is the stellar distance.  FLAM will have the same number of 
;		elements as WAVE
KURUCZ, lambda_A, vega_bbflux, kurucz_t[neart_index], logg, 0
nlambda=n_elements(lambda_A)
openw, 1, 'kuruz_flux_star9602K.txt'
printf, 1, 'Kuruz model for Star Temp 9602 K (Vega 9602+/-180 Kinman and Castelli, 2002)'
printf, 1, 'Wavelength (A)   Flux (erg/cm2/s/A)'
for n=0, nlambda-1 do printf, 1, lambda_A(n), ' ', vega_bblux(n)
close, 1

KURUCZ, 5556., bbflux5556, kurucz_t[neart_index],logg, 0




bbflux5556=17468704D0

stop


;Normalize instead by using the STIS HST value of Vega (3.44e-9 at 5556A)
norm=3.44e-9/bbflux5556
vega_bbflux=vega_bbflux*norm

;stop



;d=25.04D0*9.46e12                    ;25.04 light years (in km)
;R=mean([2.362D0,2.818D0])*6.957e5    ;radius of Vega in Km (from Rsolar = 6.957e10^5 km)

;new_flux=vega_bbflux*!pi*((R/d)^2)
;new=bbflux5556*!pi*((R/d)^2)
;vega_bbflux=new_flux

plot, lambda_A, vega_bbflux



;-----------------
;Load Charon color (from central grey part of Charon, from P_CHARON2)
;-----------------
pivot_w=[4920.00, 6240.00, 8610.00];, 8830.00] (don't use CH4)
;Charon_color=[0.0127237,0.0111456,0.00662754];,0.00676188]


restore, filename='quicklook8_prenewscrange.idlsave'
;In this external routine the flux of charon is defined as:
;calib_charon_red_flux=charon_red_counts/(red_exptime*PCHARON_RED)
nobsused=n_elements(sc_range)

obs0=0
obs1=nobsused-1

;Set up output arrays
mag_target_red=fltarr(obs1-obs0+1)
mag_target_blu=mag_target_red
mag_target_nir=mag_target_red
mag_target_ch4=mag_target_red
mag_target_pan=mag_target_red
mag_oth_target=mag_target_red
rangeout=fltarr(obs1-obs0+1)

flux_vega=fltarr(obs1-obs0+1, 4)
flux_charon=fltarr(obs1-obs0+1, 4)


;plot, lambda_A, target_spectrum, charsize=2
;oplot, lambda_A, BLU_CHI
;oplot, lambda_A,RED_CHI
;oplot, lambda_A, NIR_CHI


;Run through each of the MVIC Charon unresolved observations
set_plot, 'x'
setdecomposedstate, 0
for obsused=obs0, obs1 do begin

   ;----------------------------------------------
   ;Find target color using MVIC
   ;----------------------------------------------
   
   ;------------------------------
   ;For this observation get the range and color of Charon as seen by MVIC
   ;------------------------------
   range=sc_range(obsused)
   Charon_color=[calib_charon_blue_flux(obsused),calib_charon_red_flux(obsused), calib_charon_nir_flux(obsused)]
   rangeout(obsused)=range

   ;------------------------------
   ;Plot the color
   ;------------------------------
   loadct, 0
   !P.Multi=[0,1,1]
   plot, pivot_w, Charon_color, charsize=2, xrange=[min(lambda_A), max(lambda_A)], $
         thick=3, psym=-1

   ;------------------------------
   ;Normalize the Charon spectra at the red pivot wavelength of MVIC
   ;------------------------------
   flux_red=interpol(target_spectrum, lambda_A, pivot_w(1))
                                ;Find the normalization factor that
                                ;normalizes the Charon spectra at the
                                ;Red pivot wavelength
   normfac=Charon_color(1)/flux_red
   ;Normalize the spectra & plot it
   norm_target_spectrum=target_spectrum*normfac
   oplot, lambda_A, norm_target_spectrum
   stop

   ;------------------------------
   ;Find the total flux weighted by the filter throughputs
   ;needed to determine the target's magnitude as seen by MVIC
   ;------------------------------
   mvic_blu_target_flux=int_tabulated(lambda_A,norm_target_spectrum*BLU_CHI)
   mvic_red_target_flux=int_tabulated(lambda_A,norm_target_spectrum*RED_CHI)
   mvic_nir_target_flux=int_tabulated(lambda_A,norm_target_spectrum*NIR_CHI)


   ;----------------------------------------------
   ;Other filter
   ;----------------------------------------------
   
   ;------------------------------
   ;Find flux weighted by the filter throughputs   
   ;(These are the fluxes at the distance
   ;of New Horizons and the Solar distance at the time of the observation)
   ;------------------------------
   otherfilter_target_spec=norm_target_spectrum*filter_throughput
   oth_target_flux=int_tabulated(lambda_A, otherfilter_target_spec)

   ;------------------------------
   ;Adjust flux for the changing distance of the Sun from target
   ;and the observer distance
   ;------------------------------
   r1=30.517 ;heliocentric distance (AU) from Marc's (2010) paper
   delta1=29.521 ;observer distance (AU) from Marc's (2010) paper

   r1=39.5
   delta1=38.5

   
   ;target-observer and target-Sun distances from the headers in km
   oneAUinkm=149597870.7D0
   r2=(targ_sun(obsused))/oneAUinkm ;from SPCTSORN
   delta2=sc_range(obsused)/oneAUinkm ;from SPCTRANG

   helio_dist_norm=(r2/r1)^2
   obsvr_dist_norm=(delta2/delta1)^2


   oth_target_flux_norm=oth_target_flux*helio_dist_norm*obsvr_dist_norm
   

   ;----------------------------------------------
   ;Vega
   ;----------------------------------------------

   ;-----------------
   ;Multiple Kurucz model by throughput from MVIC and other filter
   ;-----------------
   red_filter_vega=vega_bbflux*red_chi
   blu_filter_vega=vega_bbflux*blu_chi
   nir_filter_vega=vega_bbflux*nir_chi

   oth_filter_vega=vega_bbflux*filter_throughput
   
   ;-----------------
   ;Find total flux
   ;-----------------
   red_flux_vega=int_tabulated(lambda_A, red_filter_vega)
   blu_flux_vega=int_tabulated(lambda_A, blu_filter_vega)
   nir_flux_vega=int_tabulated(lambda_A, nir_filter_vega)

   oth_flux_vega=int_tabulated(lambda_A, oth_filter_vega)


   ;----------------------------------------------
   ;Find Magnitudes
   ;----------------------------------------------


   ;Use equation: Mag(Charon in MVIC) = -2.5 log_10 (Fc/F0)

   mag_target_red(obsused)=(-2.5)*alog10(mvic_red_target_flux/red_flux_vega)
   mag_target_blu(obsused)=(-2.5)*alog10(mvic_blu_target_flux/blu_flux_vega)
   mag_target_nir(obsused)=(-2.5)*alog10(mvic_nir_target_flux/nir_flux_vega)
  
   mag_oth_target(obsused)=(-2.5)*alog10(oth_target_flux_norm/oth_flux_vega)

  
   
   print, 'Magnitudes: Blue, Red, NIR, and other filter'
   print, mag_target_blu(obsused), mag_target_red(obsused), mag_target_nir(obsused), mag_oth_target(obsused)

   stop



endfor


stop


;range*20e-6=radius

resrange=(606*2)/20e-6

setupplot, 'Charon_mag_withrange.ps'
!P.Multi=[0,1,1]
loadct, 0
plot, rangeout, mag_charon_blu, charsize=1.5, xrange=[1.2e8, 0], yrange=[-3,6], psym=-1, $
      ytitle='V Magnitudes', xtitle='Range (km)', title='Charon'
oplot, [resrange, resrange], !y.crange, linestyle=1, thick=3
XYOUTS, 5.8e7, 5.5, 'Range at which Charon is resolved', charsize=1.25
loadct, 13
oplot, rangeout, mag_charon_blu, psym=-1, color=50
oplot, rangeout, mag_charon_red, psym=-1, color=255
oplot, rangeout, mag_charon_nir, psym=-1, color=240
;oplot, rangeout, mag_charon_ch4, psym=-1, color=140


XYOUTS, 1e8, 0, 'Blu', color=50, charsize=1.5
XYOUTS, 1e8, -0.5, 'Red', color=255, charsize=1.5
XYOUTS, 1e8, -1.0, 'NIR', color=240, charsize=1.5
;XYOUTS, 1e8, -1.5, 'CH4', color=140, charsize=1.5


;plot, rangeout, flux_vega(0,*)



device,/close

end

[INFO]    :: Logging to /Users/kryanhinds/sedm_phot/grb260310a_final/grb260310a_final_log.log
---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------
[INFO]    :: For SDSS-G, there are 38 fits
[INFO]    :: For SDSS-G, there are 38 fits
[█░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 1/38 (2%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 1 OF 38 [1/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 1 OF 38 [1/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_06_50_20_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=18.2, thresholds: std>36.3, med>555.3) ~185 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.47 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (819,335)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61111.285
[INFO]    :: Observation time: 2026-03-12T06:50:20.767433
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 9 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.58", ΔDec=1.24" from 9 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-3.58", ΔDec=1.24") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 7 sources in sci
[INFO]    :: [V2] star_match_shift: 320 sources in ref
[INFO]    :: [V2] star_match_shift: 7 matches → dx=1.624±0.869 px, dy=1.524±1.259 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (6 inliers) → dx=1.736±0.638 px, dy=1.653±0.659 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.736, 1.653) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1977 clipped px (>= 9.95 nMgy), dilated by 14px -> 9487 masked (1.2% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 90.8% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.122 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_44662.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_44662.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_44662.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_44662.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.138, -0.154) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (41, 41) (FWHM=6.66px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-12T06_24620ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 45/120 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_44662.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.03634909
[INFO]    :: [V2] PSF centroid corrected: shift=(0.021, -0.283) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (41, 41) (FWHM=6.66px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-12T06_24620sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -115.566345 -0.707671 209.212513
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 209
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 8
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-12T06comb_psf_44662.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[29.31175423 29.35856439 29.31831004 29.28895358]
[0.99925786 0.99810691 0.99974746 0.99901356]
[31.66072299 31.74302969 31.69449513 31.6393128 ]
[0.9994475  0.99835816 0.9997148  0.99965293]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3183 |  31.6945 |  0.999747 |  0.999715 |
|  29.3118 |  31.6607 |  0.999258 |  0.999448 |
|  29.289  |  31.6393 |  0.999014 |  0.999653 |
|  29.3586 |  31.743  |  0.998107 |  0.998358 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 29.319
[INFO]    ::   - ZP Std: 0.025
[INFO]    ::   - ZP Range: [29.289,29.359]
[INFO]    ::   - ZP Mode: 29.312
[INFO]    ::   - ZP Median: 29.315
[INFO]    ::   - ZP 16th & 84th percentiles: 29.300, 29.339
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.684
[INFO]    ::   - ZP Std: 0.039
[INFO]    ::   - ZP Range: [31.639,31.743]
[INFO]    ::   - ZP Mode: 31.661
[INFO]    ::   - ZP Median: 31.678
[INFO]    ::   - ZP 16th & 84th percentiles: 31.650, 31.720
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=29.315 std=0.025 No. stars=4
[INFO]    :: ZP Ref=31.678 std=0.039 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.113
[INFO]    :: [V2] Blanked subtraction at 9487 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-12T06_24620bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61111.28496249998
[INFO]    :: BACKGROUND: -1.685 counts
[INFO]    :: SN FLUX 22632.917 counts
[INFO]    :: SN MAG 18.428 mag
[INFO]    :: SN MAG - BACKGROUND 18.428 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.391 0.922
[INFO]    :: xoff_arc=0.145 arcsec, yoff_arc=0.342 arcsec
[INFO]    :: [V2] Injection pool: 713582 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 102.817
[INFO]    :: S/N (std artifical sn)= 102.817
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.266 mag
[INFO]    :: Mag = 18.428+/-0.011 
[INFO]    :: 5-sig limit = 21.711
[INFO]    :: 3-sig limit = 22.266
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.155 arcsec (0.418 pixels) 
[INFO]    :: DEC error: 0.198 arcsec (0.534 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.039
[INFO]    :: Science zeropoint std: 0.025
--------------------------------------------------------------------------------
 Mag = 18.428+/-0.048 lim=22.266 MJD=61111.285
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  698.4
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-12T06_24620_photometry.txt
[INFO]    :: Total time: 8.0 seconds
✓ COMPLETED: 1/38 | 37 remaining
[██░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 2/38 (5%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 2 OF 38 [2/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 2 OF 38 [2/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260312_10_44_01_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=18.4, thresholds: std>36.8, med>562.1) ~191 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.72 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (822,346)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61111.447
[INFO]    :: Observation time: 2026-03-12T10:44:01.251439
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 6 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-10.07", ΔDec=-2.98" from 6 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-10.07", ΔDec=-2.98") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 10 sources in sci
[INFO]    :: [V2] star_match_shift: 334 sources in ref
[INFO]    :: [V2] star_match_shift: 10 matches → dx=1.299±0.363 px, dy=1.069±2.119 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (8 inliers) → dx=1.299±0.290 px, dy=0.827±0.956 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.299, 0.827) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1967 clipped px (>= 9.95 nMgy), dilated by 10px -> 6798 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 90.6% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.118 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06112.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06112.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06112.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_06112.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.0
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.068, -0.132) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.64px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-12T10_38641ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/114 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_06112.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01813793
[INFO]    :: [V2] PSF centroid corrected: shift=(0.011, -0.251) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.64px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-12T10_38641sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -117.203475 -0.753481 212.017244
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 212
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 5(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 7
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-12T10comb_psf_06112.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[29.31213323 29.35626418 29.31929677 29.26277144]
[0.99979828 0.99792747 0.99990882 0.99980776]
[31.64242147 31.7161237  31.67987705 31.63069536]
[0.99956799 0.99856407 0.99970202 0.99933108]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3193 |  31.6799 |  0.999909 |  0.999702 |
|  29.2628 |  31.6307 |  0.999808 |  0.999331 |
|  29.3121 |  31.6424 |  0.999798 |  0.999568 |
|  29.3563 |  31.7161 |  0.997927 |  0.998564 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 29.313
[INFO]    ::   - ZP Std: 0.033
[INFO]    ::   - ZP Range: [29.263,29.356]
[INFO]    ::   - ZP Mode: 29.312
[INFO]    ::   - ZP Median: 29.316
[INFO]    ::   - ZP 16th & 84th percentiles: 29.286, 29.339
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.667
[INFO]    ::   - ZP Std: 0.034
[INFO]    ::   - ZP Range: [31.631,31.716]
[INFO]    ::   - ZP Mode: 31.642
[INFO]    ::   - ZP Median: 31.661
[INFO]    ::   - ZP 16th & 84th percentiles: 31.636, 31.699
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=29.316 std=0.033 No. stars=4
[INFO]    :: ZP Ref=31.661 std=0.034 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.115
[INFO]    :: [V2] Blanked subtraction at 6798 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-12T10_38641bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61111.44723669998
[INFO]    :: BACKGROUND: -0.660 counts
[INFO]    :: SN FLUX 21335.740 counts
[INFO]    :: SN MAG 18.493 mag
[INFO]    :: SN MAG - BACKGROUND 18.493 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.273 1.617
[INFO]    :: xoff_arc=0.843 arcsec, yoff_arc=0.599 arcsec
[INFO]    :: [V2] Injection pool: 745443 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 119.489
[INFO]    :: S/N (std artifical sn)= 119.489
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.493 mag
[INFO]    :: Mag = 18.493+/-0.009 
[INFO]    :: 5-sig limit = 21.939
[INFO]    :: 3-sig limit = 22.493
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.046 arcsec (0.123 pixels) 
[INFO]    :: DEC error: 0.360 arcsec (0.971 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.034
[INFO]    :: Science zeropoint std: 0.033
--------------------------------------------------------------------------------
 Mag = 18.493+/-0.048 lim=22.493 MJD=61111.447
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  917.9
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-12T10_38641_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 2/38 | 36 remaining
[███░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 3/38 (7%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 3 OF 38 [3/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 3 OF 38 [3/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_06_38_29_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 124
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=121, XR=909, YL=0, YT=895 (interior std=24.9, thresholds: std>49.8, med>835.6) ~106 bright sources remain
Cutting out image with borders: XL=124, XR=899, YL=10, YT=889 (geometry: 124,899,10,889  noise: 121,909,0,895)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.16 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (810,350)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61112.277
[INFO]    :: Observation time: 2026-03-13T06:38:29.591224
[INFO]    :: Exposure time: 240.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=0.11", ΔDec=0.25" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 2 inliers within 3.0" of median (ΔRA=0.11", ΔDec=0.25") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 7 sources in sci
[INFO]    :: [V2] star_match_shift: 328 sources in ref
[INFO]    :: [V2] star_match_shift: 7 matches → dx=2.142±0.437 px, dy=1.187±1.176 px
[INFO]    :: [V2] [star_match] Applied fine shift (2.142, 1.187) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1935 clipped px (>= 9.95 nMgy), dilated by 12px -> 8058 masked (1.0% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 91.9% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.256 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32093.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32093.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32093.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_32093.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 37.5
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[INFO]    :: ePSF science fallback OK: FWHM=4.60 px  elong=1.14  n_stars=36  scatter=2.489 px
[INFO]    :: Using ePSF kernel for science convolution (PSFEx chi² too high)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-13T06_23909ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/117 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_32093.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.00669613
[INFO]    :: [V2] PSF centroid corrected: shift=(0.000, -0.268) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.83px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-13T06_23909sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -126.192242 -1.228302 277.208496
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 277
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 6
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-13T06comb_psf_32093.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[31.05464192 31.1171906  31.0685768  31.03522439]
[0.63126112 0.63896293 0.63337931 0.63666491]
[31.5838464  31.67066355 31.64159419 31.58561445]
[0.96975973 0.95669163 0.99946628 0.99927637]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  31.1172 |  31.6707 |  0.638963 |  0.956692 |
|  31.0352 |  31.5856 |  0.636665 |  0.999276 |
|  31.0686 |  31.6416 |  0.633379 |  0.999466 |
|  31.0546 |  31.5838 |  0.631261 |  0.96976  |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.069
[INFO]    ::   - ZP Std: 0.030
[INFO]    ::   - ZP Range: [31.035,31.117]
[INFO]    ::   - ZP Mode: 31.055
[INFO]    ::   - ZP Median: 31.062
[INFO]    ::   - ZP 16th & 84th percentiles: 31.045, 31.094
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.635
[INFO]    ::   - RSQ Std: 0.003
[INFO]    ::   - RSQ Range: [0.631,0.639]
[INFO]    ::   - RSQ Mode: 0.631
[INFO]    ::   - RSQ Median: 0.635
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.632, 0.638
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.620
[INFO]    ::   - ZP Std: 0.037
[INFO]    ::   - ZP Range: [31.584,31.671]
[INFO]    ::   - ZP Mode: 31.584
[INFO]    ::   - ZP Median: 31.614
[INFO]    ::   - ZP 16th & 84th percentiles: 31.585, 31.657
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.981
[INFO]    ::   - RSQ Std: 0.019
[INFO]    ::   - RSQ Range: [0.957,0.999]
[INFO]    ::   - RSQ Mode: 0.970
[INFO]    ::   - RSQ Median: 0.985
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.963, 0.999
[INFO]    :: ZP Sci=31.062 std=0.030 No. stars=4
[INFO]    :: ZP Ref=31.614 std=0.037 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.601
[INFO]    :: [V2] Blanked subtraction at 8058 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-13T06_23909bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61112.27673099982
[INFO]    :: BACKGROUND: -34.654 counts
[INFO]    :: SN FLUX 93270.252 counts
[INFO]    :: SN MAG 18.637 mag
[INFO]    :: SN MAG - BACKGROUND 18.638 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.887 1.113
[INFO]    :: xoff_arc=0.329 arcsec, yoff_arc=0.413 arcsec
[INFO]    :: [V2] Injection pool: 674073 valid positions in valid-data region
[INFO]    :: Background injection positions used: 435
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 36.047
[INFO]    :: S/N (std artifical sn)= 36.047
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.337 mag
[INFO]    :: Mag = 18.637+/-0.030 
[INFO]    :: 5-sig limit = 20.782
[INFO]    :: 3-sig limit = 21.337
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.091 arcsec (0.246 pixels) 
[INFO]    :: DEC error: 0.319 arcsec (0.860 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.037
[INFO]    :: Science zeropoint std: 0.030
--------------------------------------------------------------------------------
 Mag = 18.637+/-0.056 lim=21.337 MJD=61112.277
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  967.3
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-13T06_23909_photometry.txt
[INFO]    :: Total time: 22.0 seconds
✓ COMPLETED: 3/38 | 35 remaining
[████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 4/38 (10%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 4 OF 38 [4/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 4 OF 38 [4/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260313_11_58_16_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=19.1, thresholds: std>38.3, med>624.0) ~191 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.28 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (818,349)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61112.499
[INFO]    :: Observation time: 2026-03-13T11:58:16.519607
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.66", ΔDec=0.24" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-3.66", ΔDec=0.24") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 6 sources in sci
[INFO]    :: [V2] star_match_shift: 318 sources in ref
[INFO]    :: [V2] star_match_shift: 6 matches → dx=1.608±0.392 px, dy=0.890±0.408 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.608, 0.890) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1948 clipped px (>= 9.95 nMgy), dilated by 13px -> 8698 masked (1.1% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 91.0% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.130 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_56762.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_56762.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_56762.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_56762.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(0.025, 0.106) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (37, 37) (FWHM=6.15px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-13T11_43096ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 43/116 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_56762.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01844956
[INFO]    :: [V2] PSF centroid corrected: shift=(0.012, -0.262) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (37, 37) (FWHM=6.15px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-13T11_43096sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -130.874648 -0.927799 235.991328
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 235
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 8
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-13T11comb_psf_56762.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[29.36044513 29.41316751 29.37172751 29.33136065]
[0.99953301 0.9987742  0.9995447  0.99923047]
[31.65147169 31.73195675 31.68486455 31.63153522]
[0.99950368 0.99844073 0.99973015 0.99964286]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3717 |  31.6849 |  0.999545 |  0.99973  |
|  29.3604 |  31.6515 |  0.999533 |  0.999504 |
|  29.3314 |  31.6315 |  0.99923  |  0.999643 |
|  29.4132 |  31.732  |  0.998774 |  0.998441 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 29.369
[INFO]    ::   - ZP Std: 0.029
[INFO]    ::   - ZP Range: [29.331,29.413]
[INFO]    ::   - ZP Mode: 29.360
[INFO]    ::   - ZP Median: 29.366
[INFO]    ::   - ZP 16th & 84th percentiles: 29.345, 29.393
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.675
[INFO]    ::   - ZP Std: 0.038
[INFO]    ::   - ZP Range: [31.632,31.732]
[INFO]    ::   - ZP Mode: 31.651
[INFO]    ::   - ZP Median: 31.668
[INFO]    ::   - ZP 16th & 84th percentiles: 31.641, 31.709
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=29.366 std=0.029 No. stars=4
[INFO]    :: ZP Ref=31.668 std=0.038 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.120
[INFO]    :: [V2] Blanked subtraction at 8698 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-13T11_43096bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61112.49880270008
[INFO]    :: BACKGROUND: -1.720 counts
[INFO]    :: SN FLUX 18939.006 counts
[INFO]    :: SN MAG 18.673 mag
[INFO]    :: SN MAG - BACKGROUND 18.673 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.930 0.789
[INFO]    :: xoff_arc=0.345 arcsec, yoff_arc=0.292 arcsec
[INFO]    :: [V2] Injection pool: 724171 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 75.393
[INFO]    :: S/N (std artifical sn)= 75.393
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.173 mag
[INFO]    :: Mag = 18.673+/-0.014 
[INFO]    :: 5-sig limit = 21.619
[INFO]    :: 3-sig limit = 22.173
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.084 arcsec (0.228 pixels) 
[INFO]    :: DEC error: 0.425 arcsec (1.147 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.038
[INFO]    :: Science zeropoint std: 0.029
--------------------------------------------------------------------------------
 Mag = 18.673+/-0.050 lim=22.173 MJD=61112.499
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  967.3
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-13T11_43096_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 4/38 | 34 remaining
[█████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 5/38 (13%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 5 OF 38 [5/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 5 OF 38 [5/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_21_32_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=116, XR=909, YL=0, YT=895 (interior std=32.2, thresholds: std>64.4, med>1352.8) ~96 bright sources remain
Cutting out image with borders: XL=116, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 116,909,0,895)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.75 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (806,334)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61113.265
[INFO]    :: Observation time: 2026-03-14T06:21:32.925241
[INFO]    :: Exposure time: 340.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 30 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 30 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 30 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 6 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-10.03", ΔDec=-2.73" from 6 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-10.03", ΔDec=-2.73") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 7 sources in sci
[INFO]    :: [V2] star_match_shift: 308 sources in ref
[INFO]    :: [V2] star_match_shift: 7 matches → dx=2.029±0.388 px, dy=1.627±1.221 px
[INFO]    :: [V2] [star_match] Applied fine shift (2.029, 1.627) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1950 clipped px (>= 9.95 nMgy), dilated by 10px -> 6815 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.4% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.157 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_75982.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_75982.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_75982.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_75982.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 11.1
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[INFO]    :: ePSF science fallback OK: FWHM=4.22 px  elong=1.29  n_stars=40  scatter=3.188 px
[INFO]    :: Using ePSF kernel for science convolution (PSFEx chi² too high)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-14T06_22892ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/118 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_75982.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02530576
[INFO]    :: [V2] PSF centroid corrected: shift=(0.025, -0.251) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.71px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-14T06_22892sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -199.892801 -3.568025 445.396746
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 445
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 6
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-14T06comb_psf_75982.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[30.18385921 30.23464127 30.19303471 30.15107953]
[0.98326358 0.98107314 0.98376    0.9840044 ]
[31.60515655 31.68244519 31.64362026 31.59243378]
[0.99918894 0.99777805 0.9996251  0.99950913]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  30.1511 |  31.5924 |  0.984004 |  0.999509 |
|  30.193  |  31.6436 |  0.98376  |  0.999625 |
|  30.1839 |  31.6052 |  0.983264 |  0.999189 |
|  30.2346 |  31.6824 |  0.981073 |  0.997778 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 30.191
[INFO]    ::   - ZP Std: 0.030
[INFO]    ::   - ZP Range: [30.151,30.235]
[INFO]    ::   - ZP Mode: 30.184
[INFO]    ::   - ZP Median: 30.188
[INFO]    ::   - ZP 16th & 84th percentiles: 30.167, 30.215
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.983
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.981,0.984]
[INFO]    ::   - RSQ Mode: 0.983
[INFO]    ::   - RSQ Median: 0.984
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.982, 0.984
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.631
[INFO]    ::   - ZP Std: 0.035
[INFO]    ::   - ZP Range: [31.592,31.682]
[INFO]    ::   - ZP Mode: 31.605
[INFO]    ::   - ZP Median: 31.624
[INFO]    ::   - ZP 16th & 84th percentiles: 31.599, 31.664
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.998, 1.000
[INFO]    :: ZP Sci=30.188 std=0.030 No. stars=4
[INFO]    :: ZP Ref=31.624 std=0.035 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.266
[INFO]    :: [V2] Blanked subtraction at 6815 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-14T06_22892bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61113.26496409997
[INFO]    :: BACKGROUND: -11.942 counts
[INFO]    :: SN FLUX 31592.004 counts
[INFO]    :: SN MAG 18.940 mag
[INFO]    :: SN MAG - BACKGROUND 18.940 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.902 -0.512
[INFO]    :: xoff_arc=0.705 arcsec, yoff_arc=0.190 arcsec
[INFO]    :: [V2] Injection pool: 695456 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 54.495
[INFO]    :: S/N (std artifical sn)= 54.495
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.088 mag
[INFO]    :: Mag = 18.940+/-0.020 
[INFO]    :: 5-sig limit = 21.533
[INFO]    :: 3-sig limit = 22.088
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.076 arcsec (0.205 pixels) 
[INFO]    :: DEC error: 0.303 arcsec (0.817 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.035
[INFO]    :: Science zeropoint std: 0.030
--------------------------------------------------------------------------------
 Mag = 18.940+/-0.050 lim=22.088 MJD=61113.265
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  979.5
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-14T06_22892_photometry.txt
[INFO]    :: Total time: 15.0 seconds
✓ COMPLETED: 5/38 | 33 remaining
[██████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 6/38 (15%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 6 OF 38 [6/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 6 OF 38 [6/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_06_40_49_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 117
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=113, XR=909, YL=0, YT=888 (interior std=24.6, thresholds: std>49.2, med>850.1) ~54 bright sources remain
Cutting out image with borders: XL=117, XR=899, YL=10, YT=888 (geometry: 117,899,10,889  noise: 113,909,0,888)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.53 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (803,336)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61113.278
[INFO]    :: Observation time: 2026-03-14T06:40:49.785702
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  4  detections along bad columns.
[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  4  detections along bad columns.
[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  4  detections along bad columns.
[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 6 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-10.13", ΔDec=-3.00" from 6 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-10.13", ΔDec=-3.00") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 6 sources in sci
[INFO]    :: [V2] star_match_shift: 314 sources in ref
[INFO]    :: [V2] star_match_shift: 6 matches → dx=1.722±0.323 px, dy=0.819±0.494 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.722, 0.819) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1928 clipped px (>= 9.95 nMgy), dilated by 9px -> 6199 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.7% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.182 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_22480.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_22480.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_22480.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_22480.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 2.0
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.107, -0.127) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.14px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-14T06_24049ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 50/122 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_22480.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.00601697
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.030, -0.226) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.14px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-14T06_24049sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -123.379621 -1.093641 274.093421
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 274
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 4(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 5(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 5(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 5, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 5
[INFO]    :: Length of matched catalog= 5
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 5
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-14T06comb_psf_22480.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 5
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 3)
[28.91461993 28.94575077 28.91809615]
[0.99987995 0.9979362  0.99986341]
[31.77189848 31.84589189 31.80988865]
[0.99726044 0.99745233 0.99572206]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  28.9146 |  31.7719 |  0.99988  |  0.99726  |
|  28.9181 |  31.8099 |  0.999863 |  0.995722 |
|  28.9458 |  31.8459 |  0.997936 |  0.997452 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 3
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 3
[INFO]    :: Number of stars used for zeropoint calculation: 3
[INFO]    :: Science Zeropoints from (3 stars)
[INFO]    ::   - ZP Mean: 28.926
[INFO]    ::   - ZP Std: 0.014
[INFO]    ::   - ZP Range: [28.915,28.946]
[INFO]    ::   - ZP Mode: 28.915
[INFO]    ::   - ZP Median: 28.918
[INFO]    ::   - ZP 16th & 84th percentiles: 28.916, 28.937
[INFO]    ::   - # above threshold (0.35): 3
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (3 stars)
[INFO]    ::   - ZP Mean: 31.809
[INFO]    ::   - ZP Std: 0.030
[INFO]    ::   - ZP Range: [31.772,31.846]
[INFO]    ::   - ZP Mode: 31.772
[INFO]    ::   - ZP Median: 31.810
[INFO]    ::   - ZP 16th & 84th percentiles: 31.784, 31.834
[INFO]    ::   - # above threshold (0.35): 3
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.997
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.996,0.997]
[INFO]    ::   - RSQ Mode: 0.997
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.996, 0.997
[INFO]    :: ZP Sci=28.918 std=0.014 No. stars=3
[INFO]    :: ZP Ref=31.810 std=0.030 No. stars=3
[WARNING] :: [V2] QUALITY WARNING: only 3 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.070
[INFO]    :: [V2] Blanked subtraction at 6199 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-14T06_24049bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61113.27835390018
[INFO]    :: BACKGROUND: -2.165 counts
[INFO]    :: SN FLUX 10383.555 counts
[INFO]    :: SN MAG 18.877 mag
[INFO]    :: SN MAG - BACKGROUND 18.877 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.328 0.953
[INFO]    :: xoff_arc=0.863 arcsec, yoff_arc=0.353 arcsec
[INFO]    :: [V2] Injection pool: 755196 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 46.761
[INFO]    :: S/N (std artifical sn)= 46.761
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.859 mag
[INFO]    :: Mag = 18.877+/-0.023 
[INFO]    :: 5-sig limit = 21.305
[INFO]    :: 3-sig limit = 21.859
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.087 arcsec (0.235 pixels) 
[INFO]    :: DEC error: 0.577 arcsec (1.558 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.030
[INFO]    :: Science zeropoint std: 0.014
--------------------------------------------------------------------------------
 Mag = 18.877+/-0.040 lim=21.859 MJD=61113.278
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1191.9
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-14T06_24049_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 6/38 | 32 remaining
[███████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 7/38 (18%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 7 OF 38 [7/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 7 OF 38 [7/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260314_10_43_31_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=11.8, thresholds: std>23.6, med>266.5) ~135 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.74 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (817,346)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61113.447
[INFO]    :: Observation time: 2026-03-14T10:43:31.332753
[INFO]    :: Exposure time: 90.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 8 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=0.11", ΔDec=-0.37" from 8 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 2 inliers within 3.0" of median (ΔRA=0.11", ΔDec=-0.37") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 10 sources in sci
[INFO]    :: [V2] star_match_shift: 334 sources in ref
[INFO]    :: [V2] star_match_shift: 10 matches → dx=1.616±0.363 px, dy=0.949±1.362 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (9 inliers) → dx=1.544±0.221 px, dy=0.974±0.928 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.544, 0.974) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1953 clipped px (>= 9.95 nMgy), dilated by 10px -> 6810 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 91.2% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.080 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_20156.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_20156.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_20156.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_20156.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.091, -0.165) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.69px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-14T10_38611ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 49/118 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_20156.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02672941
[INFO]    :: [V2] PSF centroid corrected: shift=(0.036, -0.280) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.69px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-14T10_38611sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -52.428753 -0.33833 96.195224
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 96
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 7
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-14T10comb_psf_20156.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[28.5302824  28.5749898  28.53927356 28.49105984]
[0.99964829 0.99749901 0.99987943 0.99964143]
[31.62926877 31.70181906 31.67132558 31.62385825]
[0.9993097  0.99809699 0.99963625 0.99961713]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  28.5393 |  31.6713 |  0.999879 |  0.999636 |
|  28.5303 |  31.6293 |  0.999648 |  0.99931  |
|  28.4911 |  31.6239 |  0.999641 |  0.999617 |
|  28.575  |  31.7018 |  0.997499 |  0.998097 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 28.534
[INFO]    ::   - ZP Std: 0.030
[INFO]    ::   - ZP Range: [28.491,28.575]
[INFO]    ::   - ZP Mode: 28.530
[INFO]    ::   - ZP Median: 28.535
[INFO]    ::   - ZP 16th & 84th percentiles: 28.510, 28.558
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.997,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.657
[INFO]    ::   - ZP Std: 0.032
[INFO]    ::   - ZP Range: [31.624,31.702]
[INFO]    ::   - ZP Mode: 31.629
[INFO]    ::   - ZP Median: 31.650
[INFO]    ::   - ZP 16th & 84th percentiles: 31.626, 31.687
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=28.535 std=0.030 No. stars=4
[INFO]    :: ZP Ref=31.650 std=0.032 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.057
[INFO]    :: [V2] Blanked subtraction at 6810 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-14T10_38611bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61113.44688999979
[INFO]    :: BACKGROUND: -1.511 counts
[INFO]    :: SN FLUX 7170.941 counts
[INFO]    :: SN MAG 18.896 mag
[INFO]    :: SN MAG - BACKGROUND 18.896 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.164 -0.383
[INFO]    :: xoff_arc=0.431 arcsec, yoff_arc=0.142 arcsec
[INFO]    :: [V2] Injection pool: 745401 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 68.287
[INFO]    :: S/N (std artifical sn)= 68.287
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.289 mag
[INFO]    :: Mag = 18.896+/-0.016 
[INFO]    :: 5-sig limit = 21.734
[INFO]    :: 3-sig limit = 22.289
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.100 arcsec (0.271 pixels) 
[INFO]    :: DEC error: 0.180 arcsec (0.486 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.032
[INFO]    :: Science zeropoint std: 0.030
--------------------------------------------------------------------------------
 Mag = 18.896+/-0.046 lim=22.289 MJD=61113.447
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1214.9
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-14T10_38611_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 7/38 | 31 remaining
[████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 8/38 (21%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 8 OF 38 [8/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 8 OF 38 [8/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 150
Right X border at 950
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=136, XR=909, YL=0, YT=886 (interior std=20.0, thresholds: std>40.1, med>584.5) ~212 bright sources remain
Cutting out image with borders: XL=150, XR=899, YL=10, YT=850 (geometry: 150,899,10,850  noise: 136,909,0,886)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 3.10 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (833,339)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61114.253
[INFO]    :: Observation time: 2026-03-15T06:04:28.616002
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-15T06_21868bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[WARNING] :: Unhandled exception on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260315_06_04_28_f_a_b_ZTF26aakjzdt_g_g: ERROR 8 in wcss2p() at line 3966 of file cextern/wcslib/C/wcs.c:
One or more of the pixel coordinates were invalid.
ERROR 6 in linx2p() at line 979 of file cextern/wcslib/C/lin.c:
De-distort error.
ERROR 5 in disx2p() at line 1411 of file cextern/wcslib/C/dis.c:
Convergence not achieved after 30 iterations, residual 1.3e-10.

[█████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 9/38 (23%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 9 OF 38 [9/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 9 OF 38 [9/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260316_05_39_16_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 121
Right X border at 950
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=108, XR=909, YL=0, YT=899 (interior std=27.5, thresholds: std>54.9, med>1054.1) ~119 bright sources remain
Cutting out image with borders: XL=121, XR=899, YL=10, YT=850 (geometry: 121,899,10,850  noise: 108,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.17 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (798,343)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61115.236
[INFO]    :: Observation time: 2026-03-16T05:39:16.481552
[INFO]    :: Exposure time: 313.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 6 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-9.94", ΔDec=-2.97" from 6 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-9.94", ΔDec=-2.97") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 6 sources in sci
[INFO]    :: [V2] star_match_shift: 303 sources in ref
[INFO]    :: [V2] star_match_shift: 6 matches → dx=1.279±1.141 px, dy=0.881±0.428 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.279, 0.881) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1944 clipped px (>= 9.95 nMgy), dilated by 12px -> 8097 masked (1.0% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 93.4% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.261 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38305.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38305.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38305.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_38305.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 43.6
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[INFO]    :: ePSF science fallback OK: FWHM=5.02 px  elong=1.03  n_stars=40  scatter=2.352 px
[INFO]    :: Using ePSF kernel for science convolution (PSFEx chi² too high)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-16T05_20356ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/119 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_38305.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01374824
[INFO]    :: [V2] PSF centroid corrected: shift=(0.005, -0.267) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (37, 37) (FWHM=5.85px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-16T05_20356sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -196.364252 -2.577457 381.243478
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 381
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 6
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-16T05comb_psf_38305.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[30.87627197 30.91167468 30.87958408 30.82429861]
[0.75253806 0.73705167 0.74497946 0.73992251]
[31.62120664 31.69858555 31.65908763 31.59778894]
[0.99479352 0.9932722  0.99970751 0.99946533]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  30.8763 |  31.6212 |  0.752538 |  0.994794 |
|  30.8796 |  31.6591 |  0.744979 |  0.999708 |
|  30.8243 |  31.5978 |  0.739923 |  0.999465 |
|  30.9117 |  31.6986 |  0.737052 |  0.993272 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 30.873
[INFO]    ::   - ZP Std: 0.031
[INFO]    ::   - ZP Range: [30.824,30.912]
[INFO]    ::   - ZP Mode: 30.876
[INFO]    ::   - ZP Median: 30.878
[INFO]    ::   - ZP 16th & 84th percentiles: 30.849, 30.896
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.744
[INFO]    ::   - RSQ Std: 0.006
[INFO]    ::   - RSQ Range: [0.737,0.753]
[INFO]    ::   - RSQ Mode: 0.753
[INFO]    ::   - RSQ Median: 0.742
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.738, 0.749
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.644
[INFO]    ::   - ZP Std: 0.038
[INFO]    ::   - ZP Range: [31.598,31.699]
[INFO]    ::   - ZP Mode: 31.621
[INFO]    ::   - ZP Median: 31.640
[INFO]    ::   - ZP 16th & 84th percentiles: 31.609, 31.680
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.997
[INFO]    ::   - RSQ Std: 0.003
[INFO]    ::   - RSQ Range: [0.993,1.000]
[INFO]    ::   - RSQ Mode: 0.995
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.994, 1.000
[INFO]    :: ZP Sci=30.878 std=0.031 No. stars=4
[INFO]    :: ZP Ref=31.640 std=0.038 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.496
[INFO]    :: [V2] Blanked subtraction at 8097 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-16T05_20356bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61115.23560710019
[INFO]    :: BACKGROUND: -21.540 counts
[INFO]    :: SN FLUX 45473.878 counts
[INFO]    :: SN MAG 19.234 mag
[INFO]    :: SN MAG - BACKGROUND 19.234 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.113 -0.340
[INFO]    :: xoff_arc=0.783 arcsec, yoff_arc=0.126 arcsec
[INFO]    :: [V2] Injection pool: 675503 valid positions in valid-data region
[INFO]    :: Background injection positions used: 435
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 21.929
[INFO]    :: S/N (std artifical sn)= 21.929
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.393 mag
[INFO]    :: Mag = 19.234+/-0.048 
[INFO]    :: 5-sig limit = 20.839
[INFO]    :: 3-sig limit = 21.393
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.231 arcsec (0.622 pixels) 
[INFO]    :: DEC error: 0.162 arcsec (0.436 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.038
[INFO]    :: Science zeropoint std: 0.031
--------------------------------------------------------------------------------
 Mag = 19.234+/-0.069 lim=21.393 MJD=61115.236
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1319.7
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-16T05_20356_photometry.txt
[INFO]    :: Total time: 9.0 seconds
✓ COMPLETED: 9/38 | 29 remaining
[██████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 10/38 (26%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 10 OF 38 [10/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 10 OF 38 [10/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g.fits
[WARNING] :: Error converting SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_08_49_17_f_b_ZTF26aakjzdt_g_g.fits, error: "Keyword 'PC1_1' not found."
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 150
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=26.8, thresholds: std>53.6, med>1241.1) ~104 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=850 (geometry: 50,899,150,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 4.64 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (436,139)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61116.368
[INFO]    :: Observation time: 2026-03-17T08:49:17.062280
[INFO]    :: Exposure time: 230.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=1.46", ΔDec=22.90" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 1 inliers within 3.0" of median (ΔRA=1.46", ΔDec=22.90") — likely non-stellar science detections; falling back to star_match
[WARNING] :: [V2] star_match_shift: only 1 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=351.700 px, dy=445.000 px
[WARNING] :: [V2] Fine shift (351.700, 445.000) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 132 clipped px (>= 9.95 nMgy), dilated by 26px -> 10601 masked (1.3% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.4% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=-0.000 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-17T08_31757bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_43321.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_43321.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_43321.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_43321.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.8
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.004, -0.054) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (77, 77) (FWHM=12.52px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-17T08_31757ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 63/139 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_43321.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.0515591
[INFO]    :: [V2] PSF centroid corrected: shift=(0.063, -0.212) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (77, 77) (FWHM=12.52px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-17T08_31757sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -316.529725 -4.060247 503.77553
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 503
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 0, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 1
[INFO]    :: Length of matched catalog= 0
[WARNING] :: Few stars (0) matched!
[WARNING] :: ZTF26aakjzdt 61116.36755870003 g: Less than 2 matched calibration stars 
[WARNING] :: Exiting: not enough stars to calibrate!
[WARNING] :: Pipeline stopped after: Reference catalog
[███████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 11/38 (28%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 11 OF 38 [11/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 11 OF 38 [11/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260317_10_44_03_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 150
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=22.3, thresholds: std>44.7, med>873.0) ~127 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=850 (geometry: 50,899,150,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.69 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (824,363)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61116.447
[INFO]    :: Observation time: 2026-03-17T10:44:03.379216
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 11 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 11 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 11 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.45", ΔDec=3.13" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-3.45", ΔDec=3.13") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 6 sources in sci
[INFO]    :: [V2] star_match_shift: 318 sources in ref
[INFO]    :: [V2] star_match_shift: 6 matches → dx=0.128±3.164 px, dy=0.234±9.653 px
[WARNING] :: [V2] star_match_shift: scatter too large (MAD_x=3.16, MAD_y=9.65 px > 3.0px) — suppressing shift, relying on WCS alignment
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=41.800 px, dy=5.100 px
[INFO]    :: [V2] [pcc] Applied fine shift (41.800, 5.100) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1886 clipped px (>= 9.95 nMgy), dilated by 15px -> 10072 masked (1.2% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.9% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.209 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-17T10_38643bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_51551.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_51551.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_51551.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_51551.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.8
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.071, 0.121) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.26px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-17T10_38643ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/111 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_51551.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.03865433
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.032, -0.248) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.26px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-17T10_38643sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -218.649045 -2.053852 350.85058
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 350
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 2, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 3
[INFO]    :: Length of matched catalog= 2
[WARNING] :: Few stars (2) matched!
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 2
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-17T10comb_psf_51551.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 2
[WARNING] :: [V2] Saturation filter would reject all 2 matched calibrators — keeping the original list.  ZP_ref scatter will likely be inflated by Legacy-Survey clipping (consider using PS1 reference for this field).
[WARNING] :: [V2] chi2_shift requested (+19.3,-0.4) px — clipped to (+5.0,-0.4) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+19.3,-1.9) px — clipped to (+5.0,-1.9) px (likely latching on a saturated-star residual)
[28.90848925 28.84708491]
[0.99963522 0.99936325]
[nan nan]
[0.02845543 0.02826907]
[WARNING] :: Few stars in the field, continuing if 1 or 2 stars in science
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  28.9085 |      nan |  0.999635 | 0.0284554 |
|  28.8471 |      nan |  0.999363 | 0.0282691 |
+----------+----------+-----------+-----------+
[WARNING] :: No stars with rsq>0.35. Lowering threshold to 0.30
[WARNING] :: No stars with rsq>0.30. Lowering threshold to 0.25
[WARNING] :: No stars with rsq>0.25. Lowering threshold to 0.20
[WARNING] :: No stars with rsq>0.20. Lowering threshold to 0.15
[WARNING] :: No stars with rsq>0.15. Lowering threshold to 0.10
[WARNING] :: No stars with rsq>0.10. Lowering threshold to 0.05
[WARNING] :: No stars with rsq>0.05. Lowering threshold to 0.00
[WARNING] :: [V2] PSF-fit quality poor for all stars — using best available at floor threshold 0.05
[WARNING] :: All PSF-fit zeropoints failed for reference image (ZTF26aakjzdt g)
[WARNING] :: [V2] Attempting aperture-photometry zeropoint fallback...
[WARNING] :: [V2] Aperture ZP fallback also found 0 valid stars — exiting
[WARNING] :: Pipeline stopped after: Zeropoint calibration
[████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 12/38 (31%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 12 OF 38 [12/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 12 OF 38 [12/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_05_56_58_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=136, XR=909, YL=0, YT=899 (interior std=24.9, thresholds: std>49.8, med>854.2) ~170 bright sources remain
Cutting out image with borders: XL=136, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 136,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.69 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (839,339)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61117.248
[INFO]    :: Observation time: 2026-03-18T05:56:58.390594
[INFO]    :: Exposure time: 234.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 9 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.85", ΔDec=3.19" from 9 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-3.85", ΔDec=3.19") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 8 sources in sci
[INFO]    :: [V2] star_match_shift: 352 sources in ref
[INFO]    :: [V2] star_match_shift: 8 matches → dx=0.608±0.395 px, dy=0.342±2.228 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (6 inliers) → dx=0.537±0.280 px, dy=0.342±1.422 px
[INFO]    :: [V2] [star_match] Applied fine shift (0.537, 0.342) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1977 clipped px (>= 9.95 nMgy), dilated by 15px -> 10129 masked (1.2% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 88.7% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.141 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06922.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06922.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_06922.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_06922.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.1
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.025, -0.004) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.27px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-18T05_21418ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/121 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_06922.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01739473
[INFO]    :: [V2] PSF centroid corrected: shift=(0.015, -0.273) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.27px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-18T05_21418sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -167.384369 -1.690137 316.318225
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 316
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 2, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 4
[INFO]    :: Length of matched catalog= 2
[WARNING] :: Few stars (2) matched!
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 2
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-18T05comb_psf_06922.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 2
[WARNING] :: [V2] Saturation filter would reject all 2 matched calibrators — keeping the original list.  ZP_ref scatter will likely be inflated by Legacy-Survey clipping (consider using PS1 reference for this field).
[29.48122322 29.44188188]
[0.99937078 0.99974128]
[32.4814811  32.41052807]
[0.98377595 0.9919625 ]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.4419 |  32.4105 |  0.999741 |  0.991963 |
|  29.4812 |  32.4815 |  0.999371 |  0.983776 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 2
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 2
[INFO]    :: Number of stars used for zeropoint calculation: 2
[INFO]    :: Science Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 29.462
[INFO]    ::   - ZP Std: 0.020
[INFO]    ::   - ZP Range: [29.442,29.481]
[INFO]    ::   - ZP Mode: 29.481
[INFO]    ::   - ZP Median: 29.462
[INFO]    ::   - ZP 16th & 84th percentiles: 29.448, 29.475
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 1.000
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 32.446
[INFO]    ::   - ZP Std: 0.035
[INFO]    ::   - ZP Range: [32.411,32.481]
[INFO]    ::   - ZP Mode: 32.481
[INFO]    ::   - ZP Median: 32.446
[INFO]    ::   - ZP 16th & 84th percentiles: 32.422, 32.470
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.988
[INFO]    ::   - RSQ Std: 0.004
[INFO]    ::   - RSQ Range: [0.984,0.992]
[INFO]    ::   - RSQ Mode: 0.984
[INFO]    ::   - RSQ Median: 0.988
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.985, 0.991
[INFO]    :: ZP Sci=29.462 std=0.020 No. stars=2
[INFO]    :: ZP Ref=32.446 std=0.035 No. stars=2
[WARNING] :: [V2] QUALITY WARNING: only 2 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.064
[INFO]    :: [V2] Blanked subtraction at 10129 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-18T05_21418bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[WARNING] :: [V2] chi2_shift requested (-1.9,-7.4) px — clipped to (-1.9,-5.0) px (likely latching on a saturated-star residual)
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61117.24789769994
[INFO]    :: BACKGROUND: 3.767 counts
[INFO]    :: SN FLUX 10882.230 counts
[INFO]    :: SN MAG 19.370 mag
[INFO]    :: SN MAG - BACKGROUND 19.370 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.914 -5.000
[INFO]    :: xoff_arc=0.709 arcsec, yoff_arc=1.853 arcsec
[INFO]    :: [V2] Injection pool: 703148 valid positions in valid-data region
[WARNING] :: [V2] chi2_shift requested (-1.9,-7.4) px — clipped to (-1.9,-5.0) px (likely latching on a saturated-star residual)
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 30.628
[INFO]    :: S/N (std artifical sn)= 30.628
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.892 mag
[INFO]    :: Mag = 19.370+/-0.035 
[INFO]    :: 5-sig limit = 21.338
[INFO]    :: 3-sig limit = 21.892
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.125 arcsec (0.338 pixels) 
[INFO]    :: DEC error: 0.702 arcsec (1.894 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.035
[INFO]    :: Science zeropoint std: 0.020
--------------------------------------------------------------------------------
 Mag = 19.370+/-0.054 lim=21.892 MJD=61117.248
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-18T05_21418_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 12/38 | 26 remaining
[█████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░] 13/38 (34%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 13 OF 38 [13/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 13 OF 38 [13/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260318_12_40_50_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=18.4, thresholds: std>36.8, med>603.9) ~188 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.89 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (791,331)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61117.528
[INFO]    :: Observation time: 2026-03-18T12:40:50.002721
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 10 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=0.48", ΔDec=-0.59" from 10 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 2 inliers within 3.0" of median (ΔRA=0.48", ΔDec=-0.59") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 10 sources in sci
[INFO]    :: [V2] star_match_shift: 292 sources in ref
[INFO]    :: [V2] star_match_shift: 10 matches → dx=2.311±0.690 px, dy=0.797±1.342 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (9 inliers) → dx=2.471±0.475 px, dy=0.627±1.418 px
[INFO]    :: [V2] [star_match] Applied fine shift (2.471, 0.627) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1940 clipped px (>= 9.95 nMgy), dilated by 11px -> 7418 masked (0.9% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 94.0% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.114 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71850.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71850.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71850.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_71850.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.1
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.047, -0.153) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.10px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-18T12_45650ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/128 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_71850.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02858417
[INFO]    :: [V2] PSF centroid corrected: shift=(0.006, -0.258) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.10px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-18T12_45650sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.245316 0.136704 5.141331
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 5
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 14, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 154
[INFO]    :: Length of matched catalog= 14
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 14
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-18T12comb_psf_71850.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 14
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 13)
[29.2741749  29.18612176 29.27640357 29.31609777 29.28490475 29.38907262
 29.33044266 29.25180847 29.39818943 29.24612881 29.28753554 29.18985435]
[0.99682244 0.80442244 0.99949214 0.99743155 0.99977504 0.93140035
 0.88556271 0.99947139 0.84085414 0.9969828  0.96509695 0.97406193]
[31.68319894 31.63638938 31.65664252 31.73152538 31.69092254 31.70466469
 31.65469279 31.64187538 31.63335881 31.63644783 31.68821518 31.66242992]
[0.99906678 0.97793099 0.99948829 0.99838134 0.99967731 0.99293126
 0.99493221 0.99946028 0.97861163 0.99961734 0.996648   0.99886804]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.2849 |  31.6909 |  0.999775 |  0.999677 |
|  29.2764 |  31.6566 |  0.999492 |  0.999488 |
|  29.2518 |  31.6419 |  0.999471 |  0.99946  |
|  29.3161 |  31.7315 |  0.997432 |  0.998381 |
|  29.2461 |  31.6364 |  0.996983 |  0.999617 |
|  29.2742 |  31.6832 |  0.996822 |  0.999067 |
|  29.1899 |  31.6624 |  0.974062 |  0.998868 |
|  29.2875 |  31.6882 |  0.965097 |  0.996648 |
|  29.3891 |  31.7047 |  0.9314   |  0.992931 |
|  29.3304 |  31.6547 |  0.885563 |  0.994932 |
|  29.3982 |  31.6334 |  0.840854 |  0.978612 |
|  29.1861 |  31.6364 |  0.804422 |  0.977931 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 12
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 10
[INFO]    :: Number of stars used for zeropoint calculation: 10
[INFO]    :: Science Zeropoints from (10 stars)
[INFO]    ::   - ZP Mean: 29.285
[INFO]    ::   - ZP Std: 0.051
[INFO]    ::   - ZP Range: [29.190,29.389]
[INFO]    ::   - ZP Mode: 29.274
[INFO]    ::   - ZP Median: 29.281
[INFO]    ::   - ZP 16th & 84th percentiles: 29.249, 29.324
[INFO]    ::   - # above threshold (0.35): 12
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.949
[INFO]    ::   - RSQ Std: 0.066
[INFO]    ::   - RSQ Range: [0.804,1.000]
[INFO]    ::   - RSQ Mode: 0.997
[INFO]    ::   - RSQ Median: 0.985
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.875, 0.999
[INFO]    :: Reference Zeropoints from (10 stars)
[INFO]    ::   - ZP Mean: 31.666
[INFO]    ::   - ZP Std: 0.023
[INFO]    ::   - ZP Range: [31.636,31.705]
[INFO]    ::   - ZP Mode: 31.683
[INFO]    ::   - ZP Median: 31.660
[INFO]    ::   - ZP 16th & 84th percentiles: 31.639, 31.690
[INFO]    ::   - # above threshold (0.35): 12
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.995
[INFO]    ::   - RSQ Std: 0.008
[INFO]    ::   - RSQ Range: [0.978,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.989, 1.000
[INFO]    :: ZP Sci=29.281 std=0.051 No. stars=10
[INFO]    :: ZP Ref=31.660 std=0.023 No. stars=10
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.112
[INFO]    :: [V2] Blanked subtraction at 7418 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-18T12_45650bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61117.52835649997
[INFO]    :: BACKGROUND: -2.599 counts
[INFO]    :: SN FLUX 10159.516 counts
[INFO]    :: SN MAG 19.263 mag
[INFO]    :: SN MAG - BACKGROUND 19.264 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.410 0.176
[INFO]    :: xoff_arc=0.523 arcsec, yoff_arc=0.065 arcsec
[INFO]    :: [V2] Injection pool: 739949 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 51.814
[INFO]    :: S/N (std artifical sn)= 51.814
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.357 mag
[INFO]    :: Mag = 19.263+/-0.021 
[INFO]    :: 5-sig limit = 21.802
[INFO]    :: 3-sig limit = 22.357
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.139 arcsec (0.375 pixels) 
[INFO]    :: DEC error: 0.252 arcsec (0.681 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.023
[INFO]    :: Science zeropoint std: 0.051
--------------------------------------------------------------------------------
 Mag = 19.263+/-0.060 lim=22.357 MJD=61117.528
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-18T12_45650_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 13/38 | 25 remaining
[██████████████░░░░░░░░░░░░░░░░░░░░░░░░░░] 14/38 (36%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 14 OF 38 [14/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 14 OF 38 [14/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260319_06_20_47_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=136, XR=909, YL=0, YT=892 (interior std=21.1, thresholds: std>42.2, med>623.2) ~251 bright sources remain
Cutting out image with borders: XL=136, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 136,909,0,892)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.15 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (844,354)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61118.264
[INFO]    :: Observation time: 2026-03-19T06:20:47.575357
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 16 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 16 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 16 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 8 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-4.10", ΔDec=1.53" from 8 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-4.10", ΔDec=1.53") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 9 sources in sci
[INFO]    :: [V2] star_match_shift: 368 sources in ref
[INFO]    :: [V2] star_match_shift: 9 matches → dx=0.965±0.410 px, dy=0.310±1.219 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (6 inliers) → dx=0.964±0.154 px, dy=0.612±0.704 px
[INFO]    :: [V2] [star_match] Applied fine shift (0.964, 0.612) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1941 clipped px (>= 9.95 nMgy), dilated by 12px -> 8061 masked (1.0% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 87.7% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.130 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_09538.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_09538.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_09538.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_09538.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 5.0
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[WARNING] :: ePSF science fallback rejected (unphysical): FWHM=44.21 px (limit 17.4)  elong=3.43 (limit 1.5) — using PSFEx kernel instead
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.165, -0.122) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.79px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-19T06_22847ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 45/113 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_09538.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.00937842
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.005, -0.247) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.79px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-19T06_22847sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -119.852579 -0.858986 229.013516
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 229
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 8
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-19T06comb_psf_09538.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[29.25156877 29.30090008 29.25664762 29.21222805]
[0.99978417 0.99843887 0.9998888  0.99966857]
[31.64780955 31.72365466 31.68411678 31.62745105]
[0.99947121 0.99843594 0.99971388 0.99965783]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.2566 |  31.6841 |  0.999889 |  0.999714 |
|  29.2516 |  31.6478 |  0.999784 |  0.999471 |
|  29.2122 |  31.6275 |  0.999669 |  0.999658 |
|  29.3009 |  31.7237 |  0.998439 |  0.998436 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 29.255
[INFO]    ::   - ZP Std: 0.031
[INFO]    ::   - ZP Range: [29.212,29.301]
[INFO]    ::   - ZP Mode: 29.252
[INFO]    ::   - ZP Median: 29.254
[INFO]    ::   - ZP 16th & 84th percentiles: 29.231, 29.280
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.671
[INFO]    ::   - ZP Std: 0.037
[INFO]    ::   - ZP Range: [31.627,31.724]
[INFO]    ::   - ZP Mode: 31.648
[INFO]    ::   - ZP Median: 31.666
[INFO]    ::   - ZP 16th & 84th percentiles: 31.637, 31.705
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.998,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 1.000
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=29.254 std=0.031 No. stars=4
[INFO]    :: ZP Ref=31.666 std=0.037 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.108
[INFO]    :: [V2] Blanked subtraction at 8061 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-19T06_22847bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61118.26443939982
[INFO]    :: BACKGROUND: 1.569 counts
[INFO]    :: SN FLUX 8808.002 counts
[INFO]    :: SN MAG 19.392 mag
[INFO]    :: SN MAG - BACKGROUND 19.392 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=1.465 -2.199
[INFO]    :: xoff_arc=0.543 arcsec, yoff_arc=0.815 arcsec
[INFO]    :: [V2] Injection pool: 729822 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 38.789
[INFO]    :: S/N (std artifical sn)= 38.789
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.171 mag
[INFO]    :: Mag = 19.392+/-0.028 
[INFO]    :: 5-sig limit = 21.616
[INFO]    :: 3-sig limit = 22.171
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.396 arcsec (1.068 pixels) 
[INFO]    :: DEC error: 0.305 arcsec (0.824 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.037
[INFO]    :: Science zeropoint std: 0.031
--------------------------------------------------------------------------------
 Mag = 19.392+/-0.056 lim=22.171 MJD=61118.264
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-19T06_22847_photometry.txt
[INFO]    :: Total time: 11.0 seconds
✓ COMPLETED: 14/38 | 24 remaining
[███████████████░░░░░░░░░░░░░░░░░░░░░░░░░] 15/38 (39%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 15 OF 38 [15/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 15 OF 38 [15/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260320_06_37_16_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=17.7, thresholds: std>35.4, med>531.3) ~186 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.75 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (807,350)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61119.276
[INFO]    :: Observation time: 2026-03-20T06:37:16.872727
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 22 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 24 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-0.26", ΔDec=-0.20" from 24 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 4 inliers within 3.0" of median (ΔRA=-0.26", ΔDec=-0.20") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 30 sources in sci
[INFO]    :: [V2] star_match_shift: 309 sources in ref
[INFO]    :: [V2] star_match_shift: 27 matches → dx=1.715±0.471 px, dy=1.279±2.137 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (24 inliers) → dx=1.875±0.523 px, dy=1.217±1.934 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.875, 1.217) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1955 clipped px (>= 9.95 nMgy), dilated by 10px -> 6853 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.3% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.112 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32634.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32634.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32634.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_32634.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.5
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.098, -0.059) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.72px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-20T06_23836ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/116 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_32634.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01413926
[INFO]    :: [V2] PSF centroid corrected: shift=(0.029, -0.255) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.72px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-20T06_23836sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.313383 0.217949 2.872727
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 2
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 18, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 162
[INFO]    :: Length of matched catalog= 18
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 18
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-20T06comb_psf_32634.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 18
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 17)
[29.13472083 29.49844129 29.23414336 29.22666021 29.24446479 29.45307186
 29.289578   29.2560156  29.24345513 29.24729128 29.11417719 29.20380835
 29.15608341 29.21899309 29.19953672 29.21844207 29.23657059]
[0.78190052 0.81518209 0.99897894 0.88807659 0.99988803 0.73721424
 0.99787247 0.99994072 0.88355988 0.95987213 0.91823503 0.99986136
 0.81168904 0.99816897 0.99350129 0.98816272 0.98867693]
[31.54780963 31.72094996 31.6655561  31.61571174 31.63448155 31.70159228
 31.70701875 31.67346093 31.67668463 31.63084828 31.63407514 31.62694978
 31.62112473 31.62254741 31.59766816 31.67481829 31.64601693]
[0.96056984 0.98634699 0.99894285 0.99201419 0.99938233 0.97488753
 0.99821605 0.9996268  0.99351992 0.99457459 0.97729648 0.99961106
 0.98434878 0.99960609 0.99698868 0.99727973 0.99890626]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.256  |  31.6735 |  0.999941 |  0.999627 |
|  29.2445 |  31.6345 |  0.999888 |  0.999382 |
|  29.2038 |  31.6269 |  0.999861 |  0.999611 |
|  29.2341 |  31.6656 |  0.998979 |  0.998943 |
|  29.219  |  31.6225 |  0.998169 |  0.999606 |
|  29.2896 |  31.707  |  0.997872 |  0.998216 |
|  29.1995 |  31.5977 |  0.993501 |  0.996989 |
|  29.2366 |  31.646  |  0.988677 |  0.998906 |
|  29.2184 |  31.6748 |  0.988163 |  0.99728  |
|  29.2473 |  31.6308 |  0.959872 |  0.994575 |
|  29.1142 |  31.6341 |  0.918235 |  0.977296 |
|  29.2267 |  31.6157 |  0.888077 |  0.992014 |
|  29.2435 |  31.6767 |  0.88356  |  0.99352  |
|  29.4984 |  31.7209 |  0.815182 |  0.986347 |
|  29.1561 |  31.6211 |  0.811689 |  0.984349 |
|  29.1347 |  31.5478 |  0.781901 |  0.96057  |
|  29.4531 |  31.7016 |  0.737214 |  0.974888 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 17
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 15
[INFO]    :: Number of stars used for zeropoint calculation: 15
[INFO]    :: Science Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 29.238
[INFO]    ::   - ZP Std: 0.068
[INFO]    ::   - ZP Range: [29.135,29.453]
[INFO]    ::   - ZP Mode: 29.135
[INFO]    ::   - ZP Median: 29.234
[INFO]    ::   - ZP 16th & 84th percentiles: 29.201, 29.254
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.927
[INFO]    ::   - RSQ Std: 0.088
[INFO]    ::   - RSQ Range: [0.737,1.000]
[INFO]    ::   - RSQ Mode: 0.782
[INFO]    ::   - RSQ Median: 0.988
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.814, 0.999
[INFO]    :: Reference Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 31.649
[INFO]    ::   - ZP Std: 0.032
[INFO]    ::   - ZP Range: [31.598,31.707]
[INFO]    ::   - ZP Mode: 31.666
[INFO]    ::   - ZP Median: 31.634
[INFO]    ::   - ZP 16th & 84th percentiles: 31.621, 31.676
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.991
[INFO]    ::   - RSQ Std: 0.011
[INFO]    ::   - RSQ Range: [0.961,1.000]
[INFO]    ::   - RSQ Mode: 0.961
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.981, 0.999
[INFO]    :: ZP Sci=29.234 std=0.068 No. stars=15
[INFO]    :: ZP Ref=31.634 std=0.032 No. stars=15
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.110
[INFO]    :: [V2] Blanked subtraction at 6853 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-20T06_23836bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61119.27588969981
[INFO]    :: BACKGROUND: -0.540 counts
[INFO]    :: SN FLUX 7364.206 counts
[INFO]    :: SN MAG 19.566 mag
[INFO]    :: SN MAG - BACKGROUND 19.566 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.945 1.273
[INFO]    :: xoff_arc=0.350 arcsec, yoff_arc=0.472 arcsec
[INFO]    :: [V2] Injection pool: 745328 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 41.683
[INFO]    :: S/N (std artifical sn)= 41.683
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.423 mag
[INFO]    :: Mag = 19.566+/-0.026 
[INFO]    :: 5-sig limit = 21.869
[INFO]    :: 3-sig limit = 22.423
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.125 arcsec (0.336 pixels) 
[INFO]    :: DEC error: 0.492 arcsec (1.328 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.032
[INFO]    :: Science zeropoint std: 0.068
--------------------------------------------------------------------------------
 Mag = 19.566+/-0.080 lim=22.423 MJD=61119.276
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-20T06_23836_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 15/38 | 23 remaining
[████████████████░░░░░░░░░░░░░░░░░░░░░░░░] 16/38 (42%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 16 OF 38 [16/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 16 OF 38 [16/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_06_17_58_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=17.9, thresholds: std>35.8, med>554.3) ~173 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.52 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (803,363)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61120.262
[INFO]    :: Observation time: 2026-03-21T06:17:58.133297
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 25 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 25 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 25 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31791144999983,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 29 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-0.13", ΔDec=-0.34" from 29 nearest-neighbour pairs
[INFO]    :: [V2] ps1_catalog_shift (pass B): 5 inliers → dx=-0.145±1.205 px, dy=0.141±1.697 px
[INFO]    :: [V2] [ps1_catalog] Applied fine shift (-0.145, 0.141) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1938 clipped px (>= 9.95 nMgy), dilated by 9px -> 6199 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.8% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.128 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_39171.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_39171.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_39171.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_39171.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.4
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.100, -0.083) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.09px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-21T06_22678ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 43/112 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_39171.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.04559759
[INFO]    :: [V2] PSF centroid corrected: shift=(0.039, -0.248) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.09px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-21T06_22678sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.493236 0.286437 3.040858
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 3
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 20, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 163
[INFO]    :: Length of matched catalog= 20
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 20
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-21T06comb_psf_39171.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 20
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 18)
[29.09613086 29.1319274  29.28342453 29.24211347 29.12473896 29.24524071
 29.30071743 29.28283526 29.2485497  29.37232905 29.22493955 29.18104363
 29.19151963 29.28920109 29.2002338  29.69616404 29.16360484 29.23106033]
[0.70404318 0.36633365 0.70518905 0.99915536 0.97392768 0.99991952
 0.59356622 0.99816286 0.99993009 0.90903966 0.92419443 0.95627281
 0.99987775 0.90090124 0.99867613 0.29634708 0.98127451 0.99268337]
[31.53778579 31.57661389 31.7013099  31.6414803  31.59899527 31.61277719
 31.67961706 31.68160485 31.65353765 31.65099985 31.61642525 31.61576532
 31.60976951 31.62679013 31.60045299 31.58246378 31.65626279 31.63012993]
[0.97972937 0.96091858 0.98593398 0.99844835 0.99595182 0.99812828
 0.98085015 0.99627995 0.99928938 0.99447857 0.9954449  0.98431054
 0.99844222 0.9900944  0.9981703  0.99464968 0.99576351 0.99859531]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.2485 |  31.6535 |  0.99993  |  0.999289 |
|  29.2452 |  31.6128 |  0.99992  |  0.998128 |
|  29.1915 |  31.6098 |  0.999878 |  0.998442 |
|  29.2421 |  31.6415 |  0.999155 |  0.998448 |
|  29.2002 |  31.6005 |  0.998676 |  0.99817  |
|  29.2828 |  31.6816 |  0.998163 |  0.99628  |
|  29.2311 |  31.6301 |  0.992683 |  0.998595 |
|  29.1636 |  31.6563 |  0.981275 |  0.995764 |
|  29.1247 |  31.599  |  0.973928 |  0.995952 |
|  29.181  |  31.6158 |  0.956273 |  0.984311 |
|  29.2249 |  31.6164 |  0.924194 |  0.995445 |
|  29.3723 |  31.651  |  0.90904  |  0.994479 |
|  29.2892 |  31.6268 |  0.900901 |  0.990094 |
|  29.2834 |  31.7013 |  0.705189 |  0.985934 |
|  29.0961 |  31.5378 |  0.704043 |  0.979729 |
|  29.3007 |  31.6796 |  0.593566 |  0.98085  |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 17
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 15
[INFO]    :: Number of stars used for zeropoint calculation: 15
[INFO]    :: Science Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 29.223
[INFO]    ::   - ZP Std: 0.054
[INFO]    ::   - ZP Range: [29.125,29.301]
[INFO]    ::   - ZP Mode: 29.132
[INFO]    ::   - ZP Median: 29.231
[INFO]    ::   - ZP 16th & 84th percentiles: 29.168, 29.283
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.883
[INFO]    ::   - RSQ Std: 0.177
[INFO]    ::   - RSQ Range: [0.366,1.000]
[INFO]    ::   - RSQ Mode: 0.704
[INFO]    ::   - RSQ Median: 0.974
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.705, 0.999
[INFO]    :: Reference Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 31.630
[INFO]    ::   - ZP Std: 0.029
[INFO]    ::   - ZP Range: [31.577,31.682]
[INFO]    ::   - ZP Mode: 31.577
[INFO]    ::   - ZP Median: 31.627
[INFO]    ::   - ZP 16th & 84th percentiles: 31.603, 31.656
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.991
[INFO]    ::   - RSQ Std: 0.010
[INFO]    ::   - RSQ Range: [0.961,0.999]
[INFO]    ::   - RSQ Mode: 0.980
[INFO]    ::   - RSQ Median: 0.996
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.983, 0.998
[INFO]    :: ZP Sci=29.231 std=0.054 No. stars=15
[INFO]    :: ZP Ref=31.627 std=0.029 No. stars=15
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.110
[INFO]    :: [V2] Blanked subtraction at 6199 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-21T06_22678bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61120.26247850014
[INFO]    :: BACKGROUND: 2.264 counts
[INFO]    :: SN FLUX 6665.470 counts
[INFO]    :: SN MAG 19.671 mag
[INFO]    :: SN MAG - BACKGROUND 19.672 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.234 0.422
[INFO]    :: xoff_arc=0.828 arcsec, yoff_arc=0.156 arcsec
[INFO]    :: [V2] Injection pool: 755196 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 37.104
[INFO]    :: S/N (std artifical sn)= 37.104
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.402 mag
[INFO]    :: Mag = 19.671+/-0.029 
[INFO]    :: 5-sig limit = 21.848
[INFO]    :: 3-sig limit = 22.402
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.773 arcsec (2.086 pixels) 
[INFO]    :: DEC error: 0.316 arcsec (0.851 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.029
[INFO]    :: Science zeropoint std: 0.054
--------------------------------------------------------------------------------
 Mag = 19.671+/-0.068 lim=22.402 MJD=61120.262
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-21T06_22678_photometry.txt
[INFO]    :: Total time: 6.0 seconds
✓ COMPLETED: 16/38 | 22 remaining
[█████████████████░░░░░░░░░░░░░░░░░░░░░░░] 17/38 (44%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 17 OF 38 [17/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 17 OF 38 [17/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260321_10_39_30_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=16.3, thresholds: std>32.6, med>462.9) ~202 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.64 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (805,358)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61120.444
[INFO]    :: Observation time: 2026-03-21T10:39:30.519116
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 2 matches found between science image 1 and reference image 0
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 24 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 2 matches found between science image 1 and reference image 0
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3178244999998,71.84164445000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 24 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-1.01", ΔDec=-4.93" from 24 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-1.01", ΔDec=-4.93") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 33 sources in sci
[INFO]    :: [V2] star_match_shift: 317 sources in ref
[INFO]    :: [V2] star_match_shift: 29 matches → dx=-0.228±1.505 px, dy=0.732±1.702 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (25 inliers) → dx=-0.256±1.349 px, dy=0.663±1.589 px
[INFO]    :: [V2] [star_match] Applied fine shift (-0.256, 0.663) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 2015 clipped px (>= 9.95 nMgy), dilated by 9px -> 6331 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.9% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.100 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_86211.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_86211.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_86211.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_86211.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.1
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.093, -0.107) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (27, 27) (FWHM=4.42px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-21T10_38370ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/117 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_86211.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01712813
[INFO]    :: [V2] PSF centroid corrected: shift=(0.036, -0.253) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (27, 27) (FWHM=4.42px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-21T10_38370sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.548857 0.36505 2.673061
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 2
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 7(<=7) stars found in reference catalog, increasing search radius: 1->2
[INFO]    :: Catalog stars in PS1/SDSS found = 20, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 173
[INFO]    :: Length of matched catalog= 20
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 20
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-21T10comb_psf_86211.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 20
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 18)
[WARNING] :: [V2] chi2_shift requested (-3.5,+5.1) px — clipped to (-3.5,+5.0) px (likely latching on a saturated-star residual)
[29.1880476  29.27792463 29.294526   29.30028759 29.21439205 29.29028682
 29.57978265 29.34195607 29.30481401 29.43280891 29.27566421 29.24914753
 29.26301904 29.35002394 29.2950686  29.54091013 29.29312039 29.22361006]
[0.8301271  0.7194972  0.99891671 0.99629602 0.95956587 0.99982496
 0.76829207 0.99804267 0.99992858 0.92737317 0.95751243 0.95822972
 0.99970607 0.90487817 0.99907653 0.51071684 0.98497853 0.99141736]
[31.53237383 31.70686577 31.66822233 31.61474823 31.62842131 31.64312163
 31.68076268 31.71436713 31.67705707 31.67867964 31.63843081 31.63526741
 31.6306735  31.63501195 31.62127326 31.59961532 31.67219836 31.64482063]
[0.92110931 0.97060275 0.99730708 0.9960164  0.99350329 0.999229
 0.95174829 0.9981117  0.99920233 0.99373468 0.99494955 0.98249082
 0.99949958 0.98511709 0.99964914 0.99714414 0.99758564 0.99878893]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3048 |  31.6771 |  0.999929 |  0.999202 |
|  29.2903 |  31.6431 |  0.999825 |  0.999229 |
|  29.263  |  31.6307 |  0.999706 |  0.9995   |
|  29.2951 |  31.6213 |  0.999077 |  0.999649 |
|  29.2945 |  31.6682 |  0.998917 |  0.997307 |
|  29.342  |  31.7144 |  0.998043 |  0.998112 |
|  29.3003 |  31.6147 |  0.996296 |  0.996016 |
|  29.2236 |  31.6448 |  0.991417 |  0.998789 |
|  29.2931 |  31.6722 |  0.984979 |  0.997586 |
|  29.2144 |  31.6284 |  0.959566 |  0.993503 |
|  29.2491 |  31.6353 |  0.95823  |  0.982491 |
|  29.2757 |  31.6384 |  0.957512 |  0.99495  |
|  29.4328 |  31.6787 |  0.927373 |  0.993735 |
|  29.35   |  31.635  |  0.904878 |  0.985117 |
|  29.188  |  31.5324 |  0.830127 |  0.921109 |
|  29.5798 |  31.6808 |  0.768292 |  0.951748 |
|  29.2779 |  31.7069 |  0.719497 |  0.970603 |
|  29.5409 |  31.5996 |  0.510717 |  0.997144 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 18
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 16
[INFO]    :: Number of stars used for zeropoint calculation: 16
[INFO]    :: Science Zeropoints from (16 stars)
[INFO]    ::   - ZP Mean: 29.309
[INFO]    ::   - ZP Std: 0.078
[INFO]    ::   - ZP Range: [29.214,29.541]
[INFO]    ::   - ZP Mode: 29.278
[INFO]    ::   - ZP Median: 29.294
[INFO]    ::   - ZP 16th & 84th percentiles: 29.255, 29.347
[INFO]    ::   - # above threshold (0.35): 18
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.917
[INFO]    ::   - RSQ Std: 0.128
[INFO]    ::   - RSQ Range: [0.511,1.000]
[INFO]    ::   - RSQ Mode: 0.830
[INFO]    ::   - RSQ Median: 0.972
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.813, 0.999
[INFO]    :: Reference Zeropoints from (16 stars)
[INFO]    ::   - ZP Mean: 31.648
[INFO]    ::   - ZP Std: 0.028
[INFO]    ::   - ZP Range: [31.600,31.707]
[INFO]    ::   - ZP Mode: 31.707
[INFO]    ::   - ZP Median: 31.641
[INFO]    ::   - ZP 16th & 84th percentiles: 31.624, 31.678
[INFO]    ::   - # above threshold (0.35): 18
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.988
[INFO]    ::   - RSQ Std: 0.020
[INFO]    ::   - RSQ Range: [0.921,1.000]
[INFO]    ::   - RSQ Mode: 0.921
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.979, 0.999
[INFO]    :: ZP Sci=29.294 std=0.078 No. stars=16
[INFO]    :: ZP Ref=31.641 std=0.028 No. stars=16
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.115
[INFO]    :: [V2] Blanked subtraction at 6331 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-21T10_38370bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61120.44410309987
[INFO]    :: BACKGROUND: 0.451 counts
[INFO]    :: SN FLUX 6998.702 counts
[INFO]    :: SN MAG 19.681 mag
[INFO]    :: SN MAG - BACKGROUND 19.681 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-3.746 2.254
[INFO]    :: xoff_arc=1.388 arcsec, yoff_arc=0.835 arcsec
[INFO]    :: [V2] Injection pool: 750576 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 44.690
[INFO]    :: S/N (std artifical sn)= 44.690
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.614 mag
[INFO]    :: Mag = 19.681+/-0.024 
[INFO]    :: 5-sig limit = 22.059
[INFO]    :: 3-sig limit = 22.614
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.394 arcsec (1.062 pixels) 
[INFO]    :: DEC error: 0.347 arcsec (0.935 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.028
[INFO]    :: Science zeropoint std: 0.078
--------------------------------------------------------------------------------
 Mag = 19.681+/-0.086 lim=22.614 MJD=61120.444
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-21T10_38370_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 17/38 | 21 remaining
[██████████████████░░░░░░░░░░░░░░░░░░░░░░] 18/38 (47%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 18 OF 38 [18/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 18 OF 38 [18/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 65
Right X border at 909
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=21.4, thresholds: std>42.8, med>826.9) ~34 bright sources remain
Cutting out image with borders: XL=65, XR=899, YL=150, YT=889 (geometry: 65,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.42 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (404,369)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61121.232
[INFO]    :: Observation time: 2026-03-22T05:33:43.382055
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 14 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 14 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 14 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3178244999998,71.84164445000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[WARNING] :: [V2] ps1_catalog_shift: only 2 sources detected (need ≥5) — skipping
[WARNING] :: [V2] star_match_shift: only 1 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-13.200 px, dy=-16.800 px
[INFO]    :: [V2] [pcc] Applied fine shift (-13.200, -16.800) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 58 clipped px (>= 9.95 nMgy), dilated by 14px -> 2337 masked (0.3% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 100.0% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=0.036 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-22T05_20023bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_62291.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_62291.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_62291.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_62291.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.8
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.179, 0.289) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (41, 41) (FWHM=6.53px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-22T05_20023ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 51/121 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_62291.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 0.99299062
[INFO]    :: [V2] PSF centroid corrected: shift=(0.151, -0.276) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (41, 41) (FWHM=6.53px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-22T05_20023sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -187.130045 -1.475684 318.266481
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 318
[WARNING] :: Unhandled exception on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_33_43_f_a_b_ZTF26aakjzdt_g_g: 'NoneType' object is not subscriptable
[████████████████████░░░░░░░░░░░░░░░░░░░░] 19/38 (50%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 19 OF 38 [19/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 19 OF 38 [19/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260322_05_45_54_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 0
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=898 (interior std=33.1, thresholds: std>66.2, med>1465.9) ~39 bright sources remain
Cutting out image with borders: XL=10, XR=899, YL=10, YT=889 (geometry: 10,899,10,889  noise: 0,909,0,898)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.36 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (406,367)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61121.24
[INFO]    :: Observation time: 2026-03-22T05:45:54.541138
[INFO]    :: Exposure time: 320.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3178244999998,71.84164445000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 26 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-8.24", ΔDec=-0.39" from 26 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-8.24", ΔDec=-0.39") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 25 sources in sci
[INFO]    :: [V2] star_match_shift: 251 sources in ref
[INFO]    :: [V2] star_match_shift: 20 matches → dx=0.465±2.908 px, dy=-0.143±1.552 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (16 inliers) → dx=-0.090±2.454 px, dy=-0.251±1.085 px
[INFO]    :: [V2] [star_match] Applied fine shift (-0.090, -0.251) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 63 clipped px (>= 9.95 nMgy), dilated by 13px -> 2099 masked (0.3% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 100.0% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.243 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_68189.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_68189.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_68189.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_68189.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.0
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.097, -0.140) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (39, 39) (FWHM=6.36px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-22T05_20754ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 55/131 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_68189.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01207868
[INFO]    :: [V2] PSF centroid corrected: shift=(0.149, -0.284) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (39, 39) (FWHM=6.36px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-22T05_20754sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -2.855512 -2.899093 4.670519
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 4
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 8(<=7) stars found in reference catalog, increasing search radius: 1->2
[INFO]    :: Catalog stars in PS1/SDSS found = 14, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 130
[INFO]    :: Length of matched catalog= 14
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 14
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-22T05comb_psf_68189.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 14
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 13)
[29.40240212 29.53929077 29.54723448 29.53102395 29.3282851  29.39798429
 29.41552037 29.45080204 29.43222498 29.50040942 28.9232849  29.58445104]
[0.98501058 0.88046267 0.88540855 0.97252834 0.87350049 0.97366502
 0.97849297 0.96856832 0.90340954 0.99902519 0.36155342 0.89378633]
[31.64202383 31.69551415 31.66378165 31.70176009 31.59118179 31.65448927
 31.56763905 31.59116712 31.62093611 31.63636173 31.5167634  31.55895786]
[0.99900178 0.99420116 0.99830667 0.97972099 0.9977455  0.99869622
 0.99934798 0.99921464 0.99446594 0.99935613 0.96070076 0.99498498]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.5004 |  31.6364 |  0.999025 |  0.999356 |
|  29.4024 |  31.642  |  0.985011 |  0.999002 |
|  29.4155 |  31.5676 |  0.978493 |  0.999348 |
|  29.398  |  31.6545 |  0.973665 |  0.998696 |
|  29.531  |  31.7018 |  0.972528 |  0.979721 |
|  29.4508 |  31.5912 |  0.968568 |  0.999215 |
|  29.4322 |  31.6209 |  0.90341  |  0.994466 |
|  29.5845 |  31.559  |  0.893786 |  0.994985 |
|  29.5472 |  31.6638 |  0.885409 |  0.998307 |
|  29.5393 |  31.6955 |  0.880463 |  0.994201 |
|  29.3283 |  31.5912 |  0.8735   |  0.997745 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 12
[INFO]    :: [V2] Median-zp filter: dropped 1 outlier(s) >|0.30| mag from median (sci_med=29.442, ref_med=31.629; remaining: 11)
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 9
[INFO]    :: Number of stars used for zeropoint calculation: 9
[INFO]    :: Science Zeropoints from (9 stars)
[INFO]    ::   - ZP Mean: 29.469
[INFO]    ::   - ZP Std: 0.058
[INFO]    ::   - ZP Range: [29.398,29.547]
[INFO]    ::   - ZP Mode: 29.402
[INFO]    ::   - ZP Median: 29.451
[INFO]    ::   - ZP 16th & 84th percentiles: 29.406, 29.537
[INFO]    ::   - # above threshold (0.35): 11
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.938
[INFO]    ::   - RSQ Std: 0.047
[INFO]    ::   - RSQ Range: [0.874,0.999]
[INFO]    ::   - RSQ Mode: 0.985
[INFO]    ::   - RSQ Median: 0.969
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.883, 0.981
[INFO]    :: Reference Zeropoints from (9 stars)
[INFO]    ::   - ZP Mean: 31.629
[INFO]    ::   - ZP Std: 0.038
[INFO]    ::   - ZP Range: [31.568,31.696]
[INFO]    ::   - ZP Mode: 31.642
[INFO]    ::   - ZP Median: 31.636
[INFO]    ::   - ZP 16th & 84th percentiles: 31.591, 31.661
[INFO]    ::   - # above threshold (0.35): 11
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.996
[INFO]    ::   - RSQ Std: 0.005
[INFO]    ::   - RSQ Range: [0.980,0.999]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.998
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.994, 0.999
[INFO]    :: ZP Sci=29.451 std=0.058 No. stars=9
[INFO]    :: ZP Ref=31.636 std=0.038 No. stars=9
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.134
[INFO]    :: [V2] Blanked subtraction at 2099 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-22T05_20754bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[WARNING] :: [V2] chi2_shift requested (-8.1,+3.6) px — clipped to (-5.0,+3.6) px (likely latching on a saturated-star residual)
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61121.24021490011
[INFO]    :: BACKGROUND: -6.755 counts
[INFO]    :: SN FLUX 7138.440 counts
[INFO]    :: SN MAG 19.817 mag
[INFO]    :: SN MAG - BACKGROUND 19.818 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-5.000 3.629
[INFO]    :: xoff_arc=1.853 arcsec, yoff_arc=1.345 arcsec
[INFO]    :: [V2] Injection pool: 736556 valid positions in valid-data region
[WARNING] :: [V2] chi2_shift requested (-8.1,+3.6) px — clipped to (-5.0,+3.6) px (likely latching on a saturated-star residual)
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 14.543
[INFO]    :: S/N (std artifical sn)= 14.543
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.531 mag
[INFO]    :: Mag = 19.817+/-0.072 
[INFO]    :: 5-sig limit = 20.976
[INFO]    :: 3-sig limit = 21.531
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.584 arcsec (1.575 pixels) 
[INFO]    :: DEC error: 0.361 arcsec (0.975 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.038
[INFO]    :: Science zeropoint std: 0.058
--------------------------------------------------------------------------------
 Mag = 19.817+/-0.100 lim=21.531 MJD=61121.240
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-22T05_20754_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 19/38 | 19 remaining
[█████████████████████░░░░░░░░░░░░░░░░░░░] 20/38 (52%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 20 OF 38 [20/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 20 OF 38 [20/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_05_16_18_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 0
Right X border at 950
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=19.0, thresholds: std>38.0, med>642.0) ~34 bright sources remain
Cutting out image with borders: XL=10, XR=899, YL=10, YT=889 (geometry: 10,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.17 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (338,382)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61122.22
[INFO]    :: Observation time: 2026-03-23T05:16:18.345327
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 1 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (1found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 2 matches found between science image 1 and reference image 0
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 2 matches found between science image 1 and reference image 0
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.31747689999978,71.84179205000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 27 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=0.37", ΔDec=-5.48" from 27 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=0.37", ΔDec=-5.48") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 27 sources in sci
[INFO]    :: [V2] star_match_shift: 262 sources in ref
[INFO]    :: [V2] star_match_shift: 21 matches → dx=1.055±1.719 px, dy=0.200±0.844 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (20 inliers) → dx=0.972±1.561 px, dy=0.132±0.763 px
[INFO]    :: [V2] [star_match] Applied fine shift (0.972, 0.132) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 6 clipped px (>= 9.95 nMgy), dilated by 12px -> 727 masked (0.1% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 100.0% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.143 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_07846.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_07846.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_07846.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_07846.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.030, -0.121) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (37, 37) (FWHM=5.87px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-23T05_18978ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 57/131 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_07846.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01517317
[INFO]    :: [V2] PSF centroid corrected: shift=(0.185, -0.275) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (37, 37) (FWHM=5.87px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-23T05_18978sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.178544 0.155957 2.245182
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 2
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 13, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 260
[INFO]    :: Length of matched catalog= 13
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 13
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-23T05comb_psf_07846.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 13
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 12)
[29.25098693 29.30961491 29.27698965 29.32908522 29.24570623 29.244233
 29.18524871 29.23330989 29.22336001 29.2903631  29.32697189 29.3715972 ]
[0.99590544 0.96795802 0.97261163 0.97379672 0.95373393 0.98296362
 0.99393351 0.98871638 0.96914355 0.99956535 0.8671274  0.93771241]
[31.67326678 31.72489675 31.69355384 31.73338006 31.61485062 31.68219409
 31.59662026 31.63114686 31.65736298 31.69069744 31.56226154 31.59911429]
[0.99914467 0.99566143 0.99846554 0.97620329 0.99811713 0.99896
 0.99950706 0.99927251 0.99494395 0.99963676 0.9851262  0.99549745]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.2904 |  31.6907 |  0.999565 |  0.999637 |
|  29.251  |  31.6733 |  0.995905 |  0.999145 |
|  29.1852 |  31.5966 |  0.993934 |  0.999507 |
|  29.2333 |  31.6311 |  0.988716 |  0.999273 |
|  29.2442 |  31.6822 |  0.982964 |  0.99896  |
|  29.3291 |  31.7334 |  0.973797 |  0.976203 |
|  29.277  |  31.6936 |  0.972612 |  0.998466 |
|  29.2234 |  31.6574 |  0.969144 |  0.994944 |
|  29.3096 |  31.7249 |  0.967958 |  0.995661 |
|  29.2457 |  31.6149 |  0.953734 |  0.998117 |
|  29.3716 |  31.5991 |  0.937712 |  0.995497 |
|  29.327  |  31.5623 |  0.867127 |  0.985126 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 12
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 10
[INFO]    :: Number of stars used for zeropoint calculation: 10
[INFO]    :: Science Zeropoints from (10 stars)
[INFO]    ::   - ZP Mean: 29.273
[INFO]    ::   - ZP Std: 0.037
[INFO]    ::   - ZP Range: [29.223,29.329]
[INFO]    ::   - ZP Mode: 29.251
[INFO]    ::   - ZP Median: 29.264
[INFO]    ::   - ZP 16th & 84th percentiles: 29.238, 29.319
[INFO]    ::   - # above threshold (0.35): 12
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.967
[INFO]    ::   - RSQ Std: 0.035
[INFO]    ::   - RSQ Range: [0.867,1.000]
[INFO]    ::   - RSQ Mode: 0.996
[INFO]    ::   - RSQ Median: 0.973
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.950, 0.994
[INFO]    :: Reference Zeropoints from (10 stars)
[INFO]    ::   - ZP Mean: 31.656
[INFO]    ::   - ZP Std: 0.042
[INFO]    ::   - ZP Range: [31.597,31.725]
[INFO]    ::   - ZP Mode: 31.673
[INFO]    ::   - ZP Median: 31.665
[INFO]    ::   - ZP 16th & 84th percentiles: 31.606, 31.692
[INFO]    ::   - # above threshold (0.35): 12
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.995
[INFO]    ::   - RSQ Std: 0.007
[INFO]    ::   - RSQ Range: [0.976,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.998
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.993, 0.999
[INFO]    :: ZP Sci=29.264 std=0.037 No. stars=10
[INFO]    :: ZP Ref=31.665 std=0.042 No. stars=10
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.110
[INFO]    :: [V2] Blanked subtraction at 727 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-23T05_18978bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61122.21965650003
[INFO]    :: BACKGROUND: -1.942 counts
[INFO]    :: SN FLUX 6625.121 counts
[INFO]    :: SN MAG 19.711 mag
[INFO]    :: SN MAG - BACKGROUND 19.711 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=1.352 -0.789
[INFO]    :: xoff_arc=0.501 arcsec, yoff_arc=0.292 arcsec
[INFO]    :: [V2] Injection pool: 745708 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 28.513
[INFO]    :: S/N (std artifical sn)= 28.513
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.156 mag
[INFO]    :: Mag = 19.711+/-0.037 
[INFO]    :: 5-sig limit = 21.601
[INFO]    :: 3-sig limit = 22.156
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.511 arcsec (1.379 pixels) 
[INFO]    :: DEC error: 0.212 arcsec (0.571 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.042
[INFO]    :: Science zeropoint std: 0.037
--------------------------------------------------------------------------------
 Mag = 19.711+/-0.067 lim=22.156 MJD=61122.220
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-23T05_18978_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 20/38 | 18 remaining
[██████████████████████░░░░░░░░░░░░░░░░░░] 21/38 (55%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 21 OF 38 [21/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 21 OF 38 [21/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260323_09_49_04_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 0
Right X border at 950
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=14.8, thresholds: std>29.5, med>388.8) ~45 bright sources remain
Cutting out image with borders: XL=10, XR=899, YL=150, YT=889 (geometry: 10,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.86 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (336,370)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61122.409
[INFO]    :: Observation time: 2026-03-23T09:49:04.759557
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 26 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 2 matches found between science image 1 and reference image 0
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 26 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[INFO]    :: 3 matches found between science image 1 and reference image 0
[INFO]    :: X-Y nudge complete
[INFO]    :: Nudging RA by -1.940000001354747e-05 arcsec and Dec by -1.8699999998261774e-05 arcsec
[INFO]    :: Removing temporary files
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 6 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=5.80", ΔDec=-10.37" from 6 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=5.80", ΔDec=-10.37") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 6 sources in sci
[INFO]    :: [V2] star_match_shift: 286 sources in ref
[INFO]    :: [V2] star_match_shift: 6 matches → dx=0.238±0.628 px, dy=1.929±0.871 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (5 inliers) → dx=0.006±0.451 px, dy=2.008±0.818 px
[INFO]    :: [V2] [star_match] Applied fine shift (0.006, 2.008) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: no clipped pixels (>= 9.95 nMgy) found
[INFO]    :: [V2] SEDM reproject alignment complete: 100.0% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.267 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71150.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71150.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_71150.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_71150.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.1
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.726, 0.042) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.01px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-23T09_35344ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 56/130 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_71150.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01735881
[INFO]    :: [V2] PSF centroid corrected: shift=(0.173, -0.264) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.01px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-23T09_35344sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -66.043642 -0.156169 131.683901
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 131
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 1(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 2, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 2
[INFO]    :: Length of matched catalog= 2
[WARNING] :: Few stars (2) matched!
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 2
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-23T09comb_psf_71150.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 2
[29.36063798 29.37712178]
[0.99865251 0.99886524]
[31.66678436 31.59279459]
[0.9997605  0.58120972]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3771 |  31.5928 |  0.998865 |  0.58121  |
|  29.3606 |  31.6668 |  0.998653 |  0.999761 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 2
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 2
[INFO]    :: Number of stars used for zeropoint calculation: 2
[INFO]    :: Science Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 29.369
[INFO]    ::   - ZP Std: 0.008
[INFO]    ::   - ZP Range: [29.361,29.377]
[INFO]    ::   - ZP Mode: 29.361
[INFO]    ::   - ZP Median: 29.369
[INFO]    ::   - ZP 16th & 84th percentiles: 29.363, 29.374
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,0.999]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 0.999
[INFO]    :: Reference Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 31.630
[INFO]    ::   - ZP Std: 0.037
[INFO]    ::   - ZP Range: [31.593,31.667]
[INFO]    ::   - ZP Mode: 31.667
[INFO]    ::   - ZP Median: 31.630
[INFO]    ::   - ZP 16th & 84th percentiles: 31.605, 31.655
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.790
[INFO]    ::   - RSQ Std: 0.209
[INFO]    ::   - RSQ Range: [0.581,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 0.790
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.648, 0.933
[INFO]    :: ZP Sci=29.369 std=0.008 No. stars=2
[INFO]    :: ZP Ref=31.630 std=0.037 No. stars=2
[WARNING] :: [V2] QUALITY WARNING: reference RSQ std=0.209 > 0.10 — PSF matching quality is inconsistent across field stars. Photometry for this epoch may have systematic errors > 0.2 mag.
[WARNING] :: [V2] QUALITY WARNING: only 2 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.125
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-23T09_35344bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61122.40908259992
[INFO]    :: BACKGROUND: -1.349 counts
[INFO]    :: SN FLUX 5971.962 counts
[INFO]    :: SN MAG 19.929 mag
[INFO]    :: SN MAG - BACKGROUND 19.929 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.293 2.473
[INFO]    :: xoff_arc=0.850 arcsec, yoff_arc=0.916 arcsec
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 43.188
[INFO]    :: S/N (std artifical sn)= 43.188
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.824 mag
[INFO]    :: Mag = 19.929+/-0.025 
[INFO]    :: 5-sig limit = 22.270
[INFO]    :: 3-sig limit = 22.824
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.109 arcsec (0.294 pixels) 
[INFO]    :: DEC error: 0.648 arcsec (1.749 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.037
[INFO]    :: Science zeropoint std: 0.008
--------------------------------------------------------------------------------
 Mag = 19.929+/-0.045 lim=22.824 MJD=61122.409
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-23T09_35344_photometry.txt
[INFO]    :: Total time: 5.0 seconds
✓ COMPLETED: 21/38 | 17 remaining
[███████████████████████░░░░░░░░░░░░░░░░░] 22/38 (57%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 22 OF 38 [22/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 22 OF 38 [22/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260324_11_04_31_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=17.6, thresholds: std>35.3, med>555.0) ~137 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.99 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (768,391)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61123.461
[INFO]    :: Observation time: 2026-03-24T11:04:31.272828
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.73", ΔDec=-0.00" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-3.73", ΔDec=-0.00") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 7 sources in sci
[INFO]    :: [V2] star_match_shift: 255 sources in ref
[INFO]    :: [V2] star_match_shift: 7 matches → dx=1.489±0.797 px, dy=0.428±0.427 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (6 inliers) → dx=1.558±0.502 px, dy=0.352±0.326 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.558, 0.352) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1973 clipped px (>= 9.95 nMgy), dilated by 11px -> 7494 masked (0.9% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 97.0% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.107 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_90665.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_90665.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_90665.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_90665.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.234, -0.049) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (33, 33) (FWHM=5.36px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-24T11_39871ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/122 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_90665.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.00538998
[INFO]    :: [V2] PSF centroid corrected: shift=(0.002, -0.179) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (33, 33) (FWHM=5.36px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-24T11_39871sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -115.674059 -0.760934 207.559353
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 207
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 6(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 6, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 8
[INFO]    :: Length of matched catalog= 6
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 6
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-24T11comb_psf_90665.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 6
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 4)
[28.9155756  28.98224094 28.9511099  28.91835322]
[0.99905962 0.9972638  0.9996723  0.9995293 ]
[31.66493809 31.7430839  31.70388223 31.65281203]
[0.99951823 0.99854665 0.99962112 0.99923169]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  28.9511 |  31.7039 |  0.999672 |  0.999621 |
|  28.9184 |  31.6528 |  0.999529 |  0.999232 |
|  28.9156 |  31.6649 |  0.99906  |  0.999518 |
|  28.9822 |  31.7431 |  0.997264 |  0.998547 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 4
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 28.942
[INFO]    ::   - ZP Std: 0.027
[INFO]    ::   - ZP Range: [28.916,28.982]
[INFO]    ::   - ZP Mode: 28.916
[INFO]    ::   - ZP Median: 28.935
[INFO]    ::   - ZP 16th & 84th percentiles: 28.917, 28.967
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.001
[INFO]    ::   - RSQ Range: [0.997,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.998, 1.000
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.691
[INFO]    ::   - ZP Std: 0.035
[INFO]    ::   - ZP Range: [31.653,31.743]
[INFO]    ::   - ZP Mode: 31.665
[INFO]    ::   - ZP Median: 31.684
[INFO]    ::   - ZP 16th & 84th percentiles: 31.659, 31.724
[INFO]    ::   - # above threshold (0.35): 4
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,1.000]
[INFO]    ::   - RSQ Mode: 1.000
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 1.000
[INFO]    :: ZP Sci=28.935 std=0.027 No. stars=4
[INFO]    :: ZP Ref=31.684 std=0.035 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.079
[INFO]    :: [V2] Blanked subtraction at 7494 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-24T11_39871bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61123.46147269988
[INFO]    :: BACKGROUND: -1.186 counts
[INFO]    :: SN FLUX 3510.480 counts
[INFO]    :: SN MAG 20.071 mag
[INFO]    :: SN MAG - BACKGROUND 20.072 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.812 -0.688
[INFO]    :: xoff_arc=0.301 arcsec, yoff_arc=0.255 arcsec
[INFO]    :: [V2] Injection pool: 736689 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 21.407
[INFO]    :: S/N (std artifical sn)= 21.407
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.205 mag
[INFO]    :: Mag = 20.071+/-0.050 
[INFO]    :: 5-sig limit = 21.650
[INFO]    :: 3-sig limit = 22.205
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.080 arcsec (0.215 pixels) 
[INFO]    :: DEC error: 0.131 arcsec (0.355 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.035
[INFO]    :: Science zeropoint std: 0.027
--------------------------------------------------------------------------------
 Mag = 20.071+/-0.067 lim=22.205 MJD=61123.461
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-24T11_39871_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 22/38 | 16 remaining
[████████████████████████░░░░░░░░░░░░░░░░] 23/38 (60%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 23 OF 38 [23/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 23 OF 38 [23/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260325_06_47_25_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=20.0, thresholds: std>39.9, med>693.9) ~171 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.54 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (759,360)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61124.283
[INFO]    :: Observation time: 2026-03-25T06:47:25.399043
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 20 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 26 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=1.11", ΔDec=-2.52" from 26 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 2 inliers within 3.0" of median (ΔRA=1.11", ΔDec=-2.52") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 30 sources in sci
[INFO]    :: [V2] star_match_shift: 262 sources in ref
[INFO]    :: [V2] star_match_shift: 25 matches → dx=1.911±0.593 px, dy=0.572±1.833 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (23 inliers) → dx=1.911±0.554 px, dy=0.609±1.707 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.911, 0.609) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1992 clipped px (>= 9.95 nMgy), dilated by 9px -> 6317 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 98.0% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.110 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_85574.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_85574.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_85574.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_85574.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.3
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.189, -0.149) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.17px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-25T06_24445ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 49/124 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_85574.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.0372523
[INFO]    :: [V2] PSF centroid corrected: shift=(0.047, -0.255) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.17px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-25T06_24445sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.100075 0.038864 2.839347
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 2
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 17, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 164
[INFO]    :: Length of matched catalog= 17
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 17
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-25T06comb_psf_85574.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 17
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 16)
[29.26722052 29.34988774 29.2652942  29.27505074 29.15164431 29.32443329
 29.29011079 29.30789521 29.27624778 29.23249378 29.22993085 29.10943812
 29.24594404 29.73452818 29.2522417  29.24224441]
[0.8229819  0.75994062 0.9988411  0.99988622 0.57214827 0.99798786
 0.99987116 0.93130593 0.9192232  0.95544406 0.9999101  0.89850841
 0.99883804 0.28787208 0.98751345 0.98536903]
[31.55616721 31.70649397 31.64974828 31.62193424 31.68662648 31.69492533
 31.66101278 31.66017423 31.62408462 31.61794223 31.61653629 31.62858042
 31.60645336 31.58332649 31.65946131 31.63184025]
[0.98003045 0.98564219 0.99793677 0.99923391 0.97610923 0.99817619
 0.99954218 0.99461105 0.99482405 0.98348462 0.99940363 0.9886081
 0.99922566 0.99733996 0.99730141 0.99868835]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.2299 |  31.6165 |  0.99991  |  0.999404 |
|  29.2751 |  31.6219 |  0.999886 |  0.999234 |
|  29.2901 |  31.661  |  0.999871 |  0.999542 |
|  29.2653 |  31.6497 |  0.998841 |  0.997937 |
|  29.2459 |  31.6065 |  0.998838 |  0.999226 |
|  29.3244 |  31.6949 |  0.997988 |  0.998176 |
|  29.2522 |  31.6595 |  0.987513 |  0.997301 |
|  29.2422 |  31.6318 |  0.985369 |  0.998688 |
|  29.2325 |  31.6179 |  0.955444 |  0.983485 |
|  29.3079 |  31.6602 |  0.931306 |  0.994611 |
|  29.2762 |  31.6241 |  0.919223 |  0.994824 |
|  29.1094 |  31.6286 |  0.898508 |  0.988608 |
|  29.2672 |  31.5562 |  0.822982 |  0.98003  |
|  29.3499 |  31.7065 |  0.759941 |  0.985642 |
|  29.1516 |  31.6866 |  0.572148 |  0.976109 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 15
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 13
[INFO]    :: Number of stars used for zeropoint calculation: 13
[INFO]    :: Science Zeropoints from (13 stars)
[INFO]    ::   - ZP Mean: 29.259
[INFO]    ::   - ZP Std: 0.041
[INFO]    ::   - ZP Range: [29.152,29.324]
[INFO]    ::   - ZP Mode: 29.267
[INFO]    ::   - ZP Median: 29.265
[INFO]    ::   - ZP 16th & 84th percentiles: 29.232, 29.292
[INFO]    ::   - # above threshold (0.35): 15
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.922
[INFO]    ::   - RSQ Std: 0.117
[INFO]    ::   - RSQ Range: [0.572,1.000]
[INFO]    ::   - RSQ Mode: 0.823
[INFO]    ::   - RSQ Median: 0.985
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.841, 1.000
[INFO]    :: Reference Zeropoints from (13 stars)
[INFO]    ::   - ZP Mean: 31.643
[INFO]    ::   - ZP Std: 0.027
[INFO]    ::   - ZP Range: [31.606,31.695]
[INFO]    ::   - ZP Mode: 31.650
[INFO]    ::   - ZP Median: 31.632
[INFO]    ::   - ZP 16th & 84th percentiles: 31.618, 31.663
[INFO]    ::   - # above threshold (0.35): 15
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.993
[INFO]    ::   - RSQ Std: 0.008
[INFO]    ::   - RSQ Range: [0.976,1.000]
[INFO]    ::   - RSQ Mode: 0.980
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.984, 0.999
[INFO]    :: ZP Sci=29.265 std=0.041 No. stars=13
[INFO]    :: ZP Ref=31.632 std=0.027 No. stars=13
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.113
[INFO]    :: [V2] Blanked subtraction at 6317 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-25T06_24445bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61124.28293249989
[INFO]    :: BACKGROUND: -0.812 counts
[INFO]    :: SN FLUX 4183.132 counts
[INFO]    :: SN MAG 20.212 mag
[INFO]    :: SN MAG - BACKGROUND 20.212 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.141 1.859
[INFO]    :: xoff_arc=0.423 arcsec, yoff_arc=0.689 arcsec
[INFO]    :: [V2] Injection pool: 756160 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 22.780
[INFO]    :: S/N (std artifical sn)= 22.780
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.413 mag
[INFO]    :: Mag = 20.212+/-0.047 
[INFO]    :: 5-sig limit = 21.858
[INFO]    :: 3-sig limit = 22.413
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.164 arcsec (0.442 pixels) 
[INFO]    :: DEC error: 0.445 arcsec (1.202 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.027
[INFO]    :: Science zeropoint std: 0.041
--------------------------------------------------------------------------------
 Mag = 20.212+/-0.068 lim=22.413 MJD=61124.283
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-25T06_24445_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 23/38 | 15 remaining
[█████████████████████████░░░░░░░░░░░░░░░] 24/38 (63%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 24 OF 38 [24/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 24 OF 38 [24/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_05_51_18_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=24.7, thresholds: std>49.3, med>1060.7) ~173 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 3.55 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (782,378)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61125.244
[INFO]    :: Observation time: 2026-03-26T05:51:18.490357
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 11 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-3.37", ΔDec=0.65" from 11 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 1 inliers within 3.0" of median (ΔRA=-3.37", ΔDec=0.65") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 5 sources in sci
[INFO]    :: [V2] star_match_shift: 294 sources in ref
[INFO]    :: [V2] star_match_shift: 5 matches → dx=-0.677±3.521 px, dy=0.550±17.020 px
[WARNING] :: [V2] star_match_shift: scatter too large (MAD_x=3.52, MAD_y=17.02 px > 3.0px) — suppressing shift, relying on WCS alignment
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-43.300 px, dy=-151.200 px
[WARNING] :: [V2] Fine shift (-43.300, -151.200) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1970 clipped px (>= 9.95 nMgy), dilated by 20px -> 14095 masked (1.7% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 94.9% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.115 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88289.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88289.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88289.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_88289.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.0
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.205, -0.035) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (59, 59) (FWHM=9.57px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-26T05_21078ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 44/116 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_88289.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02255013
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.057, -0.184) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (59, 59) (FWHM=9.57px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-26T05_21078sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -1.674872 -1.203384 7.479282
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 7
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 7(<=7) stars found in reference catalog, increasing search radius: 1->2
[INFO]    :: Catalog stars in PS1/SDSS found = 9, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 153
[INFO]    :: Length of matched catalog= 9
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 9
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-26T05comb_psf_88289.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 9
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 7)
[29.31796729 29.30663864 29.29774359 29.3726991  29.3244433  29.27689586
 29.31649878]
[0.98813567 0.7957023  0.99879461 0.99763599 0.99942471 0.99862474
 0.97504305]
[31.75866606 31.68172328 31.72332214 31.81290347 31.75039458 31.681617
 31.70457568]
[0.99871199 0.80473279 0.99972348 0.99920145 0.99974164 0.9993832
 0.99954744]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3244 |  31.7504 |  0.999425 |  0.999742 |
|  29.2977 |  31.7233 |  0.998795 |  0.999723 |
|  29.2769 |  31.6816 |  0.998625 |  0.999383 |
|  29.3727 |  31.8129 |  0.997636 |  0.999201 |
|  29.318  |  31.7587 |  0.988136 |  0.998712 |
|  29.3165 |  31.7046 |  0.975043 |  0.999547 |
|  29.3066 |  31.6817 |  0.795702 |  0.804733 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 7
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 5
[INFO]    :: Number of stars used for zeropoint calculation: 5
[INFO]    :: Science Zeropoints from (5 stars)
[INFO]    ::   - ZP Mean: 29.313
[INFO]    ::   - ZP Std: 0.009
[INFO]    ::   - ZP Range: [29.298,29.324]
[INFO]    ::   - ZP Mode: 29.318
[INFO]    ::   - ZP Median: 29.316
[INFO]    ::   - ZP 16th & 84th percentiles: 29.303, 29.320
[INFO]    ::   - # above threshold (0.35): 7
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.965
[INFO]    ::   - RSQ Std: 0.070
[INFO]    ::   - RSQ Range: [0.796,0.999]
[INFO]    ::   - RSQ Mode: 0.988
[INFO]    ::   - RSQ Median: 0.998
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.968, 0.999
[INFO]    :: Reference Zeropoints from (5 stars)
[INFO]    ::   - ZP Mean: 31.724
[INFO]    ::   - ZP Std: 0.029
[INFO]    ::   - ZP Range: [31.682,31.759]
[INFO]    ::   - ZP Mode: 31.759
[INFO]    ::   - ZP Median: 31.723
[INFO]    ::   - ZP 16th & 84th percentiles: 31.696, 31.753
[INFO]    ::   - # above threshold (0.35): 7
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.972
[INFO]    ::   - RSQ Std: 0.068
[INFO]    ::   - RSQ Range: [0.805,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.991, 1.000
[INFO]    :: ZP Sci=29.316 std=0.009 No. stars=5
[INFO]    :: ZP Ref=31.723 std=0.029 No. stars=5
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.109
[INFO]    :: [V2] Blanked subtraction at 14095 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-26T05_21078bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61125.24396369979
[INFO]    :: BACKGROUND: -1.679 counts
[INFO]    :: SN FLUX 4831.783 counts
[INFO]    :: SN MAG 20.106 mag
[INFO]    :: SN MAG - BACKGROUND 20.107 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=0.113 0.113
[INFO]    :: xoff_arc=0.042 arcsec, yoff_arc=0.042 arcsec
[INFO]    :: [V2] Injection pool: 666565 valid positions in valid-data region
[INFO]    :: Background injection positions used: 438
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 10.422
[INFO]    :: S/N (std artifical sn)= 10.422
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.458 mag
[INFO]    :: Mag = 20.106+/-0.099 
[INFO]    :: 5-sig limit = 20.904
[INFO]    :: 3-sig limit = 21.458
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.866 arcsec (2.336 pixels) 
[INFO]    :: DEC error: 0.115 arcsec (0.310 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.029
[INFO]    :: Science zeropoint std: 0.009
--------------------------------------------------------------------------------
 Mag = 20.106+/-0.104 lim=21.458 MJD=61125.244
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-26T05_21078_photometry.txt
[INFO]    :: Total time: 9.0 seconds
✓ COMPLETED: 24/38 | 14 remaining
[██████████████████████████░░░░░░░░░░░░░░] 25/38 (65%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 25 OF 38 [25/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 25 OF 38 [25/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260326_07_03_27_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 141
Right X border at 909
Bottom Y border at 150
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=23.7, thresholds: std>47.3, med>966.3) ~174 bright sources remain
Cutting out image with borders: XL=141, XR=899, YL=150, YT=850 (geometry: 141,899,150,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.72 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (776,368)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61125.294
[INFO]    :: Observation time: 2026-03-26T07:03:27.713177
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 5 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-0.33", ΔDec=-5.71" from 5 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=-0.33", ΔDec=-5.71") — likely non-stellar science detections; falling back to star_match
[WARNING] :: [V2] star_match_shift: only 4 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=314.400 px, dy=361.100 px
[WARNING] :: [V2] Fine shift (314.400, 361.100) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1964 clipped px (>= 9.95 nMgy), dilated by 15px -> 10145 masked (1.2% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 95.7% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.059 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38552.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38552.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_38552.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_38552.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.8
[INFO]    :: [V2] PSF centroid corrected: shift=(0.063, 0.292) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.34px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-26T07_25407ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 43/117 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_38552.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02131113
[INFO]    :: [V2] PSF centroid corrected: shift=(0.011, -0.270) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (45, 45) (FWHM=7.34px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-26T07_25407sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -314.412562 -3.244186 419.598672
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 419
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 1(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 1(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 1(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 1, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 1
[INFO]    :: Length of matched catalog= 1
[WARNING] :: Few stars (1) matched!
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 1
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-26T07comb_psf_38552.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 1
[WARNING] :: [V2] Saturation filter would reject all 1 matched calibrators — keeping the original list.  ZP_ref scatter will likely be inflated by Legacy-Survey clipping (consider using PS1 reference for this field).
[29.35666435]
[0.9993959]
[32.39958459]
[0.99199252]
[WARNING] :: Few stars in the field, continuing if 1 or 2 stars in science
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.3567 |  32.3996 |  0.999396 |  0.991993 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 1
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 1
[INFO]    :: Number of stars used for zeropoint calculation: 1
[INFO]    :: Science Zeropoints from (1 stars)
[INFO]    ::   - ZP Mean: 29.357
[INFO]    ::   - ZP Std: 0.000
[INFO]    ::   - ZP Range: [29.357,29.357]
[INFO]    ::   - ZP Mode: 29.357
[INFO]    ::   - ZP Median: 29.357
[INFO]    ::   - ZP 16th & 84th percentiles: 29.357, 29.357
[INFO]    ::   - # above threshold (0.35): 1
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,0.999]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 0.999
[INFO]    :: Reference Zeropoints from (1 stars)
[INFO]    ::   - ZP Mean: 32.400
[INFO]    ::   - ZP Std: 0.000
[INFO]    ::   - ZP Range: [32.400,32.400]
[INFO]    ::   - ZP Mode: 32.400
[INFO]    ::   - ZP Median: 32.400
[INFO]    ::   - ZP 16th & 84th percentiles: 32.400, 32.400
[INFO]    ::   - # above threshold (0.35): 1
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.992
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.992,0.992]
[INFO]    ::   - RSQ Mode: 0.992
[INFO]    ::   - RSQ Median: 0.992
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.992, 0.992
[INFO]    :: ZP Sci=29.357 std=0.013 No. stars=1
[INFO]    :: ZP Ref=32.400 std=0.013 No. stars=1
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.061
[INFO]    :: [V2] Blanked subtraction at 10145 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-26T07_25407bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61125.29407070018
[INFO]    :: BACKGROUND: 2.127 counts
[INFO]    :: SN FLUX 5490.257 counts
[INFO]    :: SN MAG 20.008 mag
[INFO]    :: SN MAG - BACKGROUND 20.008 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-4.086 -0.430
[INFO]    :: xoff_arc=1.514 arcsec, yoff_arc=0.159 arcsec
[INFO]    :: [V2] Injection pool: 705280 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 15.235
[INFO]    :: S/N (std artifical sn)= 15.235
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.772 mag
[INFO]    :: Mag = 20.008+/-0.069 
[INFO]    :: 5-sig limit = 21.217
[INFO]    :: 3-sig limit = 21.772
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 1.023 arcsec (2.761 pixels) 
[INFO]    :: DEC error: 0.596 arcsec (1.607 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.000
[INFO]    :: Science zeropoint std: 0.000
--------------------------------------------------------------------------------
 Mag = 20.008+/-0.070 lim=21.772 MJD=61125.294
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-26T07_25407_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 25/38 | 13 remaining
[███████████████████████████░░░░░░░░░░░░░] 26/38 (68%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 26 OF 38 [26/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 26 OF 38 [26/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_05_15_47_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=29.5, thresholds: std>59.0, med>1419.3) ~122 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.09 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (777,359)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61126.219
[INFO]    :: Observation time: 2026-03-27T05:15:47.704361
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 13 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 7 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-4.18", ΔDec=0.50" from 7 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 1 inliers within 3.0" of median (ΔRA=-4.18", ΔDec=0.50") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 7 sources in sci
[INFO]    :: [V2] star_match_shift: 260 sources in ref
[INFO]    :: [V2] star_match_shift: 7 matches → dx=0.011±4.616 px, dy=0.196±1.150 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (6 inliers) → dx=0.787±3.774 px, dy=0.508±0.968 px
[WARNING] :: [V2] star_match_shift: scatter too large (MAD_x=3.77, MAD_y=0.97 px > 3.0px) — suppressing shift, relying on WCS alignment
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-343.900 px, dy=395.000 px
[WARNING] :: [V2] Fine shift (-343.900, 395.000) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 2016 clipped px (>= 9.95 nMgy), dilated by 12px -> 8281 masked (1.0% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 97.2% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.113 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_96948.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_96948.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_96948.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_96948.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.2
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.086, 0.570) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.63px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-27T05_18947ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/122 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_96948.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01253826
[INFO]    :: [V2] PSF centroid corrected: shift=(0.005, -0.229) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.63px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-27T05_18947sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -307.544106 -4.019156 547.860433
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 547
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 2(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 2, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 4
[INFO]    :: Length of matched catalog= 2
[WARNING] :: Few stars (2) matched!
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 2
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-27T05comb_psf_96948.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 2
[WARNING] :: [V2] Saturation filter would reject all 2 matched calibrators — keeping the original list.  ZP_ref scatter will likely be inflated by Legacy-Survey clipping (consider using PS1 reference for this field).
[WARNING] :: [V2] chi2_shift requested (-5.2,-0.5) px — clipped to (-5.0,-0.5) px (likely latching on a saturated-star residual)
[29.4532267  29.40796705]
[0.99857411 0.99861158]
[32.43329583 32.38479904]
[0.94152768 0.97478189]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.408  |  32.3848 |  0.998612 |  0.974782 |
|  29.4532 |  32.4333 |  0.998574 |  0.941528 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 2
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Number of stars after sigma clipping: 2
[INFO]    :: Number of stars used for zeropoint calculation: 2
[INFO]    :: Science Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 29.431
[INFO]    ::   - ZP Std: 0.023
[INFO]    ::   - ZP Range: [29.408,29.453]
[INFO]    ::   - ZP Mode: 29.453
[INFO]    ::   - ZP Median: 29.431
[INFO]    ::   - ZP 16th & 84th percentiles: 29.415, 29.446
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.999
[INFO]    ::   - RSQ Std: 0.000
[INFO]    ::   - RSQ Range: [0.999,0.999]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.999, 0.999
[INFO]    :: Reference Zeropoints from (2 stars)
[INFO]    ::   - ZP Mean: 32.409
[INFO]    ::   - ZP Std: 0.024
[INFO]    ::   - ZP Range: [32.385,32.433]
[INFO]    ::   - ZP Mode: 32.433
[INFO]    ::   - ZP Median: 32.409
[INFO]    ::   - ZP 16th & 84th percentiles: 32.393, 32.426
[INFO]    ::   - # above threshold (0.35): 2
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.958
[INFO]    ::   - RSQ Std: 0.017
[INFO]    ::   - RSQ Range: [0.942,0.975]
[INFO]    ::   - RSQ Mode: 0.942
[INFO]    ::   - RSQ Median: 0.958
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.947, 0.969
[INFO]    :: ZP Sci=29.431 std=0.023 No. stars=2
[INFO]    :: ZP Ref=32.409 std=0.024 No. stars=2
[WARNING] :: [V2] QUALITY WARNING: only 2 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.064
[INFO]    :: [V2] Blanked subtraction at 8281 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-27T05_18947bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[WARNING] :: [V2] chi2_shift requested (-11.8,-0.5) px — clipped to (-5.0,-0.5) px (likely latching on a saturated-star residual)
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61126.21930210013
[INFO]    :: BACKGROUND: -3.593 counts
[INFO]    :: SN FLUX 1893.752 counts
[INFO]    :: SN MAG 21.237 mag
[INFO]    :: SN MAG - BACKGROUND 21.239 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-5.000 -0.465
[INFO]    :: xoff_arc=1.853 arcsec, yoff_arc=0.172 arcsec
[INFO]    :: [V2] Injection pool: 731210 valid positions in valid-data region
[WARNING] :: [V2] chi2_shift requested (-11.8,-0.5) px — clipped to (-5.0,-0.5) px (likely latching on a saturated-star residual)
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 3.743
[INFO]    :: S/N (std artifical sn)= 3.743
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.478 mag
[INFO]    :: Mag = 21.237+/-0.257 
[INFO]    :: 5-sig limit = 20.923
[INFO]    :: 3-sig limit = 21.478
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.782 arcsec (2.109 pixels) 
[INFO]    :: DEC error: 0.259 arcsec (0.699 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.024
[INFO]    :: Science zeropoint std: 0.023
--------------------------------------------------------------------------------
 Mag = 21.237+/-0.259 lim=21.478 MJD=61126.219
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-27T05_18947_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 26/38 | 12 remaining
[████████████████████████████░░░░░░░░░░░░] 27/38 (71%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 27 OF 38 [27/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 27 OF 38 [27/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260327_07_19_50_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=26.1, thresholds: std>52.3, med>1140.9) ~122 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.52 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (762,367)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61126.305
[INFO]    :: Observation time: 2026-03-27T07:19:50.683435
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  39  detections along bad columns.
[INFO]    :: 31 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  39  detections along bad columns.
[INFO]    :: 31 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  39  detections along bad columns.
[INFO]    :: 31 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 27 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=0.29", ΔDec=-1.18" from 27 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 3 inliers within 3.0" of median (ΔRA=0.29", ΔDec=-1.18") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 31 sources in sci
[INFO]    :: [V2] star_match_shift: 253 sources in ref
[INFO]    :: [V2] star_match_shift: 27 matches → dx=1.507±0.602 px, dy=0.672±0.951 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (22 inliers) → dx=1.530±0.530 px, dy=0.758±0.759 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.530, 0.758) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1996 clipped px (>= 9.95 nMgy), dilated by 9px -> 6322 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 97.8% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.125 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_10921.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_10921.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_10921.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_10921.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.174, -0.158) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.10px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-27T07_26390ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 46/125 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_10921.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.01219671
[INFO]    :: [V2] PSF centroid corrected: shift=(0.036, -0.237) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (25, 25) (FWHM=4.10px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-27T07_26390sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -0.792754 -1.207846 6.171662
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 6
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 17, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 160
[INFO]    :: Length of matched catalog= 17
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 17
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-27T07comb_psf_10921.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 17
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 15)
[29.31119316 29.30731238 29.18738005 29.27742304 29.33874665 29.30497991
 29.36366438 29.33055635 29.26714446 29.25115727 29.16371439 29.25144369
 30.06466676 29.2387213  29.26738656]
[0.76012588 0.99796653 0.9674266  0.99970232 0.99772593 0.99989813
 0.86246264 0.87754459 0.93294867 0.99974481 0.75420578 0.99725153
 0.14340556 0.9679573  0.97506272]
[31.55218631 31.65402837 31.611631   31.6265483  31.69547503 31.66688795
 31.6616722  31.62593898 31.62805251 31.62036583 31.63960644 31.61027462
 31.59199189 31.66418613 31.63713962]
[0.9785965  0.99864801 0.99587582 0.99920606 0.99808371 0.99952944
 0.99356164 0.99481777 0.98422157 0.99943106 0.98806553 0.99962247
 0.99715847 0.99769144 0.99871305]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.305  |  31.6669 |  0.999898 |  0.999529 |
|  29.2512 |  31.6204 |  0.999745 |  0.999431 |
|  29.2774 |  31.6265 |  0.999702 |  0.999206 |
|  29.3073 |  31.654  |  0.997967 |  0.998648 |
|  29.3387 |  31.6955 |  0.997726 |  0.998084 |
|  29.2514 |  31.6103 |  0.997252 |  0.999622 |
|  29.2674 |  31.6371 |  0.975063 |  0.998713 |
|  29.2387 |  31.6642 |  0.967957 |  0.997691 |
|  29.1874 |  31.6116 |  0.967427 |  0.995876 |
|  29.2671 |  31.6281 |  0.932949 |  0.984222 |
|  29.3306 |  31.6259 |  0.877545 |  0.994818 |
|  29.3637 |  31.6617 |  0.862463 |  0.993562 |
|  29.3112 |  31.5522 |  0.760126 |  0.978597 |
|  29.1637 |  31.6396 |  0.754206 |  0.988066 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 14
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 12
[INFO]    :: Number of stars used for zeropoint calculation: 12
[INFO]    :: Science Zeropoints from (12 stars)
[INFO]    ::   - ZP Mean: 29.278
[INFO]    ::   - ZP Std: 0.041
[INFO]    ::   - ZP Range: [29.187,29.339]
[INFO]    ::   - ZP Mode: 29.311
[INFO]    ::   - ZP Median: 29.272
[INFO]    ::   - ZP 16th & 84th percentiles: 29.248, 29.316
[INFO]    ::   - # above threshold (0.35): 14
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.935
[INFO]    ::   - RSQ Std: 0.084
[INFO]    ::   - RSQ Range: [0.754,1.000]
[INFO]    ::   - RSQ Mode: 0.760
[INFO]    ::   - RSQ Median: 0.972
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.864, 1.000
[INFO]    :: Reference Zeropoints from (12 stars)
[INFO]    ::   - ZP Mean: 31.637
[INFO]    ::   - ZP Std: 0.019
[INFO]    ::   - ZP Range: [31.610,31.667]
[INFO]    ::   - ZP Mode: 31.654
[INFO]    ::   - ZP Median: 31.633
[INFO]    ::   - ZP 16th & 84th percentiles: 31.618, 31.662
[INFO]    ::   - # above threshold (0.35): 14
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.995
[INFO]    ::   - RSQ Std: 0.006
[INFO]    ::   - RSQ Range: [0.979,1.000]
[INFO]    ::   - RSQ Mode: 0.979
[INFO]    ::   - RSQ Median: 0.998
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.989, 0.999
[INFO]    :: ZP Sci=29.272 std=0.041 No. stars=12
[INFO]    :: ZP Ref=31.633 std=0.019 No. stars=12
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.114
[INFO]    :: [V2] Blanked subtraction at 6322 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-27T07_26390bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61126.30544769997
[INFO]    :: BACKGROUND: -3.665 counts
[INFO]    :: SN FLUX 2994.893 counts
[INFO]    :: SN MAG 20.581 mag
[INFO]    :: SN MAG - BACKGROUND 20.583 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-2.047 0.766
[INFO]    :: xoff_arc=0.759 arcsec, yoff_arc=0.284 arcsec
[INFO]    :: [V2] Injection pool: 756058 valid positions in valid-data region
[INFO]    :: Background injection positions used: 440
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 11.595
[INFO]    :: S/N (std artifical sn)= 11.595
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.049 mag
[INFO]    :: Mag = 20.581+/-0.090 
[INFO]    :: 5-sig limit = 21.495
[INFO]    :: 3-sig limit = 22.049
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.138 arcsec (0.373 pixels) 
[INFO]    :: DEC error: 0.158 arcsec (0.425 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.019
[INFO]    :: Science zeropoint std: 0.041
--------------------------------------------------------------------------------
 Mag = 20.581+/-0.101 lim=22.049 MJD=61126.305
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-27T07_26390_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 27/38 | 11 remaining
[█████████████████████████████░░░░░░░░░░░] 28/38 (73%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 28 OF 38 [28/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 28 OF 38 [28/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g.fits
[WARNING] :: Error converting SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_05_37_38_f_b_ZTF26aakjzdt_g_g.fits, error: "Keyword 'PC1_1' not found."
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 109
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=863 (interior std=118.2, thresholds: std>236.4, med>7243.8) ~6 bright sources remain
Cutting out image with borders: XL=109, XR=899, YL=10, YT=863 (geometry: 109,899,10,889  noise: 0,909,0,863)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.69 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (436,139)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61127.234
[INFO]    :: Observation time: 2026-03-28T05:37:38.424175
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  19  detections along bad columns.
[INFO]    :: 32 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  19  detections along bad columns.
[INFO]    :: 32 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  19  detections along bad columns.
[INFO]    :: 32 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[WARNING] :: [V2] ps1_catalog_shift: only 0 sources detected (need ≥5) — skipping
[WARNING] :: [V2] star_match_shift: only 0 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=205.800 px, dy=198.000 px
[WARNING] :: [V2] Fine shift (205.800, 198.000) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 134 clipped px (>= 9.95 nMgy), dilated by 15px -> 4202 masked (0.5% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.4% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=0.011 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-28T05_20258bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_65764.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_65764.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_65764.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: No source with appropriate FWHM found!!


> WARNING: No source with appropriate FWHM found!!


> WARNING: Not a positive definite matrix in homogenization solver

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_65764.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.0
[WARNING] :: Warning: PSF model is a perfect fit, may be overfitting
[WARNING] :: Trying again with only the highest signal-to-noise stars
[INFO]    :: Writing high SNR stars to /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_65764_highSNR.cat
[INFO]    :: Running PSFex with only the highest signal-to-noise stars

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: No source with appropriate FWHM found!!


> WARNING: 1st context group removed (not enough samples)


> WARNING: No source with appropriate FWHM found!!


> WARNING: Not a positive definite matrix in homogenization solver

[INFO]    :: High SNR PSFex status: 0
[INFO]    :: Reduced Chi^2 of science image PSF fit with high SNR stars: 0.0
[WARNING] :: Pipeline stopped after: PSF convolution
[██████████████████████████████░░░░░░░░░░] 29/38 (76%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 29 OF 38 [29/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 29 OF 38 [29/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260328_08_14_27_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 0
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=42.9, thresholds: std>85.9, med>3109.7) ~6 bright sources remain
Cutting out image with borders: XL=10, XR=899, YL=10, YT=889 (geometry: 10,899,10,889  noise: 0,909,0,899)
[WARNING] :: This image is not on target, please check images, passing on this image
[███████████████████████████████░░░░░░░░░] 30/38 (78%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 30 OF 38 [30/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 30 OF 38 [30/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 150
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=43.9, thresholds: std>87.7, med>3391.8) ~65 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=150, YT=889 (geometry: 50,899,150,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.94 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (774,359)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61128.218
[INFO]    :: Observation time: 2026-03-29T05:13:30.328542
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 12 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[WARNING] :: [V2] ps1_catalog_shift: only 4 sources detected (need ≥5) — skipping
[WARNING] :: [V2] star_match_shift: only 4 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=40.800 px, dy=1.800 px
[INFO]    :: [V2] [pcc] Applied fine shift (40.800, 1.800) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1966 clipped px (>= 9.95 nMgy), dilated by 11px -> 7464 masked (0.9% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 95.9% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.152 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-29T05_18810bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_53278.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_53278.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_53278.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_53278.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 3.7
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[INFO]    :: ePSF science fallback OK: FWHM=5.10 px  elong=1.08  n_stars=40  scatter=4.179 px
[INFO]    :: Using ePSF kernel for science convolution (PSFEx chi² too high)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-29T05_18810ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 42/112 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_53278.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.0278755
[INFO]    :: [V2] PSF centroid corrected: shift=(0.033, -0.270) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (33, 33) (FWHM=5.23px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-29T05_18810sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -768.080691 -18.240996 1331.10992
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 1331
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 0, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 2
[INFO]    :: Length of matched catalog= 0
[WARNING] :: Few stars (0) matched!
[WARNING] :: ZTF26aakjzdt 61128.2177122999 g: Less than 2 matched calibration stars 
[WARNING] :: Exiting: not enough stars to calibrate!
[WARNING] :: Pipeline stopped after: Reference catalog
[████████████████████████████████░░░░░░░░] 31/38 (81%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 31 OF 38 [31/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 31 OF 38 [31/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260329_06_25_09_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=45.9, thresholds: std>91.8, med>3656.5) ~82 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.69 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (751,363)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61128.267
[INFO]    :: Observation time: 2026-03-29T06:25:09.607937
[INFO]    :: Exposure time: 240.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 19 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 19 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 19 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 49 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-0.73", ΔDec=-0.83" from 49 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 4 inliers within 3.0" of median (ΔRA=-0.73", ΔDec=-0.83") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 82 sources in sci
[INFO]    :: [V2] star_match_shift: 248 sources in ref
[INFO]    :: [V2] star_match_shift: 33 matches → dx=1.277±1.957 px, dy=0.640±2.298 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (21 inliers) → dx=1.715±0.942 px, dy=0.869±1.105 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.715, 0.869) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1991 clipped px (>= 9.95 nMgy), dilated by 10px -> 6761 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 98.9% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.103 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_61945.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_61945.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_61945.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_61945.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.3
[INFO]    :: [V2] PSF centroid corrected: shift=(0.043, -0.180) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.57px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-29T06_23109ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 45/122 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_61945.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.0070047
[INFO]    :: [V2] PSF centroid corrected: shift=(0.007, -0.260) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (29, 29) (FWHM=4.57px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-29T06_23109sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -17.213512 -16.793701 9.135433
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 9
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 15, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 172
[INFO]    :: Length of matched catalog= 15
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 15
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-29T06comb_psf_61945.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 15
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 13)
[WARNING] :: [V2] chi2_shift requested (+0.1,-5.8) px — clipped to (+0.1,-5.0) px (likely latching on a saturated-star residual)
[29.31962043 29.48609372 29.4364578  29.49319169 29.4728417  30.03778668
 29.78562562 29.41725629 29.7791165  29.46500219 31.24376562 29.61409422
 29.66663667]
[0.40786361 0.99261907 0.99959647 0.99739751 0.9995389  0.71517672
 0.6179607  0.99942319 0.57267964 0.98507298 0.08411442 0.91246416
 0.86204212]
[31.55495927 31.66903523 31.63828611 31.71064474 31.67234307 31.68205757
 31.64207844 31.62276162 31.62421761 31.61963223 31.58734127 31.67063955
 31.64558987]
[0.96016884 0.99894869 0.99933733 0.99820068 0.99960633 0.9931901
 0.9938639  0.99949981 0.98367587 0.99959162 0.99677884 0.99708587
 0.99884098]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.4365 |  31.6383 |  0.999596 |  0.999337 |
|  29.4728 |  31.6723 |  0.999539 |  0.999606 |
|  29.4173 |  31.6228 |  0.999423 |  0.9995   |
|  29.4932 |  31.7106 |  0.997398 |  0.998201 |
|  29.4861 |  31.669  |  0.992619 |  0.998949 |
|  29.465  |  31.6196 |  0.985073 |  0.999592 |
|  29.6141 |  31.6706 |  0.912464 |  0.997086 |
|  29.6666 |  31.6456 |  0.862042 |  0.998841 |
|  30.0378 |  31.6821 |  0.715177 |  0.99319  |
|  29.7856 |  31.6421 |  0.617961 |  0.993864 |
|  29.7791 |  31.6242 |  0.57268  |  0.983676 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 12
[INFO]    :: [V2] Median-zp filter: dropped 1 outlier(s) >|0.30| mag from median (sci_med=29.490, ref_med=31.644; remaining: 11)
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 9
[INFO]    :: Number of stars used for zeropoint calculation: 9
[INFO]    :: Science Zeropoints from (9 stars)
[INFO]    ::   - ZP Mean: 29.537
[INFO]    ::   - ZP Std: 0.115
[INFO]    ::   - ZP Range: [29.417,29.779]
[INFO]    ::   - ZP Mode: 29.486
[INFO]    ::   - ZP Median: 29.486
[INFO]    ::   - ZP 16th & 84th percentiles: 29.444, 29.652
[INFO]    ::   - # above threshold (0.35): 11
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.850
[INFO]    ::   - RSQ Std: 0.204
[INFO]    ::   - RSQ Range: [0.408,1.000]
[INFO]    ::   - RSQ Mode: 0.408
[INFO]    ::   - RSQ Median: 0.985
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.600, 0.999
[INFO]    :: Reference Zeropoints from (9 stars)
[INFO]    ::   - ZP Mean: 31.645
[INFO]    ::   - ZP Std: 0.020
[INFO]    ::   - ZP Range: [31.620,31.672]
[INFO]    ::   - ZP Mode: 31.669
[INFO]    ::   - ZP Median: 31.642
[INFO]    ::   - ZP 16th & 84th percentiles: 31.623, 31.670
[INFO]    ::   - # above threshold (0.35): 11
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.994
[INFO]    ::   - RSQ Std: 0.011
[INFO]    ::   - RSQ Range: [0.960,1.000]
[INFO]    ::   - RSQ Mode: 0.960
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.990, 1.000
[INFO]    :: ZP Sci=29.486 std=0.115 No. stars=9
[INFO]    :: ZP Ref=31.642 std=0.020 No. stars=9
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.137
[INFO]    :: [V2] Blanked subtraction at 6761 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-29T06_23109bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[WARNING] :: [V2] chi2_shift requested (-9.4,+3.2) px — clipped to (-5.0,+3.2) px (likely latching on a saturated-star residual)
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61128.2674722001
[INFO]    :: BACKGROUND: -17.648 counts
[INFO]    :: SN FLUX 1966.977 counts
[INFO]    :: SN MAG 21.252 mag
[INFO]    :: SN MAG - BACKGROUND 21.261 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-5.000 3.164
[INFO]    :: xoff_arc=1.853 arcsec, yoff_arc=1.173 arcsec
[INFO]    :: [V2] Injection pool: 747623 valid positions in valid-data region
[WARNING] :: [V2] chi2_shift requested (-9.4,+3.2) px — clipped to (-5.0,+3.2) px (likely latching on a saturated-star residual)
[INFO]    :: Background injection positions used: 438
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 2.856
[INFO]    :: S/N (std artifical sn)= 2.856
[INFO]    :: Limiting magnitude (3-sigma scatter): 21.198 mag
[INFO]    :: Mag = 21.252+/-0.326 
[INFO]    :: 5-sig limit = 20.644
[INFO]    :: 3-sig limit = 21.198
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.219 arcsec (0.590 pixels) 
[INFO]    :: DEC error: 0.265 arcsec (0.715 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.020
[INFO]    :: Science zeropoint std: 0.115
--------------------------------------------------------------------------------
 Mag = 21.252+/-0.346 lim=21.198 MJD=61128.267
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-29T06_23109_photometry.txt
[INFO]    :: Total time: 8.0 seconds
✓ COMPLETED: 31/38 | 7 remaining
[█████████████████████████████████░░░░░░░] 32/38 (84%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 32 OF 38 [32/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 32 OF 38 [32/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g.fits
[WARNING] :: Error converting SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g.fits, error: "Keyword 'PC1_1' not found."
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 125
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=277.5, thresholds: std>555.0, med>12283.3) ~6 bright sources remain
Cutting out image with borders: XL=125, XR=899, YL=10, YT=889 (geometry: 125,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 3.99 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (436,139)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61130.269
[INFO]    :: Observation time: 2026-03-31T06:27:47.432221
[INFO]    :: Exposure time: 234.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  13  detections along bad columns.
[INFO]    :: 185 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  13  detections along bad columns.
[INFO]    :: 185 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: Removed  13  detections along bad columns.
[INFO]    :: 185 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[WARNING] :: [V2] ps1_catalog_shift: only 4 sources detected (need ≥5) — skipping
[WARNING] :: [V2] star_match_shift: only 0 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=42.000 px, dy=367.000 px
[WARNING] :: [V2] Fine shift (42.000, 367.000) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 134 clipped px (>= 9.95 nMgy), dilated by 22px -> 8032 masked (1.0% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.4% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=-0.007 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-31T06_23267bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32057.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32057.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_32057.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_32057.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 4.0
[WARNING] :: Warning: PSF model may not be accurate
[INFO]    :: Attempting ePSF fallback (build_psf) for science image …
[WARNING] :: ePSF science fallback rejected (unphysical): FWHM=54.93 px (limit 32.3)  elong=2.71 (limit 1.5) — using PSFEx kernel instead
[INFO]    :: [V2] PSF centroid corrected: shift=(0.309, 1.903) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (65, 65) (FWHM=10.77px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-31T06_23267ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 61/137 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_32057.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.04726707
[INFO]    :: [V2] PSF centroid corrected: shift=(0.053, -0.222) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (65, 65) (FWHM=10.77px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-31T06_23267sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -1988.092378 -92.319541 4065.78581
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 4065
[WARNING] :: Unhandled exception on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_06_27_47_f_b_ZTF26aakjzdt_g_g: 'NoneType' object is not subscriptable
[██████████████████████████████████░░░░░░] 33/38 (86%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 33 OF 38 [33/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 33 OF 38 [33/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260331_08_47_31_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=104.7, thresholds: std>209.5, med>4308.1) ~41 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.02 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (766,358)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61130.366
[INFO]    :: Observation time: 2026-03-31T08:47:31.110717
[INFO]    :: Exposure time: 180.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 17 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 154 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=3.57", ΔDec=4.82" from 154 nearest-neighbour pairs
[INFO]    :: [V2] ps1_catalog_shift (pass B): 5 inliers → dx=12.287±3.717 px, dy=11.188±3.169 px
[WARNING] :: [V2] ps1_catalog_shift: scatter too large (MAD_x=3.72, MAD_y=3.17 px > 2.0px) — suppressing correction
[INFO]    :: [V2] star_match_shift: 158 sources in sci
[INFO]    :: [V2] star_match_shift: 252 sources in ref
[INFO]    :: [V2] star_match_shift: 53 matches → dx=1.290±4.863 px, dy=0.852±9.858 px
[WARNING] :: [V2] star_match_shift: scatter too large (MAD_x=4.86, MAD_y=9.86 px > 3.0px) — suppressing shift, relying on WCS alignment
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-36.100 px, dy=-144.000 px
[WARNING] :: [V2] Fine shift (-36.100, -144.000) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1977 clipped px (>= 9.95 nMgy), dilated by 11px -> 7518 masked (0.9% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 97.3% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.089 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_49211.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_49211.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_49211.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_49211.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.7
[INFO]    :: [V2] PSF centroid corrected: shift=(1.984, -1.971) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (33, 33) (FWHM=5.45px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-03-31T08_31651ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/121 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_49211.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02020725
[INFO]    :: [V2] PSF centroid corrected: shift=(0.012, -0.240) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (33, 33) (FWHM=5.45px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-03-31T08_31651sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -20.835451 -20.50211 13.409265
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 13
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 11, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 192
[INFO]    :: Length of matched catalog= 11
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 11
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-03-31T08comb_psf_49211.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 11
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 9)
[WARNING] :: [V2] chi2_shift requested (-4.2,-5.1) px — clipped to (-4.2,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-3.6,-5.3) px — clipped to (-3.6,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-4.3,-5.3) px — clipped to (-4.3,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-4.1,-5.8) px — clipped to (-4.1,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-4.3,-5.6) px — clipped to (-4.3,-5.0) px (likely latching on a saturated-star residual)
[30.14800817 30.40108021 30.09839463 30.17255438 30.0585029  29.97478536
 29.99651207 30.44961016 29.97633044]
[0.41393781 0.2334429  0.43885478 0.45770651 0.43547367 0.42226305
 0.3740012  0.41637003 0.22948921]
[31.5744694  31.62581509 31.51955018 31.64694763 31.56681294 31.5217418
 31.54414222 31.6013964  31.55621688]
[0.97809476 0.26142612 0.93129626 0.98568712 0.96616376 0.98872558
 0.99443189 0.97962676 0.98898943]
+----------+----------+-----------+-----------+
| zp_sci   | zp_ref   | rsq_sci   | rsq_ref   |
|----------+----------+-----------+-----------|
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 7
[INFO]    :: [V2] Median-zp filter: dropped 1 outlier(s) >|0.30| mag from median (sci_med=30.098, ref_med=31.567; remaining: 6)
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 4
[INFO]    :: Number of stars used for zeropoint calculation: 4
[INFO]    :: Science Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 30.075
[INFO]    ::   - ZP Std: 0.055
[INFO]    ::   - ZP Range: [29.997,30.148]
[INFO]    ::   - ZP Mode: 30.148
[INFO]    ::   - ZP Median: 30.078
[INFO]    ::   - ZP 16th & 84th percentiles: 30.026, 30.124
[INFO]    ::   - # above threshold (0.35): 6
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.424
[INFO]    ::   - RSQ Std: 0.026
[INFO]    ::   - RSQ Range: [0.374,0.458]
[INFO]    ::   - RSQ Mode: 0.414
[INFO]    ::   - RSQ Median: 0.429
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.406, 0.443
[INFO]    :: Reference Zeropoints from (4 stars)
[INFO]    ::   - ZP Mean: 31.552
[INFO]    ::   - ZP Std: 0.021
[INFO]    ::   - ZP Range: [31.522,31.574]
[INFO]    ::   - ZP Mode: 31.574
[INFO]    ::   - ZP Median: 31.555
[INFO]    ::   - ZP 16th & 84th percentiles: 31.532, 31.571
[INFO]    ::   - # above threshold (0.35): 6
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.974
[INFO]    ::   - RSQ Std: 0.021
[INFO]    ::   - RSQ Range: [0.931,0.994]
[INFO]    ::   - RSQ Mode: 0.978
[INFO]    ::   - RSQ Median: 0.982
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.959, 0.990
[INFO]    :: ZP Sci=30.078 std=0.055 No. stars=4
[INFO]    :: ZP Ref=31.555 std=0.021 No. stars=4
[WARNING] :: [V2] QUALITY WARNING: only 4 ZP stars — zeropoint poorly constrained, photometry may be unreliable.
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.257
[INFO]    :: [V2] Blanked subtraction at 7518 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-03-31T08_31651bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[WARNING] :: [V2] chi2_shift requested (+10.7,+15.4) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61130.3663319
[INFO]    :: BACKGROUND: -35.265 counts
[INFO]    :: SN FLUX -1634.990 counts
[INFO]    :: SN MAG nan mag
[INFO]    :: SN MAG - BACKGROUND nan mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=5.000 5.000
[INFO]    :: xoff_arc=1.853 arcsec, yoff_arc=1.853 arcsec
[INFO]    :: [V2] Injection pool: 736855 valid positions in valid-data region
[WARNING] :: [V2] chi2_shift requested (+10.7,+15.4) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= -0.474
[INFO]    :: S/N (std artifical sn)= -0.474
[INFO]    :: Less than 2 sigma — reporting limit
[INFO]    :: 5-sig limit =19.486
[INFO]    :: 3-sig limit =20.041
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.705 arcsec (1.902 pixels) 
[INFO]    :: DEC error: 0.280 arcsec (0.756 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.021
[INFO]    :: Science zeropoint std: 0.055
[WARNING] :: Magnitude of 99 measured, SNR= -0.4738960667030816
--------------------------------------------------------------------------------
 Mag = 99.000+/-9.950 lim=20.041 MJD=61130.366
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Flux = -0.450+/-0.950
--------------------------------------------------------------------------------
[INFO]    :: Photometry time: 1.0seconds
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-03-31T08_31651_photometry.txt
[INFO]    :: Total time: 8.0 seconds
✓ COMPLETED: 33/38 | 5 remaining
[███████████████████████████████████░░░░░] 34/38 (89%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 34 OF 38 [34/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 34 OF 38 [34/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260404_04_39_27_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=101, XR=909, YL=0, YT=899 (interior std=61.7, thresholds: std>123.4, med>3176.6) ~91 bright sources remain
Cutting out image with borders: XL=101, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 101,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 3.14 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (793,349)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61134.194
[INFO]    :: Observation time: 2026-04-04T04:39:27.723213
[INFO]    :: Exposure time: 320.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 15 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 8 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=2.74", ΔDec=2.00" from 8 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=2.74", ΔDec=2.00") — likely non-stellar science detections; falling back to star_match
[WARNING] :: [V2] star_match_shift: only 3 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=41.900 px, dy=10.000 px
[INFO]    :: [V2] [pcc] Applied fine shift (41.900, 10.000) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1961 clipped px (>= 9.95 nMgy), dilated by 17px -> 11507 masked (1.4% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 93.4% of reference frame has real data
[INFO]    :: [V2] Alignment quality: NCC=0.165 ✓ (acceptable for shallow science vs deep reference)
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-04-04T04_16767bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_40230.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_40230.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_40230.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_40230.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 2.1
[INFO]    :: [V2] PSF centroid corrected: shift=(-3.142, -4.111) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (51, 51) (FWHM=8.46px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-04-04T04_16767ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 46/113 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_40230.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.02642571
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.003, -0.250) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (51, 51) (FWHM=8.46px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-04-04T04_16767sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -543.916322 -16.841062 1101.302339
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 1101
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 0, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 2
[INFO]    :: Length of matched catalog= 0
[WARNING] :: Few stars (0) matched!
[WARNING] :: ZTF26aakjzdt 61134.19407099998 g: Less than 2 matched calibration stars 
[WARNING] :: Exiting: not enough stars to calibrate!
[WARNING] :: Pipeline stopped after: Reference catalog
[████████████████████████████████████░░░░] 35/38 (92%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 35 OF 38 [35/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 35 OF 38 [35/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g.fits
[WARNING] :: Error converting SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260406_07_13_47_f_b_ZTF26aakjzdt_g_g.fits, error: "Keyword 'PC1_1' not found."
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 88
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=50.3, thresholds: std>100.7, med>3890.2) ~37 bright sources remain
Cutting out image with borders: XL=88, XR=899, YL=10, YT=889 (geometry: 88,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.89 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (436,139)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61136.301
[INFO]    :: Observation time: 2026-04-06T07:13:47.225708
[INFO]    :: Exposure time: 334.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[WARNING] :: [V2] ps1_catalog_shift: only 2 sources detected (need ≥5) — skipping
[WARNING] :: [V2] star_match_shift: only 2 sources in sci (need ≥5) — skipping
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-23.100 px, dy=415.100 px
[WARNING] :: [V2] Fine shift (-23.100, 415.100) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 134 clipped px (>= 9.95 nMgy), dilated by 11px -> 2642 masked (0.3% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.4% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=-0.002 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-04-06T07_26027bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_74069.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_74069.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_74069.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_74069.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.4
[INFO]    :: [V2] PSF centroid corrected: shift=(-1.179, 3.199) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.11px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-04-06T07_26027ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 61/137 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_74069.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.04726707
[INFO]    :: [V2] PSF centroid corrected: shift=(0.045, -0.245) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=5.11px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-04-06T07_26027sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -19.351254 -18.838459 9.795653
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 9
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 5(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 11, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 167
[INFO]    :: Length of matched catalog= 11
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 11
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-04-06T07comb_psf_74069.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 11
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 10)
[WARNING] :: [V2] chi2_shift requested (+7.3,+4.9) px — clipped to (+5.0,+4.9) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+7.3,+8.8) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+7.3,-10.9) px — clipped to (+5.0,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+7.3,-5.2) px — clipped to (+5.0,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-6.6,-5.3) px — clipped to (-5.0,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+15.3,+9.3) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-7.7,+8.7) px — clipped to (-5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+3.1,-6.3) px — clipped to (+3.1,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+1.8,+8.8) px — clipped to (+1.8,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-5.3,+7.3) px — clipped to (-5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-6.3,+10.3) px — clipped to (-5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+3.8,+6.3) px — clipped to (+3.8,+5.0) px (likely latching on a saturated-star residual)
[33.48728719 33.25726735 33.24525776 33.48186877 30.65417443         nan]
[0.07697392 0.07600422 0.07506177 0.07642572 0.24892231 0.05058673]
[28.64075487 30.39355524 30.85496707 29.97620509 29.65569829 31.17165096]
[0.00581094 0.13497364 0.28683144 0.05092913 0.02551998 0.58816659]
+----------+----------+-----------+-----------+
| zp_sci   | zp_ref   | rsq_sci   | rsq_ref   |
|----------+----------+-----------+-----------|
+----------+----------+-----------+-----------+
[WARNING] :: No stars with rsq>0.35. Lowering threshold to 0.30
[WARNING] :: No stars with rsq>0.30. Lowering threshold to 0.25
[WARNING] :: No stars with rsq>0.25. Lowering threshold to 0.20
[WARNING] :: All PSF-fit zeropoints failed for reference image (ZTF26aakjzdt g)
[WARNING] :: [V2] Attempting aperture-photometry zeropoint fallback...
[WARNING] :: [V2] Aperture ZP fallback also found 0 valid stars — exiting
[WARNING] :: Pipeline stopped after: Zeropoint calibration
[█████████████████████████████████████░░░] 36/38 (94%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 36 OF 38 [36/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 36 OF 38 [36/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260415_08_58_50_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 950
Bottom Y border at 0
Top Y border at 850
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=22.5, thresholds: std>45.0, med>763.9) ~257 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=850 (geometry: 50,899,10,850  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.29 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (793,395)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61145.374
[INFO]    :: Observation time: 2026-04-15T08:58:50.063102
[INFO]    :: Exposure time: 300.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 21 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 24 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=3.24", ΔDec=-4.83" from 24 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=3.24", ΔDec=-4.83") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 25 sources in sci
[INFO]    :: [V2] star_match_shift: 276 sources in ref
[INFO]    :: [V2] star_match_shift: 21 matches → dx=1.774±0.466 px, dy=0.475±0.411 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (19 inliers) → dx=1.784±0.388 px, dy=0.475±0.312 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.784, 0.475) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1959 clipped px (>= 9.95 nMgy), dilated by 13px -> 8723 masked (1.1% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 93.8% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.132 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_00670.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_00670.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_00670.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_00670.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.0
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.137, -0.160) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (39, 39) (FWHM=6.19px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-04-15T08_32330ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 47/116 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_00670.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.00383188
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.039, -0.198) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (39, 39) (FWHM=6.19px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-04-15T08_32330sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -0.195256 -0.112744 9.788419
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 9
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 13, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 157
[INFO]    :: Length of matched catalog= 13
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 13
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-04-15T08comb_psf_00670.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 13
[INFO]    :: [V2] Dropped 2 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 11)
[WARNING] :: [V2] chi2_shift requested (-0.4,+16.9) px — clipped to (-0.4,+5.0) px (likely latching on a saturated-star residual)
[29.81271298 29.16191396 29.80805926 29.87179472 29.83450936 29.79612257
 29.88056721 29.80300965 29.81849697 29.81862335 29.75024301]
[0.99823671 0.07708457 0.9995681  0.99831687 0.99983856 0.89496748
 0.91443399 0.99904372 0.99519    0.96991931 0.97659665]
[31.71676385 31.59158029 31.68324229 31.76235853 31.71937495 31.6711954
 31.70400961 31.66686443 31.66821332 31.72064579 31.69161169]
[0.99929434 0.49124069 0.99962971 0.99874276 0.99969158 0.99389296
 0.98129112 0.99927128 0.99968719 0.99471872 0.9988366 ]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.8345 |  31.7194 |  0.999839 |  0.999692 |
|  29.8081 |  31.6832 |  0.999568 |  0.99963  |
|  29.803  |  31.6669 |  0.999044 |  0.999271 |
|  29.8718 |  31.7624 |  0.998317 |  0.998743 |
|  29.8127 |  31.7168 |  0.998237 |  0.999294 |
|  29.8185 |  31.6682 |  0.99519  |  0.999687 |
|  29.7502 |  31.6916 |  0.976597 |  0.998837 |
|  29.8186 |  31.7206 |  0.969919 |  0.994719 |
|  29.8806 |  31.704  |  0.914434 |  0.981291 |
|  29.7961 |  31.6712 |  0.894967 |  0.993893 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 10
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 8
[INFO]    :: Number of stars used for zeropoint calculation: 8
[INFO]    :: Science Zeropoints from (8 stars)
[INFO]    ::   - ZP Mean: 29.820
[INFO]    ::   - ZP Std: 0.022
[INFO]    ::   - ZP Range: [29.796,29.872]
[INFO]    ::   - ZP Mode: 29.813
[INFO]    ::   - ZP Median: 29.816
[INFO]    ::   - ZP 16th & 84th percentiles: 29.804, 29.833
[INFO]    ::   - # above threshold (0.35): 10
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.975
[INFO]    ::   - RSQ Std: 0.037
[INFO]    ::   - RSQ Range: [0.895,1.000]
[INFO]    ::   - RSQ Mode: 0.998
[INFO]    ::   - RSQ Median: 0.997
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.939, 0.999
[INFO]    :: Reference Zeropoints from (8 stars)
[INFO]    ::   - ZP Mean: 31.697
[INFO]    ::   - ZP Std: 0.020
[INFO]    ::   - ZP Range: [31.668,31.721]
[INFO]    ::   - ZP Mode: 31.717
[INFO]    ::   - ZP Median: 31.698
[INFO]    ::   - ZP 16th & 84th percentiles: 31.673, 31.719
[INFO]    ::   - # above threshold (0.35): 10
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.997
[INFO]    ::   - RSQ Std: 0.005
[INFO]    ::   - RSQ Range: [0.981,1.000]
[INFO]    ::   - RSQ Mode: 0.999
[INFO]    ::   - RSQ Median: 0.999
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.994, 1.000
[INFO]    :: ZP Sci=29.816 std=0.022 No. stars=8
[INFO]    :: ZP Ref=31.698 std=0.020 No. stars=8
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.177
[INFO]    :: [V2] Blanked subtraction at 8723 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-04-15T08_32330bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61145.37419060012
[INFO]    :: BACKGROUND: -1.099 counts
[INFO]    :: SN FLUX 2290.819 counts
[INFO]    :: SN MAG 21.416 mag
[INFO]    :: SN MAG - BACKGROUND 21.416 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-0.074 0.629
[INFO]    :: xoff_arc=0.028 arcsec, yoff_arc=0.233 arcsec
[INFO]    :: [V2] Injection pool: 719588 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 7.825
[INFO]    :: S/N (std artifical sn)= 7.825
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.457 mag
[INFO]    :: Mag = 21.416+/-0.131 
[INFO]    :: 5-sig limit = 21.902
[INFO]    :: 3-sig limit = 22.457
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.079 arcsec (0.212 pixels) 
[INFO]    :: DEC error: 0.066 arcsec (0.177 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.020
[INFO]    :: Science zeropoint std: 0.022
--------------------------------------------------------------------------------
 Mag = 21.416+/-0.134 lim=22.457 MJD=61145.374
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-04-15T08_32330_photometry.txt
[INFO]    :: Total time: 8.0 seconds
✓ COMPLETED: 36/38 | 2 remaining
[██████████████████████████████████████░░] 37/38 (97%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 37 OF 38 [37/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 37 OF 38 [37/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g.fits
[WARNING] :: Error converting SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260419_08_40_17_f_b_ZTF26aakjzdt_g_g.fits, error: "Keyword 'PC1_1' not found."
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 76
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=26.2, thresholds: std>52.4, med>1179.9) ~40 bright sources remain
Cutting out image with borders: XL=76, XR=899, YL=10, YT=889 (geometry: 76,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 2.05 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (436,138)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61149.361
[INFO]    :: Observation time: 2026-04-19T08:40:17.784016
[INFO]    :: Exposure time: 240.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 8 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 15 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=6.59", ΔDec=18.26" from 15 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 0 inliers within 3.0" of median (ΔRA=6.59", ΔDec=18.26") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 15 sources in sci
[INFO]    :: [V2] star_match_shift: 266 sources in ref
[INFO]    :: [V2] star_match_shift: 12 matches → dx=-3.528±28.554 px, dy=7.754±34.575 px
[WARNING] :: [V2] star_match_shift: scatter too large (MAD_x=28.55, MAD_y=34.57 px > 3.0px) — suppressing shift, relying on WCS alignment
[INFO]    :: [V2] FFT phase_cross_correlation fine reg: dx=-269.200 px, dy=88.600 px
[WARNING] :: [V2] Fine shift (-269.200, 88.600) px exceeds 50-pixel safety limit — skipping
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 139 clipped px (>= 9.95 nMgy), dilated by 12px -> 3057 masked (0.4% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 89.4% of reference frame has real data
[WARNING] :: [V2] POOR ALIGNMENT: NCC=-0.000 (< 0.05) — check aligned reference image. Photometry will be unreliable.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-04-19T08_31217bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88033.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88033.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_88033.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_88033.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 1.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.801, 1.259) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.52px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-04-19T08_31217ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 61/134 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_88033.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.0487509
[INFO]    :: [V2] PSF centroid corrected: shift=(0.059, -0.234) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (35, 35) (FWHM=5.52px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-04-19T08_31217sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: -1.64268 -1.565544 3.827064
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 3
[INFO]    :: Searching for stars in reference catalog
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius: 1->2
[WARNING] :: 0(<=7) stars found in reference catalog, increasing search radius again: 2->5
[WARNING] :: 5(<=7) stars found in reference catalog, increasing search radius again: 5->7
[INFO]    :: Catalog stars in PS1/SDSS found = 11, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 286
[INFO]    :: Length of matched catalog= 11
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 11
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-04-19T08comb_psf_88033.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 11
[WARNING] :: [V2] chi2_shift requested (+7.8,+9.2) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+7.8,+1.2) px — clipped to (+5.0,+1.2) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-2.7,+14.5) px — clipped to (-2.7,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+13.3,-5.3) px — clipped to (+5.0,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+5.9,+3.5) px — clipped to (+5.0,+3.5) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+14.1,-3.2) px — clipped to (+5.0,-3.2) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-3.4,+7.7) px — clipped to (-3.4,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+9.7,-2.7) px — clipped to (+5.0,-2.7) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-11.2,+1.3) px — clipped to (-5.0,+1.3) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+5.5,+12.7) px — clipped to (+5.0,+5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (+12.3,-7.5) px — clipped to (+5.0,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-0.5,-11.7) px — clipped to (-0.5,-5.0) px (likely latching on a saturated-star residual)
[WARNING] :: [V2] chi2_shift requested (-10.5,-12.5) px — clipped to (-5.0,-5.0) px (likely latching on a saturated-star residual)
[30.95494464 31.18587216 31.04762659 29.48469875 25.75638067 26.96672959
 27.24351327 27.34999426]
[0.05671011 0.05700244 0.03615332 0.06538918 0.00091579 0.07873071
 0.16140025 0.18801295]
[        nan 31.27952164 30.77320009 30.07762589 28.8329612  29.31987727
 30.16705563         nan]
[1.27924226e-04 6.41281568e-01 2.61193138e-01 5.66866335e-02
 5.06981593e-03 1.37871521e-02 8.09056670e-02 6.15166207e-04]
+----------+----------+-----------+-----------+
| zp_sci   | zp_ref   | rsq_sci   | rsq_ref   |
|----------+----------+-----------+-----------|
+----------+----------+-----------+-----------+
[WARNING] :: No stars with rsq>0.35. Lowering threshold to 0.30
[WARNING] :: No stars with rsq>0.30. Lowering threshold to 0.25
[WARNING] :: No stars with rsq>0.25. Lowering threshold to 0.20
[WARNING] :: No stars with rsq>0.20. Lowering threshold to 0.15
[WARNING] :: All PSF-fit zeropoints failed for reference image (ZTF26aakjzdt g)
[WARNING] :: [V2] Attempting aperture-photometry zeropoint fallback...
[WARNING] :: [V2] Aperture ZP fallback also found 0 valid stars — exiting
[WARNING] :: Pipeline stopped after: Zeropoint calibration
[████████████████████████████████████████] 38/38 (100%)════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 38 OF 38 [38/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g
════════════════════════════════════════════════════════════════════════════════════════════════
  ⭐ PROCESSING FILE 38 OF 38 [38/38]
[INFO]    :: Performing image subtraction on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g
[INFO]    :: Starting reduction sequence on /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g
/Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Telescope: SEDM-P60
[INFO]    :: Trying to convert SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g.fits
[INFO]    :: Successfully converted SIP to TPV for /Users/kryanhinds/sedm_phot/data/ZTF26aakjzdt//rc20260420_07_08_19_f_a_b_ZTF26aakjzdt_g_g_tpv.fits
[INFO]    :: This is an ACQ image, trimming the image to remove the edges
[INFO]    :: Original image size: (899, 909)
[INFO]    ::Cutting out the image using the borders...
[INFO]    ::Finding borders using two methods...
[INFO]    ::Finding borders using gradient method...
X left border at 10.0, X right border at 899.0
Y bottom border at 10.0, Y top border at 889.0
[INFO]    ::Method 1 done 
Left X border at 50
Right X border at 909
Bottom Y border at 0
Top Y border at 899
[INFO]    ::Method 2 done 
[INFO]    ::[V2] Noise-based edge trim: XL=0, XR=909, YL=0, YT=899 (interior std=18.9, thresholds: std>37.8, med>579.2) ~227 bright sources remain
Cutting out image with borders: XL=50, XR=899, YL=10, YT=889 (geometry: 50,899,10,889  noise: 0,909,0,899)
[INFO]    :: Single filter: g
[INFO]    :: Seeing (observed,estimated): 1.84 0.00
[INFO]    :: Object name: ZTF26aakjzdt
[INFO]    :: Object position (RA,Dec) J2000: (14h37m16.141s +71d50m30.3176s)
[INFO]    :: Object position (RA,Dec) deg: (219.317254,71.841755)
[INFO]    :: Object position (x,y) (805,350)
[INFO]    :: Image dimensions (x,y): (909x899)
[INFO]    :: 6.108517708333334 arcmin reference width
[INFO]    :: Observation date: 61150.297
[INFO]    :: Observation time: 2026-04-20T07:08:19.870352
[INFO]    :: Exposure time: 240.0
[INFO]    :: Subtracting background
[INFO]    :: Background subtraction time: 0.0 seconds
[INFO]    :: Removing cosmic rays
[INFO]    :: Cosmic Removal time: 1.0 seconds
[INFO]    :: [DEBUG] survey_name=legacy, use_legacy_survey=True
[INFO]    :: Aligning science with reference image
[INFO]    :: Downloading Legacy Survey reference image for g-band
[INFO]    ::✓ Legacy Survey already cached: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits (shape (2000, 2000))
[INFO]    :: Reference image size verified: 2000x2000
[INFO]    :: Using Legacy Survey reference image: /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No alignment center specified!
[INFO]    :: Attempting to match science image to reference image 0
[INFO]    :: Using search radius of 100.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 1
[INFO]    :: Expanding search radius to 110.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 2
[INFO]    :: Expanding search radius to 120.0 pixels
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 23 objects detected in science image /Users/kryanhinds/sedm_phot/bkg_subtracted_science/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.fits
[INFO]    :: Running sextractor on image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: 285 objects detected in reference image /Users/kryanhinds/sedm_phot/ref_imgs/ZTF26aakjzdt_legacysurvey_g.fits
[WARNING] :: No matches.
[INFO]    :: 0 matches found between science image 1 and reference image 0
[WARNING] :: Not enough matches to nudge (0found, 2required)
[INFO]    :: Attempting to match science image to reference image 3
[INFO]    :: Expanding search radius to 130.0 pixels
[WARNING] :: Could not find enough matches to nudge
[INFO]    :: [V2] sedm_align_quick WCS update completed successfully
[INFO]    :: Original science size: (899, 909)
[INFO]    :: Original reference size: (2000, 2000)
[INFO]    :: [V2] SEDM detected: bypassing SWarp, reprojecying reference directly onto science grid
[INFO]    :: [V2] Reference WCS already in CD-matrix form (CD1_1=-7.278e-05, CRPIX=(1000.5,1000.5), CRVAL=(219.3179114499998,71.84177695000002)) — skipping normalisation
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: [V2] PS1 alignment catalog: 135 stars
[INFO]    :: [V2] ps1_catalog_shift: 28 science sources detected
[INFO]    :: [V2] ps1_catalog_shift (pass A): median offset ΔRA=-0.32", ΔDec=-0.72" from 28 nearest-neighbour pairs
[WARNING] :: [V2] ps1_catalog_shift: only 4 inliers within 3.0" of median (ΔRA=-0.32", ΔDec=-0.72") — likely non-stellar science detections; falling back to star_match
[INFO]    :: [V2] star_match_shift: 28 sources in sci
[INFO]    :: [V2] star_match_shift: 313 sources in ref
[INFO]    :: [V2] star_match_shift: 26 matches → dx=1.787±0.696 px, dy=0.982±1.798 px
[INFO]    :: [V2] star_match_shift: after sigma-clip (23 inliers) → dx=1.875±0.547 px, dy=0.879±1.573 px
[INFO]    :: [V2] [star_match] Applied fine shift (1.875, 0.879) px to reprojected reference
[INFO]    :: [V2] Pre-scaled Legacy reference by 1e+04 for PSFEx noise-model compatibility (zp_ref will absorb the offset)
[INFO]    :: [V2] Legacy Survey saturation mask: 1935 clipped px (>= 9.95 nMgy), dilated by 10px -> 6829 masked (0.8% of frame)
[INFO]    :: [V2] SEDM reproject alignment complete: 92.5% of reference frame has real data
[WARNING] :: [V2] MARGINAL ALIGNMENT: NCC=0.128 (< 0.15) — photometry may be noisy. Check aligned images.
[INFO]    :: Beginning to convolve images with SeXtractor and PSFEx
[INFO]    :: Convolving the reference with the PSF of the science image
[INFO]    :: Running SExtractor: /opt/homebrew/bin/sex /Users/kryanhinds/sedm_phot/aligned_images/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub.aa.fits -c /Users/kryanhinds/sedm_phot/config_files/prepsfex.sex -VERBOSE_TYPE QUIET -CATALOG_NAME /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_42516.cat -MAG_ZEROPOINT 25.0
[INFO]    :: Creating PSFex catalog with SExtractor

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: SExtractor status: 0
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_42516.cat
[INFO]    :: Running PSFex with SExtractor catalog: /Users/kryanhinds/sedm_phot/temp_config_files/sci_prepsfex_42516.cat -c /Users/kryanhinds/sedm_phot/config_files/psfex_conf.psfex -VERBOSE_TYPE QUIET

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFex status: 0
[INFO]    :: PSFEx science image created: /Users/kryanhinds/sedm_phot/out/proto_sci_prepsfex_42516.fits
[INFO]    :: Reduced Chi^2 of science image PSF fit: 0.9
[INFO]    :: [V2] PSF centroid corrected: shift=(-0.013, -0.213) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=4.96px)
[INFO]    :: Saving convolved reference image to /Users/kryanhinds/sedm_phot/convolved_ref/ZTF26aakjzdt_g_2026-04-20T07_25699ref_convolved.fits
[INFO]    :: Convolving the science with the PSF of the reference image

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: [V2] Ref star-locus filter: kept 45/117 sources (SNR≥20, elong≤1.3, FLUX_RADIUS within 25% of median)

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: PSFEx reference image created /Users/kryanhinds/sedm_phot/out/proto_ref_prepsfex_42516.fits
[INFO]    :: Reduced Chi^2 of reference image PSF fit: 1.03861997
[INFO]    :: [V2] PSF centroid corrected: shift=(0.037, -0.250) px
[INFO]    :: [V2] PSF trimmed+centred: (59, 59) → (31, 31) (FWHM=4.96px)
[INFO]    :: Saving convolved science image to /Users/kryanhinds/sedm_phot/convolved_sci/ZTF26aakjzdt_g_2026-04-20T07_25699sci_convolved.fits
[INFO]    :: Reference convolved image size: (899, 909)
[INFO]    :: Science convolved image size: (899, 909)
[INFO]    :: Generating reference catalogs
[INFO]    :: PS1 Catalog already downloaded to: /Users/kryanhinds/sedm_phot/ps_catalogs/ps_219.317254_71.841755_0.092553.xml
[INFO]    :: Catalog stars in PS1/SDSS found=109,333.191875arcsec search radius
[INFO]    :: Mean, Median, Std of sci_conv: 0.450861 0.349656 3.250135
[INFO]    :: Detecting stars with SExtractor
[INFO]    :: Detecting stars with IRAFStarFinder
[INFO]    :: Threshold for detecting stars: 3
[INFO]    :: Searching for stars in reference catalog
[INFO]    :: Catalog stars in PS1/SDSS found = 19, in a 450.0arcsec search radius
1
[INFO]    :: Stars detected in the image= 167
[INFO]    :: Length of matched catalog= 19
[INFO]    :: Finished the star matching process
[INFO]    :: Number of stars kept after magnitude cut: 19
[INFO]    :: Combining PSFs
[INFO]    :: Saving combined PSF to /Users/kryanhinds/sedm_phot/convolved_psf/ZTF26aakjzdt_g_2026-04-20T07comb_psf_42516.fits
[INFO]    :: Calculating zeropoints
[INFO]    :: Number of stars after removing duplicates: 19
[INFO]    :: [V2] Dropped 1 zp calibrators whose reference cutouts overlap the Legacy Survey saturation mask (remaining: 18)
[WARNING] :: [V2] chi2_shift requested (+0.6,-6.9) px — clipped to (+0.6,-5.0) px (likely latching on a saturated-star residual)
[29.47103873 29.51646804 29.54255388 29.43715917 29.52959929 29.47741787
 29.59181383 29.55886519 29.56803986 29.54416486 32.14688419 29.61796979
 29.51647098 29.40993555 29.51900714 29.54558484 29.5625652  29.53802297]
[0.81787207 0.63938495 0.9988585  0.85465035 0.99975967 0.64787121
 0.99800576 0.99990785 0.89493984 0.91610045 0.0669883  0.90713043
 0.99969838 0.89158269 0.99769837 0.98860649 0.97258173 0.99082261]
[31.54799979 31.72903219 31.67032166 31.61964992 31.63779777 31.72157021
 31.71249439 31.67805117 31.6858207  31.64025649 31.57357874 31.64418689
 31.62587111 31.61333425 31.62409553 31.59156372 31.67793069 31.64931363]
[0.94526636 0.98446406 0.9988976  0.98467839 0.99932522 0.97887729
 0.99814488 0.99963544 0.99319726 0.99422058 0.98472536 0.97974782
 0.99966477 0.97794462 0.99955908 0.99641376 0.99676338 0.99883032]
+----------+----------+-----------+-----------+
|   zp_sci |   zp_ref |   rsq_sci |   rsq_ref |
|----------+----------+-----------+-----------|
|  29.5589 |  31.6781 |  0.999908 |  0.999635 |
|  29.5296 |  31.6378 |  0.99976  |  0.999325 |
|  29.5165 |  31.6259 |  0.999698 |  0.999665 |
|  29.5426 |  31.6703 |  0.998858 |  0.998898 |
|  29.5918 |  31.7125 |  0.998006 |  0.998145 |
|  29.519  |  31.6241 |  0.997698 |  0.999559 |
|  29.538  |  31.6493 |  0.990823 |  0.99883  |
|  29.5456 |  31.5916 |  0.988606 |  0.996414 |
|  29.5626 |  31.6779 |  0.972582 |  0.996763 |
|  29.5442 |  31.6403 |  0.9161   |  0.994221 |
|  29.618  |  31.6442 |  0.90713  |  0.979748 |
|  29.568  |  31.6858 |  0.89494  |  0.993197 |
|  29.4099 |  31.6133 |  0.891583 |  0.977945 |
|  29.4372 |  31.6196 |  0.85465  |  0.984678 |
|  29.471  |  31.548  |  0.817872 |  0.945266 |
|  29.4774 |  31.7216 |  0.647871 |  0.978877 |
|  29.5165 |  31.729  |  0.639385 |  0.984464 |
+----------+----------+-----------+-----------+
[INFO]    :: Number of stars after rsq cut: 17
[WARNING] :: No stars removed after sigma clipping
[INFO]    :: Removing stars with zp values outside 5th and 95th percentiles
[INFO]    :: Number of stars after sigma clipping: 15
[INFO]    :: Number of stars used for zeropoint calculation: 15
[INFO]    :: Science Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 29.528
[INFO]    ::   - ZP Std: 0.039
[INFO]    ::   - ZP Range: [29.437,29.592]
[INFO]    ::   - ZP Mode: 29.471
[INFO]    ::   - ZP Median: 29.538
[INFO]    ::   - ZP 16th & 84th percentiles: 29.487, 29.562
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Science RSQ
[INFO]    ::   - RSQ Mean: 0.913
[INFO]    ::   - RSQ Std: 0.113
[INFO]    ::   - RSQ Range: [0.639,1.000]
[INFO]    ::   - RSQ Mode: 0.818
[INFO]    ::   - RSQ Median: 0.973
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.838, 0.999
[INFO]    :: Reference Zeropoints from (15 stars)
[INFO]    ::   - ZP Mean: 31.653
[INFO]    ::   - ZP Std: 0.036
[INFO]    ::   - ZP Range: [31.592,31.722]
[INFO]    ::   - ZP Mode: 31.670
[INFO]    ::   - ZP Median: 31.644
[INFO]    ::   - ZP 16th & 84th percentiles: 31.621, 31.684
[INFO]    ::   - # above threshold (0.35): 17
[INFO]    :: Reference RSQ
[INFO]    ::   - RSQ Mean: 0.990
[INFO]    ::   - RSQ Std: 0.014
[INFO]    ::   - RSQ Range: [0.945,1.000]
[INFO]    ::   - RSQ Mode: 0.945
[INFO]    ::   - RSQ Median: 0.996
[INFO]    ::   - RSQ 16th & 84th percentiles: 0.979, 0.999
[INFO]    :: ZP Sci=29.538 std=0.039 No. stars=15
[INFO]    :: ZP Ref=31.644 std=0.036 No. stars=15
[INFO]    :: Subtracting scaled & convolved refrence image
[INFO]    :: Scale factor: 0.144
[INFO]    :: [V2] Blanked subtraction at 6829 reference-saturated px
[INFO]    :: Saved scaled subtracted science image to: /Users/kryanhinds/sedm_phot/scaled_subtracted_imgs/ZTF26aakjzdt_g2026-04-20T07_25699bkgsub_scaled_subtraction.fits
[INFO]    :: Extracting photometry
[INFO]    :: Calculating magnitudes
[INFO]    :: Preliminary:
[INFO]    :: MJD: 61150.29745199997
[INFO]    :: BACKGROUND: -1.220 counts
[INFO]    :: SN FLUX 1200.907 counts
[INFO]    :: SN MAG 21.839 mag
[INFO]    :: SN MAG - BACKGROUND 21.840 mag
[INFO]    :: Offset of shifted PSF fit (xpos,ypos)=-1.293 1.410
[INFO]    :: xoff_arc=0.479 arcsec, yoff_arc=0.523 arcsec
[INFO]    :: [V2] Injection pool: 740890 valid positions in valid-data region
[INFO]    :: Background injection positions used: 439
[INFO]    :: Sigma clipping the background flux
[INFO]    :: Sigma clipping the artificial supernova flux
[INFO]    :: S/N (std background)= 5.796
[INFO]    :: S/N (std artifical sn)= 5.796
[INFO]    :: Limiting magnitude (3-sigma scatter): 22.554 mag
[INFO]    :: Mag = 21.839+/-0.173 
[INFO]    :: 5-sig limit = 22.000
[INFO]    :: 3-sig limit = 22.554
[INFO]    :: Checking for nearby residuals
[INFO]    :: Finding astrometric error

> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.


> WARNING: This executable has been compiled using a version of the ATLAS library without support for multithreading. Performance will be degraded.

[INFO]    :: RA error: 0.120 arcsec (0.324 pixels) 
[INFO]    :: DEC error: 0.406 arcsec (1.095 pixels) 
[INFO]    :: Systematic zeropoint error: 0.000
[INFO]    :: Reference zeropoint std: 0.036
[INFO]    :: Science zeropoint std: 0.039
--------------------------------------------------------------------------------
 Mag = 21.839+/-0.181 lim=22.554 MJD=61150.297
--------------------------------------------------------------------------------
[INFO]    :: Memory usage Mbyte:  1376.6
[INFO]    :: Photometry results written to /Users/kryanhinds/sedm_phot/grb260310a_final/ZTF26aakjzdt_g2026-04-20T07_25699_photometry.txt
[INFO]    :: Total time: 7.0 seconds
✓ COMPLETED: 38/38 | 0 remaining
---------------------------------------------------------------------------------------------
[INFO]    :: Finished SDSS-G in 4.91 minutes
+----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------+
|    | obj          | filt   |     mjd |     mag |   mag_err |   lim_mag |      ra |     dec |   exp_t |      flux |   flux_err |   seeing |        SNR |
|----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------|
|  0 | ZTF26aakjzdt | sdssg  | 61111.3 | 18.4282 | 0.0476839 |   22.2655 | 219.317 | 71.8418 |     180 |  6.23325  |  0.0606249 |    2.47  | 102.817    |
|  1 | ZTF26aakjzdt | sdssg  | 61111.4 | 18.4929 | 0.0481217 |   22.4935 | 219.317 | 71.8418 |     180 |  5.876    |  0.0491759 |    1.718 | 119.489    |
|  2 | ZTF26aakjzdt | sdssg  | 61112.3 | 18.6373 | 0.056396  |   21.3366 | 219.317 | 71.8418 |     240 | 25.6872   |  0.712612  |    2.159 |  36.0465   |
|  3 | ZTF26aakjzdt | sdssg  | 61112.5 | 18.6727 | 0.0501269 |   22.1732 | 219.317 | 71.8418 |     180 |  5.21592  |  0.0691832 |    2.281 |  75.3929   |
|  4 | ZTF26aakjzdt | sdssg  | 61113.3 | 18.9395 | 0.0501823 |   22.0876 | 219.317 | 71.8418 |     340 |  8.70063  |  0.159658  |    1.745 |  54.4953   |
|  5 | ZTF26aakjzdt | sdssg  | 61113.3 | 18.8772 | 0.0404293 |   21.8591 | 219.317 | 71.8418 |     180 |  2.8597   |  0.0611555 |    1.533 |  46.7611   |
|  6 | ZTF26aakjzdt | sdssg  | 61113.4 | 18.8958 | 0.0464834 |   22.2889 | 219.317 | 71.8418 |      90 |  1.97492  |  0.0289209 |    1.738 |  68.287    |
|  7 | ZTF26aakjzdt | sdssg  | 61115.2 | 19.2335 | 0.0692102 |   21.3933 | 219.317 | 71.8418 |     313 | 12.5238   |  0.571111  |    2.167 |  21.9288   |
|  8 | ZTF26aakjzdt | sdssg  | 61117.2 | 19.3698 | 0.053501  |   21.8922 | 219.317 | 71.8418 |     234 |  2.99703  |  0.097854  |    2.695 |  30.6276   |
|  9 | ZTF26aakjzdt | sdssg  | 61117.5 | 19.2635 | 0.0595786 |   22.3568 | 219.317 | 71.8418 |     180 |  2.79799  |  0.0540008 |    1.889 |  51.8139   |
| 10 | ZTF26aakjzdt | sdssg  | 61118.3 | 19.3919 | 0.0556415 |   22.1709 | 219.317 | 71.8418 |     180 |  2.42578  |  0.0625371 |    2.147 |  38.7894   |
| 11 | ZTF26aakjzdt | sdssg  | 61119.3 | 19.5663 | 0.0796496 |   22.4234 | 219.317 | 71.8418 |     180 |  2.02815  |  0.0486568 |    1.749 |  41.6828   |
| 12 | ZTF26aakjzdt | sdssg  | 61120.3 | 19.6715 | 0.0679295 |   22.4022 | 219.317 | 71.8418 |     180 |  1.83571  |  0.0494745 |    1.517 |  37.1042   |
| 13 | ZTF26aakjzdt | sdssg  | 61120.4 | 19.6813 | 0.0861855 |   22.614  | 219.317 | 71.8418 |     180 |  1.92749  |  0.0431299 |    1.639 |  44.6903   |
| 14 | ZTF26aakjzdt | sdssg  | 61121.2 | 19.8168 | 0.100017  |   21.5306 | 219.317 | 71.8418 |     320 |  1.96597  |  0.135181  |    2.356 |  14.5433   |
| 15 | ZTF26aakjzdt | sdssg  | 61122.2 | 19.711  | 0.0672601 |   22.1558 | 219.317 | 71.8418 |     180 |  1.8246   |  0.0639928 |    2.175 |  28.5126   |
| 16 | ZTF26aakjzdt | sdssg  | 61122.4 | 19.9286 | 0.0453235 |   22.8242 | 219.317 | 71.8418 |     180 |  1.64472  |  0.0380825 |    1.856 |  43.1882   |
| 17 | ZTF26aakjzdt | sdssg  | 61123.5 | 20.0713 | 0.0667162 |   22.2049 | 219.317 | 71.8418 |     180 |  0.966808 |  0.0451629 |    1.986 |  21.4072   |
| 18 | ZTF26aakjzdt | sdssg  | 61124.3 | 20.2115 | 0.0676455 |   22.4126 | 219.317 | 71.8418 |     180 |  1.15206  |  0.0505726 |    1.544 |  22.7803   |
| 19 | ZTF26aakjzdt | sdssg  | 61125.2 | 20.1062 | 0.103905  |   21.4583 | 219.317 | 71.8418 |     180 |  1.3307   |  0.127681  |    3.546 |  10.4221   |
| 20 | ZTF26aakjzdt | sdssg  | 61125.3 | 20.0077 | 0.0702395 |   21.772  | 219.317 | 71.8418 |     180 |  1.51205  |  0.0992478 |    2.719 |  15.2351   |
| 21 | ZTF26aakjzdt | sdssg  | 61126.2 | 21.2373 | 0.25919   |   21.4776 | 219.317 | 71.8418 |     180 |  0.521551 |  0.139326  |    2.088 |   3.74339  |
| 22 | ZTF26aakjzdt | sdssg  | 61126.3 | 20.5815 | 0.100765  |   22.0493 | 219.317 | 71.8418 |     180 |  0.824812 |  0.0711343 |    1.521 |  11.5951   |
| 23 | ZTF26aakjzdt | sdssg  | 61128.3 | 21.2516 | 0.346331  |   21.1981 | 219.317 | 71.8418 |     240 |  0.541718 |  0.189689  |    1.694 |   2.85583  |
| 24 | ZTF26aakjzdt | sdssg  | 61130.4 | 99      | 9.95005   |   20.0411 | 219.317 | 71.8418 |     180 | -0.450287 |  0.95018   |    2.021 |  -0.473896 |
| 25 | ZTF26aakjzdt | sdssg  | 61145.4 | 21.4156 | 0.133955  |   22.4566 | 219.317 | 71.8418 |     300 |  0.630906 |  0.0806234 |    2.293 |   7.82534  |
| 26 | ZTF26aakjzdt | sdssg  | 61150.3 | 21.8392 | 0.180814  |   22.5543 | 219.317 | 71.8418 |     240 |  0.330737 |  0.0570595 |    1.84  |   5.79636  |
+----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------+
---------------------------------------------------------------------------------------------
[INFO]    :: For SDSS-R, there are 0 fits
[INFO]    :: For SDSS-R, there are 0 fits
[WARNING] :: No fits files found in filter: SDSS-R
---------------------------------------------------------------------------------------------
[INFO]    :: For SDSS-I, there are 0 fits
[INFO]    :: For SDSS-I, there are 0 fits
[WARNING] :: No fits files found in filter: SDSS-I
---------------------------------------------------------------------------------------------
[INFO]    :: For SDSS-Z, there are 0 fits
[INFO]    :: For SDSS-Z, there are 0 fits
[WARNING] :: No fits files found in filter: SDSS-Z
---------------------------------------------------------------------------------------------
[INFO]    :: For SDSS-U, there are 0 fits
[INFO]    :: For SDSS-U, there are 0 fits
[WARNING] :: No fits files found in filter: SDSS-U
---------------------------------------------------------------------------------------------
[INFO]    :: For Bessell-V, there are 0 fits
[INFO]    :: For Bessell-V, there are 0 fits
[WARNING] :: No fits files found in filter: Bessell-V
---------------------------------------------------------------------------------------------
[INFO]    :: For Bessell-R, there are 0 fits
[INFO]    :: For Bessell-R, there are 0 fits
[WARNING] :: No fits files found in filter: Bessell-R
---------------------------------------------------------------------------------------------
[INFO]    :: For Bessell-I, there are 0 fits
[INFO]    :: For Bessell-I, there are 0 fits
[WARNING] :: No fits files found in filter: Bessell-I
---------------------------------------------------------------------------------------------
[INFO]    :: For Bessell-B, there are 0 fits
[INFO]    :: For Bessell-B, there are 0 fits
[WARNING] :: No fits files found in filter: Bessell-B
---------------------------------------------------------------------------------------------
[INFO]    :: Completed photometry on  in g in 294.51 seconds 
---------------------------------------------------------------------------------------------
+----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------+
|    | obj          | filt   |     mjd |     mag |   mag_err |   lim_mag |      ra |     dec |   exp_t |      flux |   flux_err |   seeing |        SNR |
|----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------|
|  0 | ZTF26aakjzdt | sdssg  | 61111.3 | 18.4282 | 0.0476839 |   22.2655 | 219.317 | 71.8418 |     180 |  6.23325  |  0.0606249 |    2.47  | 102.817    |
|  1 | ZTF26aakjzdt | sdssg  | 61111.4 | 18.4929 | 0.0481217 |   22.4935 | 219.317 | 71.8418 |     180 |  5.876    |  0.0491759 |    1.718 | 119.489    |
|  2 | ZTF26aakjzdt | sdssg  | 61112.3 | 18.6373 | 0.056396  |   21.3366 | 219.317 | 71.8418 |     240 | 25.6872   |  0.712612  |    2.159 |  36.0465   |
|  3 | ZTF26aakjzdt | sdssg  | 61112.5 | 18.6727 | 0.0501269 |   22.1732 | 219.317 | 71.8418 |     180 |  5.21592  |  0.0691832 |    2.281 |  75.3929   |
|  4 | ZTF26aakjzdt | sdssg  | 61113.3 | 18.9395 | 0.0501823 |   22.0876 | 219.317 | 71.8418 |     340 |  8.70063  |  0.159658  |    1.745 |  54.4953   |
|  5 | ZTF26aakjzdt | sdssg  | 61113.3 | 18.8772 | 0.0404293 |   21.8591 | 219.317 | 71.8418 |     180 |  2.8597   |  0.0611555 |    1.533 |  46.7611   |
|  6 | ZTF26aakjzdt | sdssg  | 61113.4 | 18.8958 | 0.0464834 |   22.2889 | 219.317 | 71.8418 |      90 |  1.97492  |  0.0289209 |    1.738 |  68.287    |
|  7 | ZTF26aakjzdt | sdssg  | 61115.2 | 19.2335 | 0.0692102 |   21.3933 | 219.317 | 71.8418 |     313 | 12.5238   |  0.571111  |    2.167 |  21.9288   |
|  8 | ZTF26aakjzdt | sdssg  | 61117.2 | 19.3698 | 0.053501  |   21.8922 | 219.317 | 71.8418 |     234 |  2.99703  |  0.097854  |    2.695 |  30.6276   |
|  9 | ZTF26aakjzdt | sdssg  | 61117.5 | 19.2635 | 0.0595786 |   22.3568 | 219.317 | 71.8418 |     180 |  2.79799  |  0.0540008 |    1.889 |  51.8139   |
| 10 | ZTF26aakjzdt | sdssg  | 61118.3 | 19.3919 | 0.0556415 |   22.1709 | 219.317 | 71.8418 |     180 |  2.42578  |  0.0625371 |    2.147 |  38.7894   |
| 11 | ZTF26aakjzdt | sdssg  | 61119.3 | 19.5663 | 0.0796496 |   22.4234 | 219.317 | 71.8418 |     180 |  2.02815  |  0.0486568 |    1.749 |  41.6828   |
| 12 | ZTF26aakjzdt | sdssg  | 61120.3 | 19.6715 | 0.0679295 |   22.4022 | 219.317 | 71.8418 |     180 |  1.83571  |  0.0494745 |    1.517 |  37.1042   |
| 13 | ZTF26aakjzdt | sdssg  | 61120.4 | 19.6813 | 0.0861855 |   22.614  | 219.317 | 71.8418 |     180 |  1.92749  |  0.0431299 |    1.639 |  44.6903   |
| 14 | ZTF26aakjzdt | sdssg  | 61121.2 | 19.8168 | 0.100017  |   21.5306 | 219.317 | 71.8418 |     320 |  1.96597  |  0.135181  |    2.356 |  14.5433   |
| 15 | ZTF26aakjzdt | sdssg  | 61122.2 | 19.711  | 0.0672601 |   22.1558 | 219.317 | 71.8418 |     180 |  1.8246   |  0.0639928 |    2.175 |  28.5126   |
| 16 | ZTF26aakjzdt | sdssg  | 61122.4 | 19.9286 | 0.0453235 |   22.8242 | 219.317 | 71.8418 |     180 |  1.64472  |  0.0380825 |    1.856 |  43.1882   |
| 17 | ZTF26aakjzdt | sdssg  | 61123.5 | 20.0713 | 0.0667162 |   22.2049 | 219.317 | 71.8418 |     180 |  0.966808 |  0.0451629 |    1.986 |  21.4072   |
| 18 | ZTF26aakjzdt | sdssg  | 61124.3 | 20.2115 | 0.0676455 |   22.4126 | 219.317 | 71.8418 |     180 |  1.15206  |  0.0505726 |    1.544 |  22.7803   |
| 19 | ZTF26aakjzdt | sdssg  | 61125.2 | 20.1062 | 0.103905  |   21.4583 | 219.317 | 71.8418 |     180 |  1.3307   |  0.127681  |    3.546 |  10.4221   |
| 20 | ZTF26aakjzdt | sdssg  | 61125.3 | 20.0077 | 0.0702395 |   21.772  | 219.317 | 71.8418 |     180 |  1.51205  |  0.0992478 |    2.719 |  15.2351   |
| 21 | ZTF26aakjzdt | sdssg  | 61126.2 | 21.2373 | 0.25919   |   21.4776 | 219.317 | 71.8418 |     180 |  0.521551 |  0.139326  |    2.088 |   3.74339  |
| 22 | ZTF26aakjzdt | sdssg  | 61126.3 | 20.5815 | 0.100765  |   22.0493 | 219.317 | 71.8418 |     180 |  0.824812 |  0.0711343 |    1.521 |  11.5951   |
| 23 | ZTF26aakjzdt | sdssg  | 61128.3 | 21.2516 | 0.346331  |   21.1981 | 219.317 | 71.8418 |     240 |  0.541718 |  0.189689  |    1.694 |   2.85583  |
| 24 | ZTF26aakjzdt | sdssg  | 61130.4 | 99      | 9.95005   |   20.0411 | 219.317 | 71.8418 |     180 | -0.450287 |  0.95018   |    2.021 |  -0.473896 |
| 25 | ZTF26aakjzdt | sdssg  | 61145.4 | 21.4156 | 0.133955  |   22.4566 | 219.317 | 71.8418 |     300 |  0.630906 |  0.0806234 |    2.293 |   7.82534  |
| 26 | ZTF26aakjzdt | sdssg  | 61150.3 | 21.8392 | 0.180814  |   22.5543 | 219.317 | 71.8418 |     240 |  0.330737 |  0.0570595 |    1.84  |   5.79636  |
+----+--------------+--------+---------+---------+-----------+-----------+---------+---------+---------+-----------+------------+----------+------------+
(sedm_phot) dhcp-vl2041-52810:sedm_phot kryanhinds$ 

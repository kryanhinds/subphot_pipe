Meeting with Christoffer

- Apply a low order shift to try and improve local astrometry
- Median background
  - Fit polynomial to whole image with clipping
  - Should be constant but sometimes needs a gradient (smooth through the whole image)
  - Increase filter size
    - Look at masking bright sources
- Where PSFex fails
  - Take cutout and sum to 1 (normalise)
    - Assume it is the PSF
- More detailed error messsages
- PANSTARRS filtering
  - List of stars in science, fit PSF to stars and actual stars return 1
    - scipy.stats.linregress
- More checks
    - Count fraction of negative pixels at the position of the transient
        - Half of the PSF size
    - Put a cut
- Error term
    - Uncertainty in alignment
        - Ofek optimal subtraction
        - Simulate
            - Take WCS uncertainty (rms in arc or pixel)
            - Create simulated images and photometry by shifting reference and get std (100)
            - Assume sigma and a gaussian and find shift
- Upload subtracted image fits/png for EP
    - Check header for specific to EP


 - Minar Webpage
 - Minar Server
from astropy.time import Time
import sys
print()
for mjd in sys.argv[1:]:
    if ':' not in mjd:
        print(mjd,'=====>',Time(float(mjd), format='mjd').isot)
    else:
        print(mjd,'=====>',Time(mjd, format='isot').mjd)
    print()
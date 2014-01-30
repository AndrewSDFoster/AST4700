#!/usr/bin/env python
'''
dms2deg    - converst from dms to deg
hms2deg    - converts from hms to deg
deg2dms    - converts from deg to dms
deg2hms    - converts from deg to hms
AngSepReal - uses spherical trig to find angle between two points
AngSepEucl - estimates angle between two points with a right triangle
AngSepPole - estimates angle between two points with a polar triangle
'''

import numpy as np

def dms2deg(DEC):
 '''DEC of form [degree, arcminute, arcsecond], returns as decimal'''
 #distribute the sign of the first nonzero entry to the others
 dec = DEC.copy()
 if dec[0] == 0:
  if dec[1] < 0:
   dec[2] *= -1
 elif dec[0] < 0:
  dec[1] *= -1
  dec[2] *= -1

 #check input
 if (not ((dec[0] <= 0 and dec[1] <= 0 and dec[2] <= 0) \
      or (dec[0] >= 0 and dec[1] >= 0 and dec[2] >= 0))\
    or dec[0] > 90 or dec[0] < -90 \
    or dec[1] > 60 or dec[1] < -60 \
    or dec[2] > 60 or dec[2] < -60):
  print("ERROR bad input for dms2deg, trying to continue anyway")
 
 return dec[0] + dec[1]/60. + dec[2]/3600.

def hms2deg(ra):
 '''ra of form [hour, minute, second], returns as decimal'''
 #check input
 if (ra[0] < 0 or ra[0] > 24 \
  or ra[1] < 0 or ra[1] > 60 \
  or ra[2] < 0 or ra[2] > 60):
  print("ERROR bad input for hms2deg, trying to compute anyway")

 return (ra[0] + ra[1]/60. + ra[2]/3600.)*15

def deg2dms(c):
 deg = int(  c)                     
 amn = int( (c-deg)*60.)            
 asc = int((((c-deg)*60.)-amn)*60.) 

 return np.array([deg, amn, asc])

def deg2hms(c):
 return deg2dms(c*15)

def AngSepReal(ra1, dec1, ra2, dec2):
 '''  ra1 and  ra2 are lists of the form [  hour, minute, second]
     dec1 and dec2 are lists of the form [degree, arcmin, arcsec]
     returns angle between them
 '''
 #get degrees of each angle
 RA1  = hms2deg( ra1)
 RA2  = hms2deg( ra2)
 DEC1 = dms2deg(dec1)
 DEC2 = dms2deg(dec2)

 #A is angle at pt 2, a is side length across from A (from pole to pt 1)
 #B is angle at pt 1, b is side length across from B (from pole to pt 2)
 #C is angle at pole, c is side length across from C (from pt 1 to pt 2)

 #Find angle C (difference of RA's)
 C = np.abs(RA1-RA2)

 #find sides a and b (90-dec)
 a = 90 - DEC1
 b = 90 - DEC2

 #correct to radians
 a *= np.pi/180.
 b *= np.pi/180.
 C *= np.pi/180.

 #spherical law of cosines
 c = np.arccos(np.cos(a)*np.cos(b) + np.sin(a)*np.sin(b)*np.cos(C))

 #back to degrees
 c *= 180./np.pi

 return deg2dms(c)

def AngSepEucl(ra1, dec1, ra2, dec2):
 '''  ra1 and  ra2 are lists of the form [  hour, minute, second]
     dec1 and dec2 are lists of the form [degree, arcmin, arcsec]
     returns angle between them estimated from euclidean right triangle
 '''

 #convert to degrees
 RA1  = hms2deg( ra1)
 RA2  = hms2deg( ra2)
 DEC1 = dms2deg(dec1)
 DEC2 = dms2deg(dec2)

 #find avg dec, and the differences in dec and ra
 aDEC = (DEC1 + DEC2) / 2.
 dDEC = np.abs(DEC1 - DEC2)
 dRA  = np.abs( RA1 -  RA2)

 #find read "dist" in RA
 rdRA = dRA*np.cos(aDEC*np.pi/180.)

 #compute length of hypotenuse
 dist = np.sqrt(rdRA*rdRA + dDEC*dDEC)

 return deg2dms(dist)

def AngSepPole(ra1, dec1, ra2, dec2):
 '''  ra1 and  ra2 are lists of the form [  hour, minute, second]      
     dec1 and dec2 are lists of the form [degree, arcmin, arcsec]      
     returns angle between them as estimated by a euclidean polar triangle
 '''                                                                   
 #get degrees of each angle                                            
 RA1  = hms2deg( ra1)                                                  
 RA2  = hms2deg( ra2)                                                  
 DEC1 = dms2deg(dec1)                                                  
 DEC2 = dms2deg(dec2)                                                  
                                                                       
 #A is angle at pt 2, a is side length across from A (from pole to pt 1)
 #B is angle at pt 1, b is side length across from B (from pole to pt 2)
 #C is angle at pole, c is side length across from C (from pt 1 to pt 2)
                                                                       
 #Find angle C (difference of RA's)                                    
 C = np.abs(RA1-RA2)                                                   
                                                                       
 #find sides a and b (90-dec)                                          
 a = 90 - DEC1                                                         
 b = 90 - DEC2                                                         

 #switch to using south pole if it is closer
 if (a + b)/2 > 90:
  a = 180 - a
  b = 180 - b
                                                                       
 #correct to radians                                                   
 a *= np.pi/180.                                                       
 b *= np.pi/180.                                                       
 C *= np.pi/180.                                                       
                                                                       
 #euclidean law of cosines                                             
 c = np.sqrt(a*a + b*b - 2*a*b*np.cos(C))    
                                                                       
 #back to degrees                                                      
 c *= 180./np.pi                                                       
                                                                       
 return deg2dms(c)                                                     

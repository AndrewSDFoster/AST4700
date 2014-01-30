#!/usr/bin/env python
'''
dms2deg    - converst from dms to deg
hms2deg    - converts from hms to deg
deg2dms    - converts from deg to dms
deg2hms    - converts from deg to hms
AngSepReal - uses spherical trig to find angle between two points
AngSepEucl - estimates angle between two points with a right triangle
AngSepPole - estimates angle between two points with a polar triangle
sphericalLawOfCosines - what the name suggests
sphericalLawOfSines   - ^
euclideanLawOfCosines - ^^
euclideanLawOfSines   - ^^^
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
 '''deg is decimal degrees, converts to array of [degrees, minutes, seconds]'''
 deg = int(  c)                     
 amn = int( (c-deg)*60.)            
 asc = int((((c-deg)*60.)-amn)*60.) 

 return np.array([deg, amn, asc])

def deg2hms(c):
 '''c is decimal degrees, converts to an array of [hours, minutes, seconds]'''
 return deg2dms(c*15)

def AngSepReal(ra1, dec1, ra2, dec2):
 '''  ra1 and  ra2 are lists of the form [  hour, minute, second]
     dec1 and dec2 are lists of the form [degree, arcmin, arcsec]
     returns angle between them
 '''
 #A is angle at pt 2, a is side length across from A (from pole to pt 1)
 #B is angle at pt 1, b is side length across from B (from pole to pt 2)
 #C is angle at pole, c is side length across from C (from pt 1 to pt 2)

 #Find angle C (difference of RA's)
 C = np.abs(ra1-ra2)

 #find sides a and b (90-dec)
 a = 90 - dec1
 b = 90 - dec2

 #spherical law of cosines
 c = sphericalLawOfCosines(a=a,b=b,c=None,C=C)

 return c

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
 C = deg2dms(np.abs(RA1-RA2))
                                                                       
 #find sides a and b (90-dec)                                          
 a = 90 - DEC1                                                         
 b = 90 - DEC2                                                         

 #switch to using south pole if it is closer
 if (a + b)/2 > 90:
  a = 180 - a
  b = 180 - b
                                                                       
 #euclidean law of cosines                                             
 c = euclideanLawOfCosines(a=a, b=b, c=None, C=C)
                                                                       
 return deg2dms(c)                                                     

def sphericalLawOfCosines(a, b, c, C=None):
 '''Law of cosines, lowercase variables are side lengths (dms arrays)
    either c or C (opposite side/angle pair) is calculated from the
    other parameters.
 '''
 #if angle must be found
 if C == None:
  #do math, notice the conversion to/from radians
  a = dms2deg(a)*np.pi/180
  b = dms2deg(b)*np.pi/180
  c = dms2deg(c)*np.pi/180
  C = np.arccos((np.cos(c) - np.cos(a)*np.cos(b))/np.sin(a)*np.sin(b))

  return deg2dms(C*180/np.pi)

 #if side must be found
 if c == None:
  #do math, notice the conversion to/from radians
  a = dms2deg(a)*np.pi/180 
  b = dms2deg(b)*np.pi/180 
  C = dms2deg(C)*np.pi/180
  c = np.arccos(np.cos(a)*np.cos(b) + np.sin(a)*np.sin(b)*np.cos(C))

  return deg2dms(c*180/np.pi)

 #tell the user they fucked up
 else:
  print("Error, law of cosines invalid parameters")
  return np.array([0,0,0])

def sphericalLawOfSines(angle1, side1, angle2, side2=None):
 '''Law of sines. 4 possible inputs, angle1 and side1 must be dms arrays and 
    either angle2 or side2 must also be dms arrays, with the other as None
    This function will determing the missing side/angle from the other three
    parameters
 '''
 #if side2 must be found
 if side2 == None:
  #do math, notice the conversion to/from radians
  angle1 = dms2deg(angle1)*np.pi/180
  angle2 = dms2deg(angle2)*np.pi/180
  side1  = dms2deg( side2)*np.pi/180
  side2  = np.arcsin(np.sin(angle2)*np.sin(side1)/np.sin(angle1))*180/np.pi

  return deg2dms(side2)

 #if angle2 must be found
 if angle2 == None:
  #do math, notice the conversion to/from radians
  angle1 = dms2deg(angle1)*np.pi/180
  side1  = dms2deg( side1)*np.pi/180
  side2  = dms2deg( side2)*np.pi/180
  angle2 = np.arcsin(np.sin(side2)*np.sin(angle1)/np.sin(side1))*180/np.pi

  return deg2dms(angle2)

 #tell the user that they fucked up
 else:
  print('Error, law of sines overconstrained, returning 0')
  return np.array([0,0,0])

def euclideanLawOfCosines(a, b, c, C=None):
 ''' euclidean law of cosines
    parameters of the sides of a triangle and one angle,
    when given 3, the fourth is found. C and c are opposite each other
 '''
 #if angle must be found
 if C == None:
  #do math (notice degree radian conversion)
  C = np.arccos((c*c - a*a - b*b)/(-2*a*b))*180/np.pi
  return deg2dms(C)

 #if side must be found
 if c == None:
  #do math (notice degree radian conversion)
  C = dms2deg(C)*np.pi/180
  c = np.sqrt(a*a + b*b - 2*a*b*np.cos(C))
  return c

 #tell the user that they fucked up
 else:
  print("invalid parameters for euclidean law of cosines")
  return 0

def euclideanLawOfSines(A, a, B, b = None):
 ''' euclidean Law of sines
    parameters of opposite side/angle pairs (A/a and B/b)
    from three values, the fourth is calculated
 '''
 #if side must be found
 if b == None:
  #do math, notice radians/degrees conversions
  A = dms2deg(A)*np.pi/180
  B = dms2deg(B)*np.pi/180
  b = np.sin(B)*a/np.sin(A)
  return b

 #if angle must be found
 if B == None:
  #do math, notice radians/degrees conversions
  A = dms2deg(A)*np.pi/180
  B = np.arcsin(b*np.sin(A)/a)*180/np.pi
  return deg2dms(B)

 #tell the user that they fucked up
 else:
  print("invalid parameters for euclidean law of sines")
  return 0

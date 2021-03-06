#!/usr/bin/env python
'''
dms2deg    - converst from dms to deg
hms2deg    - converts from hms to deg
deg2dms    - converts from deg to dms
deg2hms    - converts from deg to hms
AngSepReal - uses spherical trig to find angle between two points
AngSepEucl - estimates angle between two points with a right triangle
AngSepPole - estimates angle between two points with a polar triangle
equ2ecl    - converts equatorial to ecliptic coordinates 
ecl2gal    - converts ecliptic to galactic coordinates   
gal2equ    - converts galactic to equatorial coordinates 
equ2gal    - converts equatorial to galactic coordinates 
gal2ecl    - converts galactic to ecliptic coordinates   
ecl2equ    - converts ecliptic to equatorial coordinates 
fuck horizon coordinates. Seriously. I'm not doing that shit.
sphericalLawOfCosines - what the name suggests
sphericalLawOfSines   - ^
euclideanLawOfCosines - ^^
euclideanLawOfSines   - ^^^
EquinoxToJ2000        - Converts from any equinox to J2000
EpochWithJ2000equinox - finds a nearby epoch with J2000 equinox
B1950toJ2000          - converts from the B1950 equinox to the J2000 equinox
refractionAnglehor    - computes angle of atmospheric refraction
trueAltitudehor	      - computes true altitude from apparent
apparentAltitudehor   - computes apparent altitude from true
'''

import numpy as np

def dms2deg(DEC):
 '''DEC of form [degree, arcminute, arcsecond], returns as decimal'''
 #distribute the sign of the first nonzero entry to the others
 dec = np.float64(DEC.copy())
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
 #type conversion
 ra = np.float64(ra.copy())
 #check input
 if (ra[0] < 0 or ra[0] > 24 \
  or ra[1] < 0 or ra[1] > 60 \
  or ra[2] < 0 or ra[2] > 60):
  print("ERROR bad input for hms2deg, trying to compute anyway")

 return (ra[0] + ra[1]/60. + ra[2]/3600.)*15

def deg2dms(c):
 '''deg is decimal degrees, converts to array of [degrees, minutes, seconds]'''
 deg = int(  c)                     
 amn = np.abs(int( (c-deg)*60.))
 asc = ((np.abs((c-deg)*60.))-amn)*60.

 if c < 0 and deg == 0:
  if amn == 0:
   asc *= -1
  else:
   amn *= -1

 return np.array([deg, amn, asc])

def deg2hms(c):
 '''c is decimal degrees, converts to an array of [hours, minutes, seconds]'''
 return deg2dms(c/15)

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

def equ2ecl((alpha, delta)):
 '''accepts tuple of dms arrays, returns the same.
    converts hms RA and dms DEC into dms beta and lambda'''
 #define constants and get things into decimals/radians
 epsilon  = dms2deg([23, 26, 21]) * np.pi/180
 alpha    = hms2deg(alpha) * np.pi/180
 delta    = dms2deg(delta) * np.pi/180

 #calculate beta
 beta     = np.arcsin(np.sin(delta)*np.cos(epsilon) - \
                      np.cos(delta)*np.sin(epsilon)*np.sin(alpha))

 #calculate cos and sin of lmbda
 coslmbda = np.cos(delta)*np.cos(alpha)/np.cos(beta)
 sinlmbda = (np.sin(delta) - np.cos(epsilon)*np.sin(beta)) / \
                (np.sin(epsilon)*np.cos(beta))
 #use arctan2 to get the right quadrant
 lmbda    = np.arctan2(sinlmbda, coslmbda)

 #bring it back to dms
 beta     = deg2dms( beta * 180/np.pi)
 lmbda    = deg2dms(lmbda * 180/np.pi)

 #return a tuple
 return (beta, lmbda)

def equ2gal((alpha, delta)):
 '''accepts a tuple of dms arrays, returns the same.
    converts hms RA and dms DEC into dms b an l'''
 #define constants and get things into decimals/radians
 deltag = dms2deg([27,07,42])   * np.pi/180
 alphag = hms2deg([12,51,26.3]) * np.pi/180
 lnode  = dms2deg([32,55,55])   * np.pi/180
 delta  = dms2deg(delta)        * np.pi/180
 alpha  = hms2deg(alpha)        * np.pi/180

 #calculate b
 b      = np.arcsin(np.sin(deltag)*np.sin(delta) + \
                    np.cos(deltag)*np.cos(delta)*np.cos(alpha-alphag))

 #calculate cos and sin of l
 cosl   = np.cos(delta)*np.sin(alpha-alphag)/np.cos(b)
 sinl   = (np.sin(delta)-np.sin(deltag)*np.sin(b))/(np.cos(deltag)*np.cos(b))

 #use arctan2 to get the right quadrant
 l      = np.arctan2(sinl, cosl)

 #readd lnode back in
 l     += lnode

 #bring it back to dms
 l      = deg2dms(l * 180/np.pi)
 b      = deg2dms(b * 180/np.pi)

 #return a tuple
 return (b,l)

def ecl2equ((beta, lmbda)):
 '''accepts a tuple of dms arrays, returns the same.
    converts dms beta and lmbda to hms RA and dms DEC'''
 #define constants and get things into decimals/radians
 epsilon  = dms2deg([23, 26, 21]) * np.pi/180
 beta     = dms2deg(beta)  * np.pi/180
 lmbda    = dms2deg(lmbda) * np.pi/180

 #calculate delta
 delta    = np.arcsin(np.sin(beta)*np.cos(epsilon) + \
                      np.cos(beta)*np.sin(epsilon)*np.sin(lmbda))

 #calculate cos and sin of alpha
 cosalpha = np.cos(lmbda)*np.cos(beta)/np.cos(delta)
 sinalpha = (np.cos(epsilon)*np.sin(delta) - np.sin(beta)) / \
                    (np.sin(epsilon)*np.cos(delta))
 #use arctan2 to get the right quadrant
 alpha    = np.arctan2(sinalpha, cosalpha)

 #bring it back to dms
 alpha    = deg2hms(alpha * 180/np.pi)
 delta    = deg2dms(delta * 180/np.pi)

 #return a tuple
 #return (alpha, delta)

def gal2equ((b, l)):
 '''accepts a tuple of dms arrays, returns the same.
    converts dms b and l to hms RA and dms DEC'''
 #define constants and get things into decimals/radians
 deltag = dms2deg([27,07,42])   * np.pi/180
 alphag = hms2deg([12,51,26.3]) * np.pi/180
 lnode  = dms2deg([32,55,55])   * np.pi/180
 b      = dms2deg(b)            * np.pi/180
 l      = hms2deg(l)            * np.pi/180

 #calculate delta
 delta  = np.arcsin(np.sin(deltag)*np.sin(b) + \
                    np.cos(deltag)*np.cos(b)*np.sin(l-lnode))

 #calculate cos and sin of alpha
 cosalp = (np.sin(b)-np.sin(deltag)*np.sin(delta))/(np.cos(deltag)*np.cos(delta))
 sinalp = np.cos(l-lnode)*np.cos(b)/np.cos(delta)

 #use arctan2 to get the right quadrant
 alpha  = np.arctan2(sinalp, cosalp)

 #add alphag back in
 alpha += alphag

 #bring it back to dms
 alpha  = deg2hms(alpha * 180/np.pi)
 delta  = deg2dms(delta * 180/np.pi)

 #return a tuple
 return (alpha, delta)

def ecl2gal((beta, lmbda)):
 '''accepts a tuple of dms arrays, returns the same.
    converts dms beta and lambda to dms b and l'''
 return equ2gal(ecl2equ((beta, lmbda)))

def gal2ecl((b, l)):
 '''accepts a tuple of dms arrays, returns the same.
    converts dms b and l to dms beta and lambda'''
 return equ2ecl(gal2equ((b, l)))

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

def EquinoxToJ2000(alpha, delta, pmA, pmD, date, BJD=False):
 '''converts Ephemeri from one time to J2000.0
    alpha and delta are the location of the star at start (RA and dec)
    pmA and pmD are the proper motions in alpha and delta in arcsec/yr
    returns alpha and delta in J2000.0
    Can also be done in BJD rather than years
 '''
 #compute time for time standard
 if BJD:
  T = np.float64(date-2451545.0)/36525.
  year = (T*100)+2000.0
 else:
  T = np.float64(date-2000.0)/100.
  year = date

 #Compute precession constants for time and get everything into radians
 M = (1.2812323*T + 0.0003879*T*T + 0.0000101*T*T*T)*np.pi/180.
 N = (0.5567530*T - 0.0001185*T*T - 0.0000116*T*T*T)*np.pi/180.
 alpha  = hms2deg(alpha)*np.pi/180
 delta  = dms2deg(delta)*np.pi/180

 #find the mean epoch for each time
 alpham = alpha - 0.5*(M + N*np.sin(alpha)*np.tan(delta))
 deltam = delta - 0.5*N*np.cos(alpham)

 #find the new location of the star's old position
 alpha0 = alpha - M - N*np.sin(alpham)*np.tan(deltam)
 delta0 = delta - N*np.cos(alpham)

 #return to hms and dms
 delta0 = deg2dms(delta0*180/np.pi)
 alpha0 = deg2hms(alpha0*180/np.pi)


 #Account for proper motions
 (alphaf, deltaf) = EpochWithJ2000equinox(alpha0, delta0, pmA, pmD, year)

 #return
 return (alphaf, deltaf)

def EpochWithJ2000equinox(alpha0, delta0, pmA, pmD, date, BJD=False):
 '''Uses J2000 Equinox and proper motions to find new locations of  stars
    at a different epoch. Takes alpha0 and delta0 as the locations at epoch and
    equinox of J2000, and pmA, pmD, the proper motions in alpha and delta and
    calculates the starses positsdjolfjslifjoisjd ljals lksjdlkjk fuck you
 '''
 #years of how many of them go past since yeah
 #wow I was tired when I wrote this, making it more readable
 #get the amount of time in years since 2000.0
 if BJD:
  years = np.abs((date-2451545.0)/365.)
 else:
  years = np.abs(date-2000.0)

 delta0 = dms2deg(delta0)
 alpha0 = hms2deg(alpha0)

 #correct for proper motions, be sure to correct for cos(dec) factor in RA
 deltaf = delta0 + pmD*years/3600.
 #average delta, converted to radians inside of cosine
 alphaf = alpha0 + (pmA*years/3600.)/np.cos(((delta0+deltaf)/2.)*np.pi/180.)

 #back to dms/hms
 deltaf = deg2dms(deltaf)
 alphaf = deg2hms(alphaf)

 #return the values
 return (alphaf, deltaf)

def B1950toJ2000(alpha, delta, pmA, pmD):
 '''Converts from B1950 equinox/epoch to J2000 equinox/epoch
    takes alpha, delta, and proper motions in each for B1950
    returns alpha and delta in J2000
    NOTE: This is just a wrapper for EquinoxToJ2000() with the equinox
          coming from B1950
 '''
 return EquinoxToJ2000(alpha, delta, pmA, pmD, 2433282.423, BJD=True)

def refractionAnglehor(Aapp):
 ''' Finds the angle of refraction of an object at an apparent altitude of Aapp
     Aapp given in typical [deg, min, sec] numpy array
     NOTE: THIS IS AN APPROXIMATION THAT ONLY WORKS CLOSE TO THE HORIZON
 '''
 c0 =  35.338/60
 c1 = -13.059/60
 c2 =   2.765/60
 c3 =  -0.244/60

 #get to degs
 Aapp  = dms2deg(Aapp)

 #formula
 theta = deg2dms(c0 + c1*Aapp + c2*Aapp**2 + c3*Aapp**3)

 return theta

def trueAltitudehor(Aapp):
 ''' wrapper that finds true altitude from the apparent altitude, Aapp,
     Aapp given in typical [deg, min, sec] numpy array
     NOTE: THIS IS AN APPROXIMATION THAT ONLY WORKS CLOSE TO THE HORIZON
 '''
 #this is basically just a wrapper function
 a = deg2dms(dms2deg(Aapp) - dms2deg(refractionAnglehor(Aapp)))
 return a

def apparentAltitudehor(a):
 ''' Finds apparent altitude from true altitude, given in typical [deg,min,sec]
     numpy array
     NOTE: THIS IS AN APPROXIMATION THAT ONLY WORKS CLOSE TO THE HORIZON
 '''
 c0 =  35.338/60
 c1 = -13.059/60
 c2 =   2.765/60
 c3 =  -0.244/60

 #this is the polynomial
 roots = np.roots((-c3, -c2, 1-c1, -c0-dms2deg(a)))

 #only return the real one
 for root in roots:
  if np.imag(root) == 0:
   Aapp = deg2dms(np.real(root))

 return Aapp

def trueAltZen(a):
 '''Finds true altitude from apparent using approximation near zenith
    a in [degrees,minutes,seconds] numpy array
    NOTE: THIS IS AN APPROXIMATION THAT ONLY WORKS CLOSE TO THE HORIZON
 '''
 n = 1.0002923
 return deg2dms(np.arcsin(n*np.sin(dms2deg(a)*np.pi/180))*180/np.pi)

def apparAltZen(z):
 '''Finds apparent altitude from apparent using approximation near zenith
    z in [degrees,minutes,seconds] numpy array
    NOTE: THIS IS AN APPROXIMATION THAT ONLY WORKS CLOSE TO THE HORIZON
 '''
 n = 1.0002923
 return deg2dms(np.arcsin(np.sin(dms2deg(z)*np.pi/180)/n)*180/np.pi)

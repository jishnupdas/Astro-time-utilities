#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 15:36:48 2019

@author: jishnu
"""

#%%
from astropy.time import Time
from astropy import coordinates as coord, units as u
from astropy.coordinates import EarthLocation

#%%
class Times:
    """Class to represent and convert between different time formats
    JD,MJD,BJD,FITS,UTC etc..."""
    
    def __init__(self,UT=None,FITS=None,JD=None,RA=None,DEC=None,
                 MJD=None,target=None,longitude=None,latitude=None,
                 elevation=None,location=None,BJD=None):
        self.UT        = UT
        self.JD        = JD
        self.RA        = RA
        self.DEC       = DEC
        self.BJD       = BJD
        self.MJD       = MJD
        self.FITS      = FITS
        self.target    = target
        self.longitude = longitude
        self.latitude  = latitude
        self.elevation = elevation
        
    def set_jd(self,JD):
        """create an astropy time object from Julian Date format"""
        self.JD = Time(JD,format='jd',scale='utc')
    
    def set_ut(self,UT):
        """create an astropy time object from UT format"""
        self.UT = Time(UT,scale='utc')
    
    def set_fits(self,FITS):
        """create an astropy time object from FITS format time"""
        self.FITS = Time(FITS,format='fits',scale='utc')
    
    def set_mjd(self,MJD):
        """Create an astropy time object from MJD string"""
        self.MJD = Time(MJD,format='mjd',scale='utc')
        
    def jd_to_cal(self):
        """Convert JD into Calendar datetime"""
        return self.JD.isot
        
    def mjd_to_cal(self):
        """Convert Modified JD into Calendar datetime"""
        return self.MJD.isot
    
    def ut_to_jd(self):
        """Convert UT time into Julian Date"""
        return self.UT.jd
    
    def fits_to_jd(self):
        """Convert FITS time into Julian Date"""
        return self.FITS.jd
    
    def set_target(self,RA,DEC):
        """Create a SkyCoord object from the given RA and DEC"""
        self.RA = RA
        self.DEC = DEC
        self.target = coord.SkyCoord(self.RA,self.DEC,
#                                     unit=(u.hourangle,u.deg),
                                     frame='icrs')
    
    def set_location(self,latitude,longitude,elevation):
        """Create an EarthLocation Object from given lat,lon and elevation"""
        self.elevation = elevation
        self.longitude,self.latitude = longitude,latitude
        self.location = EarthLocation.from_geodetic(self.longitude,
                                                    self.latitude,
                                                    self.elevation)
    def jd_to_bjd(self):
        """COnvert the given Julian date into Barycentric Julian Date,
        the Function requires the location on earth from where observation 
        was made. the SkyCoord of the target and the JD to be converted"""
        try:
            self.JD=Time(self.JD,format='jd',scale='utc',
                         location=self.location)
            ltt_bary = self.JD.light_travel_time(self.target)
            self.BJD = self.JD.tdb + ltt_bary
            return self.BJD
        except:
            print('''Please make sure you have set all the necessary 
                  parameters,(RA,DEC,*Location)''')
                        
        
#Example:
    '''
RA = ['23h49m03.25s','02h22m19.83s','03h21m39.62s']
DEC = ['41d19m26.30s','-23d24m55.93s','47d27m18.79s']

#%%

'VBO Kavalur'
longitude = '78d49m18s'
latitude = '12d34m28s'
elevation = 725*u.m

#%%

x = Times()
x.set_jd(2451373.377923)
x.set_target(RA[0],DEC[0])
x.set_location(latitude,longitude,elevation)

#%%

print('JD  : {}, {}'.format(x.JD,x.JD.iso))
print('BJD : {}, {}'.format(x.jd_to_bjd().value,x.BJD.iso))
'''

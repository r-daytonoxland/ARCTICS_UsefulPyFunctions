"""For quickly converting location and UTC time to Solar elevation angle
Rowan Dayton-Oxland for ARCTICS, 6 Jun 2024
"""

# Libraries
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Angle, get_body
from astropy.time import Time
from sunpy.coordinates import frames
import pandas as pd


def get_data(filename):
    """Read csv file to pandas dataframe
    Input: filename (string)
    Output: df pandas dataframe """
    return pd.read_csv(filename)


def sun_loc(obstime):
    """Get centre of solar disk SkyCoord object at observation time
    Input: obstime (astopy.time Time object) 
    Output: astropy.SkyCoord location
    """
    return SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=obstime, observer="earth", frame=frames.Helioprojective)


def moon_loc(obstime):
    """ Get moon location from Earth at observation time as SkyCoord object
    Input: obstime (astopy.time Time object) 
    Output: astropy.SkyCoord location"""
    return get_body(body="moon", time=obstime)


def get_elevation(obstime, latitude, longitude, moon=None, both=None):
    """ Default to sun, for moon instead add keyword moon
    Input: obstime (astropy.time Time object)
           latitude (float) in degrees
           longitude (float) in degrees
    Kwarg: moon (optional argument) switches to moon elevation
    Output: astropy.coordinates Angle object"""

    # Observers location (hard-coded to be at sea level but this could be changed)
    obsloca = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=0*u.km)

    # Define observers frame
    frame_altaz = AltAz(obstime=Time(obstime), location=obsloca)

    # Get location of body in altaz
    if moon:
        m = moon_loc(obstime)
        altaz_moon = m.transform_to(frame_altaz)
        moon_elevation = Angle(altaz_moon.T.alt)
        return moon_elevation
    if both:
        s = sun_loc(obstime)
        altaz_sun = s.transform_to(frame_altaz)
        sun_elevation = Angle(altaz_sun.T.alt)
        return sun_elevation, moon_elevation
    else:
        return sun_elevation


def convert_todataframe(filename, moon=None, both=None):
    """Convert a csv file to a pandas dataframe including an angle"""
    
    df = get_data(filename)

    els=[]
    for index, row in df.iterrows():
        el = get_elevation(Time(row['time']), row['lat'], row['lon'], moon=moon, both=both)
        els.append(el)

    df['Angle'] = els

    return df


def convert_csv(filename, tofile=None, moon=None, both=None):

    df = convert_todataframe(filename, moon=moon, both=both)
    if tofile:
        df.to_csv(tofile)
    else:
        df.to_csv(filename)


""" Example use """

# Get the solar elevation angle from direct input
obstime = "2024-09-21 19:00:00" # UTC
latitude = 56
longitude = -5
get_elevation(obstime, latitude, longitude)

# Convert a csv file of time, lat, lon data to include solar elevation angle
filename = 'timelatlon.csv'
convert_csv(filename, tofile='with_sunelevation.csv')
convert_csv(filename, tofile='with_moonelevation.csv', moon=True)
"""For quickly converting location and UTC time to Solar elevation angle"""

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Angle
from astropy.time import Time
from sunpy.coordinates import frames
import pandas as pd


def get_data(filename):

    return pd.read_csv(filename)


def sun_loc(obstime):

    return SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=obstime, observer="earth", frame=frames.Helioprojective)


def moon_loc(obstime):

    return get_body(body="moon", obstime=obstime, observer="earth")


def get_elevation(obstime, latitude, longitude, moon=None):
    """Default to sun, for moon add keyword moon"""

    obsloca = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=0*u.km)
    
    # Define centre of the solar disk
    if moon:
        c = moon_loc(obstime)
    else:
        c = sun_loc(obstime)

    # Define observers frame
    frame_altaz = AltAz(obstime=Time(obstime), location=obsloca)
    
    # Transform sun location to observers altaz frame
    altaz = c.transform_to(frame_altaz)
    
    # Return solar elevation angle
    return Angle(altaz.T.alt)




def convert_todataframe(filename):
    
    df = get_data(filename)

    els=[]
    for index, row in df.iterrows():
        el = get_elevation(Time(row['time']), row['lat'], row['lon'])
        els.append(el)

    df['Angle'] = els

    return df


def convert_csv(filename, tofile=None):

    df = convert_todataframe(filename)
    if tofile:
        df.to_csv(tofile)
    else:
        df.to_csv(filename)


""" Example use """

# Get the solar elevation angle from direct input
obstime = "2024-09-21 19:00:00" # UTC
latitude = 56
longitude = -5
get_solel(obstime, latitude, longitude)

# Convert a csv file of time, lat, lon data to include solar elevation angle
filename = 'timelatlon.csv'
convert_csv(filename, tofile='with_sunelevation.csv')
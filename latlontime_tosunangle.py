import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from sunpy.coordinates import frames, sun

def get_solel(obstime, latitude, longitude):

    obsloca = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=0*u.km)
    
    # Define centre of the solar disk
    c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=obstime,
             observer="earth", frame=frames.Helioprojective)

    # Define observers frame
    frame_altaz = AltAz(obstime=Time(obstime), location=obsloca)
    
    # Transform sun location to observers altaz frame
    sun_altaz = c.transform_to(frame_altaz)
    
    # Return solar elevation angle
    return sun_altaz.T.alt


""" Example conversion """
obstime = "2024-09-21 19:00:00"
latitude = 56
longitude = -5

get_solel(obstime, latitude, longitude)
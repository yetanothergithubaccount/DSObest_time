#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Collection of Solveighs astro calculation helpers
#
#

import datetime
from datetime import date
import pytz
from skyfield.api import load, wgs84, N, W
from astropy.coordinates import AltAz
from astropy.time import Time
import ephem
import config
from skyfield.framelib import ecliptic_frame
import decimal

dec = decimal.Decimal
debug = False

eph = load('de421.bsp') # will be downloaded at first load

def compass_direction(azimuth):
  direction = ""
  '''
  N: 0
  NE: 45
  E: 90
  ES: 135
  S: 180
  SW: 225
  W: 270
  WN: 315
  '''
  if azimuth >= 0 and azimuth < 15:
    direction = "N"
  if azimuth >= 15 and azimuth < 30:
    direction = "NNE"
  if azimuth >= 30 and azimuth < 60:
    direction = "NE"
  if azimuth >= 60 and azimuth < 75:
    direction = "ENE"
  if azimuth >= 75 and azimuth < 105:
    direction = "E"
  if azimuth >= 105 and azimuth < 135:
    direction = "ESE"
  if azimuth >= 135 and azimuth < 150:
    direction = "SE"
  if azimuth >= 150 and azimuth < 165:
    direction = "SSE"
  if azimuth >= 165 and azimuth < 195:
    direction = "S"
  if azimuth >= 195 and azimuth < 225:
    direction = "SSW"
  if azimuth >= 225 and azimuth < 240:
    direction = "SW"
  if azimuth >= 240 and azimuth < 255:
    direction = "WSW"
  if azimuth >= 255 and azimuth < 285:
    direction = "W"
  if azimuth >= 285 and azimuth < 300:
    direction = "WNW"
  if azimuth >= 300 and azimuth < 330:
    direction = "NW"
  if azimuth >= 330 and azimuth < 345:
    direction = "NWN"
  if azimuth >= 345 and azimuth <= 360:
    direction = "N"
  return direction

def observation_night_directions(the_object, the_object_name, today, tomorrow, utcoffset, the_location):
  try:
    # observation directions 20 pm .. 4 am
    theDate_today = today.strftime("%Y-%m-%d")
    theDate_tomorrow = tomorrow.strftime("%Y-%m-%d")

    # 20 pm
    time = Time(str(theDate_today) + " 18:59:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_20 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_20))

    # 22 pm
    time = Time(str(theDate_today) + " 20:59:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_22 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_22))

    # 24 pm
    time = Time(str(theDate_today) + " 21:59:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_0 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_0))

    # 2 am
    time = Time(str(theDate_tomorrow) + " 00:00:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_2 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_2))

    # 4 am
    time = Time(str(theDate_tomorrow) + " 01:59:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_4 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_4))

    # 6 am
    time = Time(str(theDate_tomorrow) + " 03:59:00") + utcoffset
    if debug:
      print(time)
    the_object_altaz = the_object.transform_to(AltAz(obstime=time, location=the_location))
    to_alt = the_object_altaz.alt
    to_az = the_object_altaz.az
    if debug:
      print(str(the_object_name) + "'s altitude = " + str(to_alt) + ", azimut = " + str(to_az))
    direction_6 = compass_direction(to_az.value)
    if debug:
      print(str(time) + ": " + str(direction_6))

    return direction_20, direction_22, direction_0, direction_2, direction_4, direction_6
  except Exception as e:
    print(str(e))


def astro_night_times(theDate, latitude, longitude, debug):
  civil_night_start = None
  civil_night_end = None
  nautical_night_start = None
  nautical_night_end = None
  astronomical_night_start = None
  astronomical_night_end = None

  today = datetime.date.today()
  date_today = datetime.datetime.strptime(theDate, "%d.%m.%Y")
  theDay = date_today.strftime("%d")
  theMonth = date_today.strftime("%m")
  theYear = date_today.strftime("%Y")
  today = today.replace(day=int(theDay), month=int(theMonth), year=int(theYear))
  tomorrow = today + datetime.timedelta(days=1)
  date_tomorrow = datetime.datetime.strptime(tomorrow.strftime("%d.%m.%Y"), "%d.%m.%Y")

  earth = ephem.Observer()
  earth.lat = str(latitude)
  earth.lon = str(longitude)
  earth.date = date_today

  sun = ephem.Sun()
  sun.compute()

  try:
    earth.horizon = "0"
    #sunset = ephem.localtime(earth.next_setting(sun))
    #sunrise = ephem.localtime(earth.next_rising(sun))

    earth.horizon = "-6"
    earth.date = date_today
    sun.compute()
    civil_night_start = ephem.localtime(earth.next_setting(sun))
    earth.date = date_tomorrow  # make sure to hit the next day's rising
    sun.compute()
    civil_night_end = ephem.localtime(earth.next_rising(sun))

    earth.horizon = "-12"
    earth.date = date_today
    sun.compute()
    nautical_night_start = ephem.localtime(earth.next_setting(sun))
    earth.date = date_tomorrow  # make sure to hit the next day's rising
    sun.compute()
    nautical_night_end = ephem.localtime(earth.next_rising(sun))

    earth.horizon = "-18"
    earth.date = date_today
    sun.compute()
    astronomical_night_start = ephem.localtime(earth.next_setting(sun))
    earth.date = date_tomorrow  # make sure to hit the next day's rising
    sun.compute()
    astronomical_night_end = ephem.localtime(earth.next_rising(sun))

  # ephem throws an "AlwaysUpError" when there is no astronomical twilight (which occurs in summer in nordic countries)
  except ephem.AlwaysUpError:
    if debug:
      print("No astronomical night at the moment: " + str(theDate))

  if debug:
    print("Civil night start: " + str(civil_night_start))
    print("Civil night end: " + str(civil_night_end))
    print("Nautical night start: " + str(nautical_night_start))
    print("Nautical night end: " + str(nautical_night_end))
    print("Astronomical night start: " + str(astronomical_night_start))
    print("Astronomical night end: " + str(astronomical_night_end))

  return civil_night_start, civil_night_end, nautical_night_start, nautical_night_end, astronomical_night_start, astronomical_night_end

def moon_data(theDate, theTime):
  for_date = theDate.split(".")
  for_time = theTime.split(":")
  for_date = datetime.datetime(int(for_date[2]), int(for_date[1]), int(for_date[0]), int(for_time[0]), int(for_time[1]))
  #print(for_date)

  home = ephem.Observer()
  home.lat, home.lon = str(config.coordinates['latitude']), str(config.coordinates['longitude'])
  home.date = for_date

  moon = ephem.Moon()
  moon.compute(home)

  tz_germany = pytz.timezone(config.coordinates['timezone'])
  moon_rise = ephem.localtime(home.next_rising(moon)).astimezone(tz_germany).strftime("%d.%m.%Y %H:%M")
  moon_set  = ephem.localtime(home.next_setting(moon)).astimezone(tz_germany).strftime("%d.%m.%Y %H:%M")
  full_moon = ephem.localtime(ephem.next_full_moon(home.date)).strftime("%d.%m.%Y")

  ts = load.timescale()
  t = ts.utc(int(theDate.split(".")[2]), int(theDate.split(".")[1]), int(theDate.split(".")[0]), int(for_time[0]), int(for_time[1]))

  sun, moon, earth = eph['sun'], eph['moon'], eph['earth']
  #e = earth.at(t)
  mylocation = earth + wgs84.latlon(49.878708* N, 8.646927*W)
  l = mylocation.at(t)
  s = l.observe(sun).apparent()
  m = l.observe(moon).apparent()

  _, slon, _ = s.frame_latlon(ecliptic_frame)
  _, mlon, _ = m.frame_latlon(ecliptic_frame)
  moon_phase = (mlon.degrees - slon.degrees) % 360.0
  moon_phase_percent = 100.0 * m.fraction_illuminated(sun)

  mlat, mlon, mdist = m.frame_latlon(ecliptic_frame)
  #observer = wgs84.latlon(home.lat, home.lon)
  #apparent = (earth + observer).at(t).observe(moon).apparent()
  alt, az, distance = m.altaz()

  if debug:
    print("Moonrise: " + moon_rise)
    print("Moonset: " + moon_set)
    print("Next full moon: " + str(full_moon))
    print("Phase (0°–360°): " + str(moon_phase))
    print("Percent illuminated: " + str(moon_phase_percent))
    print("Moon lat " + str(mlat.degrees) + " lon " + str(mlon.degrees) + " dist " + str(mdist))

  return moon_rise, moon_set, full_moon, round(moon_phase,0), round(moon_phase_percent,2), round(alt.degrees,0), round(az.degrees,0), distance

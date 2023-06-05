# file for code

# Import and Style Statements
import matplotlib
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.formatter.useoffset"] = False
import numpy as np
import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from astropy.visualization import astropy_mpl_style, quantity_support
import datetime
import pandas as pd
import os
from IPython.display import display, Markdown
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
plt.style.use(astropy_mpl_style)
quantity_support()
import lightkurve as lk
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Spacer, Table
from reportlab.platypus import TableStyle, Paragraph, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
import io
from reportlab.lib.units import inch
from matplotlib.backends.backend_pdf import PdfPages
PAGE_SIZE = landscape(letter)
from fractions import Fraction

# Initialize Location and Data
path = 'tessdata2.csv'
def load_tess_data(filename=path):
    '''Loads data on all TOIs from a csv file
    
    Args:
        filename (str): the path of the csv file with the data
    
    Returns:
        data (DataFrame): read in dataframe
    '''
    data = pd.read_csv(filename, header=4)
    return data

# Define Stone Edge Observatory as 'SEO' for use as default location 
SEO = EarthLocation(lat=38.29*u.deg, lon=-122*u.deg, height=56*u.m)

# Load in the catalogue of known TOIs and their information
DATA = load_tess_data(path)

# Define TESS time offset
TESST = 2457000

class TOI():
    '''Class for a specific TOI, this holds the information for a given star
    and is worked with by the important functions later on. 
    '''
    
    def __init__(self, toi, datatable=None):
        '''Pull information on a given TOI from the database

        Args:
            toi (float, str): TOI ID for the target star
            datatable (DataFrame): the full TESS df, or None and load it in
        '''
        if datatable == None:
            datatable = load_tess_data()
            
        toirow = datatable.iloc[(datatable[datatable['Full TOI ID'] == toi].index)[0]]

        self.name = str(toi)
        self.transit0 = toirow['Epoch Value']
        self.ra = toirow['TIC Right Ascension']
        self.dec = toirow['TIC Declination']
        self.period = toirow['Orbital Period Value']
        self.duration = toirow['Transit Duration Value']
        self.depth = toirow['Transit Depth Value']
        self.tmag = toirow['TMag Value']
        self.comment = toirow['Public Comment']
        self.tic = toirow['TIC']

def calculate_transit_times(toi, fraction=1, timerange=(2457000, 2461000), mode='mid'):
    '''A function that calculates all the possible transit times in the 
    BJD scale for the entire range of TESS data for the given fraction for
    a specific TOI planet.
    
    Args:
        toi (TOI): the TOI to calculate for
        fraction (int): highest fraction to divide period by
        timerange (tuple(int, int)): JD day range
        mode ('mid' or 'ends'): if mode is set to mid, the function returns
        the TESS time mid transit times, if set to ends, each time is a tuple
        of Time objects of start and end of transit.
    
    Returns:
        midtranstimes (list(float)): if mode == mid
        transtimes (list(tuple(float, float))): list of possible transit 
        ingress and egress times in JD 
    '''
    period = toi.period
    epoch = toi.transit0
    
    # find first predicted transit time before timerange start
    prior = epoch + TESST
    while prior > timerange[0]:
        prior -= period
   
    all_transits = []
    
    
    # find all even fraction periods using fraction of period
    if fraction != 1:
        if (fraction % 2) == 0:
            even_max = fraction
        else:
            even_max = fraction - 1
        evens = []
        for num in range(2, even_max + 1):
            if (num % 2 == 0):
                evens.append(num)
        
        for frac in evens:
            even_frac_period = period / frac
            trans_time = prior
            while trans_time < timerange[1]:
                if trans_time > timerange[0]:
                    all_transits.append(trans_time)
                trans_time += even_frac_period
    
    # find all odd fraction periods
    if ((fraction % 2) == 1):
        odd_max = fraction
    else:
        odd_max = fraction - 1
    odds = []
    for num in range(1, odd_max + 1):
        if (num % 2 == 1):
            odds.append(num)
            
    for frac in odds:
        odd_frac_period = period / frac  
        trans_time = prior
        while trans_time < timerange[1]:
            if trans_time > timerange[0]:
                all_transits.append(trans_time)
            trans_time += odd_frac_period
    
    # if there are no transits in the range
    if len(all_transits) == 0:
        return []
    
    # order and remove duplicates
    all_transits = np.sort(all_transits)
    all_transits_u = [all_transits[0]]
    for i, t in enumerate(all_transits):
        prev = all_transits_u[-1]
        nex = t
        if (t - prev) > 1:
            all_transits_u.append(t)
    
    
    if mode == 'mid':
        midtranstimes = [(x - TESST) for x in all_transits_u]
        return midtranstimes
    
    transdays = [x * u.day for x in all_transits_u]
    halfdur = (toi.duration / 2) * u.hour
    transtimesM = [Time(x, format='jd') for x in transdays]
    transtimes = [(x - halfdur, x + halfdur) for x in transtimesM]
    
    return transtimes

def altaz_at_time(toi, time, location):
    '''Gives the altaz object of the TOI at the given time
    
    Args:
        toi (TOI): the TOI object
        time (Time): Time object
        location (EarthLocation): location of observatory
    
    Returns:
        altaz (AltAz): AltAz object
    '''
    coords = SkyCoord(toi.ra, toi.dec, unit='deg')
    altaz = coords.transform_to(AltAz(obstime=time, location=location))
    return altaz

def calculate_transit_alts(toi, location, times):
    '''Calculates the altitude and azimuth at ingress and egress for a TOI
    
    Args: 
        toi (TOI): the TOI to calculate for 
        times (list(tuple(float, float))): list of ingress and egress times
        location (EarthLocation): location of observatory 
        
    Returns: 
        altazs (list(tuple(float, float), tuple(float, float))): alt az list
    '''
    altazs = []
    for trans in times: 
        altazs.append((altaz_at_time(toi, trans[0], location), altaz_at_time(toi, trans[1], location)))
    return altazs

def decimal_to_fraction(decimal):
    '''Takes a decimal and converts it to the best matching mixed
    fraction. This is useful for labelling which period fraction a given
    possible transit is, and is used in multiple places. Works down to 
    1/20th fraction. 
    
    Args:
        decimal (float): decimal 
    
    Returns: 
        fraction (str): mixed number representation
    '''
    neg = ''
    if decimal < 0:
        neg = '-'
    absdec = abs(decimal)
    fraction = Fraction(absdec).limit_denominator(20)
    mixed_fraction = f'{neg}{fraction.numerator // fraction.denominator} {fraction.numerator % fraction.denominator}/{fraction.denominator}'
    return mixed_fraction

def generate_all_transits(toi_obj, location, timerange, utcoffset, minfrac, minalt):
    '''Function used by the find_all_transits shell in order to generate a
    dataframe of all the transists for a specific toi in the parameters.
    
    Args: 
        toi_obj (TOI): the TOI object to calculate for
        location (EarthLocation): location of observatory
        timerange (Tuple(int, int)): tuple of start and end day w.r.t current time
        utcoffset (int): hours difference from UTC time, default to PST
        minfrac (int): fraction of given period to look for, defualt to 1
        minalt (float): minimum altitude of ingress or egrees above horizon
    
    Returns:
        df (DataFrame): dataframe containing transits for given TOI
    '''
    transit_times = calculate_transit_times(toi_obj, minfrac, timerange, mode='ends')
    transit_altaz = calculate_transit_alts(toi_obj, location, transit_times)
    transit_count = len(transit_times)
    
    df = pd.DataFrame()
    df['TOI ID'] = [toi_obj.name] * transit_count
    df['Time Ingress JD'] = [f'{float(x[0].value):.4f}' for x in transit_times]
    
    itimes = [Time((x[0] + (utcoffset*u.hour)), format='fits').value for x in transit_times]
    itimes2 = [x.split('T') for x in itimes]
    ingress_times = [" ".join(x) for x in itimes2]
    df['Time Ingress Local'] = ingress_times
    
    etimes = [Time((x[1] + (utcoffset*u.hour)), format='fits').value for x in transit_times]
    etimes2 = [x.split('T') for x in etimes]
    egress_times = [" ".join(x) for x in etimes2]
    df['Time Egress Local'] = egress_times
    
    df['Azimuth Ingress'] = [float(f'{(x[0].az.degree):.4f}') for x in transit_altaz]
    df['Altitude Ingress'] = [float(f'{(x[0].alt.degree):.4f}') for x in transit_altaz]
    df['Altitude Egress'] = [float(f'{(x[1].alt.degree):.4f}') for x in transit_altaz]
    
    def sun18(time):
        '''Returns True if the sun is below 18 degrees, making it astro dark
        at the specified time which must be given as a Time object
        
        Args:
            time (Time): time object of when to calculate for
        
        Returns: 
            bool: True if dark, False if day
        '''
        altaz = AltAz(obstime=time, location=SEO)
        alt = float(get_sun(Time.now()).transform_to(altaz).alt.deg)
        return alt < -18
    
    nights1 = []
    for x in transit_times:
        nights1.append(sun18(x[0]))
    df['Nighttime Ingress?'] = nights1

    nights2 = []
    for x in transit_times:
        nights2.append(sun18(x[1])) 
    df['Nighttime Egress?'] = nights2
    
    df['Depth (ppt)'] = [toi_obj.depth / 1000] * transit_count
    df['TMag'] = [toi_obj.tmag] * transit_count
    df['Duration'] = [round(toi_obj.duration, 3)] * transit_count
    midtrans = [x[0].value + ((toi_obj.duration) / 24) for x in transit_times]
    df['Period Fraction'] = [decimal_to_fraction((x - (toi_obj.transit0 + TESST)) / toi_obj.period) for x in midtrans]
    df['Period'] = [toi_obj.period] * transit_count
    df['Comments'] = [toi_obj.comment] * transit_count
    
    # sort by altitude of ingress and egress
    df1 = df[df['Altitude Egress'] > 30]
    df2 = df[df['Altitude Ingress'] > 30]
    df3 = pd.concat([df1, df2], ignore_index=True).drop_duplicates()

    # sort by at least one gress during night
    df4 = df3[df3['Nighttime Ingress?']]
    df5 = df3[df3['Nighttime Egress?']]
    df6 = pd.concat([df4, df5], ignore_index=True).drop_duplicates()
    
    return df6

def find_all_transits(location=SEO, 
                      data = DATA,
                      timerange = (2460000, 2460100),  
                      utcoffset = -8,
                      minfrac = 1,
                      periods = (80, 1000),
                      mindepth = 0, 
                      minalt = 25):
    '''Finds all transits of TOIs that are in the datasheet with 
    the given parameters and returns and displays a datatable of them. 
        
    Args:
        location (EarthLocation): location of observatory
        timerange (tuple(int, int)): tuple of start and end in JD
        data (DataFrame): full table loaded from csv
        utcoffset (int): hours difference from UTC time, default to PST
        minfrac (int): fraction of given period to look for, defualt to 1
        periods (tuple(float, float)): range of periods in days to look for
        mindepth (float): minimum transit depth to include in result
        minalt (float): minimum altitude of ingress or egrees above horizon
    
    Returns:
        mainframe (DataFrame): table with info on transits in the parameters
    '''
    
    # restrict TOIs to search for by transit depth and orbital period
    pmin, pmax = periods
    dmin = 1000 * mindepth
    selected = data[(data['Orbital Period Value'] > pmin) 
                    & (data['Orbital Period Value'] < pmax) 
                    & (data['Transit Depth Value'] > dmin)]
    selected.reset_index(drop=True, inplace=True)
    
    # create mainframe that will be returned 
    mainframe = pd.DataFrame(columns=['TOI ID', 'Time Ingress JD', 'Time Ingress Local', 'Time Egress Local',
                                    'Azimuth Ingress', 'Altitude Ingress', 'Altitude Egress', 
                                    'Nighttime Ingress?', 'Nighttime Egress?', 'Depth (ppt)', 'TMag', 'Duration', 
                                    'Period Fraction', 'Period', 'Comments'])
    
    # for each TOI fitting the parameters, find its transits and add them to mainframe
    for row in range(len(selected)):
        toi_obj = TOI(selected.iloc[row]['Full TOI ID'], selected)
        toidata = generate_all_transits(toi_obj, location, timerange, utcoffset, minfrac, minalt)
        mainframe = pd.concat([mainframe, toidata], ignore_index=True)
        
    mainframe.sort_values('Time Ingress JD', inplace=True)
    mainframe.reset_index(drop=True, inplace=True)
    
    return mainframe


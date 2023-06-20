# file for code

# Import and Style Statements
import matplotlib
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import math
import warnings
import datetime
import time
import pandas as pd
import os
from IPython.display import display, Markdown
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
import lightkurve as lk
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Spacer, Table
from reportlab.platypus import TableStyle, Paragraph, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from astropy.visualization import astropy_mpl_style, quantity_support
import io
from reportlab.lib.units import inch
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction
plt.style.use(astropy_mpl_style)
matplotlib.rcParams["axes.formatter.useoffset"] = False
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="More than 20 figures have been opened.*")

quantity_support()
PAGE_SIZE = landscape(letter)

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
    
    def __init__(self, toi, datatable='none'):
        '''Pull information on a given TOI from the database

        Args:
            toi (float, str): TOI ID for the target star
            datatable (DataFrame): the full TESS df, or None and load it in
        '''
        if isinstance(datatable, str):
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
        self.sectors = toirow['Sectors']
        self.infodict = {'transit0':self.transit0, 'RA':self.ra, 'DEC': self.dec, 
                         'period':self.period, 'duration': self.duration, 'TIC':self.tic,
                         'depth':self.depth, 'tmag':self.tmag, 'sectors':self.sectors,
                         'comment':self.comment}

    def __repr__(self):
        return f'TOI {self.name}: {self.infodict}'
    
    def __str__(self):
        return f'TOI {self.name}'
    
def toi_info(toi):
    '''
    '''
    toi = TOI(toi)
    print(f"Recorded info in TOI database for TOI {toi.name}: ")
    for x in toi.infodict:
        print(f'  {x} : {toi.infodict[x]}')
    return

# functions used for finding all upcoming TESS transits
def earth_location(loc):
    '''Takes an inputted location coordinate set from exec
    and returns an actual EarthLocation object

    Args:
        loc (list(float, float, float)): inputted location
    
    Returns:
        eloc (EarthLocation): corresponding EarthLocation object
    '''
    return EarthLocation(lat=loc[0]*u.deg, lon=loc[1]*u.deg, height=loc[2]*u.m)

def calculate_transit_times(toi, 
                            fraction=1, 
                            timerange=(2457000, 2461000), 
                            mode='mid'):
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

def altaz_at_time(toi, 
                  time, 
                  location):
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

def calculate_transit_alts(toi, 
                           location, 
                           times):
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
        altazs.append((altaz_at_time(toi, trans[0], location), 
                       altaz_at_time(toi, trans[1], location)))
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
    mixed_fraction = f'{neg}{fraction.numerator // fraction.denominator} \
    {fraction.numerator % fraction.denominator}/{fraction.denominator}'
    if (mixed_fraction[-1]=='1') and (mixed_fraction[-2]=='/'):
        return mixed_fraction[:-4]
    return mixed_fraction

def generate_all_transits(toi_obj, 
                          location, 
                          timerange, 
                          utcoffset, minfrac, minalt):
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
    
    itimes = [Time((x[0] + (utcoffset*u.hour)), 
                   format='fits').value for x in transit_times]
    itimes2 = [x.split('T') for x in itimes]
    ingress_times = [" ".join(x) for x in itimes2]
    df['Time Ingress Local'] = ingress_times
    
    etimes = [Time((x[1] + (utcoffset*u.hour)), 
                   format='fits').value for x in transit_times]
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
    df['Period Fraction'] = [decimal_to_fraction((x - (toi_obj.transit0 + TESST)) /\
                                                  toi_obj.period) for x in midtrans]
    df['Period'] = [toi_obj.period] * transit_count
    df['Comments'] = [toi_obj.comment] * transit_count
    
    # sort by altitude of ingress and egress
    df1 = df[df['Altitude Egress'] > minalt]
    df2 = df[df['Altitude Ingress'] > minalt]
    df3 = pd.concat([df1, df2], ignore_index=True).drop_duplicates()

    # sort by at least one gress during night
    df4 = df3[df3['Nighttime Ingress?']]
    df5 = df3[df3['Nighttime Egress?']]
    df6 = pd.concat([df4, df5], ignore_index=True).drop_duplicates()
    
    return df6

def save_df_to_pdf(df):
    '''Saves the large dataframe with all the possible transits to a 
    nicer looking pdf 
    
    Args:
        df (DataFrame): the df of all the transits (or any df really)
    
    Returns:
        writes pdf "transits.pdf" containing the styled dataframe
    '''
    pagesz = landscape((2000, 600))
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=pagesz)
    doc.topMargin = .1 * inch
    doc.bottomMargin = .1 * inch
    
    # Create a list to hold all the elements in the PDF
    elements = []
    
    # Add title to PDF
    title_style = getSampleStyleSheet()['Title']
    title_paragraph = Paragraph('Predicted Transits', title_style)
    elements.append(title_paragraph)
    elements.append(Spacer(1, inch * 0.25))
    
    # Add pandas dataframe to PDF
    data = [df.columns[:,].tolist()] + df.values.tolist()
    table = Table(data)
    table.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.grey), 
                               ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
                               ('ALIGN', (0,0), (-1,-1), 'CENTER'), 
                               ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                               ('FONTSIZE', (0,0), (-1,0), 14), 
                               ('BOTTOMPADDING', (0,0), (-1,0), 12),
                               ('BACKGROUND', (0,1), (-1,-1), colors.beige), 
                               ('TEXTCOLOR', (0,1), (-1,-1), colors.black),
                               ('FONTNAME', (0,1), (-1,-1), 'Helvetica'), 
                               ('FONTSIZE', (0,1), (-1,-1), 10),
                               ('ALIGN', (0,1), (-1,-1), 'LEFT'), 
                               ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
                               ('GRID', (0,0), (-1,-1), 1, colors.black)]))
    
    elements.append(table)

    # Build PDF
    doc.build(elements)

    # Save PDF from buffer
    pdf = buffer.getvalue()
    buffer.close()
    with open("transits.pdf", 'wb') as f:
        f.write(pdf)
    return

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
    mainframe = pd.DataFrame(columns=['TOI ID', 'Time Ingress JD', 'Time Ingress Local', 
                                      'Time Egress Local','Azimuth Ingress', 
                                      'Altitude Ingress', 'Altitude Egress', 
                                      'Nighttime Ingress?', 'Nighttime Egress?', 
                                      'Depth (ppt)', 'TMag', 'Duration', 
                                      'Period Fraction', 'Period', 'Comments'])
    
    # for each TOI fitting the params, find its transits and add them to mainframe
    for row in range(len(selected)):
        toi_obj = TOI(selected.iloc[row]['Full TOI ID'], selected)
        toidata = generate_all_transits(toi_obj, location, timerange,
                                         utcoffset, minfrac, minalt)
        mainframe = pd.concat([mainframe, toidata], ignore_index=True)
        
    mainframe.sort_values('Time Ingress JD', inplace=True)
    mainframe.reset_index(drop=True, inplace=True)
    
    display(mainframe)
    mainframe.to_csv('transits.csv', index=False)
    save_df_to_pdf(mainframe)
    return mainframe

# function used for plotting a specific transit
def toi_plot_transit(toi, 
                     ingress_time, 
                     utcoffset=-7, 
                     location=SEO):
    '''Makes a pretty plot of the target location and the sun and moon
    locations over the course of a given tranist.
    
    Args:
        toi (TOI, float): either a TOI object or a float TOI ID
        ingress_time (float): JD time of ingress
        location (EarthLocation): location of observatory
        utcoffset (int): hour difference from UTC
        
    Returns:
        displays a graph and saves it to 'transplot.png'
    '''
    if isinstance(toi, float):
        toi = TOI(toi)
    
    # create star line
    duration = toi.duration
    ra = toi.ra
    dec = toi.dec
    tjd = Time(ingress_time, format='jd')
    coords = SkyCoord(ra, dec, unit="deg")
    delta = np.linspace(-5, duration+5, 1000)*u.hour
    times = delta + tjd
    toi_altazs = coords.transform_to(AltAz(obstime=times, location=location))
    plt.scatter(delta, toi_altazs.alt, label=str(toi), c=toi_altazs.az, lw=0, s=8, cmap='viridis')
    
    # create line for moon
    frame = AltAz(obstime=times, location=location)
    moon_locs = get_moon(times)
    moon_altaz = moon_locs.transform_to(frame)
    plt.plot(delta, moon_altaz.alt, color=[0.75]*3, ls='--', label='Moon')
    
    # create line for sun
    frame = AltAz(obstime=times, location=location)
    sun_locs = get_sun(times)
    sun_altaz = sun_locs.transform_to(frame)
    plt.plot(delta, sun_altaz.alt, color='r', ls='-', label='Sun')
    
    # format plot attributes and save and show plot
    plt.colorbar().set_label('Target Azimuth [deg]')
    plt.axvspan(0, duration, alpha=.25, color='grey')
    local = tjd + (utcoffset*u.hour)
    local.format = 'fits'
    local2 = str(local).split('T')
    local3 = ' at '.join(local2)
    plt.xlabel(f'Hours since {local3} Ingress')
    plt.ylabel('Altitude [deg]')
    plt.title(f'Altitude of {toi}')
    plt.legend(loc='upper left')
    plt.savefig('transplot.png')
    plt.show()
    return 

# functions for retrieiving and displaying light curves
def save_figures_to_pdf(figures, 
                        toi):
    '''Takes a list of figures and saves them to a pdf on top of the other
    
    Args:
        figures (list(Figure)): list of figures
        toi (str): toi name for titling purposes
    
    Returns:
        creates an output pdf with the figure pngs in it
    '''
    # Create buffer to hold PDF
    buffer = io.BytesIO()

    # Create PDF canvas using buffer
    doc = SimpleDocTemplate(buffer, pagesize=PAGE_SIZE)
    
    doc.topMargin = .1 * inch
    doc.bottomMargin = .1 * inch
    
    # Create a list to hold all the elements in the PDF
    elements = []

    # Add title to PDF
    title_style = getSampleStyleSheet()['Title']
    title_paragraph = Paragraph(f'TESS Lightcurves for TOI {toi}', title_style)
    elements.append(title_paragraph)
    elements.append(Spacer(1, inch * 0.25))
    
    # Add figures to PDF
    for i, fig in enumerate(figures):
        # Save figure to buffer as PNG
        img_buffer = io.BytesIO()
        fig.savefig(img_buffer, format='png')
        plt.close(fig)

        # Add PNG to PDF
        img = Image(img_buffer)
        img.drawHeight = 100
        img.drawWidth = 500
        elements.append(img)
        elements.append(Spacer(1, inch * 0.25))
        
    # Build PDF
    doc.build(elements)

    # Save PDF from buffer
    pdf = buffer.getvalue()
    buffer.close()
    with open("lightcurves.pdf", 'wb') as f:
        f.write(pdf)
    return

def toi_extract_all_curves(toi, 
                           fraction=1):
    '''Find all the light curves in the TESS database for the given TOI
    that correspond to a transit time, or possible transit time based on
    the fraction given. Download the products and arrange into a dict 
    that contains the necessary info for plotting.
    
    Args:
        toi (TOI, float): a TOI object or a float TOI ID
        fraction (int): largest fraction to divide period into
    
    Returns:
        tess_data_dict (dict): dictionary containing 
    '''
    if isinstance(toi, float):
        toi = TOI(toi)
    
    # initialize variables
    sectors = [int(x) for x in toi.sectors[1:-1].split(',')]
    TIC = toi.tic
    epoch = toi.transit0
    period = toi.period
    duration = toi.duration / 24
    depth = toi.depth / 10000
    
    # find available light curve data
    dataprods = lk.search_lightcurve(f'TIC {TIC}')
    
    # find all transit times
    all_transits = calculate_transit_times(toi, fraction, timerange=(TESST, TESST + 4000), mode='mid')
    
    # extract applicable products into dictionary with entry for each transit time
    tess_data_dict = {}
    for i, prod in enumerate(dataprods):
        try:
            product = prod.download()
        except:
            continue
        times = product['time'].value
        start = times[0]
        end = times[-1]
        fluxes = product['flux'].value
        author = prod.author.data[0]
        sector = prod.mission[0]
        info2 = {'author': author, 'sector': sector, 'duration':duration, 
                 'depth':depth, 'epoch':epoch, 'period':period}
        for trans in all_transits:
            if start < trans < end:
                iden = (trans)
                ingr = trans - (duration * .5) - .5
                egr = trans + (duration * .5) + .5
                
                ingr_ind = 0
                for time in times:
                    if time > ingr:
                        break
                    else:
                        ingr_ind += 1   
                egr_ind = 0
                for time in times:
                    if time > egr:
                        break
                    else:
                        egr_ind += 1
                
                time_range = times[ingr_ind:egr_ind]
                flux_range = fluxes[ingr_ind:egr_ind]
                if iden not in tess_data_dict:
                    tess_data_dict[iden] = []
                tess_data_dict[iden].append((time_range, flux_range, info2))
    return tess_data_dict

def toi_plot_curves(toi, 
                    fraction=1):
    '''Plots all the light curves for a given TOI at the given fraction
    by downloading them with toi_extract_all_curves. Also saves each fig
    to a pdf containing all transits for this TOI.
    
    Args:
        toi (TOI, float): a TOI object or a float TOI ID
        fraction (int): largest fraction to divide period into
    
    Returns:
        creates an output pdf with each figure using auxiliary function
    '''
    if isinstance(toi, float):
        toi = TOI(toi)
    epoch = toi.transit0
    period = toi.period
    duration = toi.duration / 24
    depth = toi.depth / 10000
    
    tess_data_dict = toi_extract_all_curves(toi, fraction)
    figures = []
    for key in tess_data_dict:
        
        datasets = tess_data_dict[key]
        curve_count = len(datasets)
        
        fig = plt.figure(figsize=(10, 2))
        rows = 1
        columns = curve_count
        
        for i, plot in enumerate(datasets):
            ax = fig.add_subplot(rows, columns, i + 1)
            y = np.array(plot[1]) / np.median(plot[1])
            #plt.plot(plot[0], plot[1], "Black", linewidth=.75, alpha=.85)
            plt.plot(plot[0], y, "Black", linewidth=.75, alpha=.85)
            plt.title(plot[2]['author'], size=12)
            #ax.set_yticklabels([])
            plt.tight_layout()
            plt.yticks(fontsize=8)
            plt.xticks(fontsize=8)
            plt.axvline(x=float(key), ls='--', color='red', alpha=0.5)
            ax.axvspan((float(key) - (duration / 2)), (float(key) + (duration / 2)), 
                       alpha=0.25, color='grey')
            if i == 0:
                plt.ylabel('Flux', size=12)
       
        # generate title and show plots
        sector_id = datasets[0][2]['sector']
        n1 = ((key - epoch) / period)
        n = decimal_to_fraction(n1)
        title = f'{sector_id} (E = {n})'
        fig.suptitle(title, size=10, y=1, fontweight='bold')
        plt.text(.5, -.1, 'Time (BJD-2457000)', transform=fig.transFigure, 
                 horizontalalignment='center')
        figures.append(fig)
        plt.close('all')
        #plt.show()
        #print('')
        
    save_figures_to_pdf(figures, toi.name)
    return

# functions for analyzing a TOI and all its curves
def least_square(fluxes_r, 
                 model_fluxes):
    '''Calculate the chi squared value between the points and the model points
    
    Args:
        fluxes_r (array): array of data fluxes
        model_fluxes (array): array of model fluxes
        
    Returns
        float: chi squared test statistic
    '''
    fluxes_r = np.nan_to_num(fluxes_r, nan=1)
    return float(sum((fluxes_r - model_fluxes)**2))

def test_least_square(midtransit, 
                      args): 
    '''Helper function to calculate the model flux points
    
    Args:
        midtransit (float): middle of transit time
        args (list): list of transit data time, flux, duration, depth
        
    Returns:
        calls least square function, returns float chi square statistic
    '''
    # args = [times, fluxes, duration, depth]
    
    # extract transit data constants
    times, fluxes, duration, depth = args[0], args[1], args[2], args[3]
    
    # calculate ingress and egress times for the given midtransit point
    ing, egr = (midtransit - (duration / 2)), (midtransit + (duration / 2))
    
    # initialize an array of model fluxes the same size as actual fluxes
    model_fluxes = np.zeros_like(fluxes)
    
    # loop through times and create model 
    for i, t in enumerate(times):
        
        # if time is during transit, value is (1 - depth), 1 otherwise
        if t < ing:
            model_fluxes[i] = 1
        elif egr > t > ing:
            drop = 1 - (depth / 100) 
            model_fluxes[i] = drop
        elif t > egr:
            model_fluxes[i] = 1
    
    # check the whole model got filled in 
    assert np.min(model_fluxes) > 0 
    
    return least_square(fluxes, model_fluxes)

def transit_time_fit(times, 
                     fluxes, 
                     duration, 
                     depth, 
                     n=1000):
    '''Uses statistic tests to find the best fit transit midtime 
    
    Args:
        times (array): time data
        fluxes (array): flux data
        duration (float): time length of the transit
        depth (float): depth of the transit
        n (int): maximum precision of fit
    
    Returns:
        minimum_time (float): the calculated transit midpoint
    '''
    # create constant args [times, norm_flux, duration, depth]
    norm_flux = fluxes / np.median(fluxes)
    args = [times, norm_flux, duration, depth]
    
    # find start and end of data range
    if len(times) == 0:
        return np.nan
    start, end = times[0], times[-1]

    # range of possible mid transit times in the data range
    possible_transits = np.linspace(start, end, n)
    
    chi2s = []
    
    # initialize minimum values
    minimum_time = 0 
    minimum_chi2 = math.inf
    
    # try each possible choice
    for mid in possible_transits:
        
        # calculate the chi2 statistic for given midtransit time
        chi2 = test_least_square(mid, args)
        
        # if this value fits better than previous best, update the best fit
        if chi2 < minimum_chi2:
            minimum_time = mid
            minimum_chi2 = chi2
        
        #chi2s.append(chi2)
    
    return minimum_time#, chi2s

def transit_depth_fit(times, 
                      fluxes, 
                      midtransit, 
                      duration, 
                      depth, 
                      n=1000, 
                      sideperc=.3):
    '''Calculate the transit depth with the best fit curve
    
    Args:
        times (array): time data
        fluxes (array): flux data
        midtransit (float): calculated mid transit time
        duration (float): transit duration
        depth (float): given transit depth
        n (int): fineness of fit, not used... 
        
    Return:
        float calculated depth value
        sigmas: i forget what this is atm...
    '''
    # create and select in transit array from data
    #norm_flux = fluxes / np.median(fluxes)
    
    ing, egr = (midtransit - (duration / 2)), (midtransit + (duration / 2))
    baseline = []
    for i, time in enumerate(times):
        if (time < ing) or (time > egr):
            baseline.append(fluxes[i])
    median = np.median(np.array(baseline))
    norm_flux = fluxes / median
    
    intransit_full = []
    for i, time in enumerate(times):
        if ing < time < egr:
            intransit_full.append(norm_flux[i])
    
    # estimate # of points in egress and ingress to cut off, for now, just cut max 20% per side
    intr_count = len(intransit_full)
    perc_10 = int(intr_count * sideperc)
    intransit = intransit_full[perc_10:intr_count - perc_10]
    
    # determine if depth is significant 
    norm_base = baseline / np.median(baseline)
    base_mean = np.mean(norm_base)
    base_std = np.std(norm_base)
    diff = base_mean - np.mean(intransit)
    sigmas = diff / base_std
              #((1 - minimum_depth) * 100)
    return (1 - np.median(intransit)) * 100, sigmas

def plot_fit(times, 
             fluxes, 
             midtransit, 
             depth, 
             duration):
    '''
    '''
    #norm_fluxes = fluxes / np.median(fluxes)
    
    # calculate ingress and egress times for the given midtransit point
    ing, egr = (midtransit - (duration / 2)), (midtransit + (duration / 2))
    
    
    baseline = []
    for i, time in enumerate(times):
        if (time < ing) or (time > egr):
            baseline.append(fluxes[i])
    median = np.median(np.array(baseline))
    norm_flux = fluxes / median
    
    
    # initialize an array of model fluxes the same size as actual fluxes
    model_fluxes = np.zeros_like(fluxes)
    
    # loop through times and create model 
    for i, t in enumerate(times):
        
        # if time is during transit, value is (1 - depth), 1 otherwise
        if t < ing:
            model_fluxes[i] = 1
        elif egr > t > ing:
            drop = 1 - (depth / 100) 
            model_fluxes[i] = drop
        elif t > egr:
            model_fluxes[i] = 1
    
    
    # graph data times and fluxes
    plt.scatter(times, norm_flux, s=6)
    plt.plot(times, model_fluxes, c='r', lw=10, alpha=.5)
    plt.xlabel("TESS Time BJD-2457000")
    plt.ylabel("Median Normalized Flux")
    plt.title("Transit Fit")
    return

def transit_analysis(times, 
                     fluxes, 
                     duration, 
                     depth, 
                     show=False):
    '''Function that does initial analysis of a single light curve and 
    tries to determine if it contains a transit by fitting the depth and
    the time, and doing a few extra checks. 
    
    Args:
        times (array): time data
        fluxes (array): flux data
        duration (float): transit duration
        depth (float): given transit depth
        show (bool): obsolete now I believe
        
    Returns:
        information on the transit fit
    '''
    # replace any nans in the data
    fluxes_1 = fluxes
    median_val = np.nanmedian(fluxes)
    fluxes_1[np.isnan(fluxes)] = median_val
    broken = False
    
    ttime = transit_time_fit(times, fluxes_1, duration, depth)
    tdepth, sigma = transit_depth_fit(times, fluxes_1, ttime, duration, depth)
    
    transit = True
    # ADD MORE CONSTRAINTS HERE

    if math.isnan(ttime) or math.isnan(tdepth):
        if show:
            print("Not likely to be a transit!")
        transit = False
    
    if tdepth < .05:
        if show:
            print("Not likely to be a transit!")
        transit = False
        
    if abs(depth - tdepth) > .2:
        if show:
            print("Not likely to be a transit!")
        transit = False
    
    if len(times) > 0:
        if (abs(times[0] - ttime) < .1) or (abs(times[-1] - ttime) < .1):
            if show:
                print("Not likely to be a transit!")
            transit = False
    
    if show:
        print(f"Best fit transit with {sigma:.3f} sigma identified at \
              {round(ttime, 3)} days with depth of {round(tdepth, 3)}%")
        plot_fit(times, fluxes, ttime, tdepth, duration)

    return ttime, tdepth, sigma, transit

def plot_toi_fits(datalist): 
    '''Given the list of light curves and info for a TOI, does the analysis
    of each transit and prepares a datastructure to be reduced to a pdf.
    
    Args:
        datalist (dict): dictionary contianing all the information to be processed
    
    Returns:
        figure_list (list(Fig)): list of generated figures for each light curve
        df (DataFrame): dataframe containing information on each transit
    '''
    sectors = ['----']
    fractions = ['----']
    depths = []
    sigmas = ['----']
    transited = ['----']
    ttv = ['----']
    
    figure_list = []
    
    for transtime in datalist.keys():
        dataset = datalist[transtime]
        if len(depths) == 0:
            depths.append(dataset[0][2]['depth'])
        
        fig = plt.figure(figsize=(10, 2))
        rows = 1
        columns = len(dataset)
        
        period = 0
        epoch = 0
        ttimed = []
        tdepths = []
        tsigmas = []
        transits = []
        ttvs = []
        
        
        expected_mid = transtime
        single = len(dataset)
        for i, trans in enumerate(dataset):
            epoch = trans[2]['epoch']
            period = trans[2]['period']

            ax = fig.add_subplot(rows, columns, i + 1)
            ttime, tdepth, sigma, transit = transit_analysis(trans[0], trans[1], 
                                                             trans[2]['duration'], 
                                                             trans[2]['depth'], 
                                                             False)
            
            ttimed.append(ttime)
            tsigmas.append(sigma)
            tdepths.append(tdepth)
            ttvs.append(ttime - expected_mid)
            if transit:
                transits.append(1)
            else:
                transits.append(0)
            
            
            ing, egr = (ttime - (trans[2]['duration'] / 2)), \
                       (ttime + (trans[2]['duration'] / 2))
    
            baseline = []
            for j, time in enumerate(trans[0]):
                if (time < ing) or (time > egr):
                    baseline.append(trans[1][j])
            median = np.median(np.array(baseline))
            norm_flux = trans[1] / median


            # initialize an array of model fluxes the same size as actual fluxes
            model_fluxes = np.zeros_like(trans[1])

            # loop through times and create model 
            for k, t in enumerate(trans[0]):

                # if time is during transit, value is (1 - depth), 1 otherwise
                if t < ing:
                    model_fluxes[k] = 1
                elif egr > t > ing:
                    drop = 1 - (tdepth / 100) 
                    model_fluxes[k] = drop
                elif t > egr:
                    model_fluxes[k] = 1


            # graph data times and fluxes
            plt.scatter(trans[0], norm_flux, s=6)
            if transit:
                lc = 'g'
            else:
                lc = 'r'
            plt.plot(trans[0], model_fluxes, c=lc, lw=5, alpha=.5)
            plt.ylim(((1 - (depths[0] / 100) * 2.5), (1 +((depths[0] / 100) * 1.5))))
            plt.axvline(x=float(expected_mid), ls='--', color='red', alpha=0.5)
            if i != 0:
                ax.set_yticklabels([])
            plt.tight_layout()
            plt.yticks(fontsize=8)
            plt.xticks(fontsize=8)
            if i == 0:
                plt.ylabel('Flux', size=12)
            if columns >= 6:
                fsize = 5
            elif columns >= 4:
                fsize = 6
            else:
                fsize = 8
            if single != 1:
                plt.title(f"{trans[2]['author']} - T: {ttime:.3f}, D: {tdepth:.3f}%, \u03C3: {sigma:.3f}", fontsize=fsize)
        
        # ADD CODE HERE TO THROW OUT BAD FITS from lists
        # remove from: ttimed, tsigmas, ttvs, tdepths
        def sigma_calc(*args):
            csigmas = []
            for lst in args:
                csigmas.append(np.std(lst))
            return csigmas
        
        # protect against one bad data analysis missing a transit
        if (len(transits) - sum(transits)) == 1:
            idx = np.argmin(transits)
            ttimed.remove(ttimed[idx])
            tsigmas.remove(tsigmas[idx])
            ttvs.remove(ttvs[idx])
            tdepths.remove(tdepths[idx])
        
        # protect against a false positive transit skewing good transit results
        if len(ttvs) > 1:
            if ((max(ttvs) > .3) and (min(ttvs) < .3)) or (np.std(ttvs) > .1):
                idx = np.argmax(np.abs(ttvs))
                ttimed.remove(ttimed[idx])
                tsigmas.remove(tsigmas[idx])
                ttvs.remove(ttvs[idx])
                tdepths.remove(tdepths[idx])
                transits[idx] = False
        
        # protect against a nan result somehow evaluating as True
        for idx, dep in enumerate(tdepths): 
            if math.isnan(dep):
                idx = np.argmax(np.abs(ttvs))
                ttimed.remove(ttimed[idx])
                tsigmas.remove(tsigmas[idx])
                ttvs.remove(ttvs[idx])
                tdepths.remove(tdepths[idx])
                transits[idx] = False
                


        ntswitch = False
        positive = sum(transits) / len(transits)
        if positive >= .75:
            yn = ' : Likely Transit : '           
        elif positive >= .5:
            yn = ' : Possible Transit : '
        elif positive >= .25:
            yn = ' : Unlikely Transit : '           
        else:
            yn = ' : No Transit : '
            ntswitch = True
        
        tstr = str(len(transits)) + yn + str(round(positive, 2))
        transited.append(tstr)
        # for adding the fraction to the graph title
        n1time = expected_mid
        n1 = ((n1time - epoch) / period)
        n = decimal_to_fraction(n1)
        
        sectors.append(dataset[0][2]['sector'])
        fractions.append(n)
        
        if not ntswitch:
            if len(tdepths) > 0:
                mean_depth = np.mean(tdepths)
            else:
                mean_depth = 0
            depths.append(mean_depth)

            if len(tsigmas) > 0:
                mean_sigma = np.mean(tsigmas)
            else:
                mean_sigma = 0
            sigmas.append(mean_sigma)

            if len(ttvs) > 0:
                mean_ttv = np.mean(ttvs)
            else:
                mean_ttv = 0                 
            ttv.append(mean_ttv)
        else:
            depths.append(0)
            sigmas.append('----')
            ttv.append('----')
        
        title = str(dataset[0][2]['sector']) + yn + f'(E = {n})'
        if single != 1:
            fig.suptitle(title, size=10, y=1, fontweight='bold')
        if single == 1:
            title2 = (f"  {trans[2]['author']} - T: {ttime:.3f}, D: {tdepth:.3f}%, \u03C3: {sigma:.3f}")
            fig.suptitle(title + title2, size=10, y=1, fontweight='bold')
            
        plt.text(.5, .01, 'Time (BJD-2457000)', transform=fig.transFigure, 
                 horizontalalignment='center')
        #plt.show()
        #print('')
        
        if ((yn == ' : Likely Transit : ') or 
            (yn == ' : Possible Transit : ') or 
            ('/' not in n)):
            figure_list.append(fig)

    display_data = {"TESS Sector": sectors, 
    "Epoch": fractions,
    "# : Transit? : %": transited,
    "Depth (%)": list(np.array(depths) * 1),
    "Sigma (\u03C3)": sigmas,
    "TTV (Days)": ttv}

    df = pd.DataFrame(display_data)
    TOI_ID = "Data Summary"
    #display(Markdown(f"<h2 style='text-align:center';>{TOI_ID}</h2>"))
    #display(df)
    #print('')
    plt.close('all')

    return figure_list, df

def create_pdf(figures, 
               df, 
               title, 
               subtitle):
    '''Creates the pdf output for main TOI analysis given a list of 
    light curve graphs, a dataframe with information about them, and
    two title strings.
    
    Args:
        figures (list(Fig)): list of graphs to be shown
        df (DataFrame): dataframe to be shown
        title (str): main title
        subtitle (str): subtitle
    
    Returns: 
        pdf (Buffer): contains the pdf
    '''
    # Create buffer to hold PDF
    buffer = io.BytesIO()

    # Create PDF canvas using buffer
    pgsize = landscape(letter)
    doc = SimpleDocTemplate(buffer, pagesize=pgsize)
    
    doc.topMargin = .1 * inch
    doc.bottomMargin = .1 * inch
    
    # Create a list to hold all the elements in the PDF
    elements = []

    # Add title to PDF
    title_style = getSampleStyleSheet()['Title']
    title_paragraph = Paragraph(title, title_style)
    elements.append(title_paragraph)
    
    stitle = ParagraphStyle(name='stitle', fontSize=10, leading=12, 
                            alignment=1, spaceAfter=12,)
    elements.append(Paragraph(subtitle, stitle))
    
    elements.append(Spacer(1, inch * 0.25))
    
    
    # Add figures to PDF
    for i, fig in enumerate(figures):
        # Save figure to buffer as PNG
        img_buffer = io.BytesIO()
        fig.savefig(img_buffer, format='png')
        plt.close(fig)

        # Add PNG to PDF
        img = Image(img_buffer)
        img.drawHeight = 100
        img.drawWidth = 500
        elements.append(img)
        elements.append(Spacer(1, inch * 0.25))
        
        

    # Add pandas dataframe to PDF
    data = [df.columns[:,].tolist()] + df.values.tolist()
    table = Table(data)
    table.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.grey), 
                               ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
                               ('ALIGN', (0,0), (-1,-1), 'CENTER'), 
                               ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                               ('FONTSIZE', (0,0), (-1,0), 14), 
                               ('BOTTOMPADDING', (0,0), (-1,0), 12),
                               ('BACKGROUND', (0,1), (-1,-1), colors.beige), 
                               ('TEXTCOLOR', (0,1), (-1,-1), colors.black),
                               ('FONTNAME', (0,1), (-1,-1), 'Helvetica'), 
                               ('FONTSIZE', (0,1), (-1,-1), 10),
                               ('ALIGN', (0,1), (-1,-1), 'LEFT'), 
                               ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
                               ('GRID', (0,0), (-1,-1), 1, colors.black)]))
    elements.append(table)

    # Build PDF
    doc.build(elements)

    # Save PDF from buffer
    pdf = buffer.getvalue()
    buffer.close()
    return pdf

def toi_analysis(toi, 
                 fract):
    '''Main shell function that starts the analysis functions and writes
    the results to a pdf. This is what is called from the execution file.
    
    Args:
        toi (float): TOI ID
        fract (int): smallest fraction to divide period into
    
    Returns:
        creates a pdf containing the analyzed information for the TOI
    '''
    #display(Markdown(f"<h2 style='text-align:center';>TOI {toi} at Fraction {fract}</h2>"))
    print('Downloading and extracting light curve data...')
    x = toi_extract_all_curves(toi, fract)
    period = round(x[list(x.keys())[0]][0][2]['period'], 3)
    exptrans = [round(t, 3) for t in list(x.keys())]
    #display(Markdown(f"<p style='text-align:center';>Given period of {period} days, 
    # with transits at: {exptrans} BJD-2457000</p>"))
    
    subtitle = str(f"Given period of {period} days, with transits at: {exptrans} BJD-2457000")
    title = str(f"TOI {toi} at Fraction {fract}")  
    print(f'Analyzing {len(exptrans)} light curves...')
    st = time.time()
    figures, df = plot_toi_fits(x)
    et = time.time()
    print(f'Analyzed {len(exptrans)} light curves in {et - st} seconds.')
    pdf = create_pdf(figures, df, title, subtitle)
    with open(f"{toi}.pdf", 'wb') as f:
        f.write(pdf)
    return

def analyze_all_transits(data=DATA):
    st = time.time()
    count = 0
    skip = 0
    missed = []
    start = 0
    for toi in data['Full TOI ID']:
        
        # allow resuming at certain point in list
        if toi == 1246.04:
            start = 1
        
        # analyze rest of list
        if start == 1:
            try:
                obj = TOI(toi)
                period = obj.period
                if period < 28:
                    print(f'Skipping TOI {toi} with {period} day period.')
                    skip += 1
                    continue
                fract = int(period / 14)
                if fract > 20:
                    fract = 20
            except:
                print(f'Failed (step 1) for TOI {toi}, skipping...')
                missed.append(toi)
                continue

            try:
                print(f'Analyzing TOI {toi} with {period} day period at {fract}th fraction...')
                toi_analysis(toi, fract)
                count += 1
            except:
                print(f'Failed (step 2) for TOI {toi}, skipping...')
                missed.append(toi)

    et = time.time()
    rt = (et - st) / 60
    print(f'Succesfully analyzed {count} TOIs.')
    print(f'Skipped {skip} TOIs with short periods.') 
    print(f'Program failed for TOIs: {missed}.')
    print(f'Total runtime: {rt} minutes.')
    return
    

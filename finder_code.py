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

def generate_all_transits(toi_obj, location, timerange, 
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
    table.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.grey), ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
                               ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                               ('FONTSIZE', (0,0), (-1,0), 14), ('BOTTOMPADDING', (0,0), (-1,0), 12),
                               ('BACKGROUND', (0,1), (-1,-1), colors.beige), ('TEXTCOLOR', (0,1), (-1,-1), colors.black),
                               ('FONTNAME', (0,1), (-1,-1), 'Helvetica'), ('FONTSIZE', (0,1), (-1,-1), 10),
                               ('ALIGN', (0,1), (-1,-1), 'LEFT'), ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
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
def toi_plot_transit(toi, ingress_time, utcoffset=-7, location=SEO):
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

def save_figures_to_pdf(figures, toi):
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

def toi_extract_all_curves(toi, fraction=1):
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

def toi_plot_curves(toi, fraction=1):
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
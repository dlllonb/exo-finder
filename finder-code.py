# file for code

# Import and Style Statements
%matplotlib inline
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


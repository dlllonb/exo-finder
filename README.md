# exo-finder
Version 2.1.0 (June 2023) <br>
Exoplanet analysis code for TESS data

#### Requirements:
* Python environment with 'finder_code.py' and 'finder_exec.py' from this repo, and supporting packages.
  - Requires: astropy, numpy, pandas, matplotlib, reportlab, lightkurve. Use 'pip install *package name*'.
* Internet connection for some functions.
  - If it uses lightkurve, which pulls from online data portal.
* TESS exoplanet database csv.
  - Download from this repo, or get latest release from the [source](https://tev.mit.edu/data/collection/193/).

#### Starting:
* Open a terminal and run 'python finder_exec.py' from the folder containing code files and data file.
  - Further direction on specific use of functions is given at runtime

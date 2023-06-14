# file to execute code from terminal

import finder_code as code

# intro
print('')
print("Welcome to the TESS exoplanet tool kit. Capabilities include predicting and " +
       "plotting specific exoplanet transits for observation from the ground as " +
       "well as retrieving plotting and analyzing the archival TESS data for a " +
       "specific TOI. Additionally, all calculations have the option to include " +
       "fractional pieces of the given orbital periods to aid in possible " +
       "observation of missed transits, and reanalysis of past data in order to " +
       "refine predictions for possible real periods especially with respect to " +
       "long period planets. (V2.1 Dillon Bass, 6/14/23)")

def help():
    print("")
    print("Here is the list of commands:")
    print("  help : shows this page")
    print("  quit : quits the program")
    print("  info: give info on specific TOI")
    print("  predictor: starts the transit predictor")
    print("  plotter: starts the transit plotter")
    print("  curves: starts curve finder for a TOI")
    print("  analyze: starts TOI curve analysis")
    print("")

# show commands initially
help()

# functions called from code file
def info():
    print('-----------------------------------')
    print('Prints out info on a given TOI')
    print('-----------------------------------')

    tid = -1
    while tid == -1:
        tid1 = input("TOI ID: ")
        if tid1 == 'q' or tid1 == 'quit':
            return
        try: 
            tid2 = float(tid1)
            tid = tid2
        except:
            continue
    
    try:
        code.toi_info(tid)
    except:
        print("That is an invalid TOI number (probably).")
    return

def transit_predictor():
    print('-----------------------------------')
    print('Transit predictor gives possibly observabletransits with the given' +
          'parameters. Default settings show a few past transits. Make sure you' +
          'know the correct JD date range. Type q at anytime to quit the program.')
    print('-----------------------------------')
    while True:
        default = input("Default all parameters? (y/n): ")
        if default == 'y':
            print('Running...')
            result = code.find_all_transits()
            print("File 'transtis.csv' and 'transits.pdf' have been created.")
            return
        elif default == 'n': 
            break 
        elif default == 'q' or default == 'quit':
            return
        else:
            continue
    
    # chose not to use defaults, now set them all
    location = [38.29, -122, 56]
    while True:
        new_loc = input("Default location to SEO? (y/n): ")
        if new_loc == 'y':
            break
        elif new_loc == 'n':
            while True:
                nlat = input("Latitude: ")
                if nlat == 'q' or nlat == 'quit':
                    return
                try:
                    nlat1 = float(nlat)
                    location[0] = nlat1
                    break
                except:
                    continue
            while True:
                nlong = input("Longitude: ")
                if nlong == 'q' or nlong == 'quit':
                    return
                try:
                    nlong1 = float(nlong)
                    location[1] = nlong1
                    break
                except:
                    continue
            while True:
                nh = input("Height: ")
                if nh == 'q' or nh == 'quit':
                    return
                try:
                    nh1 = float(nh)
                    location[2] = nh1
                    break
                except:
                    continue 
            break 
        elif new_loc == 'q' or new_loc == 'quit':
            return
        else:
            continue
    elocation = code.earth_location(location)
    
    timerange = [-1, -1]
    while timerange[0] == -1:
        start = input("JD start date: ")
        if start == 'q' or start == 'quit':
            return
        try:
            flt = float(start)
            timerange[0] = flt
        except:
            continue
    while timerange[1] == -1:
        end = input("JD end date: ")
        if end == 'q' or end == 'quit':
            return
        try:
            flt = float(end)
            assert(flt > timerange[0])
            timerange[1] = flt
        except:
            continue
    
    utcoffset = 'x'
    while utcoffset == 'x':
        offset = input("Hours offset from UTC: ")
        if offset == 'q' or offset == 'quit':
            return
        try: 
            hours = int(offset)
            utcoffset = hours
        except:
            continue

    minfrac = -1
    while minfrac == -1:
        frac = input("Smallest fraction of period: ")
        if frac == 'q' or frac == 'quit':
            return
        try: 
            fract = int(frac)
            minfrac = fract
        except:
            continue

    periods = [-1, -1]
    while periods[0] == -1:
        minper = input("Minimum period in days: ")
        if minper == 'q' or minper == 'quit':
            return
        try:
            flt = float(minper)
            periods[0] = flt
        except:
            continue
    while periods[1] == -1:
        maxper = input("Maximum period in days: ")
        if maxper == 'q' or maxper == 'quit':
            return
        try:
            flt = float(maxper)
            periods[1] = flt
        except:
            continue 

    mindepth = -1
    while mindepth == -1:
        mindep = input("Minimum transit depth in ppt: ")
        if mindep == 'q' or offset == 'quit':
            return
        try: 
            mindep1 = float(mindep)
            mindepth = mindep1
        except:
            continue

    minalt = -1
    while minalt == -1:
        minal = input("Minimum transit altitude in degrees: ")
        if minal == 'q' or minal == 'quit':
            return
        try: 
            minal1 = float(minal)
            minalt= minal1
        except:
            continue
    
    print('')
    print('Currently selected parameters are:')
    print(f'  Location : {location} (deg, deg, m)')
    print(f'  Time Range : {timerange} day JD')
    print(f'  UTC offset : {utcoffset} hours')
    print(f'  Dividing Fraction : {minfrac}')
    print(f'  Period Range : {periods} days')
    print(f'  Minimum Transit Depth : {mindepth} ppt')
    print(f'  Minimum Transit Altitude : {minalt} deg')

    while True:
        confirm = input("Confirm running with these parameters? (y/n): ")
        if confirm == 'n':
            print("Process Aborted.")
            return
        elif confirm == 'y': 
            break 
        else:
            continue
    
    try:
        result = code.find_all_transits(location=elocation, timerange=timerange, 
                                        utcoffset=utcoffset, minfrac=minfrac, 
                                        periods=periods, mindepth=mindepth, 
                                        minalt=minalt)
        print("File 'transtis.csv' and 'transits.pdf' have been created.")
    except:
        print("Program failed with given parameters.")
    return 

def transit_plotter():
    print('-----------------------------------')
    print('Transit plotter shows a nice chart of a given transit.')
    print('-----------------------------------')
    location = [38.29, -122, 56]
    while True:
        new_loc = input("Default location to SEO? (y/n): ")
        if new_loc == 'y':
            break
        elif new_loc == 'n':
            while True:
                nlat = input("Latitude: ")
                if nlat == 'q' or nlat == 'quit':
                    return
                try:
                    nlat1 = float(nlat)
                    location[0] = nlat1
                    break
                except:
                    continue
            while True:
                nlong = input("Longitude: ")
                if nlong == 'q' or nlong == 'quit':
                    return
                try:
                    nlong1 = float(nlong)
                    location[1] = nlong1
                    break
                except:
                    continue
            while True:
                nh = input("Height: ")
                if nh == 'q' or nh == 'quit':
                    return
                try:
                    nh1 = float(nh)
                    location[2] = nh1
                    break
                except:
                    continue 
            break 
        elif new_loc == 'q' or new_loc == 'quit':
            return
        else:
            continue
    elocation = code.earth_location(location)

    tid = -1
    while tid == -1:
        tid1 = input("TOI ID: ")
        if tid1 == 'q' or tid1 == 'quit':
            return
        try: 
            tid2 = float(tid1)
            tid = tid2
        except:
            continue
    
    jdt = -1
    while jdt == -1:
        jdt1 = input("Ingress Time JD: ")
        if jdt1 == 'q' or jdt1 == 'quit':
            return
        try: 
            jdt2 = float(jdt1)
            jdt = jdt2
        except:
            continue

    utcoffset = 'x'
    while utcoffset == 'x':
        offset = input("Hours offset from UTC: ")
        if offset == 'q' or offset == 'quit':
            return
        try: 
            hours = int(offset)
            utcoffset = hours
        except:
            continue

    print(f'Plotting transit of {tid} at {jdt} as seen from {location}')
    try:
        code.toi_plot_transit(tid, jdt, utcoffset, elocation)
    except:
        print('Parameters are incompatible with program, please try again.')
        return
    
    print("Output png 'transplot.png' generated.")
    return

def transit_curves():
    print('-----------------------------------')
    print('Produces a output pdf with light curves from the TESS database' +
          'for the given TOI, down to the given fractional period minimum.')
    print('-----------------------------------')

    tid = -1
    while tid == -1:
        tid1 = input("TOI ID: ")
        if tid1 == 'q' or tid1 == 'quit':
            return
        try: 
            tid2 = float(tid1)
            tid = tid2
        except:
            continue
    
    frac = -1
    while frac == -1:
        frac1 = input("Smallest fraction of given period to search: ")
        if frac1 == 'q' or frac1 == 'quit':
            return
        try: 
            frac2 = int(frac1)
            frac = frac2
        except:
            continue
    
    try:
        print("Searching...")
        code.toi_plot_curves(tid, frac)
        print("Output file 'lightcurves.pdf' generated.")
    except:
        print("Program failed with given parameters.")
    return
     
def transit_analysis():
    print('-----------------------------------')
    print('Produces analysis pdf of a given TOI looking specifically at the '+
          'fractional periods down to specified fractional limit. Warning, ' +
          'this process can take some time to run, and generate a lot of warnings ' +
          'depending on given parameters.')
    print('-----------------------------------')

    tid = -1
    while tid == -1:
        tid1 = input("TOI ID: ")
        if tid1 == 'q' or tid1 == 'quit':
            return
        try: 
            tid2 = float(tid1)
            tid = tid2
        except:
            continue
    
    frac = -1
    while frac == -1:
        frac1 = input("Smallest fraction of given period to search: ")
        if frac1 == 'q' or frac1 == 'quit':
            return
        try: 
            frac2 = int(frac1)
            frac = frac2
        except:
            continue
    
    #try:
    print("Analyzing...")
    code.toi_analysis(tid, frac)
    print(f"Output file {tid}.pdf generated.")
    #except:
        #print("Program failed with given parameters.")
    return

# main loop for the program overall 
while True:
    command = input("c: ")
    if command == "q" or command == "quit":
        break
    elif command == "help":
        help()
    elif command == 'info':
        info()
    elif command == "predictor":
        transit_predictor()
    elif command == "plotter":
        transit_plotter()
    elif command == "curves":
        transit_curves()
    elif command == "analyze":
        transit_analysis()
    else:
        print("Unknown command... type help to see list.")

    

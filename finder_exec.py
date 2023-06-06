# file to execute code from terminal

import finder_code as code

def help():
    print("")
    print("Here is the list of commands:")
    print("  help : shows this page")
    print("  quit : quits the program")
    print("  predictor: starts the transit predictor")
    print("")

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
            print("File 'transtis.csv' has been created.")
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
    
    result = code.find_all_transits(location=elocation, timerange=timerange, 
                                    utcoffset=utcoffset, minfrac=minfrac, 
                                    periods=periods, mindepth=mindepth, 
                                    minalt=minalt)
    print("File 'transtis.csv' has been created.")
    return 

# main loop for the program overall 
while True:
    command = input("c: ")
    if command == "q" or command == "quit":
        break
    elif command == "help":
        help()
    elif command == "predictor":
        transit_predictor()
    elif command == "plotter":
        pass
    else:
        print("Unknown command... type help to see list.")
    

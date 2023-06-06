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
    print('Transit predictor gives possibly observable transits with the given parameters')
    while True:
        default = input("Default all parameters? (y/n): ")
        if default == 'y':
            print('Running...')
            return code.find_all_transits()
        elif default == 'n': 
            break 
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
                try:
                    nlat1 = float(nlat)
                    location[0] = nlat1
                    break
                except:
                    continue
            while True:
                nlong = input("Longitude: ")
                try:
                    nlong1 = float(nlong)
                    location[1] = nlong1
                    break
                except:
                    continue
            while True:
                nh = input("Height: ")
                try:
                    nh1 = float(nh)
                    location[2] = nh1
                    break
                except:
                    continue 
            break 
        else:
            continue
    elocation = code.earth_location(location)
    
    timerange = [-1, -1]
    while timerange[0] == -1:
        start = input("JD start date: ")
        try:
            flt = float(start)
            timerange[0] = flt
        except:
            continue
    while timerange[1] == -1:
        end = input("JD end date: ")
        try:
            flt = float(end)
            assert(flt > timerange[0])
            timerange[1] = flt
        except:
            continue
    
    utcoffset = 'x'
    while utcoffset == 'x':
        offset = input("Hours offset from UTC: ")
        try: 
            hours = int(offset)
            utcoffset = hours
        except:
            continue

    minfrac = -1
    while minfrac == -1:
        frac = input("Smallest fraction of period: ")
        try: 
            fract = int(frac)
            minfrac = fract
        except:
            continue

    periods = [-1, -1]
    while periods[0] == -1:
        minper = input("Minimum period: ")
        try:
            flt = float(minper)
            periods[0] = flt
        except:
            continue
    while periods[1] == -1:
        maxper = input("Maximum period: ")
        try:
            flt = float(maxper)
            periods[1] = flt
        except:
            continue 

    mindepth = -1
    while mindepth == -1:
        mindep = input("Minimum transit depth: ")
        try: 
            mindep1 = float(mindep)
            mindepth = mindep1
        except:
            continue

    minalt = -1
    while minalt == -1:
        minal = input("Minimum transit altitude: ")
        try: 
            minal1 = float(minal)
            minalt= minal1
        except:
            continue
    
    print('')
    print('Currently selected parameters are:')
    print(f'  Location : {location}')
    print(f'  Time Range : {timerange} day JD')
    print(f'  UTC offset : {utcoffset} hours')
    print(f'  Dividing Fraction : {minfrac}')
    print(f'  Period Range : {periods} days')
    print(f'  Minimum Transit Depth : {mindepth}')
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
                                    periods=periods, mindepth=mindepth, minalt=minalt)
    print("File 'transtis.csv' has been created.")
    return 

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
    

import ConfigParser
import numpy as np
import ephem
import math

# Create random array for source
def create_source(N_src, duration, priority):
    src_id = range(1,N_src+1)
    src_ra = np.random.randint(0, 2400, N_src)/100.
    src_dec = 90.0 - np.random.randint(0, 18000, N_src)/100.
    src_duration = np.random.randint(duration[0], duration[1], N_src)
    src_priority = np.sort(np.random.randint(10*priority[0], 10*priority[1], N_src)/10.)[::-1]
    src_input = np.transpose([src_id,src_ra,src_dec,src_duration,src_priority]) 
    return src_input

# Writing source input to a file
def write_data_file(src_input, data_file):
      np.savetxt(data_file, src_input, delimiter='\t', fmt='%s', header='ID\tRA\tDec\tDur\tPrty\tRise time\t\tSet time', comments='#')

# Writing result to a file
def write_result_file(src_result, result_file):
      np.savetxt(result_file, src_result, delimiter='\t', fmt='%s', header='ID\tStart time\t\tEnd time\t\tSlew\tDur\tPrty', comments='#')


# Reading source input from a file
def read_data_file(data_file):
    src_data = np.loadtxt(data_file, delimiter='\t', comments='#', dtype='S') 
    return src_data

# Reading config file
def read_config_file():
    config={}
    config_file='obsim_config.txt'
    config_obj = ConfigParser.RawConfigParser()
    config_obj.read(config_file)
    config['new_input_data'] = config_obj.get('Sources', 'generate_new_input_data')
    config['N_src'] = config_obj.getint('Sources', 'number_of_sources')
    config['min_duration'] = config_obj.getfloat('Sources', 'min_duration')
    config['max_duration'] = config_obj.getfloat('Sources', 'max_duration')
    config['min_priority'] = config_obj.getfloat('Sources', 'min_priority')
    config['max_priority'] = config_obj.getfloat('Sources', 'max_priority')
    config['data_file'] = config_obj.get('File', 'data_file')
    config['result_file'] = config_obj.get('File', 'result_file')
    config['date'] = config_obj.get('Date', 'start_date')
    config['start_time'] = config_obj.get('Date', 'start_time')
    config['end_time'] = config_obj.get('Date', 'end_time')
    config['time_zone'] = config_obj.getfloat('Date', 'time_zone')
    config['lat'] = config_obj.getfloat('Location', 'latitude')
    config['lon'] = config_obj.getfloat('Location', 'longitude')
    config['alt'] = config_obj.getfloat('Location', 'altitude')
    config['telescope_azimuth'] = config_obj.getfloat('Telescope_position', 'telescope_azimuth')
    config['telescope_altitude'] = config_obj.getfloat('Telescope_position', 'telescope_altitude')
    config['slew_rate_az'] = config_obj.getfloat('Telescope_position', 'slew_rate_az')
    config['slew_rate_alt'] = config_obj.getfloat('Telescope_position', 'slew_rate_alt')
    config['priority_wait_time'] = config_obj.getfloat('schedule_by_priority', 'wait_time')
    return config


# Initializing ephem location varaible
def init_loc(config):
    loc = ephem.Observer()
    loc.lon = str(config['lon'])
    loc.lat = str(config['lat'])
    loc.elevation = config['alt']
    loc.epoch="2000/1/1 12:00:00"
    loc.date = (config['date'] + " " + config['start_time'] )
    loc.pressure = 0
    return loc

# Initializing ephem fixed body source variable
def init_src():
    src = ephem.FixedBody()
    src._epoch = "2000/1/1 12:00:00"   
    return src


# Compute source rise time
def src_rise_time(src, loc, start_date_time, end_date_time):
    loc.date = start_date_time
    src.compute(loc)
    if (src.neverup or src.rise_time > end_date_time):
        raise Exception("Notup")
    if (src.rise_time == None):
        rise_time = start_date_time
    elif (src.rise_time < start_date_time and src.set_time > start_date_time):
        rise_time = start_date_time
    elif (src.rise_time > start_date_time and src.rise_time < end_date_time ):
        rise_time = src.rise_time
    else: 
        loc.date = end_date_time
        src.compute(loc)
        if (src.rise_time > start_date_time and src.rise_time < end_date_time):
            rise_time = src.rise_time
        elif (src.rise_time == None or src.rise_time < start_date_time):
            rise_time = start_date_time
        else:
            raise Exception("Notup")
    return rise_time


#Compute source set time
def src_set_time(src, loc, start_date_time, end_date_time):
    loc.date = start_date_time
    src.compute(loc)
    if (src.neverup):
        raise Exception("Notup")
    if (src.set_time == None or src.set_time > end_date_time):
        set_time = end_date_time
    elif (src.set_time > start_date_time and src.set_time < end_date_time ):
        set_time = src.set_time
    else: 
        loc.date = end_date_time
        src.compute(loc)
        if (src.set_time > start_date_time and src.set_time < end_date_time):
            set_time = src.set_time
        elif (src.set_time == None or src.set_time > end_date_time):
            set_time = end_date_time
        else:
            raise Exception("Notup")
    return set_time
 

# Get start and end date for the observation in UT time
def get_start_end_dates(config):
    start_date_time = ephem.date(config['date'] + " " + config['start_time'])
    end_date_time = ephem.date(config['date'] + " " + config['end_time'])
    if ( start_date_time > end_date_time):
            end_date_time = ephem.date(end_date_time + 1 )
    tz_offset = config['time_zone'] * ephem.hour
    start_date_time = ephem.date(start_date_time - tz_offset)
    end_date_time = ephem.date(end_date_time - tz_offset)
    return start_date_time, end_date_time


# Calculate shortest angular difference in azimuth and altitude
def calc_angle_difference(az1, alt1, az2, alt2):
    az_diff = ephem.degrees(math.fabs(az2 - az1)) 
    alt_diff = ephem.degrees(math.fabs(alt2 - alt1))
    if (az_diff > math.pi):
        az_diff = ephem.degrees(2*math.pi - az_diff)
    return az_diff, alt_diff


# Compute source azimuth and altitude at a given time
def get_src_az_alt(src, loc, ra, dec, curr_time):
    src._ra = str(ra)
    src._dec = str(dec)
    loc.date = curr_time    
    src.compute(loc)
    return src.az, src.alt 


def print_summary(schedule, src_data):
    N_src = np.shape(src_data)[0]
    x = np.array(schedule)
    N_obs_id = x[:,0]
    N_obs = np.shape(N_obs_id[N_obs_id !='Wait'])[0]
    src_priority = src_data[:,4].astype(np.float)
    obs_priority_tmp = x[:,5]
    obs_priority = obs_priority_tmp[obs_priority_tmp != '0'].astype(np.float)
    total_slew = (x[:,3]).astype(np.float).sum()
    total_duration = x[:,4].astype(np.float).sum()
    print "From total of %s sources,  %s sources is scheduled for observation" %(N_src, N_obs)
    print "Total slew time is %.2f hour" %(total_slew/3600.)
    print "Total observing duration is %.2f hour" %(total_duration/3600.)
    print "Median value of priority for the entire source is %.2f" % np.median(src_priority)
    print "Median value of priority for the scheduled source is %.2f" % np.median(obs_priority)


if __name__ == "__main__":
    print "This is a funtion module supporting the obsim.py main program"

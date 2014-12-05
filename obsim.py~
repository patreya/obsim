import numpy as np
import ephem
import os.path, math
import obsim_util as obsim

# Main routine to create source input, filter them if they are seen from the location, and then save them in the text file
def create_new_src_input(config):
    N_src = config['N_src']
    duration = [ config['min_duration'], config['max_duration'] ]
    priority =  [ config['min_priority'], config['max_priority'] ]
    src_input = obsim.create_source(N_src, duration , priority )
    loc = obsim.init_loc(config)
    start_date_time, end_date_time = obsim.get_start_end_dates(config)
    src = obsim.init_src()
    src_data=[]
    removed_src = 0
    j = 0
    tz_offset = config['time_zone'] * ephem.hour
    for i in range(0, N_src):
        src.name = str(src_input[i][0])
        src._ra = str(src_input[i][1])
        src._dec = str(src_input[i][2])
        loc.date = start_date_time
        src.compute(loc)
        try:
            rise_time = obsim.src_rise_time(src, loc, start_date_time, end_date_time)
            rise_time_local = ephem.date(rise_time + tz_offset)
            try:
                set_time = obsim.src_set_time(src, loc, start_date_time, end_date_time)
                set_time_local = ephem.date(set_time + tz_offset)
                src_data.append ([str(j), str(src_input[i][1]), str(src_input[i][2]), str(int(src_input[i][3])), str(src_input[i][4]), str(rise_time_local), str(set_time_local)])
                j += 1
            except Exception, Notup:
                removed_src += 1
        except Exception, Notup:
             removed_src += 1
    obsim.write_data_file(src_data, config['data_file'])
    print "%s source stored in the data file %s" % (j-1, config['data_file'])


# get the next visible source with highest priority
def get_next_src(src_data, N_src, src, loc, curr_time, tz_offset, curr_tel_pos, slew_rate, schedule):
    for i in range(0, N_src):
        src_obs_status = src_data[i][7]
        curr_time_local =  ephem.date(curr_time + tz_offset)
        if (src_obs_status == '0'):
            src_rise_time_local = ephem.Date(src_data[i][5])
            src_rise_time = ephem.Date(src_rise_time_local  - tz_offset)
            src_az, src_alt = obsim.get_src_az_alt(src, loc, src_data[i][1], src_data[i][2], curr_time)
            az_diff, alt_diff = obsim.calc_angle_difference(src_az, src_alt, curr_tel_pos[0], curr_tel_pos[1])
            slew_time = (math.degrees(az_diff)/slew_rate[0] + math.degrees(alt_diff)/slew_rate[0])*ephem.second
            src_obs_start_time = ephem.Date(curr_time + slew_time)
            if (src_obs_start_time >= src_rise_time):
                src_set_time_local = ephem.Date(src_data[i][6])
                src_set_time = ephem.Date(src_set_time_local  - tz_offset)
                src_obs_duration = float(src_data[i][3])*ephem.second
                src_obs_end_time = ephem.Date(src_obs_start_time + src_obs_duration)
                if ( src_obs_end_time < src_set_time):
                    src_az_end, src_alt_end = obsim.get_src_az_alt(src, loc, src_data[i][1], src_data[i][2], src_obs_end_time)
                    curr_time_local =  ephem.date(curr_time + tz_offset)
                    src_obs_end_time_local =  ephem.date(src_obs_end_time + tz_offset)
                    schedule.append( [ int(src_data[i][0]), str(curr_time_local), str(src_obs_end_time_local), int(round(slew_time*24.0*3600.0)), int(round(src_obs_duration*24.0*3600.0)), (src_data[i][4]) ] )
                    status=[1, i, src_obs_end_time, src_az_end, src_alt_end]
                    return schedule, status
    status=[0] 
    return schedule, status
                          

# wrapper for get_next_src and called by main 
def get_schedule_by_priority(src_data, config):
    N_src = np.shape(src_data)[0]
    print "Using schedule by priority algorithm for %s sources" % N_src
    src_data = np.column_stack((src_data, np.zeros(N_src, dtype=int)))
    start_date_time, end_date_time = obsim.get_start_end_dates(config)
    tz_offset = config['time_zone'] * ephem.hour
    wait_time = config['priority_wait_time'] * ephem.second
    loc = obsim.init_loc(config)
    src = obsim.init_src()
    # Initializing current position and time
    curr_time = start_date_time
    curr_tel_pos = [ ephem.degrees(str(config["telescope_azimuth"])), ephem.degrees(str(config["telescope_altitude"])) ]
    slew_rate = [ config["slew_rate_az"], config["slew_rate_alt"] ]
    schedule = []
    while (curr_time < end_date_time):
        schedule, status = get_next_src(src_data, N_src, src, loc, curr_time, tz_offset, curr_tel_pos, slew_rate, schedule)
        if (status[0]):
            src_data[status[1]][7] = 1
            curr_time = status[2]
            curr_tel_pos = [status[3], status[4]]
        else:
            from_time_local =  ephem.Date(curr_time + tz_offset)
            curr_time = ephem.Date(curr_time + wait_time)
            if (curr_time > end_date_time):
                curr_time = end_date_time
            curr_time_local =  ephem.date(curr_time + tz_offset)
            time_diff = (curr_time_local - from_time_local)*24.*3600
            schedule.append( [ 'Wait', str(from_time_local), str(curr_time_local), int(round(time_diff)), '0', '0' ] )
    return schedule               


def main():
    # Reading config settings
    config = obsim.read_config_file()
    # Generate new source data or read existing data from file
    if ( ( not os.path.isfile(config['data_file']) ) or ( config['new_input_data'] != 'no' ) ) :
        print "Generating a new input parameters"
        create_new_src_input(config)
    # Reading filtered source data from the file
    src_data = obsim.read_data_file(config['data_file'])
    # get schedule of observation using priority of the source
    schedule = get_schedule_by_priority(src_data, config)
    # Write result in to a file
    obsim.write_result_file(schedule, config['result_file'])
    # print summary of the schedule
    obsim.print_summary(schedule, src_data)


if __name__ == "__main__":
    main()



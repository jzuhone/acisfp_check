# ----------------------------------------------------------------
#
# who_in_fp
#
# ----------------------------------------------------------------
def who_in_fp(simpos=80655):
    """
    Returns a string telling you which instrument is in
    the Focal Plane. "launchlock" is returned because that's a
    position we never expect to see the sim in - it's an indicator
    to the user that there's a problem.

    Also, The ranges for detector sections use the max and min hard
    stop locations, and they also split the difference between "I"
    and "S" for each instrument.

          input: - TSC position (simpos) - INTEGER

          output - String indicating what is in the focal plane
                   "launchlock" - default
                   "ACIS-I"
                   "ACIS-S"
                   "HRC-I"
                   "HRC-S"
    """
    is_in_the_fp = 'launchlock'

    #  Set the value of is_in_the_fp to the appropriate value. It will default
    #  to "launchlock" if no value matches
    if 104839 >= simpos >= 82109:
        is_in_the_fp = 'ACIS-I'
    elif 82108 >= simpos >= 70736:
        is_in_the_fp = 'ACIS-S'
    elif -20000 >= simpos >= -86147:
        is_in_the_fp = 'HRC-I'
    elif -86148 >= simpos >= -104362:
        is_in_the_fp = 'HRC-S'

    #  return the string indicating which instrument is in the Focal Plane
    return is_in_the_fp


# ---------------------------------------------------------------------
#
# find_obsid_intervals
#
# ---------------------------------------------------------------------

def find_obsid_intervals(cmd_states):
    """
    User reads the SKA commanded states archive, via
    a call to the SKA get_cmd_states, between the
    user specified START and STOP times. 

    An example of the read would be:

        start_time = '2005:001'
        stop_time = '2005:031'  # i.e. January

        cmd_states = cmd_statesFetch(start_time, stop_time)

    Problem is, ALL commanded states that were stored
    in the archive will be returned. So then you call:

        find_obsid_intervals(cmd_states)

    And this will find the obsid intervals.
    What this program does is to extract the time interval for
    each OBSID. Said interval start is defined by a
    WSPOW00000/WSVIDALLDN, and the interval end is
    defined by the first AA000000 that follows.

    When the interval has been found,
    a list element is created from the value of
    states data at the time point of the first NPNT
    line seen - *minus* the trans_keys, tstart and tstop
    times. The values of datestart and datestop are
    the WSPOW00000/WSVIDALLDN and AA000000 times. The
    exposure time of the interval is also tacked on to
    the end of the list. This list
    is appended to a Master list of all obsid intervals
    and this list is returned. Users

    Notes: The obsid filtering method includes the
           configuration from the last OBSID, through
           a setup for the present OBSID, through the
           XTZ - AA000, down to the power down.

            - This might show a cooling from the
              last config, temp changes due to some
              possible maneuvering, past shutdown
    """
    #
    # Some inits
    #

    # a little initialization
    firstpow = False
    xtztime = None
    simpos = None

    # EXTRACTING THE OBSERVATIONS
    #
    # Find the first line with a WSPOW00000 in it. This is the start of
    # the interval. Then get the first XTZ line, the NPNT line, the
    # AA000000 line, and lastly the next WSPOW00000 line.
    # This constitutes one observation.

    obsid_interval_list = []

    for eachstate in cmd_states:

        # Make sure we skip maneuver obsids explicitly
        if 50000 > eachstate['obsid'] >= 38001:
            continue

        pow_cmd = eachstate['power_cmd']

        # is this the first WSPOW of the interval?
        if pow_cmd in ['WSPOW00000', 'WSVIDALLDN'] and not firstpow:
            firstpow = True
            datestart = eachstate['datestart']
            tstart = eachstate['tstart']

        # Process the first XTZ0000005 line you see
        if pow_cmd in ['XTZ0000005', 'XCZ0000005'] and \
                (xtztime is None and firstpow):
            xtztime = eachstate['tstart']
            # MUST fix the instrument now
            instrument = who_in_fp(eachstate['simpos'])

        # Process the first AA00000000 line you see
        if pow_cmd == 'AA00000000' and firstpow:
            datestop = eachstate['datestop']
            tstop = eachstate['tstop']

            # now calculate the exposure time
            if xtztime is not None:

                # Having found the startScience and stopScience, you have an
                # OBSID interval. Now form the element and append it to
                # the Master List. We add the text version of who is in
                # the focal plane

                obsid_dict = {"datestart": datestart,
                              "datestop": datestop,
                              "tstart": tstart,
                              "tstop": tstop,
                              "obsid": eachstate['obsid'],
                              "instrument": instrument}
                obsid_interval_list.append(obsid_dict)

            # now clear out the data values
            firstpow = False
            xtztime = None
            simpos = None

    # End of LOOP for eachstate in cmd_states:

    return obsid_interval_list


# --------------------------------------------------------------------------
#
#   hrc_science_obs_filter - filter *OUT* any HRC science observations
#
# --------------------------------------------------------------------------

def hrc_science_obs_filter(obsidinterval_list):
    """
    This method will filter *OUT* any HRC science observations from the
    input obsid interval list. Filtered are obs that have either
    HRC-I" or HRC-S" as the science instrument, AND an obsid LESS THAN
    50,000
    """
    acis_and_ecs_only = []
    for eachobservation in obsidinterval_list:
        if eachobservation["instrument"].startswith("ACIS-") or \
                eachobservation["obsid"] >= 50000:
            acis_and_ecs_only.append(eachobservation)
    return acis_and_ecs_only


# --------------------------------------------------------------------------
#
#   ecs_only_filter
#
# --------------------------------------------------------------------------

def ecs_only_filter(obsidinterval_list):
    """
    This method will filter out any science observation from the
    input obsid interval list. It keeps any observation that has an
    obsid of 50,000 or greater
    """
    ecs_only = []
    for eachobservation in obsidinterval_list:
        if eachobservation["obsid"] >= 60000 and \
                eachobservation["instrument"] == "HRC-S":
            ecs_only.append(eachobservation)
    return ecs_only


# ----------------------------------------------------------------------
#
#  get_all_specific_instrument
#
# ---------------------------------------------------------------------

def get_all_specific_instrument(observations, instrument):
    """
    Given a list  of obsid intervals extracted by this class
    class, return the list all those obsids with the specified instrument
    """
    same_inst = []
    for eachobs in observations:
        if eachobs["instrument"] == instrument:
            same_inst.append(eachobs)
    return same_inst

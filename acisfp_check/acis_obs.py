###############################################################################
#
#   ObsidList - Class that will extract CHANDRA ACIS OBSIDs using
#               the commanded states database. Also provided are
#               a series of filters the user can use to select
#               Observations of a particular configuration.
#
#               Users must supply a start and stop time for extraction
#               of states from the Commanded States Data base
#
###############################################################################

#----------------------------------------------------------------
#
# who_in_fp
#
#----------------------------------------------------------------
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


class ObsidList:

    def __init__(self):
        #
        # Define the indexes to be used in an element of an Intervals list
        #
        self.datestart  = 0
        self.datestop   = 1
        self.tstart     = 2
        self.tstop      = 3
        self.obsid      = 4
        self.in_focal_plane = 5
        # internally maintained results data structures. We do not keep
        # every result at the moment (e.g. observations filtered on pitch)
        # But that may change in the future.
        self.cmd_states = None
        self.obsid_interval_list = []

    # ---------------------------------------------------------------------
    #
    # find_obsid_intervals
    #
    # ---------------------------------------------------------------------
    def find_obsid_intervals(self, cmd_states):
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
        firstpow = None
        obsid = None
        xtztime = None
        aa0time = None

        # EXTRACTING THE OBSERVATIONS
        #
        # Find the first line with a WSPOW00000 in it. This is the start of
        # the interval. Then get the first XTZ line, the NPNT line, the
        # AA000000 line, and lastly the next WSPOW00000 line.
        # This constitutes one observation.

        for eachstate in cmd_states:

            # Make sure we skip maneuver obsids explicitly
            if 50000 > eachstate['obsid'] >= 38001:
                continue

            # is this the first WSPOW of the interval?
            if eachstate['power_cmd'] in ['WSPOW00000', 'WSVIDALLDN'] and \
               firstpow is None:
                firstpow = eachstate
                DOYfetchstart = eachstate['datestart']
                secsfetchstart = eachstate['tstart']

            # Process the first XTZ0000005 line you see
            if eachstate['power_cmd'] in ['XTZ0000005', 'XCZ0000005'] and \
               (xtztime is None and firstpow is not None):
                xtztime = eachstate['tstart']

            # Process the first NPNT line you see
            if obsid is None and firstpow is not None:
                obsid = eachstate['obsid']

            # Process the first AA00000000 line you see
            if eachstate['power_cmd'] == 'AA00000000' and aa0time is None and firstpow is not None:
                aa0time = eachstate['tstop']
                DOYfetchstop = eachstate['datestop']
                secsfetchstop = eachstate['tstop']

                # now calculate the exposure time
                if xtztime is not None:

                    # Having found the startScience and stopScience, you have an
                    # OBSID interval. Now form the list element and append it to
                    # the Master List. We add on the exposure time and the text
                    # version of who is in the focal plane
                    science_instrument = who_in_fp(eachstate['simpos'])

                    self.obsid_interval_list.append([DOYfetchstart,
                                                     DOYfetchstop,
                                                     secsfetchstart,
                                                     secsfetchstop,
                                                     obsid,
                                                     science_instrument])

                # now clear out the data values
                firstpow = None
                obsid = None
                xtztime = None
                aa0time = None

        # End of LOOP for eachstate in cmd_states:

        return self.obsid_interval_list

    ######################################################################
    #
    #   FILTERS
    #
    ######################################################################


    #--------------------------------------------------------------------------
    #
    #   hrc_science_obs_filter - filter *OUT* any HRC science observations
    #
    #--------------------------------------------------------------------------
    def hrc_science_obs_filter(self, obsidinterval_list):
        """
        This method will filter *OUT* any HRC science observations from the
        input obsid interval list. Filtered are obs that have either
        HRC-I" or HRC-S" as the science instrument, AND an obsid LESS THAN
        50,000
        """
        acis_and_ecs_only = []
        for eachobservation in obsidinterval_list:
            if eachobservation[self.in_focal_plane].startswith("ACIS-") or \
               eachobservation[self.obsid] >= 50000:
                acis_and_ecs_only.append(eachobservation)
        return acis_and_ecs_only


    #--------------------------------------------------------------------------
    #
    #   ecs_only_filter
    #
    #--------------------------------------------------------------------------
    def ecs_only_filter(self, obsidinterval_list):
        """
        This method will filter out any science observation from the
        input obsid interval list.It keeps any observation that has an
        obsid of 50,000 or greater
        """
        ecs_only = []
        for eachobservation in obsidinterval_list:
            if eachobservation[self.obsid] >= 60000 and \
               eachobservation[self.in_focal_plane] == "HRC-S":
                ecs_only.append(eachobservation)
        return ecs_only


    ######################################################################
    #
    #   GETS
    #
    ######################################################################

    #----------------------------------------------------------------------
    #
    #  get_all_specific_instrument
    #
    #---------------------------------------------------------------------
    def get_all_specific_instrument(self, observations, instrument):
        """
        Given a list  of obsid intervals extracted by this class
        class, return the list all those obsids with the specified instrument
        """
        same_inst = []
        for eachobs in observations:
            if eachobs[self.in_focal_plane] == instrument:
                same_inst.append(eachobs)
        return same_inst


    #----------------------------------------------------------------------
    #
    # get_obsid
    #
    #---------------------------------------------------------------------
    def get_obsid(self, observation):
        """
        Given a list element from the list of obsids extracted by this
        class, extract the obsid and return it.
        NOTE: type(obsid) = int!!!
        """
        return observation[self.obsid] if len(observation) > 0 else None


    #----------------------------------------------------------------------
    #
    # get_tstart
    #
    #---------------------------------------------------------------------
    def get_tstart(self, observation):
        """
        Given a list element from the list of obsids extracted by this
        class, extract the tstart and return it
        """
        return observation[self.tstart]

    #----------------------------------------------------------------------
    #
    # get_tstop
    #
    #---------------------------------------------------------------------
    def get_tstop(self, observation):
        """
        Given a list element from the list of obsids extracted by this
        class, extract the tstop and return it
        """
        return observation[self.tstop]

###############################################################################
#
#   ObsidFindFilter - Class that will extract CHANDRA ACIS OBSIDs using
#                        the commanded states database. Also provided are
#                        a series of filters the user can use to select
#                        Observations of a particular configuration.
#
#                        These filters are: Exposure time range
#                                           CCD Count range
#                                           CTI observation removal
#                                           Pitch range
#
#
#                        Users must supply a start and stop time for extraction
#                        of states from the Commanded States Data base
#
###############################################################################
from Chandra.Time import DateTime
#----------------------------------------------------------------
#
# who_in_fp.py
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

class ObsidFindFilter():

    def __init__(self):
        #
        # Define the indexes to be used in an element of an Intervals list
        #
        self.datestart  = 0
        self.datestop   = 1
        self.tstart     = 2
        self.tstop      = 3
        self.obsid      = 4
        self.power_cmd  = 5
        self.si_mode    = 6
        self.pcad_mode  = 7
        self.vid_board  = 8
        self.clocking   = 9
        self.fep_count  = 10
        self.ccd_count  = 11
        self.simpos     = 12
        self.simfa_pos  = 13
        self.pitch      = 14
        self.ra         = 15
        self.dec        = 16
        self.roll       = 17
        self.q1         = 18
        self.q2         = 19
        self.q3         = 20
        self.q4         = 21
        self.trans_keys = 22
        self.hetg       = 23
        self.letg       = 24
        self.dither     = 25
        self.exptime    = 26
        self.in_focal_plane = 27
        self.is_fp_sensitive = 28
        # internally maintained results data structures. We do not keep
        # every result at the moment (e.g. observations filtered on pitch)
        # But that may change in the future.
        self.cmd_states = None
        self.obsid_interval_list = None
        self.non_CTI_obs = None

        self.list_of_sensitive_obs = []

    # ---------------------------------------------------------------------
    #
    # find_obsid_intervals
    #
    # ---------------------------------------------------------------------
    def find_obsid_intervals(self, cmd_states, outfilespec=None):
        """
        User reads the SKA commanded states archive, via 
        a call to the SKA get_cmd_states, between the
        user specified START and STOP times. The following
        items are returned as a numpy array:

             Data item and type
             ------------------
             ('datestart', '|S21'), 
             ('datestop', '|S21'), 
             ('tstart', '<f8'), 
             ('tstop', '<f8'), 
             ('obsid', '<i8'), 
             ('power_cmd', '|S10'), 
             ('si_mode', '|S8'), 
             ('pcad_mode', '|S4'), 
             ('vid_board', '<i8'), 
             ('clocking', '<i8'), 
             ('fep_count', '<i8'), 
             ('ccd_count', '<i8'), 
             ('simpos', '<i8'), 
             ('simfa_pos', '<i8'), 
             ('pitch', '<f8'), 
             ('ra', '<f8'), 
             ('dec', '<f8'), 
             ('roll', '<f8'), 
             ('q1', '<f8'), 
             ('q2', '<f8'), 
             ('q3', '<f8'), 
             ('q4', '<f8'),
             ('trans_keys', '|S48'),
             ('hetg', '|S4'), 
             ('letg', '|S4'), 
             ('dither', '|S4'),
             ('exptime', '<f8')
      
     
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
        min_exptime = 30000.0

        # Open the output file, and write out a header
        if outfilespec != None:
            outfile = open(outfilespec, 'w')
            outfile.write('DOYstart DOYstop TSTART TSTOP OBSID PWR_CMD '
                          'SI_MODE PCAD_MODE VID_BOARD CLOCKING FEP_COUNT '
                          'CCD_CNT SIMPOS SIMFA_POS PITCH RA DEC ROLL Q1 '
                          'Q2 Q3 Q4 TRANS-KEYS HETG LETG Dither EXPTIME')

        # a little initialization
        firstpow = None
        DOYfetchstart = ' '
        obsid = None
        xtztime = None
        aa0time = None
        exptime = '-1'
        pitch = ''
        ccdcnt = ''
        DOYfetchstop = ' '

        self.obsid_interval_list = []

        # EXTRACTING THE OBSERVATIONS
        #
        # Find the first line with a WSPOW00000 in it. This is the start of 
        # the interval. Then get the first XTZ line, the NPNT line, the 
        # AA000000 line, and lastly the next WSPOW00000 line. 
        #
        #This constitutes one observation.
 
        for eachstate in cmd_states:

            # is this the first WSPOW of the interval?
            if (eachstate.power_cmd == 'WSPOW00000' or eachstate.power_cmd == 'WSVIDALLDN') and \
               firstpow is None:
                firstpow = eachstate
                DOYfetchstart = eachstate.datestart
                secsfetchstart = DateTime(DOYfetchstart).secs

            # Process the first XTZ0000005 line you see
            if (eachstate.power_cmd == 'XTZ0000005' or eachstate.power_cmd == 'XCZ0000005') and \
               (xtztime is None and firstpow is not None):
                xtztime = DateTime(eachstate.datestart).secs

            # Process the first NPNT line you see
            if obsid is None and firstpow is not None:
                obsid = eachstate.obsid
                power_cmd = eachstate.power_cmd 
                si_mode = eachstate.si_mode
                pcad_mode = eachstate.pcad_mode
                vid_board = eachstate.vid_board
                clocking = eachstate.clocking
                fep_count = eachstate.fep_count
                ccd_cnt = eachstate.ccd_count
                simpos = eachstate.simpos
                simfa_pos = eachstate.simfa_pos
                pitch = eachstate.pitch
                ra = eachstate.ra
                dec = eachstate.dec
                roll = eachstate.roll
                q1 = eachstate.q1
                q2 = eachstate.q2
                q3 = eachstate.q3
                q4 = eachstate.q4
                trans_keys = eachstate.trans_keys
                hetg = eachstate.hetg
                letg = eachstate.letg
                dither = eachstate.dither

            # Process the first AA00000000 line you see
            if eachstate.power_cmd == 'AA00000000' and aa0time is None and firstpow is not None:
                aa0time = DateTime(eachstate.datestart).secs
                DOYfetchstop = eachstate.datestop
                secsfetchstop = DateTime(DOYfetchstop).secs

                # now calculate the exposure time
                if xtztime is not None:
                    exptime = round(float(aa0time)) - round(float(xtztime))
                else:
                    exptime = -1

                # Having found the AA0000000, you have an OBSID interval. Now
                # form the list element and append it to the Master List. We
                # add on the exposure time and the text version of who is in
                # the focal plane
                science_instrument = who_in_fp(simpos)

                self.obsid_interval_list.append([DOYfetchstart, 
                                                 DOYfetchstop, 
                                                 secsfetchstart, 
                                                 secsfetchstop, 
                                                 obsid,
                                                 power_cmd,
                                                 si_mode,
                                                 pcad_mode,
                                                 vid_board,
                                                 clocking,
                                                 fep_count,
                                                 ccd_cnt, 
                                                 simpos, 
                                                 simfa_pos, 
                                                 pitch, 
                                                 ra, 
                                                 dec, 
                                                 roll, 
                                                 q1, q2, q3, q4, 
                                                 trans_keys,
                                                 hetg, 
                                                 letg, 
                                                 dither,
                                                 exptime,
                                                 science_instrument])

                # now  clear out the data values
                firstpow = None
                DOYfetchstart = ' '
                obsid = None
                xtztime = None
                aa0time = None
                exptime = '-1'
                pitch = ''
                DOYfetchstop= ' '

        # End of LOOP for eachstate in cmd_states:

        # Write out the OBSID interval list and return it upon exit

        if outfilespec != None:
            outfile.write("\n\n")
            outfile.write(str(self.obsid_interval_list) )
            outfile.close()

        return self.obsid_interval_list

    ######################################################################
    #
    #   FILTERS
    #
    ######################################################################

    #---------------------------------------------------------------------
    #
    #   exp_time_filter
    #
    #---------------------------------------------------------------------
    def exp_time_filter(self, obsidinterval_list, startstoplist):
        """
        Given the output of FindObsidIntervals(cmd_states, ' ')
        (which is a list of commanded states with the exposure
        time tacked onto the end, and the OBSID start/stop
        times at the beginning...one state per obsid.), find
        all obsid intervals whose exposure is within the range 
        specified and return those intervals as a list.
     
     
        input: Interval list to be searched
     
               List containing the min and max times
               If the list contains only one time, return
               all intervals whose exposure time is greater than or
               equal to that time
     
        output: A list of all intervals whose exposure times are
                within the start/stop time inclusive.
     
        usage: start_time = '2010:001'
               stop_time  = '2010:014' 
               cmd_states = cmd_statesFetch(start_time, stop_time)
               OBSIDIntervals = FindObsidIntervals(cmd_states, '')
     
        expintervals = ExpTimeFilter(OBSIDIntervals, [min_exp_length, <max_exp_length>])
        """
        exptimelist = []

        # Capture the start and stop time for use in the filter. If the
        # length of the stop
        min_exp_length = startstoplist[0]

        if len(startstoplist) == 1:
            max_exp_length = min_exp_length
        else:
            max_exp_length = startstoplist[1]

        for eachinterval in obsidinterval_list:
            if max_exp_length >= eachinterval[self.exptime] >= min_exp_length:
                exptimelist.append(eachinterval)

        return exptimelist

    #---------------------------------------------------------------------
    #
    #   CTI_filter
    #
    #---------------------------------------------------------------------
    def cti_filter(self, obsidinterval_list):
        """
        Given a list of obsid intervals, remove any obsid
                     interval which is a CTI observation (i.e. has 
                     an OBSID greater than 50k).
                   - return the filtered list
                     specified and return those intervals as a list.
     
                input: OBSIDIntervals - Interval list to be searched
     
               output: A list of all intervals whose OBSIDs are less 
                       than 50k
                       
               The program loops through the obsid intervals on the list 
               and if the value of the OBSID is less than 50k it appends
               that obsid interval to the output list
     
                usage: start_time = '2010:001'
                       stop_time  = '2010:014' 
                       filespec = <some file path> or " "
                       cmd_states = cmd_statesFetch(start_time, stop_time)
                       obsid_intervals = FindObsidIntervals(cmd_states, filespec)
     
                       Non_cti_intervals = CTIFilter(OBSIDIntervals)
                      
        """
        self.non_CTI_obs = []

        # check the OBSID of the interval. If it is 50k or greater, it's
        # a CTI observation and we want to remove those from our list
        for eachinterval in obsidinterval_list:
            if eachinterval[self.obsid] < 50000:
                self.non_CTI_obs.append(eachinterval)

        return self.non_CTI_obs
 
    #---------------------------------------------------------------------
    #
    #   pitch_filter 
    #
    #---------------------------------------------------------------------
    def pitch_filter(self, obsidinterval_list, pitchrangelist):
        """
        Given the output of FindObsidIntervals(cmd_states, ' ')
                      (which is a list of commanded states with the exposure
                      time tacked onto the end, and the OBSID start/stop
                      times at the beginning...one state per obsid.), find
                      all obsid intervals whose pitch is within the range 
                      specified and return those intervals as a list.
     
        input: expintervals = PitchFilter(OBSIDIntervals, [min_pitch, <max_pitch>])
     
                    Interval list to be searched
     
                    List containing the min and max pitchs
                    If the list contains only one pitch, return
                    all intervals whose pitch is greater than or
                    equal to that pitch
     
               output: A list of all intervals whose pitches are
                       within the min.max pitch inclusive.
     
                usage: start_time = '2010:001'
                       stop_time  = '2010:014' 
                       cmd_states = cmd_statesFetch(start_time, stop_time)
                       OBSIDIntervals = FindObsidIntervals(cmd_states, 'junk.dat')

          pitchintervals = PitchFilter(OBSIDIntervals, [min_pitch, <max_pitch>])
        """
        pitchlist = []

        # Capture the pitches for use in the filter.
        minpitch = pitchrangelist[0]

        if len(pitchrangelist) == 1:
            maxpitch = 180
        else:
            maxpitch = pitchrangelist[1]

        for eachinterval in obsidinterval_list:
            if maxpitch > eachinterval[self.pitch] >= minpitch:
                pitchlist.append(eachinterval)

        return pitchlist


    #--------------------------------------------------------------------------
    #
    #   CcdCountFilter
    #
    #--------------------------------------------------------------------------

    def ccdcountfilter(self,obsidinterval_list, ccdcountrangelist):
        """
        given the output of FindObsidIntervals(cmd_states, ' ')
        (which is a list of commanded states with the exposure
        time tacked onto the end, and the OBSID start/stop
        times at the beginning...one state per obsid.), find
        all obsid intervals whose CCD Count is within the range 
        specified and return those intervals as a list.

         input: expintervals = CcdCountFilter(OBSIDIntervals, 
                                              [min_count, <max_count>])

             OBSIDIntervals - Interval list to be searched

             [min_count, <max_count>] - List containing the min and max counts
                                        If the list contains only one count,
                                        return all intervals whose count
                                        exactly that. If there are two
                                        counts then return all the OBSIDS that
                                        have those counts (inclusive)

         output: A list of all intervals whose counts are
                 within the min.and max count inclusive.

         usage:

            CcdCountintervals = CcdCountFilter(OBSIDIntervals, [min_count, <max_count>])
        """
        ccdcountlist = []

        # Capture the counts for use in the filter.
        mincount = ccdcountrangelist[0]

        if len(ccdcountrangelist) == 1:
           maxcount = mincount
        else:
           maxcount = ccdcountrangelist[1]

        for eachinterval in obsidinterval_list:
            if maxcount >= eachinterval[self.ccd_count] >= mincount:
                ccdcountlist.append(eachinterval)

        return ccdcountlist


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
        acis_and_cti_only = []
        for eachobservation in obsidinterval_list:
            if eachobservation[self.in_focal_plane].startswith("ACIS-") or \
               eachobservation[self.obsid] >= 50000:
                acis_and_cti_only.append(eachobservation)
        return acis_and_cti_only


    #--------------------------------------------------------------------------
    #
    #   cti_only_filter
    #
    #--------------------------------------------------------------------------
    def cti_only_filter(self, obsidinterval_list):
        """
        This method will filter out any science observation from the 
        input obsid interval list.It keeps any observation that has an
        obsid of 50,000 or greater
        """
        cti_only = []
        for eachobservation in obsidinterval_list:
            if eachobservation[self.obsid] >= 50000 and \
               eachobservation[self.in_focal_plane] == "HRC-S":
                cti_only.append(eachobservation)
        return cti_only


    #--------------------------------------------------------------------------
    #
    #   fp_sens_filter
    #
    #--------------------------------------------------------------------------
    def fp_sens_filter(self, obsidinterval_list):
        """
        This method will return a list of science observations where 
        the value of the FP TEMP sensitive boolean is True
        It will check to be sure the length of the elements in the list
        is long enough to contain the boolean
        """
        fp_only = []
        for eachobservation in obsidinterval_list:
            if len(eachobservation) >= self.is_fp_sensitive and \
               eachobservation[self.is_fp_sensitive]:
                fp_only.append(eachobservation)
        return fp_only

    ######################################################################
    #
    #   GETS
    #
    ######################################################################

    #----------------------------------------------------------------------
    #
    # get_ccd_count 
    #
    #---------------------------------------------------------------------
    def get_ccd_count(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.ccd_count]

    #----------------------------------------------------------------------
    #
    # get_datestart 
    #
    #---------------------------------------------------------------------
    def get_datestart(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.datestart]

    #----------------------------------------------------------------------
    #
    # get_datestop 
    #
    #---------------------------------------------------------------------
    def get_datestop(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.datestop]

    #----------------------------------------------------------------------
    #
    # get_exptime 
    #
    #---------------------------------------------------------------------
    def get_exptime(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.exptime]

    #----------------------------------------------------------------------
    #
    # get_instrument 
    #
    #---------------------------------------------------------------------
    def get_instrument(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.in_focal_plane]

    #----------------------------------------------------------------------
    #
    # get_all_instruments
    #
    #---------------------------------------------------------------------
    def get_all_instruments(self, observations):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the instruments, append to a list  and return it
        """
        return [self.get_instrument(eachobs) for eachobs in observations]

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
    # get_obsid_list 
    #
    #---------------------------------------------------------------------
    def get_obsid_list(self, observation_list):
        """
        Given a list of obsids extracted by this class, loop through the 
        list and return a list of obsids.  The list contains ints
        """
        return [each_observation[self.obsid] 
                for each_observation in observation_list]


    #----------------------------------------------------------------------
    #
    # get_pitch 
    #
    #---------------------------------------------------------------------
    def get_pitch(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the obsid and return it
        """
        return observation[self.pitch]


    #----------------------------------------------------------------------
    #
    # get_sensitive 
    #
    #---------------------------------------------------------------------
    def get_sensitive(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the sim position and return it
        """
        return observation[self.is_fp_sensitive]


    #----------------------------------------------------------------------
    #
    # get_simpos 
    #
    #---------------------------------------------------------------------
    def get_simpos(self, observation):
        """
        Given a list element from the list of obsids extracted by this 
        class, extract the sim position and return it
        """
        return observation[self.simpos]

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


    #----------------------------------------------------------------------
    #
    # get_SI_mode_list 
    #
    #---------------------------------------------------------------------
    def get_si_mode_list(self, observation_list):
        """
        Given a list of obsids extracted by this class, loop through the 
        list and return a list of obsids.  The list contains ints
        """
        return [each_observation[self.si_mode] 
                for each_observation in observation_list]

    ######################################################################
    #
    #   SETS
    #
    ######################################################################
    #----------------------------------------------------------------------
    #
    # set_instrument
    #
    #---------------------------------------------------------------------
    def set_instrument(self, observation, instrument):
        """
        Given a list element from the list of obsids extracted by this 
        class, and one of the instruments from this list:

            ACIS-I
            ACIS-S
            HRC-I
            HRC-S

        Set the observation[self.in_focal_plane] value to the input 
        instrument
        """
        if instrument in ["ACIS-I", "ACIS-S", "HRC-I", "HRC-S"]:
            observation[self.in_focal_plane] = instrument
        return

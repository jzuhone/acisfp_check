#!/usr/bin/env python

"""
========================
dpa_check
========================

This code generates backstop load review outputs for checking the ACIS
focal plane temperature: FP_TEMP11. It also generates FP_TEMP11 model 
validation plots comparing predicted values to telemetry for the 
previous three weeks.
"""
from __future__ import print_function

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
matplotlib.use('Agg')

import glob
from Ska.Matplotlib import pointpair, \
    cxctime2plotdate
import Ska.engarchive.fetch_sci as fetch
from Chandra.Time import DateTime
from collections import defaultdict
import numpy as np
import xija
from acis_thermal_check import \
    ACISThermalCheck, \
    calc_off_nom_rolls, \
    get_options, \
    mylog, get_acis_limits
from acis_thermal_check.utils import \
    plot_two
import os
import sys
from kadi import events

# Import ACIS-specific observation extraction, filtering
# and attribute support routines.
from .acis_obs import ObsidFindFilter

model_path = os.path.abspath(os.path.dirname(__file__))
default_nopref_list = os.path.join(model_path, "FPS_NoPref.txt")

validation_limits = {'PITCH': [(1, 3.0), (99, 3.0)],
                     'TSCPOS': [(1, 2.5), (99, 2.5)]
                     }

hist_limit = [(-120.0, -112.0)]

def calc_model(model_spec, states, start, stop, T_acisfp=None,
               T_acisfp_times=None, dh_heater=None, dh_heater_times=None):
    """
    Create and run the Thermal Model for the Focal Plane temperature.

    Given: Model name (some string)
           Commanded States collected by make_week_predict
           start time
           Stop Time
           T_acisfp
           T_acisfp_times
    """
    # Start by creating the basic modeling framework for a XIJA Thermal Model
    # Give it some name, start and stop time and the name of the JSON file
    # --------------
    # THERMAL_MODEL
    # --------------
    model = xija.ThermalModel('acisfp', start=start, stop=stop, model_spec=model_spec)

    # create a numpy array of the start and stop times in the 
    # commanded states array fetched by make_week_predict
    times = np.array([states['tstart'], states['tstop']])

    # Now set any data values for the components of your model
    # What you have to push in manually are:
    #       any states information like vid_board or ccd count
    #       any pseudo-MSID's such as 1cbat (because the node does not reflect the MSID)
    #       any single value initializations you think ought to be made.
    #         - e.g. fptemp in this case since it's what you are looking for.

    # For each item in the Commanded States data structure which matters to us,
    # insert the values in the commanded states data structure into the model:
    # 
    # Telemetry doesn't have to be pushed in - the model handles that.But items in the states
    # array have to be manually shoved in.
    #
    # pitch comes from the telemetry
    model.comp['eclipse'].set_data(False)
    model.comp['sim_z'].set_data(states['simpos'], times)
    model.comp['fptemp'].set_data(T_acisfp, T_acisfp_times)

    if "roll" in model.comp:
        model.comp['roll'].set_data(calc_off_nom_rolls(states), times)

    for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
        model.comp[name].set_data(states[name], times)

    # model.comp Not in xija documentation
    model.comp['dh_heater'].set_data(dh_heater, dh_heater_times)

    #  "orbitephem0_x","orbitephem0_y","orbitephem0_z" are not in Commanded
    #  states  but they are in telemetry
    # We have to manually insert the aoattqt<x> valued because some items 
    # are sampled on 5 minute intervals and some are not.
    for i in range(1, 5):
        name = 'aoattqt{}'.format(i)
        state_name = 'q{}'.format(i)
        model.comp[name].set_data(states[state_name], times)

    # Get ephemeris from eng archive
    for axis in "xyz":
        name = 'orbitephem0_{}'.format(axis)
        msid = fetch.Msid(name, model.tstart - 2000, model.tstop + 2000)
        model.comp[name].set_data(msid.vals, msid.times)

    # Set some initial values. Some of these are superfluous. You do this because some
    # of these values may not be set at the actual start time. The telemetry might not have
    # reached that point
    model.comp['dpa_power'].set_data(0.0)
    model.comp['1cbat'].set_data(-53.0)
    model.comp['sim_px'].set_data(-120.0)

    # Create the model
    model.make()
    model.calc()

    return model

class ACISFPCheck(ACISThermalCheck):

    def __init__(self, msid, name, validation_limits,
                 hist_limit, calc_model, args, other_telem=None,
                 other_map=None):
        super(ACISFPCheck, self).__init__(msid, name, validation_limits,
                                          hist_limit, calc_model, args, 
                                          other_telem=other_telem, other_map=other_map)
        # Set specific limits for the focal plane model
        self.fp_sens_limit, self.acis_i_limit, self.acis_s_limit = get_acis_limits("fptemp")
        # Read in the FP Sensitive Nopref file and form nopref array from it.
        self.nopref_array = process_nopref_list(self.args.fps_nopref)
        # Create an empty observation list which will hold the results. This
        # list contains all ACIS and all CTI observations and will have the 
        # sensitivity boolean added.
        self.obs_with_sensitivity = []

    def _gather_perigee(self, run_start, load_start):
        # The first step is to build a list of all the perigee passages.
        perigee_passages = []

        # Gather the perigee passages that occur from the
        # beginning of the model run up to the start of the load
        # from kadi
        rzs = events.rad_zones.filter(run_start, load_start)
        for rz in rzs:
            perigee_passages.append([rz.start, rz.perigee])
        for rz in rzs:
            perigee_passages.append([rz.stop, rz.perigee])

        # We will get the load passages from the relevant CRM pad time file
        # (e.g. DO12143_CRM_Pad.txt) inside the bsdir directory
        # Each line is either an inbound  or outbound CTI
        #
        # The reason we are doing this is because we want to draw vertical
        # lines denoting each perigee passage on the plots
        #
        # Open the file
        crm_file_path = glob.glob(self.bsdir + "/*CRM*")[0]
        crm_file = open(crm_file_path, 'r')

        alines = crm_file.readlines()

        idx = None
        # Keep reading until you hit the last header line which is all "*"'s
        for i, aline in enumerate(alines):
            if len(aline) > 0 and aline[0] == "*":
                idx = i+1
                break

        if idx is None:
            raise RuntimeError("Couldn't find the end of the CRM Pad Time file header!")

        # Found the last line of the header. Start processing Perigee Passages

        # While there are still lines to be read
        for aline in alines[idx:]:
            # create an empty Peri. Passage instance location
            passage = []

            # split the CRM Pad Time file line read in and extract the
            # relevant information
            splitline = aline.split()
            passage.append(splitline[6])  # Radzone entry/exit
            passage.append(splitline[9])  # Perigee Passage time

            # append this passage to the passages list
            perigee_passages.append(passage)

        # Done with the CRM Pad Time file - close it
        crm_file.close()

        return perigee_passages

    def make_prediction_plots(self, outdir, states, times, temps, tstart):
        """
        Make output plots.

        :param outdir: the directory to which the products are written
        :param states: commanded states
        :param times: time stamps (sec) for temperature arrays
        :param temps: dict of temperatures
        :param tstart: load start time
        :rtype: dict of review information including plot file names

        This function assumes that ACIS Ops LR has been run and that the directory
        is populated with
        """

        # Gather perigee passages
        perigee_passages = self._gather_perigee(times[0], tstart)

        # Next we need to find all the ACIS-S observations within the start/stop
        # times so that we can paint those on the plots as well. We will get
        # those from the commanded states data structure called "states" 
        # 
        # Create an instance of the ObsidFindFilter class. This class provides
        # methods to extract obsid intervals from the commanded states based 
        # upon ACIS definitions and considerations. It also provides
        # various methods to filter the interval set based upon pitch range, 
        # number of ccd's, filter out CTI observations, and a range of exposure 
        # times.
        extract_and_filter = ObsidFindFilter()

        # extract the OBSID's from the commanded states. NOTE: this contains all
        # observations including CTI runs and HRC observations
        observation_intervals = extract_and_filter.find_obsid_intervals(states, None)

        # Filter out any HRC science observations BUT keep ACIS CTI observations
        acis_and_cti_obs = extract_and_filter.hrc_science_obs_filter(observation_intervals)

        # Ok so now you have all the ACIS observations collected. Also,
        # they have been identified by ObsidFindFilter as to who is in the focal plane.
        # Some apps, like this one, care about FP_TEMP sensitivity. Some do not. 
        # Since we do, then checking that and assigning a sensitivity must be done
        # 
        # Open the sensitive observation list file, which is found in the LR 
        # directory,
        # read each line, extract the OBSID and add that to a list.
        sensefile = open(os.path.join(self.bsdir, 'fp_sensitive.txt'), 'r')

        # The list_of_sensitive_obs is the list of all FP TEMP sensitive 
        # observations extracted from the file in the load review directory
        list_of_sensitive_obs = []

        # Get the list of FP_TEMP sensitive observations
        for eachline in sensefile.readlines()[1:]:
            # Extract the OBSID from each line; the obsid is in the second
            # column of this line. Append it to the list of FP_TEMP sensitive
            # observations
            #
            # NOTE: The obsid here is a STRING
            list_of_sensitive_obs.append(eachline.split()[1])
        # Done with the file - close it
        sensefile.close()

        # Now that you have the latest list of temperature sensitive OBSID's,
        # run through each observation and append either "*FP SENS*" or
        # "NOT FP SENS" to the end of each observation.

        # Now run through the observation list attribute of the ObsidFindFilter class
        for eachobservation in acis_and_cti_obs:
            # Pull the obsid from the observation and turn it into a string

            obsid = str(extract_and_filter.get_obsid(eachobservation))
            # See if it's in the sensitive list. If so, indicate whether or
            # not this observation is FP Senstive in the new list. This will be
            # used later in make_prediction_viols to catch violations.
            if obsid in list_of_sensitive_obs:
                eachobservation.append(True)
            else:
                eachobservation.append(False)

            self.obs_with_sensitivity.append(eachobservation)

        # create an empty dictionary called plots to contain the returned
        # figures, axes 1  and axes 2 of the plot_two call
        plots = {}

        # Start time of loads being reviewed expressed in units for plotdate()
        load_start = cxctime2plotdate([tstart])[0]
        # Value for left side of plots
        plot_start = max(load_start-2.0, cxctime2plotdate([times[0]])[0])

        w1 = None
        # Make plots of FPTEMP and pitch vs time, looping over
        # three different temperature ranges
        ylim = [(-120, -90), (-120, -119), (-120.0, -111.5)]
        ypos = [-110.0, -119.35, -116]
        capwidth = [2.0, 0.1, 0.4]
        textypos = [-108.0, -119.3, -115.7]
        fontsize = [12, 9, 9]
        for i in range(3):
            name = "%s_%d" % (self.name, i+1)
            plots[name] = plot_two(fig_id=i+1, x=times, y=temps[self.name],
                                   x2=pointpair(states['tstart'], states['tstop']),
                                   y2=pointpair(states['pitch']),
                                   title=self.msid.upper() + " (ACIS-I obs. in red; ACIS-S in green)",
                                   xlabel='Date', ylabel='Temperature (C)',
                                   ylabel2='Pitch (deg)', xmin=plot_start,
                                   ylim=ylim[i], ylim2=(40, 180),
                                   figsize=(12, 6), width=w1, load_start=load_start)
            # Draw a horizontal line indicating the FP Sensitive Observation Cut off
            plots[name]['ax'].axhline(self.fp_sens_limit, linestyle='--', color='red', linewidth=2.0)
            # Draw a horizontal line showing the ACIS-I -114 deg. C cutoff
            plots[name]['ax'].axhline(self.acis_i_limit, linestyle='--', color='purple', linewidth=1.0)
            # Draw a horizontal line showing the ACIS-S -112 deg. C cutoff
            plots[name]['ax'].axhline(self.acis_s_limit, linestyle='--', color='blue', linewidth=1.0)

            # Get the width of this plot to make the widths of all the
            # prediction plots the same
            if i == 0:
                w1, _ = plots[name]['fig'].get_size_inches()

            # Now plot any perigee passages that occur between xmin and xmax
            # for eachpassage in perigee_passages:
            paint_perigee(perigee_passages, states, plots, name)

            # Now draw horizontal lines on the plot running from start to stop
            # and label them with the Obsid
            draw_obsids(extract_and_filter, self.obs_with_sensitivity, self.nopref_array,
                        plots, name, ypos[i], ypos[i]-0.5*capwidth[i], ypos[i]+0.5*capwidth[i],
                        textypos[i], fontsize[i], plot_start)

            # Build the file name and output the plot to a file
            filename = self.msid.lower() + 'M%dtoM%d.png' % (-ylim[i][0], -ylim[i][1])
            outfile = os.path.join(outdir, filename)
            mylog.info('Writing plot file %s' % outfile)
            plots[name]['fig'].savefig(outfile)
            plots[name]['filename'] = filename

        self._make_state_plots(plots, 3, w1, plot_start,
                               outdir, states, load_start, figsize=(12, 6))

        return plots

    def make_prediction_viols(self, times, temps, load_start):
        """
        Find limit violations where predicted temperature is above the
        red minus margin.

        MSID is a global

        obs_with_sensitivity contains all ACIS and CTI observations
        and they have had FP sensitivity boolean added. In other words it's
        All ACIS and ECS runs.

        We will create a list of CTI-ONLY runs, and a list of all
        ACIS science runs without CTI runs. These two lists will
        be used to assess the categories of violations:

            1) Any ACIS-I observation that violates the -114 red limit
               is a violation and a load killer
                 - science_viols

            2) Any ACIS-S observation that violates the -112 red limit
               is a violation and a load killer
                 - science_viols

            3) Any ACIS FP TEMP sensitive obs that gets warmer than -118.7
               results in a "Preferences Not Met" indicator.
                 - fp_sense_viols

            3) Any CTI run that violates the -114 RED limit needs to be
               tracked and is NOT a load killer
                 - cti_viols

        """
        mylog.info('\nMAKE VIOLS Checking for limit violations in ' +
                   str(len(self.obs_with_sensitivity)) +
                   " total science observations")

        viols = {}

        # create an instance of ObsidFindFilter()
        eandf = ObsidFindFilter()

        # ------------------------------------------------------
        #   Create subsets of all the observations
        # ------------------------------------------------------
        # Just the CTI runs
        cti_only_obs = eandf.cti_only_filter(self.obs_with_sensitivity)
        # Now divide out observations by ACIS-S and ACIS-I
        ACIS_S_obs = eandf.get_all_specific_instrument(self.obs_with_sensitivity, "ACIS-S")
        ACIS_I_obs = eandf.get_all_specific_instrument(self.obs_with_sensitivity, "ACIS-I")

        # ACIS SCIENCE observations only  - no HRC; no CTI
        non_cti_obs = eandf.cti_filter(self.obs_with_sensitivity)

        # ACIS SCIENCE OBS which are sensitive to FP TEMP
        fp_sens_only_obs = eandf.fp_sens_filter(non_cti_obs)

        # Now if there is an Obsid in the FP Sense list that is ALSO in the
        # No Pref list - remove that Obsid from the FP Sense list:
        fp_sense_without_noprefs = []

        nopref_list = self.nopref_array['obsid'].tolist()

        for eachobs in fp_sens_only_obs:
            if str(eandf.get_obsid(eachobs)) not in nopref_list:
                fp_sense_without_noprefs.append(eachobs)

        temp = temps[self.name]

        # ------------------------------------
        #  CTI-ONLY, -118.7 violation check
        # ------------------------------------
        mylog.info('\n\nFP SENSITIVE -118.7 CTI ONLY violations')
        # Collect any -118.7C violations of CTI runs. These are not
        # load killers but need to be reported

        viols["cti"] = self.search_obsids_for_viols("CTI", self.fp_sens_limit, 
                                                    cti_only_obs, temp, times, 
                                                    load_start)

        # ------------------------------------------------------------
        #  FP TEMP sensitive observations; -118.7 violation check
        #     These are not load killers
        # ------------------------------------------------------------
        mylog.info('\n\nFP SENSITIVE -118.7 SCIENCE ONLY violations')

        viols["fp_sens"] = self.search_obsids_for_viols("FP-sensitive", self.fp_sens_limit,
                                                        fp_sense_without_noprefs, 
                                                        temp, times, load_start)
        # --------------------------------------------------------------
        #  ACIS-S - Collect any -112C violations of any non-CTI ACIS-S science run.
        #  These are load killers
        # --------------------------------------------------------------
        #
        mylog.info('\n\n ACIS-S -112 SCIENCE ONLY violations')

        viols["ACIS_S"] = self.search_obsids_for_viols("ACIS-S", self.acis_s_limit, 
                                                       ACIS_S_obs, temp, times, load_start)

        # --------------------------------------------------------------
        #  ACIS-I - Collect any -114C violations of any non-CTI ACIS science run.
        #  These are load killers
        # --------------------------------------------------------------
        #
        mylog.info('\n\n ACIS-I -114 SCIENCE ONLY violations')

        # Create the violation data structure.
        viols["ACIS_I"] = self.search_obsids_for_viols("ACIS-I", self.acis_i_limit, 
                                                       ACIS_I_obs, temp, times, load_start)

        return viols

    def search_obsids_for_viols(self, limit_name, limit, observations, temp, times,
                                load_start):
        """
        Given a planning limit and a list of observations, find those time intervals
        where the temp gets warmer than the planning limit and identify which 
        observations (if any) include part or all of those intervals.
        """
        # create an instance of ObsidFindFilter()
        eandf = ObsidFindFilter()

        viols_list = defaultdict(list)

        # Run through all observations
        for eachobs in observations:
            # Get the observation tstart and tstop times, and obsid
            obs_tstart = eandf.get_tstart(eachobs)
            obs_tstop = eandf.get_tstop(eachobs)
            # If the observation is in this load, let's look at it
            if obs_tstart > load_start:
                idxs = (times >= obs_tstart) & (times <= obs_tstop)
                viols = self._make_prediction_viols(times[idxs], temp[idxs],
                                                    load_start, limit, limit_name,
                                                    "max")
                # If we have flagged any violations, record the obsid for each
                # and add them to the list
                if len(viols) > 0:
                    for viol in viols:
                        viol["obsid"] = str(eandf.get_obsid(eachobs))
                    viols_list[self.msid] += viols

        # Finished - return the violations list
        return viols_list

#----------------------------------------------------------------------
#
#   paint_perigee
#
#----------------------------------------------------------------------
def paint_perigee(perigee_passages, states, plots, msid):
    """
    This function draws vertical dahsed lines for EEF, Perigee and XEF
    events in the load.EEF and XEF lines are black; Perigee is red.

    You supply the list of perigee passage events which are:
        Radzone Start/Stop time
        Perigee Passage time

        The states you created in main

        The dictionary of plots you created

        The MSID (in this case FP_TEMP) used to access the dictionary
    """
    #
    # Now plot any perigee passages that occur between xmin and xmax
    for eachpassage in perigee_passages:
        # The index [1] item is always the Perigee Passage time. Draw that line in red
        # If this line is between tstart and tstop then it needs to be drawn 
        # on the plot. otherwise ignore
        if states['tstop'][-1] >= DateTime(eachpassage[0]).secs >= states['tstart'][0]:
            # Have to convert this time into the new x axis time scale necessitated by SKA
            xpos = cxctime2plotdate([DateTime(eachpassage[0]).secs])

            # now plot the line.
            plots[msid]['ax'].vlines(xpos, -120, 20, linestyle=':', color='red', linewidth=2.0)

            # Plot the perigee passage time so long as it was specified in the CTI_report file
            if eachpassage[1] != "Not-within-load":
                perigee_time = cxctime2plotdate([DateTime(eachpassage[1]).secs])
                plots[msid]['ax'].vlines(perigee_time, -120, 20, linestyle=':', 
                                         color='black', linewidth=2.0)

def draw_obsids(extract_and_filter, 
                obs_with_sensitivity, 
                nopref_array,
                plots,
                msid, 
                ypos, 
                endcapstart, 
                endcapstop, 
                textypos, 
                fontsize,
                plot_start):
    """
    This functiion draws visual indicators across the top of the plot showing
    which observations are ACIS; whether they are ACIS-I (red) or ACIS-S (green)
    when they start and stop, and whether or not any observation is sensitive to the
    focal plane temperature.  The list of observations sensitive to the focal plane
    is found by reading the fp_sensitive.dat file that is located in each LR
    directory and is created by the LR script.

    No CTI measurements are indicated - only science runs.

    The caller supplies:
               Options from the Command line supplied by the user at runtime
               The instance of the ObsidFindFilter() class created 
               nopref rec array
               The plot dictionary
               The MSID used to index into the plot dictinary (superfluous but required)
               The position on the Y axis you'd like these indicators to appear
               The Y position of the bottom of the end caps
               The Y position of the top of the end caps
               The starting position of the OBSID number text
               The font size
               The left time of the plot in plot_date units
    """
    # Now run through the observation list attribute of the ObsidFindFilter class
    for eachobservation in obs_with_sensitivity:
        # extract the obsid

        obsid = str(extract_and_filter.get_obsid(eachobservation))

        # Color all ACIS-S observations green; all ACIS-I 
        # observations red
        if eachobservation[extract_and_filter.in_focal_plane] == "ACIS-I":
            color = 'red'
        else:
            color = 'green'

        # Add the sensitivity text if the observation was found to be FP TEMP
        # sensitive
        #
        # If the observation is FP sensitive in the first place............
        if eachobservation[extract_and_filter.is_fp_sensitive]:
            # extract the obsid for convenience
            this_obsid = extract_and_filter.get_obsid(eachobservation)

            # But if it's also in the nopref list AND the upcased CandS_status entry is "NO PREF"
            where_words = np.where(nopref_array['obsid'] == str(this_obsid))
            if str(this_obsid) in nopref_array['obsid'] and \
               nopref_array['CandS_status'][where_words[0][0]].upper()[:7] == 'NO_PREF':
                color = 'purple'
                obsid = obsid + ' * NO PREF *'
            else:
                obsid = obsid + ' * FP SENS *'

        # Convert the start and stop times into the Ska-required format
        obs_start = cxctime2plotdate([extract_and_filter.get_tstart(eachobservation)])
        obs_stop = cxctime2plotdate([extract_and_filter.get_tstop(eachobservation)])

        if eachobservation[extract_and_filter.in_focal_plane].startswith("ACIS-"):
            # For each ACIS Obsid, draw a horizontal line to show 
            # its start and stop
            plots[msid]['ax'].hlines(ypos, 
                                     obs_start, 
                                     obs_stop, 
                                     linestyle='-', 
                                     color=color, 
                                     linewidth=2.0)

            # Plot vertical end caps for each obsid to visually show start/stop
            plots[msid]['ax'].vlines(obs_start, 
                                     endcapstart, 
                                     endcapstop, 
                                     color=color, 
                                     linewidth=2.0)
            plots[msid]['ax'].vlines(obs_stop, 
                                     endcapstart, 
                                     endcapstop, 
                                     color=color, 
                                     linewidth=2.0)

            # Now print the obsid in the middle of the time span, 
            # above the line, and rotate 90 degrees. 

            obs_time = obs_start + (obs_stop - obs_start)/2
            if obs_time > plot_start:
                # Now plot the obsid.
                plots[msid]['ax'].text(obs_time, 
                                       textypos, 
                                       obsid,  
                                       color = color, 
                                       va='bottom', 
                                       ma='left', 
                                       rotation = 90, 
                                       fontsize = fontsize)


#------------------------------------------------------------------------------
#
#   process_nopref_list - read and store the list of observations which
#                         prefer a Cold and Stable focal plane, but
#                         which have had that desire waived.
#
#                          Input:  nopref file specification
#                         Output:  nopref_array (numpy array)
#
#------------------------------------------------------------------------------
def process_nopref_list(filespec=default_nopref_list):
    # Create the dtype for the nopref list rec array
    nopref_dtype = [('obsid', '|S10'), ('Seq_no', '|S10'), 
                    ('prop', '|S10'), ('CandS_status', '|S15')]

    # Create an empty nopref array
    nopref_array = np.array([], dtype=nopref_dtype)

    # Open the nopref list file
    nopreflist = open(filespec, "r")

    # For each line in the file......
    for line in nopreflist.readlines():

        # split the line out into it's constituent parts
        splitline = line.encode().split()

        # If the line is NOT a comment (i.e. does not start with a "#")
        if splitline[0][0] != '#':

            # Create the row
            one_entry = np.array(splitline, dtype=nopref_dtype)
            # Append this row onto the array
            nopref_array = np.append(nopref_array, one_entry, axis=0)

    # Done with the nopref list file
    nopreflist.close()

    # Now return the nopref array
    return nopref_array 

def main():
    opts = [("fps_nopref", {"default": default_nopref_list,
             "help": "Full path to the FP sensitive nopref file"})]
    args = get_options("acisfp", model_path, opts=opts)
    acisfp_check = ACISFPCheck("fptemp", "acisfp", validation_limits,
                               hist_limit, calc_model, args,
                               other_telem=['1dahtbon'],
                               other_map={'1dahtbon': 'dh_heater',
                                          "fptemp_11": "fptemp"})
    try:
        acisfp_check.run()
    except Exception as msg:
        if args.traceback:
            raise
        else:
            print("ERROR:", msg)
            sys.exit(1)

if __name__ == '__main__':
    main()

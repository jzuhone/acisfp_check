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

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
matplotlib.use('Agg')

from Ska.Matplotlib import pointpair, \
    cxctime2plotdate
from Chandra.Time import DateTime, secs2date
from collections import defaultdict
from acis_thermal_check import \
    ACISThermalCheck, \
    get_options, \
    mylog, get_acis_limits
from acis_thermal_check.utils import \
    plot_two, paint_perigee
import os
import sys
from kadi import events
from astropy.table import Table

#
# Import ACIS-specific observation extraction, filtering
# and attribute support routines.
#
from .acis_obs import ObsidFindFilter

model_path = os.path.abspath(os.path.dirname(__file__))


class ACISFPCheck(ACISThermalCheck):
    def __init__(self):
        valid_limits = {'PITCH': [(1, 3.0), (99, 3.0)],
                        'TSCPOS': [(1, 2.5), (99, 2.5)]
                        }
        hist_limit = [(-120.0, -109.0)]
        super(ACISFPCheck, self).__init__("fptemp", "acisfp", valid_limits,
                                          hist_limit,
                                          other_telem=['1dahtbon'],
                                          other_map={'1dahtbon': 'dh_heater',
                                                     "fptemp_11": "fptemp"})
        # Create an empty observation list which will hold the results. This
        # list contains all ACIS and all ECS observations and will have the
        # sensitivity boolean added.
        self.obs_with_sensitivity = []

    def _calc_model_supp(self, model, state_times, states, ephem, state0):
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

        # Input quaternions explicitly for calculating Earth heating
        for i in range(1, 5):
            name = 'aoattqt{}'.format(i)
            state_name = 'q{}'.format(i)
            model.comp[name].set_data(states[state_name], state_times)

        # Input ephemeris explicitly for calculating Earth heating
        for axis in "xyz":
            name = 'orbitephem0_{}'.format(axis)
            model.comp[name].set_data(ephem[name], model.times)

        # Set some initial values. You do this because some
        # of these values may not be set at the actual start time.
        model.comp['dpa_power'].set_data(0.0)
        model.comp['1cbat'].set_data(-53.0)
        model.comp['sim_px'].set_data(-120.0)

    def _make_state_plots(self, plots, num_figs, w1, plot_start,
                          states, load_start, figsize=(12, 6)):
        # Make a plot of ACIS CCDs and SIM-Z position
        plots['pow_sim'] = plot_two(
            fig_id=num_figs+1,
            title='ACIS CCDs and SIM-Z position',
            xlabel='Date',
            x=pointpair(states['tstart'], states['tstop']),
            y=pointpair(states['ccd_count']),
            yy=pointpair(states['fep_count']),
            ylabel='CCD/FEP Count',
            ylim=(-0.1, 6.1),
            xmin=plot_start,
            x2=pointpair(states['tstart'], states['tstop']),
            y2=pointpair(states['simpos']),
            ylabel2='SIM-Z (steps)',
            ylim2=(-105000, 105000),
            figsize=figsize, width=w1, load_start=load_start)
        plots['pow_sim']['ax'].lines[0].set_label('CCDs')
        plots['pow_sim']['ax'].lines[1].set_label('FEPs')
        plots['pow_sim']['ax'].legend(fancybox=True, framealpha=0.5, loc=2)
        paint_perigee(self.perigee_passages, states, plots, "pow_sim")
        plots['pow_sim']['filename'] = 'pow_sim.png'

        # Make a plot of off-nominal roll
        plots['roll_taco'] = plot_two(
            fig_id=num_figs+2,
            title='Off-Nominal Roll and Earth Solid Angle in Rad FOV',
            xlabel='Date',
            x=self.predict_model.times,
            y=self.predict_model.comp["roll"].dvals,
            xmin=plot_start,
            ylabel='Roll Angle (deg)',
            ylim=(-20.0, 20.0),
            x2=self.predict_model.times,
            y2=self.predict_model.comp['earthheat__fptemp'].dvals,
            ylabel2='Earth Solid Angle (sr)',
            ylim2=(1.0e-3, 1.0),
            figsize=figsize, width=w1, load_start=load_start)
        plots['roll_taco']['ax2'].set_yscale("log")
        paint_perigee(self.perigee_passages, states, plots, "roll_taco")
        plots['roll_taco']['filename'] = 'roll_taco.png'

    def make_prediction_plots(self, outdir, states, temps, load_start):
        """
        Make plots of the thermal prediction as well as associated 
        commanded states.

        Parameters
        ----------
        outdir : string
            The path to the output directory.
        states : NumPy record array
            Commanded states
        temps : dict of NumPy arrays
            Dictionary of temperature arrays
        load_start : float
            The start time of the load in seconds from the beginning of the
            mission.
        """

        times = self.predict_model.times

        # Gather perigee passages
        self._gather_perigee(times[0], load_start+86400.0)

        # Next we need to find all the ACIS-S observations within the start/stop
        # times so that we can paint those on the plots as well. We will get
        # those from the commanded states data structure called "states" 
        # 
        # Create an instance of the ObsidFindFilter class. This class provides
        # methods to extract obsid intervals from the commanded states based 
        # upon ACIS definitions and considerations. It also provides
        # various methods to filter the interval set based upon pitch range, 
        # number of ccd's, filter out ECS observations, and a range of exposure 
        # times.
        extract_and_filter = ObsidFindFilter()

        # extract the OBSID's from the commanded states. NOTE: this contains all
        # observations including ECS runs and HRC observations
        observation_intervals = extract_and_filter.find_obsid_intervals(states, None)

        # Filter out any HRC science observations BUT keep ACIS ECS observations
        acis_and_ecs_obs = extract_and_filter.hrc_science_obs_filter(observation_intervals)

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
        for eachobservation in acis_and_ecs_obs:
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
        load_start = cxctime2plotdate([load_start])[0]
        # Value for left side of plots
        plot_start = max(load_start-2.0, cxctime2plotdate([times[0]])[0])

        w1 = None
        # Make plots of FPTEMP and pitch vs time, looping over
        # three different temperature ranges
        ylim = [(-120, -90), (-120, -119), (-120.0, -107.5)]
        ypos = [-110.0, -119.35, -116]
        capwidth = [2.0, 0.1, 0.4]
        textypos = [-108.0, -119.3, -115.7]
        fontsize = [12, 9, 9]
        for i in range(3):
            name = f"{self.name}_{i+1}"
            plots[name] = plot_two(fig_id=i+1, x=times, y=temps[self.name],
                                   x2=self.predict_model.times,
                                   y2=self.predict_model.comp["pitch"].mvals,
                                   title=f"{self.msid.upper()} (ACIS-I in red; ACIS-S in green; ECS in blue)",
                                   xlabel='Date', ylabel='Temperature (C)',
                                   ylabel2='Pitch (deg)', xmin=plot_start,
                                   ylim=ylim[i], ylim2=(40, 180), 
                                   figsize=(12, 7.142857142857142),
                                   width=w1, load_start=load_start)
            # Draw a horizontal line indicating the FP Sensitive Observation Cut off
            plots[name]['ax'].axhline(self.cold_ecs_limit, linestyle='--', 
                                      color='dodgerblue', linewidth=2.0)
            # Draw a horizontal line showing the ACIS-I -114 deg. C cutoff
            plots[name]['ax'].axhline(self.acis_i_limit, linestyle='--', 
                                      color='purple', linewidth=2.0)
            # Draw a horizontal line showing the ACIS-S -112 deg. C cutoff
            plots[name]['ax'].axhline(self.acis_s_limit, linestyle='--', 
                                      color='blue', linewidth=2.0)
            # Draw a horizontal line showing the ACIS-S -109 deg. C cutoff
            plots[name]['ax'].axhline(self.acis_hot_limit, linestyle='--', 
                                      color='red', linewidth=2.0)
            # Get the width of this plot to make the widths of all the
            # prediction plots the same
            if i == 0:
                w1, _ = plots[name]['fig'].get_size_inches()

            # Now plot any perigee passages that occur between xmin and xmax
            # for eachpassage in perigee_passages:
            paint_perigee(self.perigee_passages, states, plots, name)

            # Now draw horizontal lines on the plot running from start to stop
            # and label them with the Obsid
            draw_obsids(extract_and_filter, self.obs_with_sensitivity, 
                        plots, name, ypos[i], ypos[i]-0.5*capwidth[i], 
                        ypos[i]+0.5*capwidth[i], textypos[i], 
                        fontsize[i], plot_start)

            # Build the file name
            filename = self.msid.lower() + 'M%dtoM%d.png' % (-ylim[i][0], -ylim[i][1])
            plots[name]['filename'] = filename

        self._make_state_plots(plots, 3, w1, plot_start,
                               states, load_start, 
                               figsize=(12, 6))

        plots['default'] = plots[f"{self.name}_3"]

        # This call allows the specific check tool
        # to customize plots after the fact
        self.custom_prediction_plots(plots)

        # Now write all of the plots after possible
        # customizations have been made
        for key in plots:
            if key != self.msid:
                outfile = os.path.join(outdir, plots[key]['filename'])
                mylog.info('Writing plot file %s' % outfile)
                plots[key]['fig'].savefig(outfile)

        return plots

    def make_prediction_viols(self, temps, load_start):
        """
        Find limit violations where predicted temperature is above the
        red minus margin.

        MSID is a global

        obs_with_sensitivity contains all ACIS and ECS observations
        and they have had FP sensitivity boolean added. In other words it's
        All ACIS and ECS runs.

        We will create a list of ECS-ONLY runs, and a list of all
        ACIS science runs without ECS runs. These two lists will
        be used to assess the categories of violations:

            1) Any ACIS-I observation that violates the -114 red limit
               is a violation and a load killer
                 - science_viols

            2) Any ACIS-S observation that violates the -112 red limit
               is a violation and a load killer
                 - science_viols

        """
        times = self.predict_model.times

        mylog.info('\nMAKE VIOLS Checking for limit violations in ' +
                   str(len(self.obs_with_sensitivity)) +
                   " total science observations")

        viols = {}

        # create an instance of ObsidFindFilter()
        eandf = ObsidFindFilter()

        # ------------------------------------------------------
        #   Create subsets of all the observations
        # ------------------------------------------------------
        # Now divide out observations by ACIS-S and ACIS-I
        ACIS_S_obs = eandf.get_all_specific_instrument(
            self.obs_with_sensitivity, "ACIS-S")
        ACIS_I_obs = eandf.get_all_specific_instrument(
            self.obs_with_sensitivity, "ACIS-I")

        # ACIS SCIENCE observations only  - no HRC; no ECS
        #non_ecs_obs = eandf.ecs_filter(self.obs_with_sensitivity)

        sci_ecs_obs = eandf.ecs_only_filter(self.obs_with_sensitivity)

        # ACIS SCIENCE OBS which are sensitive to FP TEMP
        #fp_sens_only_obs = eandf.fp_sens_filter(non_ecs_obs)

        temp = temps[self.name]

        # ------------------------------------------------------------
        # Science Orbit ECS -119.5 violations; -119.5 violation check
        # ------------------------------------------------------------
        mylog.info(f'\n\nScience Orbit ECS ({self.cold_ecs_limit} C) violations')

        ecs_viols = self.search_obsids_for_viols("Science Orbit ECS",
            self.cold_ecs_limit, sci_ecs_obs, temp, times, load_start)

        viols["ecs"] = {"name": f"Science Orbit ECS ({self.cold_ecs_limit} C)",
                        "type": "Min",
                        "values": ecs_viols}

        # ------------------------------------------------------------
        # ACIS-S - Collect any -111 C violations of any non-ECS ACIS-S
        # science run. These are load killers
        # ------------------------------------------------------------
        #
        mylog.info(f'\n\n ACIS-S Science ({self.acis_s_limit} C) violations')

        acis_s_viols = self.search_obsids_for_viols("ACIS-S",
            self.acis_s_limit, ACIS_S_obs, temp, times, load_start)
        viols["ACIS_S"] = {"name": f"ACIS-S ({self.acis_s_limit} C)",
                           "type": "Max",
                           "values": acis_s_viols}

        # ------------------------------------------------------------
        # ACIS-I - Collect any -112 C violations of any non-ECS ACIS-I
        # science run. These are load killers
        # ------------------------------------------------------------
        #
        mylog.info(f'\n\n ACIS-I Science ({self.acis_i_limit} C) violations')

        # Create the violation data structure.
        acis_i_viols = self.search_obsids_for_viols("ACIS-I",
            self.acis_i_limit, ACIS_I_obs, temp, times, load_start)

        viols["ACIS_I"] = {"name": f"ACIS-I ({self.acis_i_limit} C)",
                           "type": "Max",
                           "values": acis_i_viols}

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

        viols_list = []

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
                    viols_list += viols

        # Finished - return the violations list
        return viols_list

    def write_temps(self, outdir, times, temps):
        """
        Write the states record array to the file "temperatures.dat"
        and the Earth solid angles to "earth_solid_angles.dat".

        Parameters
        ----------
        outdir : string
            The directory the file will be written to.
        times : NumPy array
            Times in seconds from the start of the mission
        temps : NumPy array
            Temperatures in Celsius
        """
        super(ACISFPCheck, self).write_temps(outdir, times, temps)
        outfile = os.path.join(outdir, 'earth_solid_angles.dat')
        mylog.info('Writing Earth solid angles to %s' % outfile)
        e = self.predict_model.comp['earthheat__fptemp'].dvals
        efov_table = Table([times, secs2date(times), e],
                           names=['time', 'date', 'earth_solid_angle'],
                           copy=False)
        efov_table['time'].format = '%.2f'
        efov_table['earth_solid_angle'].format = '%.3e'
        efov_table.write(outfile, format='ascii', delimiter='\t', overwrite=True)


def draw_obsids(extract_and_filter, 
                obs_with_sensitivity, 
                plots,
                msid, 
                ypos, 
                endcapstart, 
                endcapstop, 
                textypos, 
                fontsize,
                plot_start):
    """
    This function draws visual indicators across the top of the plot showing
    which observations are ACIS; whether they are ACIS-I (red), ACIS-S (green),
    or ECS (blue); when they start and stop; and whether or not any observation 
    is sensitive to the focal plane temperature. The list of observations sensitive 
    to the focal plane is found by reading the fp_sensitive.dat file that is 
    located in each LR directory and is created by the LR script.

    The caller supplies:
               Options from the Command line supplied by the user at runtime
               The instance of the ObsidFindFilter() class created 
               The plot dictionary
               The MSID used to index into the plot dictionary (superfluous but required)
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

        obsid = extract_and_filter.get_obsid(eachobservation)
        in_fp = eachobservation[extract_and_filter.in_focal_plane]

        if obsid > 60000:
            # ECS observations during the science orbit are colored blue
            color = 'blue'
        else:
            # Color all ACIS-S observations green; all ACIS-I
            # observations red
            if in_fp == "ACIS-I":
                color = 'red'
            else:
                color = 'green'

        obsid_txt = str(obsid)
        # If this is an ECS measurement in the science orbit mark
        # it as such
        if obsid > 60000:
            obsid_txt += " (ECS)"

        # Convert the start and stop times into the Ska-required format
        obs_start = cxctime2plotdate([extract_and_filter.get_tstart(eachobservation)])
        obs_stop = cxctime2plotdate([extract_and_filter.get_tstop(eachobservation)])

        if in_fp.startswith("ACIS-") or obsid > 60000:
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
                                       obsid_txt,  
                                       color=color, 
                                       va='bottom', 
                                       ma='left', 
                                       rotation=90, 
                                       fontsize=fontsize)


def main():
    args = get_options("acisfp", model_path)
    acisfp_check = ACISFPCheck()
    try:
        acisfp_check.run(args)
    except Exception as msg:
        if args.traceback:
            raise
        else:
            print("ERROR:", msg)
            sys.exit(1)


if __name__ == '__main__':
    main()

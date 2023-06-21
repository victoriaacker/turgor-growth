# -*- coding: latin-1 -*-

import os
import logging
import datetime

import pandas as pd

from turgorgrowth import simulation as turgorgrowth_simulation, converter as turgorgrowth_converter, \
    postprocessing as turgorgrowth_postprocessing, tools as turgorgrowth_tools

"""
    main
    ~~~~

    An example to show how to run the model Turgor-Growth, compute the post-processing, and generate the plots for validation.

    Before running this script, you must first install model Turgor-Growth (see `README.md` at the root directory of the project).
    Then, you can run this script with the command `python`.

    :license: CeCILL-C, see LICENSE for details.
    
"""

# -----------------------------------------------------
#           CONFIGURATION OF THE SIMULATION           -
# -----------------------------------------------------

# ---------- INPUTS CONFIGURATION ----------

# Path of the directory which contains the inputs of the model
INPUTS_DIRPATH = 'inputs_simpleplant'
##INPUTS_DIRPATH = 'inputs'


# Name of the CSV files which describe the initial state of the system
HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
##ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state-growth.csv'
ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state-growth - Copie.csv'
ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
SOIL_INITIAL_STATE_FILENAME = 'soil_initial_state.csv'

# Name of the CSV files which contain the photosynthesis and senescence forcings
ELEMENTS_FORCINGS_FILENAME = 'elements_forcings-t=100-Ete.csv'
##ELEMENTS_FORCINGS_FILENAME = 'elements_forcings-t=999-Ete.csv'
##ELEMENTS_FORCINGS_FILENAME = 'elements_forcings.csv'
HIDDENZONES_FORCINGS_FILENAME = 'hiddenzones_forcings.csv'
ORGANS_FORCINGS_FILENAME = 'organs_forcings.csv'
SOIL_FORCINGS_FILENAME = 'soil_forcings.csv'

# Name of the CSV files which contains the meteo data
METEO_INPUTS_FILENAME = 'soil_forcings.csv'

# ---------- OUTPUTS CONFIGURATION ----------

# Path of the directory where to write the outputs of the model
OUTPUTS_DIRPATH = 'outputs'

# Name of the CSV files which will contain the outputs of the model
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
SOIL_OUTPUTS_FILENAME = 'soil_outputs.csv'

# ---------- POSTPROCESSING CONFIGURATION ----------

# Path of the directory where to write the postprocessing of the model
POSTPROCESSING_DIRPATH = 'postprocessing'

# Name of the CSV files which will contain the postprocessing of the model
AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
SOIL_POSTPROCESSING_FILENAME = 'soil_postprocessing.csv'

# ---------- GRAPHS CONFIGURATION ----------

# Path of the directory where to save the generated graphs
GRAPHS_DIRPATH = 'graphs'

# ---------- SIMULATION PARAMETERS ----------

# Start time of the simulation
START_TIME = 0

# Length of the simulation (in hours)
SIMULATION_LENGTH = 100

# Time step of the simulation (in hours)
TIME_STEP = 1

# Do run the simulation?
RUN_SIMU = True

# Do run the postprocessing?
RUN_POSTPROCESSING = True

# Do generate the graphs?
GENERATE_GRAPHS = True

# Do log the execution?
LOG_EXECUTION = False

# Config file path for logging
LOGGING_CONFIG_FILEPATH = 'logging.json'

# Do interpolate the forcings?
INTERPOLATE_FORCINGS = False

# -----------------------------------------------------
#           RUN OF THE SIMULATION                     -
# -----------------------------------------------------

# Warning: you should not modify the following code excepting you really know what you are doing

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# Precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 6


def force_inputs(t, population, elements_data_grouped, hiddenzones_data_grouped, soil_data_grouped):
    """Force input data of the population at `t` from input grouped dataframes"""
    for plant in population.plants:
        for axis in plant.axes:
            group_soil = soil_data_grouped.get_group((t, plant.index, axis.label))
            soil_data_to_use = group_soil.loc[group_soil.first_valid_index(),
                                                                  group_soil.columns.intersection(turgorgrowth_simulation.Simulation.ORGANS_STATE)].dropna().to_dict()
            axis.xylem.__dict__.update(soil_data_to_use)
            for phytomer in axis.phytomers:
                if phytomer.hiddenzone:
                    group_hiddenzone = hiddenzones_data_grouped.get_group((t, plant.index, axis.label, phytomer.index))
                    hiddenzone_data_to_use = group_hiddenzone.loc[group_hiddenzone.first_valid_index(),
                                                                  group_hiddenzone.columns.intersection(turgorgrowth_simulation.Simulation.HIDDENZONE_STATE)].dropna().to_dict()
                    phytomer.hiddenzone.__dict__.update(hiddenzone_data_to_use)
                for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    for element in (organ.exposed_element, organ.enclosed_element):
                        if element is None:
                            continue
                        # Element
                        group_element = elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                        element_data_to_use = group_element.loc[group_element.first_valid_index(),
                                                                group_element.columns.intersection(turgorgrowth_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
                        element.__dict__.update(element_data_to_use)


if RUN_SIMU:

    print('Prepare the simulation...')

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    if LOG_EXECUTION:
        # Setup the logging
        turgorgrowth_tools.setup_logging(config_filepath=LOGGING_CONFIG_FILEPATH, level=logging.DEBUG, log_model=True, log_compartments=True, log_derivatives=True, remove_old_logs=True)

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME, ORGANS_INITIAL_STATE_FILENAME, SOIL_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, mapping_topology = turgorgrowth_converter.from_dataframes(inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                                                inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                                                                                inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME],
                                                                                inputs_dataframes[SOIL_INITIAL_STATE_FILENAME])
    ## + soils

    # Set simulation parameters for interpolation
    if INTERPOLATE_FORCINGS:
        elements_forcings_delta_t = time_step_seconds
        hiddenzone_forcings_delta_t = time_step_seconds
        soil_forcings_delta_t = time_step_seconds
    else:
        elements_forcings_delta_t = None
        hiddenzone_forcings_delta_t = None
        soil_forcings_delta_t = None

    # Create the simulation
    simulation_ = turgorgrowth_simulation.Simulation(delta_t=time_step_seconds, interpolate_forcings=INTERPOLATE_FORCINGS, elements_forcings_delta_t=elements_forcings_delta_t,
                                                     hiddenzone_forcings_delta_t=hiddenzone_forcings_delta_t, soil_forcings_delta_t=soil_forcings_delta_t)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, mapping_topology)

    # Read forcings from CSV files, create dataframes, and group the dataframes by object index
    # Element
    elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_FORCINGS_FILENAME)
    elements_data_df = pd.read_csv(elements_data_filepath)
    elements_data_grouped = elements_data_df.groupby(turgorgrowth_simulation.Simulation.ELEMENTS_T_INDEXES)
    # HiddenZone
    hiddenzones_data_filepath = os.path.join(INPUTS_DIRPATH, HIDDENZONES_FORCINGS_FILENAME)
    hiddenzones_data_df = pd.read_csv(hiddenzones_data_filepath)
    hiddenzones_data_grouped = hiddenzones_data_df.groupby(turgorgrowth_simulation.Simulation.HIDDENZONE_T_INDEXES)
    # Soil
    soil_data_filepath = os.path.join(INPUTS_DIRPATH, SOIL_FORCINGS_FILENAME)
    soil_data_df = pd.read_csv(soil_data_filepath)
    soil_data_grouped = soil_data_df.groupby(turgorgrowth_simulation.Simulation.AXES_T_INDEXES)
    force_inputs(0, population, elements_data_grouped, hiddenzones_data_grouped, soil_data_grouped)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Reinitialize the simulation from forced population
    meteo = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_INPUTS_FILENAME), index_col='t')
    SRWC = meteo.loc[time_grid[0], ['SRWC']]
    simulation_.initialize(population, mapping_topology, SRWC=SRWC)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soil_outputs_df_list = []

    print('Prepare the simulation... DONE!')

    print('Run the simulation...')
    current_time_of_the_system = datetime.datetime.now()

    for t in time_grid:

        if t > 0:
            # Run the model ; the population is internally updated by the model
            print('t = {}'.format(t))
            simulation_.run()

        # Convert the model outputs to dataframes
        ##_, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soil_outputs_df = turgorgrowth_converter.to_dataframes(simulation_.population, simulation_.soil)
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df= turgorgrowth_converter.to_dataframes(simulation_.population)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list),(hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list), (organs_outputs_df, organs_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if 0 < t < SIMULATION_LENGTH:
            # Force the population
            force_inputs(t, population, elements_data_grouped, hiddenzones_data_grouped, soil_data_grouped)
            # Reinitialize the simulation from forced population
            ##simulation_.initialize(population, mapping_topology)
            SRWC = meteo.loc[t, ['SRWC']]
            ##simulation_.initialize(population, mapping_topology, soils, SRWC= SRWC)
            simulation_.initialize(population, mapping_topology, SRWC= SRWC)


    print('Run the simulation... DONE!')

    execution_time = datetime.datetime.now() - current_time_of_the_system
    print('Simulation run in ', execution_time)

    print('Total RHS evaluations: {}'.format(simulation_.nfev_total))

    print('Write the outputs to CSV files...')

    outputs_df_dict = {}
    for outputs_df_list, outputs_filename in ((axes_outputs_df_list, AXES_OUTPUTS_FILENAME),
                                              (hiddenzones_outputs_df_list, HIDDENZONES_OUTPUTS_FILENAME),
                                              (elements_outputs_df_list, ELEMENTS_OUTPUTS_FILENAME),
                                              (organs_outputs_df_list, ORGANS_OUTPUTS_FILENAME)):
                                              ##,(soil_outputs_df_list, SOIL_OUTPUTS_FILENAME)):
        outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df.to_csv(outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
        outputs_file_basename = outputs_filename.split('.')[0]
        outputs_df_dict[outputs_file_basename] = outputs_df

    print('Write the outputs to CSV files... DONE!')

if RUN_POSTPROCESSING:

    if not RUN_SIMU:

        print('Retrieve outputs dataframes from precedent simulation run...')

        outputs_df_dict = {}

        for outputs_filename in (AXES_OUTPUTS_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, ELEMENTS_OUTPUTS_FILENAME, ORGANS_OUTPUTS_FILENAME):
            outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
            outputs_df = pd.read_csv(outputs_filepath)
            outputs_file_basename = outputs_filename.split('.')[0]
            outputs_df_dict[outputs_file_basename] = outputs_df

        time_grid = list(outputs_df_dict.values())[0].t
        delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        print('Retrieve outputs dataframes from precedent simulation run... DONE!')
    else:
        delta_t = simulation_.delta_t

    print('Compute the post-processing...')

    axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
    hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
    elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
    organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
    # soil_postprocessing_file_basename = SOIL_POSTPROCESSING_FILENAME.split('.')[0]

    postprocessing_df_dict = {}

    (postprocessing_df_dict[axes_postprocessing_file_basename],
     postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
     postprocessing_df_dict[elements_postprocessing_file_basename],
     postprocessing_df_dict[organs_postprocessing_file_basename]) \
            = turgorgrowth_postprocessing.postprocessing(axes_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                        hiddenzones_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                        elements_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                        organs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                        delta_t=delta_t)

    print('Compute the post-processing... DONE!')

    print('Write the postprocessing to CSV files...')


    for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILENAME), (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILENAME), (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILENAME), (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILENAME)):
        postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
        postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

    print('Write the postprocessing to CSV files... DONE!')

if GENERATE_GRAPHS:

    if not RUN_POSTPROCESSING:

        print('Retrieve last computed post-processing dataframes...')

        postprocessing_df_dict = {}

        for postprocessing_filename in (AXES_POSTPROCESSING_FILENAME, HIDDENZONES_POSTPROCESSING_FILENAME, ELEMENTS_POSTPROCESSING_FILENAME, ORGANS_POSTPROCESSING_FILENAME):
            postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_df = pd.read_csv(postprocessing_filepath)
            postprocessing_file_basename = postprocessing_filename.split('.')[0]
            postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

        print('Retrieve last computed post-processing dataframes... DONE!')

    print('Generate graphs for validation...')

    turgorgrowth_postprocessing.generate_graphs(axes_df=postprocessing_df_dict[AXES_POSTPROCESSING_FILENAME.split('.')[0]],
                                                hiddenzones_df=postprocessing_df_dict[HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]],
                                                elements_df=postprocessing_df_dict[ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]],
                                                organs_df=postprocessing_df_dict[ORGANS_POSTPROCESSING_FILENAME.split('.')[0]],
                                                graphs_dirpath=GRAPHS_DIRPATH)

    print('Generate graphs for validation... DONE!')

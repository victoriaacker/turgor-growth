"""
    test_turgorgrowth
    ~~~~~~~~~~~~

    Test:

        * the run of a simulation with forcings,

    You must first install model Turgor-Growth before running this script with the command `python`. See `README.md` at the
    root directory of the project.

    CSV files must contain only ASCII characters and ',' as separator.

    :license: CeCILL-C, see LICENSE for details.
    
"""

import os

import pandas as pd

from turgorgrowth import simulation as turgorgrowth_simulation, converter as turgorgrowth_converter, \
    tools as turgorgrowth_tools


# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# Precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 6

# the precision to use for quantitative comparison test
PRECISION = 4


def force_inputs(t, population, elements_data_grouped, hiddenzones_data_grouped):
    """Force input data of the population at `t` from input grouped dataframes"""
    for plant in population.plants:
        for axis in plant.axes:
            for phytomer in axis.phytomers:
                if phytomer.hiddenzone:
                    group_hiddenzone = hiddenzones_data_grouped.get_group((t, plant.index, axis.label, phytomer.index))
                    hiddenzone_data_to_use = group_hiddenzone.loc[group_hiddenzone.first_valid_index(),
                                                                  group_hiddenzone.columns.intersection(turgorgrowth_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
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


def test_simulation_run(overwrite_desired_data=False):
    """Test the run of a simulation, without interpolation of the forcings."""

    TEST_DIR_PATH = 'simulation_run'

    # Inputs of the test
    INPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'inputs')
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    ELEMENTS_FORCINGS_FILENAME = 'elements_forcings.csv'
    HIDDENZONES_FORCINGS_FILENAME = 'hiddenzones_forcings.csv'

    # Outputs of the test
    OUTPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'outputs')
    DESIRED_HIDDENZONES_OUTPUTS_FILENAME = 'desired_hiddenzones_outputs.csv'
    DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
    ACTUAL_HIDDENZONES_OUTPUTS_FILENAME = 'actual_hiddenzones_outputs.csv'
    ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'

    # Simulation parameters
    START_TIME = 0
    SIMULATION_LENGTH = 50
    TIME_STEP = 1

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, mapping_topology = turgorgrowth_converter.from_dataframes(inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                                          inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME])

    # Create the simulation
    simulation_ = turgorgrowth_simulation.Simulation(delta_t=time_step_seconds)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, mapping_topology)

    # Read forcings from CSV files, create dataframes, and group the dataframes by object index
    elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_FORCINGS_FILENAME)
    elements_data_df = pd.read_csv(elements_data_filepath)
    elements_data_grouped = elements_data_df.groupby(turgorgrowth_simulation.Simulation.ELEMENTS_T_INDEXES)
    hiddenzones_data_filepath = os.path.join(INPUTS_DIRPATH, HIDDENZONES_FORCINGS_FILENAME)
    hiddenzones_data_df = pd.read_csv(hiddenzones_data_filepath)
    hiddenzones_data_grouped = hiddenzones_data_df.groupby(turgorgrowth_simulation.Simulation.HIDDENZONE_T_INDEXES)
    force_inputs(0, population, elements_data_grouped, hiddenzones_data_grouped)

    # Reinitialize the simulation from forced population and soils
    simulation_.initialize(population, mapping_topology)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []

    for t in time_grid:

        if t > 0:
            print(t)
            # Run the model of CN exchanges ; the population is internally updated by the model
            simulation_.run()

        # Convert the model outputs to dataframes
        hiddenzones_outputs_df, elements_outputs_df = turgorgrowth_converter.to_dataframes(simulation_.population)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if 0 < t < SIMULATION_LENGTH:
            # Force the population
            force_inputs(t, population, elements_data_grouped, hiddenzones_data_grouped)
            # Reinitialize the simulation from forced population
            simulation_.initialize(population, mapping_topology)

    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)
    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename,
         state_variables_names) \
            in ((hiddenzones_outputs_df_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME, ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                 turgorgrowth_simulation.Simulation.HIDDENZONE_T_INDEXES + turgorgrowth_simulation.Simulation.HIDDENZONE_STATE),
                (elements_outputs_df_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                 turgorgrowth_simulation.Simulation.ELEMENTS_T_INDEXES + turgorgrowth_simulation.Simulation.ELEMENTS_STATE)):
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df = outputs_df.loc[:, state_variables_names]  # compare only the values of the compartments
        turgorgrowth_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename,
                                                     actual_outputs_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


if __name__ == '__main__':
    test_simulation_run(overwrite_desired_data=False)

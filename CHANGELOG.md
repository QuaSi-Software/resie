# Changes made, features added and bugs fixed
In general the development follows the [semantic versioning](https://semver.org/) scheme. As the simulation is parameterized by the project file, the structure and meaning of this file represents the API of the simulation engine. Evaluating compatability then follows required changes in the project file. This is complicated by the implementation of the individual energy system components, which happens within the framework but is also reasonably independant of the simulation model. A change in the implementation may have big effects on simulation results without requiring any change in the project file. It is up to the developers to find a workable middle ground.

## Pre-1.0-releases
As per the definition of semantic versioning and the reality of early development, in versions prior to 1.0.0 any release might break compatability. To alleviate this somewhat, the meaning of major-minor-patch is "downshifted" to zero-major-minor. However some breaking changes may slip beneath notice.

### Version 0.8.5
* Fix a problem with general bounded supply and sink implementations not considering temperatures correctly. This was mostly an issue with direct 1-to-1 connections between components as busses do consider this already.

### Version 0.8.4
* Add geothermal probe as new component. It acts like a storage and there are many precalculated input parameters for a variety of probe field configurations. Both the simplified and the extended version were successfully validated.
* Add new szenario with geothermal probe and heat pump.
* Add geothermal collector files. This component is not finished yet and needs some rework.
* Add standalone julia code to create figures from profile data
* Add julia packages Roots and Plots
* Add function plot_optional_figures(), available to include in every component, that is called after initialize!() to create additional plots if the component offers them. Added optional flags in input file to control the auxiliary plots: auxiliary_plots, auxiliary_plots_path, auxiliary_plots_formats
* Bugfix in Resie.jl to call set_time_step ahead of the initialize!() of each component

### Version 0.8.3
* Tweak examples to match the description in the documentation
* Update profiles used in examples and move the files into subfolders
* Update installation and quick-use instructions

### Version 0.8.2
* Add testing and development feature "scenarios": Project files with certain expected outputs, which are used to test if the calculations have changed and thus produce different outputs than in the last verified version. This differs from unit tests in that the tests are not fully automated and require a certain expert knowledge to verify. The advantage over unit tests is that these outputs are less brittle and test cases can be added more easily.
* Add CLI for generating the output of scenarios, set reference outputs, compare the generated order of operations against the reference and generate an overview page over the scenario outputs
* Move most examples over to scenarios, as they were already used in such a manner. The remaining examples are intended to be developed as actual examples for new users of ReSiE.

### Version 0.8.1
* Rename output parameters `dump_info` and `dump_info_file` to `auxiliary_info` and `auxiliary_info_file`
* Rename output parameter `output_file` to `csv_output_file`
* Add output parameters `sankey_plot_file` and `output_plot_file` to allow control of where the plots will be saved
* Move input parameter `weather_file_path` from the `io_settings` section of the project config to the `simulation_parameters` section

### Version 0.8.0
* Standardise storage un-/loading control parameters for all components and make them customisable to the input/output media
* Remove operational strategy "extended_storage_loading_control"
* Update line plot & sankey outputs to better show missing demands
* Rework internal energy distribution calculations on busses:
  * If two or more busses are connected in a chain (necessarily of the same medium), a so-called proxy bus is created, which handles the calculations for any energy transfer from components on any of the principal/original busses
  * Remove the distinction of maximum potential energy utilisation (`max_energy`) and storage un-/loading potential, as the calculation on a bus considers storages according to the energy flow matrix as well as storage loading flags and thus does not require a distinction anymore. Direct 1-to-1 connections do not require this either, hence the removal.
  * Add automatically generated output channels for busses in a chain, where each bus tracks how much energy is transfered to other busses in the input->output direction
  * Some slight changes the generated order of operations, as this now makes use of the input and output order on proxy busses. This should not have an impact on results, but might change the order of operations compared to versions before v0.8.0.

### Version 0.7.1
* bugfix in sankey diagram if only one medium is present
* added possibility to not plot the sankey by adding "sankey_plot" to "io_settings"
* added possibility to define custom colors in sankey for each medium
* renamed input variable "output_keys" --> "csv_output_keys" in input file

### Version 0.7.0
#### Input and output
* Rename parameters and output variables across several components:
  * Load --> Demand                     in output channel of FixedSink
  * static_load --> constant_demand     for parameters of sinks
  * static_* --> constant_*             for parameters where * in (power, temperature, demand, supply)
  * fixed_cop --> constant_cop          for HeatPump
  * draw sum --> output_sum             in output channel of GridConnection
  * load sum --> input_sum              in output channel of GridConnection
  * power --> power_*                   for parameters of all transformers, where * in (el, th)
  * medium --> consider_medium          for control strategy parameters of all transformers
* Add output channel "Losses" to all components and to Sankey output. "Losses" are total losses, while "Losses_XX" are medium-wise break downs
* Constant values for fixed sources and sinks are now given as power instead of work/energy to be consistent with bounded sources and sinks
* Add global logging functionalities with the following categories: Debug, Info, BalanceWarn, Warn, Error and redirected all println() to logger (console and/or logging files, separately for general logs and Balancewarn)
* Update configuration options for busses:
  * Rename "connection_matrix" to "connections"
  * Rename "storage_loading" to "energy_flow"
  * Make "connections" a required part of the config for a bus, including items "input_order" and "output_order", however "energy_flow" remains optional
  * Remove "output_refs" item as "output_order" contains the same information and is now required

#### Functionality
* Sankey diagrams now display the difference of requested and delivered energy in fixed sinks and sources
* Remove condition "would overfill thermal buffer" in storage_driven strategy as this is now handled implicitly
* Add profile aggregation and segmentation with testcases
* Add import of weather files in EPW format and .dat format (DWD)
* Add functionality to map profiles from weather file to component profiles, like ambient temperature from the weather file to a geothermal collector

#### Fixes
* Fix generic storage implementation not being available due to the module not being included
* Fix the profile scaling factor of some components being required despite profiles being optional
* Add missing output channels to Electrolyser

#### Refactorings
* Change the input and output interfaces of busses such that the order matches the input and output priorities
* Rename helper function highest_temperature to highest and add types to inputs
* Provide docstrings for some structs and functions that were missing them
* Rename internal variables to match changes in the input and output variables mentioned above
* Remove last potential() step of transformer chains as this is not needed
* Add required Julia packages: Colors, Interpolations, Dates, Logging
* Rename argument "parameters" for simulation parameters to "sim_params" and add them to all components and profiles
* Rework communication of balances, energy/storage potentials and temperatures via busses. This is an extensive rework that touches almost all components and how the `process` and `potential` simulation steps work. Please note that the rework is not finished with v0.7.0 and will continue to support more energy systems and component configurations that might be of interest to users. However no compatability is knowingly broken with examples that worked in previous versions.

### Version 0.6.5
* "output_keys" and "output_plot" in the input file can now be "all", "nothing" or a list of entries for custom outputs of the CSV file and lineplot (backwards compatibility is given)
* output_values() of all components was changed to not only return the available channels, but also the corresponding media
* restructured "run_simulation()" for better readability

### Version 0.6.4
* Generalise implementation of gas boiler to that of a fuel boiler
  * Although this implicates the use of a chemical fuel to generate heat, the current implementation works with any kind of input including electricity

### Version 0.6.3
* Remove abstract subtype ControlledComponent and use abstract top-level type Component in its stead
* Add a generic implementation for storage components
* Standardise the generic implementations of fixed/bounded sources/sinks by:
  * adding options for static power/energy and temperature values instead of mandatory profiles
  * making Demand an alias to generic implementation FixedSink
* Add "Load%" output for storage components making outputting the state of charge more convenient

### Version 0.6.2
* Added part-load-ratio efficiency calculation for gas boiler component

### Version 0.6.1
* corrected wrong units given for output_plot in example files

### Version 0.6.0
* Refactor terms and names used to refer to various concepts in the simulation model:
  * energy system => component
  * energy network / system topology => energy system
  * production => processing, and variations thereof, including the simulation step produce => process
  * Module EnergySystems remains named the same as the functionality still concerns the modelling of energy systems. In this regard the old name was too limited on the components of an energy system
* Change keys "energy_systems" and "production_refs" in project file to "components" and "output_refs" 

### Version 0.5.7
* Refactor function watt_to_wh as a global function in module EnergySystems

### Version 0.5.6
* Calculation of COP of heat pump changed from dummy implementation to Carnot COP with 40% efficiency and adjusted tests accordingly.
* OoO for control steps changed: The OoO of the controls now corresponds to the base order without any rearrangement! This was necessary to ensure that the transformers are controlled before the storages. Tests were adapted accordingly.
* Added dummy implementation and example file for geothermal probes and geothermal heat collectors.
  
### Version 0.5.5
* minor changes in the output line plot and exemplary input files

### Version 0.5.4
* Fix link to documentation

### Version 0.5.3
* Choose license: MIT
* Update link to documentation

### Version 0.5.2
* updated example input files and their profile data (now given in specific values related to the net floor area or net PV area)

### Version 0.5.1
* Added additional rule to the calculation of OoO that handles multiple buses with storages and transformers that feed energy into leaf busses. The added function adjusts the position of the load() step of the storages in the OoO to ensure that the energy from the leaf transformer is not fed into the wrong storage.

### Version 0.5.0
* Split library and CLI part of Resie into two files. This allows including Resie as a library without immediately executing anything, which is especially useful for tests. At the moment the CLI does nothing but pipe arguments to the library part, but can be extended in the future, for example for additional commands other than running a simulation.

### Version 0.4.1
* balance_on() now always returns the temperature if it is given
* added control() on all energy systems with temperatures to write the required temperatures to the interfaces
* control() of sinks and sources reworked to do the same thing
* adapted tests and added new tests to check the temperatures returned by balance_on()
* added optional input and output temperatures of transformers as user input
* changed example_WP_Kascade to dynamic COP calculation depending on dynamic demand temperatures

### Version 0.4.0
* Rename dispatchable sinks and sources to bounded sinks and sources

### Version 0.3.9
* added optional input of custom OoO to the input file using the "order_of_operation" section and changed the format of the output of "dump_info" to allow copy and paste of OoO to the input file
* added tests for loading user-defined OoOs
* activated tests for transformers-in-a-row
* reworked balance_on() of the bus to allow transformers to differ between storage_potential and energy_potential in their input_interfaces
* added tests to check the new balance_on() of the bus and the functionality of the "unload_storages" flag in the demand-driven control strategy

### Version 0.3.8
* Fix issue with the signs of field max_energy being set for interfaces by energy systems. It is now always positive and method balance_on returns the energy potential with the correct sign depending on whether the caller is input or output in respect to the interface.
* Rework the implementations of energy systems GasBoiler, HeatPump and CHPP to use the new simulation step potential so that their potential and production calculations work correctly for any chain of transformers
* Minor readability improvements
* Remove potential calculations for transformer chains of length 1, as these calculations are superfluous by definition
* Remove reordering of order of operations for control dependencies, as these reorderings have become superfluous with the other reorderings covering the same cases

### Version 0.3.7
* adapted OoO for storages according the following rules:
  1. produce() and load() from leaves to trunk
  2. produce() and load() of leaves are in order of output interfaces of trunk busses
  3. produce() and load() of loading-limited storages are calculated after loading-unlimited storages with respect to the loading_matrix in busses
* adapted and added test for new calculation order of storages

### Version 0.3.6
* Refactor algorithm for determining order of operations for a more modular approach with reusable helper functions and add tests for new modular functions.
* Add helper functions for improving output of tests checking the order of operations
* Add simulation step for calculating energy potentials. The default method does nothing and transformer implementations ought to calculate their potentials in the step and write the information to the interfaces.
* Add handling of chains of transformers to order of operations algorithm. The algorithm iterates over each chain in one pass sink-to-source, adding potential steps along the way, then another pass source-to-sink adding the produce steps.
* Add deactivated tests highlighting problems with transformer chains - the tests fail until transformer implementations for potentials are added
* Add handling of storage loading to order of operations algorithm. The algorithm considers the output priorities on chains of busses for the loading of storages connected to all busses of a chain.

### Version 0.3.5
* added extended_storage_control strategy for all storages to allow users to decide whether a storage can load ANY other interconnected storage of the system framework or not (default = storages are not allowed to load other storages). This only works if the OoO of produce() and load() is correct which will be implemented in upcoming versions.
* changed definition in balance_on() of busses: Storages are never allowed to load themselves!
* adjusted testcases to meet new definition above
* added testcases to test extended_storage_control

### Version 0.3.4
* added feature to control transformers from user-given control profile with values within [0,1]. This adds to the already existing limitations in the control strategy of demand_driven, supply_driven and storage_driven strategy.
* bugfix in index assignement of connectivity matrix to output_interfaces in balance_on() of bus
* changed calculation of COP of heatpump: If fixed_cop is given in the input file, it will be used and not the temperature-dependend COP. Corrected tests to meet this definition.

### Version 0.3.3
* bugfix in calculation of IN and OUT energy in output_value() for chases with balance /= 0. 

### Version 0.3.2
* adapted control strategies for all transformers: All in- and outputs need to be sattisfied by default
* added optional user-defined input in control strategie of each transformer to ignore certain in- or outputs within control strategy
* added optional user-defined input in control strategie of each transformer to allow or deny system-wide storage loading or unloading
* changed max_power to max_energy for all energy systems in outputs and unit structs
* added max_energy to interfaces and control() of grids, sources and sinks in order to provide information on how many energy is available or can be taken, calculated in balance_on(), to meet new control strategies of transformers
* adapted balance_on() to also return energy_potential (if sum_abs_change /= 0) that represents the potential maximum energy in- or output at unit (former potential is now named storage_potential, balance_on() returns now a named tuple instead of tuples)
* added and corrected tests to meet new control strategies and test new functionalities 
* adapted examples to meet new functionalities

### Version 0.3.1
* Added tests for demand_driven strategy with and without bus between transformer and demand
* Adapted default control strategy of gasboiler and heatpump to enable storage-filling

### Version 0.3.0
* Implement connectivity matrix on bus systems for controlling which producer can load which storage on the same bus
* Remove operational strategy "use_surplus_in_cycle" as it is no longer needed with this new feature
* Move definitions of input and output priorities into connectivity matrix

### Version 0.2.11
* Implemented user-definable media for all energy systems
* Adaption of tests to the new user-definable media names
* Small bug fixes

### Version 0.2.10
* Refactor distribution of bus systems to better match expected behaviour

### Version 0.2.9
* Apply auto-format to existing code files
* Rework how medium categories are handled
* Adapt example energy system implementation (electrolyser) to use new functionality for enabling custom-defined medium categories

### Version 0.2.8
* Bugfix in calculation of energy flow between two interconnected busses
* Added reordering of Load() functions for storages 
* Bugfix in determination of order of Distribute() of busses
* Added user-definable medium names in input file (only for electrolyser for now)
* Bugfix in SeasonalThermalStorage (missing declaration)
* Revice of exemplary system frameworks
* Adapted testcases to meet the changes

### Version 0.2.7:
* Restructure test cases so that each case can be debugged
* Update existing and add new tests to highlight several problems that occur when bus systems are connected

### Version 0.2.6:
* Refactor operational strategies for GasBoiler, CHPP, Electrolyser
* Implement way to track the energy flow between connected busses
* Improve sankey output (hide null-flows and hide O2)
* Refactor order of operations to include distribution between connected busses
* Bugfix in distribution calculation of bus systems
* Change examples to illstrutrate current problems
* Change how temperatures are read from profile files
* Make feature of load-dependant temperatures of heat storages optional

### Version 0.2.5:
* Refactor and fix how the storage loading potential and temperatures are transfered across busses, in partiular when multiple busses are connected in series

### Version 0.2.4:
* added example of system topology Esslingen in different variantes
* added profiles for PV, heat demand
* removed deprecated code in pv_plant.jl
* bugfix in Resie.jl in collecting data of all interfaces (added div 2 for sum_abs_change and correct results for non-balanced interfaces)

### Version 0.2.3:
* Fix energy balance and distribution calculations of busses not working correctly for system topologies in which two bus systems are connected directly
* Fix distribution overwriting temperatures on system interfaces

### Version 0.2.2:
* fixed bug with coloring of the medium in Sankey occuring in file_output.jl when only one medium is present in energy topology
* added example input file with one medium and two busses in paralell (currently not working!)

### Version 0.2.1:
* added automated generation of sankey diagram in output. Requires additional package "ColorSchemes" from Plotly.

### Version 0.2.0:
* changed PV-System input from power (sin-function) to given power/energy profile scaled by scale factor

### Version 0.1.0:
* Initial release that encompasses the entire development before versioning
* Add a changelog
* Start of versioning oriented by semantic versioning
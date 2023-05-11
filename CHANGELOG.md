# Changes made, features added and bugs fixed
In general the development follows the [semantic versioning](https://semver.org/) scheme.
As the simulation is parameterized by the project file, the structure and meaning of this
file represents the API of the simulation engine. Evaluating compatability then follows
required changes in the project file. This is complicated by the implementation of the
individual energy systems, which happens within the framework but also reasonably
independant of the simulation model. A change in the implementation may have big effects
on simulation results without requiring any change in the project file. It is up to the
developers to find a workable middle ground.

## Pre-1.0-releases
As per the definition of semantic versioning and the reality of early development, in
versions prior to 1.0.0 any release might break compatability. To alleviate this somewhat,
the meaning of major-minor-patch is "downshifted" to zero-major-minor. However some
breaking changes may slip beneath notice.

## Version 0.6.0
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
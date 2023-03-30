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
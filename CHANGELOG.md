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

### Version 0.2.0:
* changed PV-System input from power (sin-function) to given power/energy profile scaled by scale factor

### Version 0.1.0:
* Initial release that encompasses the entire development before versioning
* Add a changelog
* Start of versioning oriented by semantic versioning
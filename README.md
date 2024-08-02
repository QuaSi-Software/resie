# Resie - simulation engine for district-scale networks of energy systems

Use of Resie is described in more detail in the accompanying documentation. You can find a rendered version online at [the official readthedocs page](https://quasi-software.readthedocs.io). This document describes installation and contains useful information for developers, who wish to work with ReSiE.

Resie is released under the MIT license. You can find a copy of the license in file `LICENSE.md`. For information on how you can contribute please check the documentation.

## Installation

### **Requirements**

* Julia, minimum v1.8.5 and tested up to v1.10.4. You can find installation instructions [here](https://julialang.org/downloads/).
* (Optional) Python3, tested with v3.9.13. Only required for developing ReSiE.

### Instructions

1. Get a copy: `git clone git@github.com:QuaSi-Software/resie.git`
1. Switch into the ReSiE root directory: `cd /path/to/resie`
1. Start the julia REPL with `julia`
1. Switch to the package REPL with `]` (no enter necessary)
1. Activate the project environment: `activate .`
1. Install and precompile required packages: `instantiate`. This should create a file `Manifest.toml` in the ReSiE root directory
1. Exit out of the package REPL with shortcut `Ctrl+c`
1. Exit out of the julia REPL with `exit()` or shortcut `Ctrl+d`

If you wish to develop with this installation of ReSiE you should also perform the following:

1. Install the pre-commit framework: `pip install pre-commit`
1. Install the pre-commit hooks into the ReSiE installation: `pre-commit install`

## Usage

A full description of how to use ReSiE on the examples it ships with can be found [in this chapter](https://quasi-software.readthedocs.io/en/latest/resie_exemplary_energy_systems/). In the following an abbreviated version:

1. Switch into project directory: `cd /path/to/resie`
1. Run the simulation with `julia --project=. src/resie-cli.jl examples/simple_heat_pump.json`
1. The outputs as well as log files can be found in the `output` folder. The simulation should run without errors and produce a file called `output/output_plot.html` which, when opened in a browser, shows an interactive plot of simulation results.

## Development

### Testing

There are two ways to run tests, one using just the command line and one using the julia (and PKG) REPL. The latter is prefered by the julia community, but does not suit the workflows of this project. The former is also useful for running tests in IDEs or text editors as it can be more easily automated.

**Note: At the moment the REPL method does not work, as the relative paths in the test project configs do not work with the CWD of the test runner, which is the `test` subdirectory for the REPL, but the base directory for the command-line test runner.**

#### Using the command-line

1. Switch into project directory: `cd /path/to/resie`
1. Run the tests with `julia --project=. test/runtests.jl`

#### Using the REPL

1. Switch into project directory: `cd /path/to/resie`
1. Start the julia REPL with `julia --project=.`
1. Run the test suite with `test Resie`

#### Enable the debugger in VS Code when running tests

Open `launch.json` (via the gear wheel in the run config dropdown in the debug tab) and add entry:
```json
{
    "type": "julia",
    "request": "launch",
    "name": "Run Resie tests",
    "program": "path/to/resie/test/runtests.jl",
    "stopOnEntry": false,
    "cwd": "path/to/resie/",
    "juliaEnv": ".",
    "args": []
}
```

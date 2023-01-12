# Bran - simulation engine for district-scale networks of energy systems

Use of Bran is described in more detail in the accompanying documentation. You can find a rendered version online at [TBD](http://example.com). This document describes installation and basic instructions, particularly for developers.

## Installation

### **Requirements**

* Julia v1.7.2 or later

### Instructions

1. Get a copy: `git clone git@bran.example.com`
1. Switch into project directory: `cd /path/to/repo`
1. Start the julia REPL with `julia`
1. Switch to the package REPL with `]` (no enter necessary)
1. Activate the environment and precompile Bran with `activate .`
1. Exit out of the package REPL with shortcut `Ctrl+c`
1. Exit out of the julia REPL with `exit()`

## Usage

1. Switch into project directory: `cd /path/to/bran`
1. Run the simulation with `julia src/Bran.jl examples/example_two_sector.json`
1. Outputs of the example projects can be found in `output/out.csv` and `output/info_dump.md`

## Testing

There are two ways to run tests, one using just the command line and one using the julia (and PKG) REPL. The latter is prefered by the julia community, but does not suit the workflows of this project. The former is also useful for running tests in IDEs or text editors as it can be more easily automated.

**Note: At the moment the REPL method does not work, as the relative paths in the test project configs do not work with the CWD of the test runner, which is the `test` subdirectory for the REPL, but the base directory for the command-line test runner.**

### Using the command-line

1. Switch into project directory: `cd /path/to/bran`
1. Run the tests with `julia --project=. test/runtests.jl`

**Note: At the moment this produces an erroneous output of "No project config file given", but as this does not affect test results, this can be ignored.**

### Using the REPL

1. Switch into project directory: `cd /path/to/bran`
1. Start the julia REPL with `julia`
1. Switch to the package REPL with `]` (no enter necessary)
1. Activate the environment and precompile Bran with `activate .`
1. Run the test suite with `test Bran`

### Enable testing in the VS Code debugger

Open `launch.json` (via the gear wheel in the run config dropdown in the debug tab) and add entry:
```json
{
    "type": "julia",
    "request": "launch",
    "name": "Run Bran tests",
    "program": "path/to/bran/test/runtests.jl",
    "stopOnEntry": false,
    "cwd": "path/to/bran/",
    "juliaEnv": ".",
    "args": []
}
```
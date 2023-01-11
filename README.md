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

1. Switch into project directory: `cd /path/to/bran`
1. Start the julia REPL with `julia`
1. Switch to the package REPL with `]` (no enter necessary)
1. Activate the environment and precompile Bran with `activate .`
1. Run the test suite with `test Bran`

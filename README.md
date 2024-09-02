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

To use profiles with a time-zone-aware datetime index beyond 2038, follow these steps:

1. Copy the zip file found at `additional_files/move_to_user_julia_scratchspaces_and_unpack.zip` from the repository to your local Julia scratchspaces directory, typically located at `C:/users/user_name/.julia/scratchspaces`.
2. Extract the content of the zip file directly into the `scratchspaces` folder. Ensure that the folder `f269a46b-ccf7-5d73-abea-4c690281aa53` is placed directly within `scratchspaces`, without any intermediate directories.

If you wish to develop with this installation of ReSiE you should also perform the following inside the ReSiE root directory:

1. Install the pre-commit framework: `pip install pre-commit`
2. Install the pre-commit hooks into the ReSiE installation: `pre-commit install`

## Usage

A full description of how to use ReSiE on the examples it ships with can be found [in this chapter](https://quasi-software.readthedocs.io/en/latest/resie_exemplary_energy_systems/). In the following an abbreviated version:

1. Switch into project directory: `cd /path/to/resie`
1. Run the simulation with `julia --project=. src/resie-cli.jl examples/simple_heat_pump.json`
1. The outputs as well as log files can be found in the `output` folder. The simulation should run without errors and produce a file called `output/output_plot.html` which, when opened in a browser, shows an interactive plot of simulation results.

## Development
The following sections contain information useful for developers.

### Code style
The prefered code style for ReSiE is the [YAS style](https://domluna.github.io/JuliaFormatter.jl/stable/yas_style/) available for [JuliaFormatter](https://domluna.github.io/JuliaFormatter.jl/stable/), with some custom config settings. The rules are defined in file `.JuliaFormatter.toml`, which is used automatically used by the formatter for all files and subdirectories in the ReSiE folder. Formatting happens in one of three cases:

1. Running `using JuliaFormatter; format("/file/to/format.jl")` in the REPL
2. Using `Format document` in the context menu in a file in VS Code
3. As a pre-commit hook before a commit is performed. This is described in more detail below.

The YAS style guide and the general guidelines for Julia also define rules for things the formatter cannot automatically adjust. It is up the developers to try to follow these guidelines when writing code. The hard line limit for the formatter is set to 120 characters, however 92 is preferred as a soft limit for things that can be easily distributed across lines, such as comments and docstrings.

#### Pre-commit hook for formatting
As an automatic form of formatting Julia code files a pre-commit hook is employed. If the additional installation instructions were followed, this hook will be executed when the git `commit` command is executed. The hook will apply the formatter on `.jl` files that have been changed with the commit. If this formatting caused additional changes, the commit is aborted and the changes applied to the workspace (but not the staging area). The additional changes need to be checked if they are okay, then added to the staging area and committed again. As the formatting changes have been added, this commit should go through.

This pre-commit hook was introduced with version 0.9.2, but at the time the formatting was **not** performed on all code files. We decided on this course of action as the work required to check all formatting changes at once was deemed too much and committing the changes without checking them was deemed inadvisable. The automatic formatting is supposed to increase readability by providing a unified code style across the whole repository, however the amount of unchecked formatting changes actually decreased readability due to the formatting resulting in some questionable choices.

Because of this decision, not all files are in the new automatic format. When developing ReSiE it is advised to first run the formatter on any files touched by changes and commit the changes, after checking they seem reasonable, as a separate commit. Over time this will result in all files fitting in the defined code style.

#### Downsides of the code style
At time of writing some downsides to the chosen code style are known, that do not sit well with any of the developers. In general we believe it is possible to work within the bounds of the code style to produce readable code. Some considerations are described in the following.

**Alignment of multi-line conditions**

A boolean condition spanning multiple lines might be written like this:
```julia
if (
    check_something(a, b) ||
    check_something_else(b, c) &&
    some_important_condition(a, b, c)
)
    do_something()
end
```
If the condition is longer than the line limit, the formatter will rewrite it like this:
```julia
if (check_something(a, b) || check_something_else(b, c)
    &&
    some_important_condition(a, b, c))
    do_something()
end
```
Now the body `do_something()` of the conditional looks, at first glance, to be part of the condition as the closing parenthesis is not on its own line. To increase readability it is recommended to add a comment that visually breaks up the condition from the body:
```julia
if (check_something(a, b) || check_something_else(b, c)
    &&
    some_important_condition(a, b, c))
    # end of condition
    do_something()
end
```

**Vanishing comments in parentheses**

At time of writing there is [a bug with JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl/issues/690) that will delete inline comments that occur before the first and after the last token of a list (which may be a condition or argument list). For example:
```julia
instance = SomeStruct( # first inline
    first_argument, # regular inline
    # stand-alone
    second_argument  # last inline
)
```
This will be reformatted as:
```julia
instance = SomeStruct(first_argument, # regular inline
                      # stand-alone
                      second_argument)
```

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

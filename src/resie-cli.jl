using Resie

"""
    main()

Entry point of the CLI.

# Command line arguments
## Positional arguments
- `String`: Filepath to the project config file (see documentation on file format). Can be
a path relative to the CWD of the caller.
"""
function main()
    if length(ARGS) > 0
        filepath = ARGS[1]
        if filepath !== nothing && filepath != ""
            Resie.load_and_run(filepath)
            return
        end

        println("Could not find or access project config file")
        return
    end

    println("No project config file argument given")
end

main()

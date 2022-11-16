module Profiles

export Profile

"""
Holds values from a file so they can be retrieved later and indexed by time.

For now the implementation reads values as they are in the file, but this will later change,
so that values are interpolated between the time steps of the simulation, making the profile
values somewhat independant from it (the overall timespan still needs to be the same).
"""
mutable struct Profile
    """Time step, in seconds, of the profile."""
    time_step :: Int

    """Holds the profile values indexed by the time step number"""
    data :: Vector{Float64}

    """Construct a profile from the file at the given path."""
    function Profile(file_path :: String)
        parsed = Vector{Float64}()
        time_step = 900

        open(file_path, "r") do file_handle
            line = strip(readline(file_handle))

            if line[1] == '#'
                splitted = split(strip(line, '#'), ';')
                if strip(splitted[1]) == "time_step"
                    time_step = parse(Int, splitted[2])
                end

            else
                splitted = split(line, ';')
                time = parse(Float64, splitted[1])
                value = parse(Float64, splitted[2])

                step_nr = Int(round(time / time_step))
                push!(parsed, (step_nr, value))
            end
        end

        return new(
            time_step, # time_step
            parsed # data
        )
    end
end

end
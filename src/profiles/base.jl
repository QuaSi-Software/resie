module Profiles

export Profile, power_at_time, work_at_time

"""
Holds values from a file so they can be retrieved later and indexed by time.

For now the implementation reads values as they are in the file, but this will later change,
so that values are interpolated between the time steps of the simulation, making the profile
values somewhat independant from it (the overall timespan still needs to be the same).
"""
mutable struct Profile
    """Time step, in seconds, of the profile."""
    time_step :: Int

    """Indicates whether the profile values are power or work."""
    is_power :: Bool

    """Holds the profile values indexed by the time step number"""
    data :: Vector{Float64}

    """Construct a profile from the file at the given path."""
    function Profile(file_path :: String)
        parsed = Vector{Float64}()
        time_step = 900
        is_power = false

        open(abspath(file_path), "r") do file_handle
            for line in readlines(file_handle)
                line = strip(line)

                if line[1] == '#'
                    splitted = split(strip(line, '#'), ';')
                    if strip(splitted[1]) == "time_step"
                        time_step = parse(Int, splitted[2])
                    elseif strip(splitted[1]) == "is_power"
                        is_power = parse(Bool, splitted[2])
                    end

                else
                    splitted = split(line, ';')
                    value = parse(Float64, splitted[2])
                    push!(parsed, value)
                end
            end
        end

        return new(
            time_step, # time_step
            is_power, # is_power
            parsed # data
        )
    end
end

"""
    power_at_time(profile, time)

Get the power value of the profile at the given time.
"""
function power_at_time(profile :: Profile, time :: Int)
    step_nr = Int(round(time / profile.time_step) + 1)
    if profile.is_power
        return profile.data[step_nr]
    else
        return profile.data[step_nr] * (3600 / profile.time_step)
    end
end

"""
    work_at_time(profile, time)

Get the work value of the profile at the given time.
"""
function work_at_time(profile :: Profile, time :: Int)
    step_nr = Int(round(time / profile.time_step) + 1)
    if profile.is_power
        return profile.data[step_nr] * (profile.time_step / 3600)
    else
        return profile.data[step_nr]
    end
end

end
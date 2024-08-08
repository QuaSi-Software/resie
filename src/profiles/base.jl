module Profiles
using Interpolations
using Dates

export Profile, power_at_time, work_at_time, value_at_time

"""
Holds values from a file so they can be retrieved later and indexed by time.

Profiles are automatically aggregated or segmentated to fit the simulation time step.
Currtently, this only works, if the time step of the profile is a multiple or a divisor 
of the requested simulation timestep. Otherwise, an error will arise.

This function can handle either a path to a .prf file, or the profile data can be
handed over as vectors. Then, the profile with the corresponding timestamps, the time step
and the flag is_power has to be provided as optional arguments. 

The flag is_power is essential to make sure that the profile data is converted correcty.
It has to be set true for all intensive values like temperatures, power, wind speed or
for profiles containting states. Only for energies, is_power has to be set to false.
"""
mutable struct Profile
    """Time step, in seconds, of the profile."""
    time_step::Int

    """Indicates whether the profile values are power or work."""
    is_power::Bool

    """Holds the profile values indexed by the time step number"""
    data::Vector{Float64}

    """Construct a profile from the file at the given path or convert given data to profile (optional)."""
    function Profile(file_path::String,                            # file path to a .prf file
                     sim_params::Dict{String,Any};                 # general simulation parameters
                     given_profile_values::Vector{Float64}=[0.0],  # optional: Vector{Float64} that holds values of the profile
                     given_timestamps::Vector{Int64}=[0],          # optional: Vector{Int64} that holds the time step in seconds
                     given_time_step::Int=0,                       # optional: Int that indicates the timestep in seconds of the given data
                     given_is_power::Bool=false)                   # optional: Bool that indicates if the data is intensive or extensive

        if given_profile_values == [0.0]  # read data from file_path
            profile_values = Vector{Float64}()
            profile_timestamps = Vector{Int64}()
            profile_time_step = nothing
            is_power = nothing
            time_format = nothing
            first_time_step = nothing

            file_handle = open(abspath(file_path), "r")
            for line in readlines(file_handle)
                line = strip(line)

                if isempty(line) || length(line) < 2 # handle empty lines
                    continue
                elseif line[1] == '#'
                    splitted = split(strip(line, '#'), ';')
                    if strip(splitted[1]) == "time_step"
                        profile_time_step = parse(Int, String(strip(splitted[2])))
                    elseif strip(splitted[1]) == "is_power"
                        is_power = parse(Bool, String(strip(splitted[2])))
                    elseif strip(splitted[1]) == "time_format"
                        time_format = String(strip(splitted[2]))
                    end
                else
                    if time_format === nothing
                        time_format = "seconds"
                        @warn "For the profile at $(file_path) no time format ('time_format') is given! A time stamp in seconds is assumed."
                    end
                    splitted = split(line, ';')
                    parsed_timestep_seconds, first_time_step = parse_timestamp(String(strip(splitted[1])), time_format, file_path, first_time_step)
                    push!(profile_timestamps, parsed_timestep_seconds)
                    push!(profile_values, parse(Float64, splitted[2]))
                end
            end
            close(file_handle)

            if profile_time_step === nothing
                profile_time_step = profile_timestamps[2] - profile_timestamps[1]
                @warn "For the profile at $(file_path) no timestep is given! A timestep of $(string(profile_time_step)) seconds was detected."
            end
            if is_power === nothing
                @error "For the profile at $(file_path) no profile type ('is_power') is given!"
                exit()
            end
        else # data was read in from somewhere else and is provided as vector
            profile_values = given_profile_values
            profile_time_step = given_time_step
            profile_timestamps = given_timestamps
            is_power = given_is_power
        end

        simulation_time_step = sim_params["time_step_seconds"]  # seconds 

        if !(simulation_time_step % profile_time_step == 0) && !(profile_time_step % simulation_time_step == 0) 
            @error ("The timestep of the profile " * file_path * ", which is " * string(profile_time_step) * " s,\n" *
                    "is not a multiple or a divisor of the requested simulation timestep of " * string(simulation_time_step) * " s!")
            throw(InputError)
        end

        if is_power # meaning intensive values (e.g., temperature, power)
            values_converted = convert_intensive_profile(profile_values, 
                                                         profile_timestamps, 
                                                         Int(profile_time_step), 
                                                         Int(simulation_time_step),
                                                         file_path)
        else # meaning extensive values (e.g., energy demand)
            values_converted = convert_extensive_profile(profile_values, 
                                                         profile_timestamps,
                                                         Int(profile_time_step), 
                                                         Int(simulation_time_step),
                                                         file_path)
        end

        return new(
            simulation_time_step,   # Int: time_step, equals simulation time step after conversion, in seconds
            is_power,               # Bool: to indicate if profil is intensive or extensive (true means is intensive profile)
            values_converted        # Vector{Float64} data of profiles in simulation time step
        )
    end
end

"""
    parse_timestamp()

function to pase a timestamp to seconds, counting from the begin of the year.
    datetime_str::String        Timestamp read from the profile file. Can be a String(Int) or a specified DateTime format as string.
    time_format::String         Format of datetime_str. Can be either "seconds", "hours" or a DateTime format, like "dd-mm-yyyy HH:MM:SS"
    file_path::String           File path of the profile for error handling.

Returns
    parsed_timestep_seconds::Float64                the parsed timestep in seconds with respect to `first_time_step`
    first_time_step::Union{Nothing, DateTime}       if `time_format` is a custom format, `first_time_step` holds the first timestep as DateTime
                                                    if `time_format` is `seconds` or `hours`, `first_time_step` remains nothing.
"""
function parse_timestamp(datetime_str::String, time_format::String, file_path::String, first_time_step::Union{Nothing, DateTime})
    if time_format == "seconds"
        return parse(Int64, datetime_str), first_time_step
    elseif time_format == "hours"
        return parse(Int64, datetime_str) * 60 , first_time_step
    else
        try
            parsed_datetime = Dates.DateTime(datetime_str, time_format)
            if first_time_step === nothing
                first_time_step = parsed_datetime
            end
            return Int64(Dates.value(parsed_datetime  - first_time_step) / 1_000), first_time_step # seconds
        catch e
            @error("Time given 'time_format' of profile at $(file_path) does not fit to the data.\n" *
                   "'time_format' has to be 'seconds', 'hours' or a specific daytime format, e.g. 'dd-mm-yyyy HH:MM:SS'.\n" *
                   "'time_format' is `$(time_format)` which does not fit to the timestamp given `$(datetime_str)`.\n" *
                   "The following error occured: $e")
            exit()  
        end
    end
end

"""
    power_at_time(profile, time)

Get the power value of the profile at the given time.
"""
function power_at_time(profile::Profile, time::Int)
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
function work_at_time(profile::Profile, time::Int)
    step_nr = Int(round(time / profile.time_step) + 1)
    if profile.is_power
        return profile.data[step_nr] * (profile.time_step / 3600)
    else
        return profile.data[step_nr]
    end
end

"""
    value_at_time(profile, time)

Get the value of the profile at the given time without any conversion.
The flag is_power will be ignored.
"""
function value_at_time(profile::Profile, time::Int)
    step_nr = Int(round(time / profile.time_step) + 1)
    return profile.data[step_nr]
end

"""
    minimum(profile)

The smallest value of the profile.

# Arguments
- `profile::Profile`: The profile
# Returns
- `Float64`: The smallest value
"""
function minimum(profile::Profile)::Float64
    return Base.minimum(profile.data)
end

"""
    maximum(profile)

The largest value of the profile.

# Arguments
- `profile::Profile`: The profile
# Returns
- `Float64`: The largest value
"""
function maximum(profile::Profile)::Float64
    return Base.maximum(profile.data)
end
"""
    convert_extensive_profile

Function to convert extensive profiles(e.g., energy demand) from the profile time step to the simulation time step.

Inputs: 
    values::Vector{Float64}         values of the profile to convert
    timestamps::Vector{Int64}       timestamps in seconds corresponding to values
    original_time_step::Int         the timestep in seconds of the values and the timestamps
    new_time_step::UInt64           the simulation time step to wich the profile should be converted
    file_path::String               the file path of the profile, for error messages only

Outputs:
    converted_profile::Vector{Float64}      values of the converted profile in the `new_time_step`

If `original_time_step` equals the `new_time_step`, the original `values` will be returned.
"""
 function convert_extensive_profile(values::Vector{Float64},
                                    timestamps::Vector{Int64},
                                    original_time_step::Int64,
                                    new_time_step::Int64,
                                    file_path::String)

    if original_time_step == new_time_step
        return values
    end

    # handle segmentation and aggregation at once
    new_max_timestep = ceil(Int, Base.maximum(timestamps) / new_time_step) * new_time_step
    if original_time_step > new_time_step
        new_length = ceil(Int, new_max_timestep / new_time_step + original_time_step / new_time_step)
    else
        new_length = ceil(Int, new_max_timestep / new_time_step)
    end
    converted_profile = zeros(new_length)

    for (timestamp, value) in zip(timestamps, values)
        new_time_index_start = floor(Int, timestamp / new_time_step) + 1
        new_time_index_end = floor(Int, (timestamp + original_time_step - 1) / new_time_step) + 1

        for j in max(1, new_time_index_start):min(new_time_index_end, new_length)
            overlap_timestamp_start = max((j - 1) * new_time_step, timestamp)
            overlap_timestamp_end = min(j * new_time_step, timestamp + original_time_step)
            overlap_time_duration = max(0, overlap_timestamp_end - overlap_timestamp_start)

            if overlap_time_duration > 0
                overlap_fraction = overlap_time_duration / original_time_step
                converted_profile[j] += value * overlap_fraction
            end
        end
    end

    @info "The profile at $(file_path) (extensive profile) was converted from the profile timestep $(original_time_step) s to the simulation timestep of $(new_time_step) s."

    return converted_profile
end


"""
    convert_intensive_profile

Function to convert intensive profiles (e.g., temperature, power) from the profile time step to the simulation time step.

Inputs: 
    values::Vector{Float64}         values of the profile to convert
    timestamps::Vector{Int64}       timestamps in seconds corresponding to values
    original_time_step::Int         the timestep in seconds of the values and the timestamps
    new_time_step::UInt64           the simulation time step to wich the profile should be converted
    file_path::String               the file path of the profile, for error messages only

Outputs:
    converted_profile::Vector{Float64}      values of the converted profile in the `new_time_step`

If `original_time_step` equals the `new_time_step`, the original `values` will be returned.
"""
function convert_intensive_profile(values::Vector{Float64},
                                   timestamps::Vector{Int64},
                                   original_time_step::Int64,
                                   new_time_step::Int64,
                                   file_path::String)

    if new_time_step == original_time_step # no change
        return values

    elseif new_time_step < original_time_step  # segmentation
        interp = interpolate((timestamps,), values, Gridded(Linear()))
        new_timestamps = Base.minimum(timestamps):new_time_step:Base.maximum(timestamps)
        converted_profile = [interp(t) for t in new_timestamps]

        # add missing entries by coping the last entry
        for _ in 1:ceil(Int,(original_time_step-new_time_step)/new_time_step)
            append!(converted_profile, converted_profile[end])
        end

        @info "The profile at $(file_path) (intensive profile) was converted from the profile timestep $(original_time_step) s to the simulation timestep of $(new_time_step) s."
        
        return converted_profile

    else # aggregation
        aggregation_factor = new_time_step / original_time_step   # is always > 1 and of type {Int} as only full dividers are allowed
        old_length = length(timestamps)
        new_length = ceil(Int, old_length / aggregation_factor)
        converted_profile = zeros(new_length)
        
        for n in 1:new_length
            old_index_start = Int((n-1) * aggregation_factor + 1)
            old_index_end = Int(min(n * aggregation_factor, old_length))
            old_number_of_steps = Int(old_index_end - old_index_start + 1)
            converted_profile[n] = sum(values[old_index_start:old_index_end])/old_number_of_steps
        end

        @info "The profile at $(file_path) (intensive profile) was converted from the profile timestep $(original_time_step) s to the simulation timestep of $(new_time_step) s."
        
        return converted_profile
    end
end

end
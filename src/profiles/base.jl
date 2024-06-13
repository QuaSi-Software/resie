module Profiles
using Interpolations

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
    function Profile(file_path::String, 
                     sim_params::Dict{String,Any}; 
                     given_profile_values::Vector{Any}=[],
                     given_timestamps::Vector{Int64}=[0],
                     given_time_step::Int=0,
                     given_is_power::Bool=false)

        if given_profile_values == []  # read data from file_path
            profile_values = Vector{Float64}()
            profile_timestamps = Vector{Float64}()
            profile_time_step = 900
            is_power = false
            given_timestep = false
            given_type = false

            open(abspath(file_path), "r") do file_handle
                for line in readlines(file_handle)
                    line = strip(line)

                    if isempty(line) || length(line) < 2 # handle empty lines
                        continue
                    elseif line[1] == '#'
                        splitted = split(strip(line, '#'), ';')
                        if strip(splitted[1]) == "time_step"
                            profile_time_step = parse(Int, splitted[2])
                            given_timestep = true
                        elseif strip(splitted[1]) == "is_power"
                            is_power = parse(Bool, splitted[2])
                            given_type = true
                        end
                    else
                        splitted = split(line, ';')
                        timestamp = parse(Float64, splitted[1])
                        push!(profile_timestamps, timestamp)
                        value = parse(Float64, splitted[2])
                        push!(profile_values, value)
                    end
                end
            end

            if !given_timestep
                profile_time_step = profile_timestamps[2] - profile_timestamps[1]
                @warn "For the profile at " * file_path * " no timestep is given! A timestep of " * string(profile_time_step) * " seconds was detected."
            end
            if !given_type
                @warn "For the profile at " * file_path * " no profile type ('is_power') is given! An energy-profile (extensive) is assumed!"
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
                                                         profile_time_step, 
                                                         simulation_time_step)
        else # meaning extensive values (e.g., energy demand)
            values_converted = convert_extensive_profile(profile_values, 
                                                         profile_timestamps,
                                                         profile_time_step,
                                                         simulation_time_step)
        end

        return new(
            simulation_time_step, # time_step, equals simulation time step after conversion
            is_power, # is_power (is intensive profile)
            values_converted # data of profiles in simulation time step
        )
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

# Function to handle extensive profiles to the simulation time step  (e.g., energy demand)
 function convert_extensive_profile(values, timestamps, original_time_step, new_time_step)
    if original_time_step == new_time_step
        return values
    end

    # handle segmentation and aggregation at once
    new_max_timestep = ceil(Int, maximum(timestamps) / new_time_step) * new_time_step
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

    return converted_profile
end

# Function to convert intensive profiles to the simulation time step (e.g., temperature, power)
function convert_intensive_profile(values, timestamps, original_time_step, new_time_step)
    if new_time_step == original_time_step # no change
        return values
    elseif new_time_step < original_time_step  # segmentation
        interp = interpolate((timestamps,), values, Gridded(Linear()))
        new_timestamps = minimum(timestamps):new_time_step:maximum(timestamps)
        converted_profile = [interp(t) for t in new_timestamps]
        # add missing entries by coping the last entry
        for _ in 1:ceil(Int,(original_time_step-new_time_step)/new_time_step)
            append!(converted_profile, converted_profile[end])
        end
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
        return converted_profile
    end
end

end
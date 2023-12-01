module Profiles
using Interpolations

export Profile, power_at_time, work_at_time, value_at_time

"""
Holds values from a file so they can be retrieved later and indexed by time.

For now the implementation reads values as they are in the file, but this will later change,
so that values are interpolated between the time steps of the simulation, making the profile
values somewhat independant from it (the overall timespan still needs to be the same).
"""
mutable struct Profile
    """Time step, in seconds, of the profile."""
    time_step::Int

    """Indicates whether the profile values are power or work."""
    is_power::Bool

    """Holds the profile values indexed by the time step number"""
    data::Vector{Float64}

    """Construct a profile from the file at the given path."""
    function Profile(file_path::String, parameters::Dict{String,Any})
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
            println("Warning: For the profile at " * file_path * " no timestep is given! A timestep of " * string(profile_time_step) * " seconds was detected.")
        end
        if !given_type
            println("Warning: For the profile at " * file_path * " no profile type ('is_power') is given! An energy-profile (extensive) is assumed!")
        end

        simulation_time_step = parameters["time_step_seconds"]  # seconds 

        if !(simulation_time_step % profile_time_step == 0) && !(profile_time_step % simulation_time_step == 0) 
            println("Error: The timestep of the profile " * file_path * 
                    ", which is " * string(profile_time_step) * 
                    " s, is not a multiple or a divisor of the requested simulation timestep of " * 
                    string(simulation_time_step) * " s!")
            exit()
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


# Function to handle extensive profiles to the simulation time step  (e.g., energy demand)
function convert_extensive_profile(values, timestamps, original_time_step, new_time_step)
    if original_time_step == new_time_step  # no change
        return values
    else # segmentation or aggragation 
        # Initialize an array to store the converted profile
        total_duration = ceil(Int, maximum(timestamps) / new_time_step) * new_time_step
        old_length = length(values)
        if original_time_step > new_time_step
            new_length = ceil(Int, total_duration / new_time_step + original_time_step / new_time_step)
        else
            new_length = ceil(Int, total_duration / new_time_step)
        end
        converted_profile = zeros(new_length)
    
        begin_search = 1
        for i in 1:new_length
            new_start_time = (i - 1) * new_time_step
            new_end_time = i * new_time_step
        
            for m in begin_search:old_length
                timestamp = timestamps[m]
                value = values[m]
                old_start_time = timestamp
                old_end_time = timestamp + original_time_step
        
                # Calculate the overlapping duration
                overlap_start = max(new_start_time, old_start_time)
                overlap_end = min(new_end_time, old_end_time)
                overlap_duration = max(0, overlap_end - overlap_start)
        
                # Calculate the proportion of overlap and distribute the value
                if overlap_duration > 0
                    overlap_fraction = overlap_duration / original_time_step
                    converted_profile[i] += value * overlap_fraction
                end

                if new_end_time <= old_start_time
                    begin_search = m-1
                    break
                end

            end     
        end
        return converted_profile
    end
 end

# Function to convert intensive profiles to the simulation time step (e.g., temperature, power)
function convert_intensive_profile(values, timestamps, original_time_step, new_time_step)
    if new_time_step == original_time_step # no change
        return values
    elseif new_time_step < original_time_step  # Segmentation
        interp = interpolate((timestamps,), values, Gridded(Linear()))
        new_timestamps = minimum(timestamps):new_time_step:maximum(timestamps)
        converted_profile = [interp(t) for t in new_timestamps]
        # add missing entries by coping the last entry
        for _ in 1:ceil(Int,(original_time_step-new_time_step)/new_time_step)
            append!(converted_profile, converted_profile[end])
        end
        return converted_profile
    else # Aggregation
        # Note: If the original_time_step is not a divider of the new_time_step, there may be small errors
        #       as there is no consideration of overlapping time steps here 
        # Determine the length of the new profile
        new_length = ceil(Int, maximum(timestamps) / new_time_step + original_time_step / new_time_step)
        converted_profile = zeros(new_length)

        sum_values = 0.0
        count_values = 0
        j = 1

        for i in eachindex(values)
            # Aggregate values within the new timestep
            if timestamps[i] < j * new_time_step
                sum_values += values[i]
                count_values += 1
            else
                # Calculate average for the previous timestep
                if count_values != 0
                    converted_profile[j] = sum_values / count_values
                end
                j += 1
                sum_values = values[i]
                count_values = 1
            end
        end

        # Handle the last set of values
        if count_values != 0
            converted_profile[j] = sum_values / count_values
        end
        return converted_profile
    end
end

end
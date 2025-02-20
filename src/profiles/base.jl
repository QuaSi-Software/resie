module Profiles
using Interpolations
using Dates, TimeZones

export Profile, power_at_time, work_at_time, value_at_time, remove_leap_days,
       add_ignoring_leap_days, sub_ignoring_leap_days

"""
Holds values from a file so they can be retrieved later and indexed by time.

Profiles are automatically aggregated or segmentated to fit the simulation time step.
Currtently, this only works, if the time step of the profile is a multiple or a divisor 
of the requested simulation timestep. Otherwise, an error will arise.

This function can handle either a path to a .prf file, or the profile data can be
handed over as vectors. Then, the profile with the corresponding timestamps, the time step
and the data_type (intensive/extensive) has to be provided as optional arguments. 

The data_type is essential to make sure that the profile data is converted correctly.
It has to be set to "intensive" for values like temperatures, power, wind speed or
for profiles containting states. Only for energies, data_type has to be set to "extensive".

Data should represent the time step following the time indicated (mean/sum). The shift parameter
can be used to shift the timestamp to the correct time. A positive shift will add to the timestamp,
which means a value is later in time. 

The interpolation type can be one of "stepwise", "linear_classic", "linear_time_preserving" or 
"linear_solar_radiation".
- "stepwise" means, the value given at a timestamp is spread/copied over the new finer timestamps. 
- "linear_classic" interpolates the given values from the original timestamps linearly to the new timestamps.
- "linear_time_preserving" interpolates the data by first shifting the data by half the original 
time step to make the values be measured at the time indicated. After the interpolation to a 
finer time step, the data is shifted back by 1/2 a time step to meet the required definition 
of the values representing the following time step. This should be used for time-critic data 
like solar radiation. But, this method will cut peaks and valleys in the data more than the 
classic interpolation.
- "linear_solar_radiation" interpolation uses a method described in the paper "A new method for
interpolating hourly solar radiation data" by M. A. S. Mohandes, A. Halawani, and A. Rehman.
It is a linear interpolation with a correction factor to keep the sum of the interpolated values
equal to the sum of the original values. This is also used in TRNSYS 18, but shows "wavy" curves
as result.

"""
mutable struct Profile
    """Time step, in seconds, of the profile."""
    time_step::Int

    """Indicates whether the profile values are power/temperature ("intensive") or work ("extensive")."""
    data_type::Symbol

    """Holds the profile values indexed by the time step number"""
    data::Dict{DateTime,Float64}

    """Construct a profile from the file at the given path or convert given data to profile (optional)."""
    function Profile(file_path::String,                               # file path to a .prf file
                     sim_params::Dict{String,Any};                    # general simulation parameters
                     given_profile_values::Vector{Float64}=Float64[], # optional: Vector{Float64} that holds values of 
                     #                                                            the profile
                     given_timestamps::Vector{DateTime}=DateTime[],   # optional: Vector{DateTime} that holds the 
                     #                                                            time step as DateTime
                     given_time_step::Dates.Second=Second(0),         # optional: Dates.Second that indicates the 
                     #                                                            timestep in seconds of the given data
                     given_data_type::Union{String,Nothing}=nothing,  # optional: datatype, shoule be "intensive" or
                     #                                                            "extensive"
                     shift::Dates.Second=Second(0),                   # optional: timeshift for data. A positive shift 
                     #                                                            adds to the timestamp = value is later
                     interpolation_type::String="stepwise")           # optional: interpolation type. Can be one of: "stepwise",
        #                                                                         "linear_classic", "linear_time_preserving"
        if given_profile_values == []  # read data from file_path
            profile_values = Vector{Float64}()
            profile_timestamps = Vector{String}()

            profile_time_step = nothing
            data_type = nothing
            timestamp_format = nothing
            time_definition = nothing
            profile_start_date = nothing
            profile_start_date_format = nothing
            time_zone = nothing
            time_shift_seconds = nothing

            file_handle = nothing

            try
                file_handle = open(abspath(file_path), "r")
                for (line_idx, line) in enumerate(readlines(file_handle))
                    line = strip(line)

                    if isempty(line) || length(line) < 1 # handle empty lines
                        continue
                    elseif line[1] == '#'
                        splitted = split(strip(line, '#'), ':'; limit=2)
                        if strip(splitted[1]) == "profile_time_step_seconds"
                            profile_time_step = parse(Int, String(strip(splitted[2])))
                        elseif strip(splitted[1]) == "data_type"
                            data_type = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "timestamp_format"
                            timestamp_format = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "time_definition"
                            time_definition = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "profile_start_date"
                            profile_start_date = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "profile_start_date_format"
                            profile_start_date_format = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "time_zone"
                            time_zone = String(strip(splitted[2]))
                        elseif strip(splitted[1]) == "time_shift_seconds"
                            time_shift_seconds = parse(Int, String(strip(splitted[2])))
                        elseif strip(splitted[1]) == "interpolation_type"
                            interpolation_type = String(strip(splitted[2]))
                        end
                    else
                        splitted = split(line, ';')
                        if length(splitted) == 1
                            push!(profile_values, parse(Float64, String(strip(splitted[1]))))
                        elseif length(splitted) == 2
                            push!(profile_timestamps, String(strip(splitted[1])))
                            push!(profile_values, parse(Float64, String(strip(splitted[2]))))
                        else
                            @error "The profile at $(file_path) has some corrupt data in line $(line_idx). Please check!"
                            close(file_handle)
                            throw(InputError)
                        end
                    end
                end
            catch e
                @error "While reading the data from the profile at $(file_path) the following error occured: $e\n" *
                       "Check if the file exists at the specified path!"
                throw(InputError)
            finally
                if file_handle !== nothing
                    close(file_handle)
                end
            end

            profile_timestamps_date = Vector{DateTime}(undef, length(profile_values))
            # convert profile_timestamps to DateTime
            if time_definition == "startdate_timestepsize"
                if isnothing(profile_start_date) || isnothing(profile_start_date_format) || isnothing(profile_time_step)
                    @error "For the profile at $(file_path) a time definition via startdate_timestepsize is chosen. " *
                           "Specify 'profile_start_date', 'profile_start_date_format' and 'profile_time_step_seconds' " *
                           "in the profile header!"
                    throw(InputError)
                else
                    # create time step from startdate and timestep
                    if length(profile_timestamps) > 0
                        @error "For the profile at $(file_path) a time stamp is given, but 'startdate_timestepsize' " *
                               "is chosen as 'time_definition'. This expects to not have a time stamp! Chose another " *
                               "'time_definition'."
                        throw(InputError)
                    end
                    start_date = parse_datestamp(profile_start_date,
                                                 convert_String_to_DateFormat(profile_start_date_format, file_path),
                                                 file_path,
                                                 "profile_start_date_format")
                    for idx in eachindex(profile_values)
                        profile_timestamps_date[idx] = add_ignoring_leap_days(start_date,
                                                                              Second((idx - 1) * profile_time_step))
                    end
                    if time_shift_seconds !== nothing
                        shift = Second(time_shift_seconds)
                    end
                end
            elseif time_definition == "startdate_timestamp"
                if isnothing(profile_start_date) || isnothing(profile_start_date_format) || isnothing(timestamp_format)
                    @error "For the profile at $(file_path) a time definition via startdate_timestamp is chosen. " *
                           "Specify 'profile_start_date', 'profile_start_date_format' and 'timestamp_format' in the " *
                           "profile header!"
                    throw(InputError)
                else
                    # calculate time step from startdate and continuos timestamp
                    start_date = parse_datestamp(profile_start_date,
                                                 convert_String_to_DateFormat(profile_start_date_format, file_path),
                                                 file_path,
                                                 "profile_start_date_format")
                    if timestamp_format == "seconds"
                        for (idx, entry) in enumerate(profile_timestamps)
                            profile_timestamps_date[idx] = add_ignoring_leap_days(start_date,
                                                                                  Second(parse(Int64, entry)))
                        end
                    elseif timestamp_format == "minutes"
                        for (idx, entry) in enumerate(profile_timestamps)
                            profile_timestamps_date[idx] = add_ignoring_leap_days(start_date,
                                                                                  Minute(parse(Int64, entry)))
                        end
                    elseif timestamp_format == "hours"
                        for (idx, entry) in enumerate(profile_timestamps)
                            profile_timestamps_date[idx] = add_ignoring_leap_days(start_date, Hour(parse(Int64, entry)))
                        end
                    else
                        @error "For the profile at $(file_path) the 'timestamp_format' has to be one of 'seconds', " *
                               "'minutes' or 'hours'!"
                        throw(InputError)
                    end
                    if time_shift_seconds !== nothing
                        shift = Second(time_shift_seconds)
                    end
                end
            elseif time_definition == "datestamp"
                if isnothing(timestamp_format)
                    @error "For the profile at $(file_path) a time definition via datestamp is chosen. " *
                           "Specify 'timestamp_format' in the profile header!"
                    throw(InputError)
                else
                    # get time step from datestamp
                    df = convert_String_to_DateFormat(timestamp_format, file_path)
                    for (idx, entry) in enumerate(profile_timestamps)
                        profile_timestamps_date[idx] = parse_datestamp(entry,
                                                                       df,
                                                                       file_path,
                                                                       "timestamp_format")
                    end
                end
                if time_shift_seconds !== nothing
                    shift = Second(time_shift_seconds)
                end
            else
                @error "For the profile at $(file_path) no 'time_definition' is chosen. " *
                       "It has to be one of 'startdate_timestepsize', 'startdate_timestamp' or 'datestamp'."
                throw(InputError)
            end

            if isnothing(profile_time_step) && length(profile_timestamps_date) > 1
                profile_time_step = Dates.value(Second(profile_timestamps_date[2] - profile_timestamps_date[1]))
            elseif isnothing(profile_time_step)
                @error "For the profile at $(file_path) no profile time step could be determined. " *
                       "Please specify one with the parameter 'profile_time_step_seconds'."
                throw(InputError)
            end
        else # data was read in from somewhere else and is provided as vector
            profile_values = given_profile_values
            profile_time_step = Dates.value(given_time_step)
            profile_timestamps_date = given_timestamps
            data_type = given_data_type
            time_zone = nothing
        end

        # check interpolation_type
        if !(interpolation_type in ["stepwise", "linear_classic", "linear_time_preserving", "linear_solar_radiation"])
            @error "The interpolation type of the profile at $(file_path) has to be one of " *
                   "'stepwise', 'linear_classic', 'linear_time_preserving' or 'linear_solar_radiation'! " *
                   "The given interpolation type is '$interpolation_type'."
            throw(InputError)
        end

        # remove leap days and daiylight savings from profile
        profile_timestamps_date, profile_values = remove_leap_days(profile_timestamps_date, profile_values, file_path)
        profile_timestamps_date = remove_day_light_saving(profile_timestamps_date, time_zone, file_path)

        # shift the profile timestep according to the given shift to get the correct definition:
        # Values are given as the mean/sum over the upcoming time step.
        if shift > Second(0)
            profile_timestamps_date = add_ignoring_leap_days.(profile_timestamps_date, shift)
            @info "The timestamp of the profile at $(file_path) was shifted by $shift."
        elseif shift < Second(0)
            profile_timestamps_date = sub_ignoring_leap_days.(profile_timestamps_date, -shift)
            @info "The timestamp of the profile at $(file_path) was shifted by $shift."
        end

        if interpolation_type == "linear_time_preserving" && profile_time_step > sim_params["time_step_seconds"]
            # shift data by half the original time step to make the values be measured at the time indicated. 
            # This is required for correct linear interpolation during segmentation and will be reversed afterwards.
            profile_timestamps_date = add_ignoring_leap_days.(profile_timestamps_date, Second(profile_time_step / 2))
        end

        # add first or last entry if only one time step is missing. copy first or last value.
        if profile_timestamps_date[1] > sim_params["start_date"] &&
           profile_timestamps_date[1] - Second(profile_time_step) <= sim_params["start_date"]
            profile_timestamps_date = vcat(profile_timestamps_date[1] - Second(profile_time_step),
                                           profile_timestamps_date)
            profile_values = vcat(profile_values[1], profile_values)
            @info "The profile at $(file_path) has been extended by one timestep at the begin by doubling the fist value."
        end
        temp_diff = Second(max(profile_time_step, sim_params["time_step_seconds"]))
        time_diff_timesteps = Second(max(0, Int64(sim_params["time_step_seconds"]) - profile_time_step))
        if profile_timestamps_date[end] < sim_params["end_date"] + time_diff_timesteps &&
           profile_timestamps_date[end] + temp_diff >= sim_params["end_date"]
            while profile_timestamps_date[end] < sim_params["end_date"] + time_diff_timesteps
                profile_timestamps_date = vcat(profile_timestamps_date,
                                               profile_timestamps_date[end] + Second(profile_time_step))
                profile_values = vcat(profile_values, profile_values[end])
            end
            @info "The profile at $(file_path) has been extended by one timestep at the end by doubling the last value."
        end

        # check timestamp for order, duplicates and consistency
        if !issorted(profile_timestamps_date)
            @error "The timestamp of the profile at $(file_path) is not sorted in ascending order!"
            throw(InputError)
        elseif length(profile_timestamps_date) !== length(unique(profile_timestamps_date))
            @error "The timestamp of the profile at $(file_path) has duplicates!"
            throw(InputError)
        elseif length(profile_timestamps_date) > 1 &&
               length(unique(diff_ignore_leap_days(profile_timestamps_date))) !== 1
            @error "The timestamp of the profile at $(file_path) has an inconsistent time step width! " *
                   "If the profile is defined by a datestamp and has daylight savings, please specify a 'time_zone'!" *
                   "If the profile has no daylight savings, a 'time_zone' must not be given!"
            throw(InputError)
        end

        # check time step compatibility
        if !(sim_params["time_step_seconds"] % profile_time_step == 0) &&
           !(profile_time_step % sim_params["time_step_seconds"] == 0)
            @error ("The timestep of the profile at $(file_path), which is $(string(profile_time_step)),\n" *
                    "is not a multiple or a divisor of the requested simulation timestep of " *
                    "$(string(sim_params["time_step_seconds"]))!")
            throw(InputError)
        end

        # check coverage of simulation time
        if profile_timestamps_date[1] > sim_params["start_date"] ||
           profile_timestamps_date[end] < sim_params["end_date"]
            @error "In the profile at $(file_path), the simulation time is not fully covered!\n" *
                   "Provide a profile that covers the simulation time of:\n" *
                   "Begin: $(Dates.format(sim_params["start_date"], "yyyy-mm-dd HH:MM:SS"))\n" *
                   "End: $(Dates.format(sim_params["end_date"], "yyyy-mm-dd HH:MM:SS"))\n" *
                   "The current profile only covers (in local standard time):\n" *
                   "Begin: $(Dates.format(profile_timestamps_date[1], "yyyy-mm-dd HH:MM:SS"))\n" *
                   "End: $(Dates.format(profile_timestamps_date[end], "yyyy-mm-dd HH:MM:SS"))"
            throw(InputError)
        end

        if data_type == "intensive"         # e.g., temperature, power
            data_type = :intensive
        elseif data_type == "extensive"     # e.g., energy demand
            data_type = :extensive
        else
            @error "For the profile at $(file_path) no valid 'data_type' is given! " *
                   "Has to be either 'intensive' or 'extensive'."
            throw(InputError)
        end

        values_converted,
        profile_timestamps_date_converted = convert_profile(profile_values,
                                                            profile_timestamps_date,
                                                            Millisecond(Int(1000 * profile_time_step)),
                                                            Millisecond(1000 * sim_params["time_step_seconds"]),
                                                            data_type,
                                                            file_path,
                                                            sim_params,
                                                            interpolation_type)

        profile_dict = Dict(zip(profile_timestamps_date_converted, values_converted))

        return new(sim_params["time_step_seconds"],   # Period [s]: time_step, equals simulation time step after conversion
                   data_type,                         # String: intensive or extensive profile data
                   profile_dict)                      # Dict{DateTime, Float64}() dict with timestamp as key and data of 
        #                                               profile, in simulation time step
    end
end

"""
    convert_String_to_DateFormat()

function to pase a string to a a DateFormat.
    datetime_str::String      A DateFormat as string
    file_path::String         File path of the profile for error handling, for error message.

Returns
    datetime_dt::DateFormat   The parsed DateFormat
"""
function convert_String_to_DateFormat(datetime_str::String, file_path::String)::DateFormat
    try
        return DateFormat(datetime_str)
    catch e
        @error "Time given date specifier '$(datetime_str)' of profile at $(file_path) could not be converted " *
               "to a DateFormat. Make shure, is matches the format specification of the Julia Dates package. " *
               "The following error occured: $e"
        throw(InputError)
    end
end

"""
    parse_datestamp()

function to pase a datestamp from a string to a given date format.
    datetime_str::String        Timestamp read from the profile file. Should be a specified DateTime format as string.
    time_format::DateFormat     Format of datetime_str. Shoule be a DateTime format, like "dd-mm-yyyy HH:MM:SS"
    file_path::String           File path of the profile for error handling, for error message.
    time_format_variable_name::String   The variable name of the time format, for error message.

Returns
    parsed_timestep::DateTime   The parsed timestep as DateTime
"""
function parse_datestamp(datetime_str::String,
                         time_format::DateFormat,
                         file_path::String,
                         time_format_variable_name::String)
    try
        return Dates.DateTime(datetime_str, time_format)
    catch e
        @error "Time given date specifier of profile at $(file_path) does not fit to the data. " *
               "'$(time_format_variable_name)' has to match the dataformat, e.g. 'dd-mm-yyyy HH:MM:SS'. " *
               "'$(time_format_variable_name)' is `$(time_format)` which does not fit to the " *
               "timestamp given `$(datetime_str)`.\n The following error occured: $e"
        throw(InputError)
    end
end

"""
    remove_leap_days(timestamp::Vector{DateTime}, values::Vector{Float64}, file_path::String)

Detect and remove a leap day from the timestep and the values.
"""
function remove_leap_days(timestamp::Vector{DateTime}, values::Vector{Float64}, file_path::String)
    filtered_pairs = [(dt, val) for (dt, val) in zip(timestamp, values) if !(month(dt) == 2 && day(dt) == 29)]
    filtered_timestamp = [dt for (dt, val) in filtered_pairs]
    filtered_values = [val for (dt, val) in filtered_pairs]

    if length(filtered_timestamp) < length(timestamp)
        @info "In the profile at $file_path, a leap day has been removed."
    end

    return filtered_timestamp, filtered_values
end

"""
    remove_leap_days(timestamp::Vector{DateTime})

Detect and remove a leap day from a timestep.
"""
function remove_leap_days(timestamp::Vector{DateTime})
    return [dt for dt in timestamp if !(month(dt) == 2 && day(dt) == 29)]
end

"""
    add_ignoring_leap_days(timestamp::DateTime, diff::DateTime.Period)

Adds diff to the timestamp. If the result is a leap day, skip to the next day.
If a leap day is in between, do not count it!
"""
function add_ignoring_leap_days(timestamp::DateTime, diff::Period)
    new_time = timestamp + diff
    nr_of_leap_days = 0
    for year in year(timestamp):year(new_time)
        if isleapyear(year)
            leap_day = Date(year, 2, 29)
            if leap_day >= Date(timestamp) && leap_day <= Date(new_time)
                nr_of_leap_days += 1
            end
        end
    end

    return new_time + nr_of_leap_days * Day(1)
end

"""
    sub_ignoring_leap_days(timestamp::DateTime, diff::DateTime.Period)

Subs diff of the timestamp. If the result is a leap day, skip to the previous day.
If a leap day is in between, do not count it!
"""
function sub_ignoring_leap_days(timestamp::DateTime, diff::Period)
    new_time = timestamp - diff
    nr_of_leap_days = 0
    for year in year(new_time):year(timestamp)
        if isleapyear(year)
            leap_day = Date(year, 2, 29)
            if leap_day <= Date(timestamp) && leap_day >= Date(new_time)
                nr_of_leap_days += 1
            end
        end
    end

    return new_time - nr_of_leap_days * Day(1)
end

"""
    sub_ignoring_leap_days(end_date::DateTime, begin_date::DateTime)

Substract: DateTime.Period = end_date - begin_date
If the result is a leap day, skip to the next day.
If a leap day is in between, do not count it!
"""
function sub_ignoring_leap_days(end_date::DateTime, begin_date::DateTime)
    time_span = end_date - begin_date
    nr_of_leap_days = 0
    for year in year(begin_date):year(end_date)
        if isleapyear(year)
            leap_day = Date(year, 2, 29)
            if leap_day >= Date(begin_date) && leap_day <= Date(end_date)
                nr_of_leap_days += 1
            end
        end
    end

    return time_span - Millisecond(nr_of_leap_days * 24 * 60 * 60 * 1000)
end

"""
    diff_ignore_leap_days(timestamps::Vector{DateTime})

Calculates all differences of consecutive elements of the timestep vector
while considering removed leap days. Returns the differences as vector.
"""
function diff_ignore_leap_days(timestamps::Vector{DateTime})
    number_of_timesteps = length(timestamps)
    diffs = Vector{Period}(undef, number_of_timesteps - 1)

    for i in 2:number_of_timesteps
        diff = timestamps[i] - timestamps[i - 1]

        # Check if the range includes a leap day and adjust the difference
        if isleapyear(timestamps[i]) &&
           month(timestamps[i - 1]) <= 2 && day(timestamps[i - 1]) <= 28 &&
           month(timestamps[i]) >= 3 && day(timestamps[i]) >= 1
            diff = diff - Millisecond(86400000)  # equals 1 day
        end
        diffs[i - 1] = diff
    end

    return diffs
end

"""
    remove_day_light_saving(timestamp::Vector{DateTime}, time_zone::Stringing, file_path::String)

If DST is present in the timestamp, shift it to local standard time.

This may require a recompilation of the TZData of the Julia package TimeZones. From first installation,
if only covers DST shifts until the year 2037. 

If the recompilation fails, try to follow these steps:
- create the following folder, if it doesn't exist:
  "user_home_directory"/.julia/scratchspaces/"TimeZone package UUID"/build/tzsource/2024a
  The Universally Unique Identifier (UUID) of your TimeZone package installation can be found in the Project.toml file. 
  Currently, for the TimeZone package, it is f269a46b-ccf7-5d73-abea-4c690281aa53.
- download the latest tz data from https://www.iana.org/time-zones (Data only)
- unzip the data
- delete all files except the ones named as continent, e.g. "europe", "asia"
- create a destination folder for the recompilation, if it doesn't exist:
  "user_home_directory"/.julia/scratchspaces/"TimeZone package UUID"/build/compiled/tzjf/v1/2024a
- re-run the compilation with the required max_year, e.g. TimeZones.TZData.compile(max_year=2200)

"""
function remove_day_light_saving(timestamp::Vector{DateTime}, time_zone::Union{Nothing,String}, file_path::String)
    if time_zone === nothing
        return timestamp
    else
        tz = ""
        try
            tz = TimeZone(time_zone)
        catch e
            @error "In the profile at $file_path, the given time zone $(time_zone) could not be detected. " *
                   "It has to be an IANA (also known as TZ or TZDB) time zone identifier."
            throw(InputError)
        end

        timestamp_max_year = year(timestamp[end])
        if year(tz.transitions[end].utc_datetime) < timestamp_max_year
            try
                TimeZones.TZData.compile(; max_year=Int(timestamp_max_year))
            catch e
                timezones_uuid = Base.identify_package("TimeZones").uuid
                @error "The Julia package TimeZones could not be recompiled to cover the DST for the requested year " *
                       "$timestamp_max_year. The currently supported latest year is " *
                       "$(year(tz.transitions[end].utc_datetime)). The following error occured: $e\n" *
                       "If the recompilation fails, try to follow these steps: \n" *
                       "- create the following folder, if it doesn't exist:" *
                       "'user_home_directory'/.julia/scratchspaces/$(timezones_uuid)/build/tzsource/2024a \n" *
                       "- download the latest tz data from https://www.iana.org/time-zones (Data only)\n" *
                       "- unzip the data\n" *
                       "- delete all files except the ones named as continent, e.g. 'europe', 'asia' \n" *
                       "- create a destination folder for the recompilation, if it doesn't exist:" *
                       "'user_home_directory'/.julia/scratchspaces/$(timezones_uuid)/build/compiled/tzjf/v1/2024a \n" *
                       "- re-run the compilation by running ReSiE or with TimeZones.TZData.compile(max_year=2200) in Julia."
                throw(InputError)
            end
            # re-creation of the TimeZone is required to apply the recompilation
            tz = TimeZone(time_zone)
        end

        corrected_timestamps = DateTime[]
        has_corrected = false
        zoned_time = 0

        for ts in timestamp
            # Convert the timestamp to ZonedDateTime and back to DateTime to eliminate any DST effects.
            # This results in local standard time.
            try
                zoned_time = ZonedDateTime(ts, tz)
                push!(corrected_timestamps, sub_ignoring_leap_days(DateTime(zoned_time), zoned_time.zone.offset.dst))
            catch e
                @error "In the profile at $file_path, the datestamp $ts is probably invalid for the specified time " *
                       "zone $time_zone. The following error occured: $e"
                throw(InputError)
            end

            if zoned_time.zone.offset.dst !== Second(0) && !has_corrected
                has_corrected = true
            end
        end

        if has_corrected
            @info "In the profile at $file_path, the daylight saving has been converted to local standard time."
        end

        return corrected_timestamps
    end
end

"""
    power_at_time(profile, time)

Get the power value of the profile at the given time.
"""
function power_at_time(profile::Profile, sim_params::Dict{String,Any})
    if profile.data_type == :intensive
        return profile.data[sim_params["current_date"]]
    elseif profile.data_type == :extensive
        return profile.data[sim_params["current_date"]] * (3600 / profile.time_step)
    else
        @error "For current profile has no information on the 'data_type'! " *
               "Has to be either 'intensive' or 'extensive'."
    end
end

"""
    work_at_time(profile, time)

Get the work value of the profile at the given time.
"""
function work_at_time(profile::Profile, sim_params::Dict{String,Any})
    if profile.data_type == :intensive
        return profile.data[sim_params["current_date"]] * (profile.time_step / 3600)
    elseif profile.data_type == :extensive
        return profile.data[sim_params["current_date"]]
    else
        @error "For current profile has no information on the 'data_type'! " *
               "Has to be either 'intensive' or 'extensive'."
    end
end

"""
    value_at_time(profile, time)

Get the value of the profile at the given time without any conversion.
The data_type will be ignored.
"""
function value_at_time(profile::Profile, sim_params::Dict{String,Any})
    return profile.data[sim_params["current_date"]]
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
    return Base.minimum((values(profile.data)))
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
    return Base.maximum((values(profile.data)))
end

"""
    convert_profile

Function to convert intensive and extensive profiles from the profile time step to the simulation time step.
Handles all cases for aggregation and segmentation (linear and stepwise) as well as required time shifts.

Inputs: 
    values::Vector{Float64}         values of the profile to convert
    timestamps::Vector{Int64}       timestamps in seconds corresponding to values
    original_time_step::Int         the timestep in seconds of the values and the timestamps
    new_time_step::UInt64           the simulation time step to wich the profile should be converted
    profile_type::Symbol            the profile type, can be :extensive (e.g. energies) or :intensive (e.g. temperatures)        
    file_path::String               the file path of the profile, for error messages only
    sim_params::Dict{String,Any}    simulation parameters
    interpolation_type::String      interpolation type. Can be one of "stepwise", "linear_classic" or "linear_time_preserving"
                                    For time-critic data like solar radiation, use "linear_time_preserving".

Outputs:
    values::Vector{Float64}         values of the converted profile in the `new_time_step`
    timestamp::Vector{DateTime}     corresponding new timestep

If `original_time_step` equals the `new_time_step`, the original `values` will be returned.
"""
function convert_profile(values::Vector{Float64},
                         timestamps::Vector{DateTime},
                         original_time_step::Period,
                         new_time_step::Period,
                         profile_type::Symbol,
                         file_path::String,
                         sim_params::Dict{String,Any},
                         interpolation_type::String)

    # construct info message
    info_message = "The profile at $(file_path) "

    time_is_aligned = sim_params["start_date"] in timestamps

    if original_time_step == new_time_step && time_is_aligned   # do nothing
        return values, timestamps

    elseif original_time_step == new_time_step                  # time shift only
        values, timestamps = profile_linear_interpolation(values,
                                                          timestamps,
                                                          original_time_step,
                                                          new_time_step,
                                                          sim_params)
        info_message *= "was shifted using linear interpolation to fit the simulation start time. "

    elseif original_time_step < new_time_step                   # aggregation
        if profile_type == :extensive
            values .*= new_time_step / original_time_step
            info_message *= "(extensive profile) "
        else
            info_message *= "(intensive profile) "
        end

        if !time_is_aligned  # do time shift
            values, timestamps = profile_linear_interpolation(values,
                                                              timestamps,
                                                              original_time_step,
                                                              original_time_step,
                                                              sim_params)
            info_message *= "was shifted using linear interpolation to fit the simulation start time and "
        end

        values, timestamps = profile_calculate_means(values,
                                                     timestamps,
                                                     original_time_step,
                                                     new_time_step,
                                                     sim_params)
        info_message *= "was converted from the profile timestep $(Second(original_time_step)) to the simulation " *
                        "timestep of $(new_time_step)."

    elseif original_time_step > new_time_step              # segmentation
        if profile_type == :extensive
            values .*= new_time_step / original_time_step
            info_message *= "(extensive profile) "
        else
            info_message *= "(intensive profile) "
        end

        if interpolation_type == "linear_solar_radiation"   # segmentation of solar radiation data
            if original_time_step !== Millisecond(1000 * 60 * 60)
                @error "For the profile at $(file_path) with solar radiation data, the original time step has to be " *
                       "1 hour (3600 seconds) to use the linear_solar_radiation method for interpolation!"
                throw(InputError)
            end
            if !time_is_aligned
                @error "For the profile at $(file_path) with solar radiation data, the original time step has to be " *
                       "aligned to the simulation start time to use the linear_solar_radiation method for interpolation!"
                throw(InputError)
            end

            values, timestamps = segment_profile(timestamps[1],
                                                 timestamps[end],
                                                 values,
                                                 new_time_step,
                                                 sim_params)

        elseif interpolation_type == "linear_classic"
            values, timestamps = profile_linear_interpolation(values,
                                                              timestamps,
                                                              original_time_step,
                                                              new_time_step,
                                                              sim_params)

        elseif interpolation_type == "linear_time_preserving"
            values, timestamps = profile_linear_interpolation(values,
                                                              timestamps,
                                                              original_time_step,
                                                              new_time_step / 2,
                                                              sim_params)
            # remove every second timestamp and value to match the definition of values representing
            # the mean/sum of the following time step
            deleteat!(values, 1:2:(length(values) - 1))
            deleteat!(timestamps, 2:2:length(timestamps))

        else # stepwise conversion
            values, timestamps = profile_spread_to_segments(values,
                                                            timestamps,
                                                            original_time_step,
                                                            new_time_step,
                                                            sim_params)
            if !time_is_aligned   # do time shift if not aligned
                values, timestamps = profile_linear_interpolation(values,
                                                                  timestamps,
                                                                  new_time_step,
                                                                  new_time_step,
                                                                  sim_params)
            end
        end

        if !time_is_aligned
            info_message *= "was shifted to fit the simulation start time and "
        end

        info_message *= "was converted from the profile timestep $(Second(original_time_step)) to the simulation " *
                        "timestep of $(Second(new_time_step)) using $(interpolation_type) interpolation."
    end

    @info info_message

    return values, timestamps
end

"""
    segment_interval() 
segments one hourly interval given as fractional times (relative to midnight) using the 
algorithm of TRNSYS 18, described in:
    Letellier-Duchesne, Samuel & McDowell, Timothy & Kummert, Michael. (2018). 
    A New Method for Determining Sub-hourly Solar Radiation from Hourly Data. 

Parameters:
  • Hn       : hourly integrated irradiation for the current hour [Wh/m²]
  • Hnp1     : hourly integrated irradiation for the next hour [Wh/m²]
  • interval_start, interval_end : start and end of the current hourly interval
      (in fractional hours, e.g., 7.24 or 8.0)
  • dt_h     : required simulation timestep (in hours)
  • G_start  : carry-over instantaneous radiation from the previous interval (W/m²)
  • sunrise, sunset : daylight bounds (fractional hours)

Returns (t_sub, irr_sub, G_end) where:
  • t_sub is a vector of fractional hour values (currently computed as startpoints)
  • irr_sub is the average irradiance (W/m²) for each subinterval,
  • G_end is the endpoint to carry over.
"""
function segment_interval(Hn::Float64, Hnp1::Float64,
                          interval_start::Float64, interval_end::Float64,
                          dt_h::Float64, G_start::Float64,
                          sunrise::Float64, sunset::Float64)
    # 1. Divide the full hourly interval into dt blocks (using ceil to ensure full coverage)
    N_total = max(Int(ceil((interval_end - interval_start) / dt_h)), 1)
    t_blocks = [interval_start + (j - 1) * dt_h for j in 1:N_total]

    # 2. Define the effective daylight region.
    eff_start = max(interval_start, sunrise)
    eff_end = min(interval_end, sunset)
    if eff_end <= eff_start
        # No daylight: all dt blocks get zero irradiance.
        return t_blocks, zeros(Float64, N_total), 0.0
    end
    dt_eff = eff_end - eff_start
    is_sunrise = eff_start == sunrise
    is_sunset = eff_end == sunset

    # 3. Subdivide the effective region into dt blocks (again using ceil)
    N_eff = max(Int(ceil(dt_eff / dt_h)), 1)

    # 4. Compute the effective endpoint.
    # If the hour includes sunset, force the endpoint to zero.
    if (sunset > interval_start && sunset < interval_end)
        effective_G_end = 0.0
    else
        if Hnp1 > Hn && is_sunrise
            effective_G_end = (2 * Hn) / dt_eff
        elseif Hnp1 < Hn && is_sunset
            effective_G_end = (0.25 * Hn + 0.75 * Hnp1) / dt_eff
        else
            effective_G_end = (0.5 * Hn + 0.5 * Hnp1) / dt_eff
        end
    end

    # 5. Compute the midpoint value (G_mid) over the effective region.
    if N_eff == 1
        G_mid = Hn / dt_h
    else
        N_mid = Int(floor(N_eff / 2))
        G_mid = (2 * Hn / (N_eff * dt_h)) - (N_mid * G_start) / N_eff - ((N_eff - N_mid) * effective_G_end) / N_eff
        if G_mid < 0
            G_mid = (2 * Hn / dt_h - (G_start + effective_G_end)) / (2 * (N_eff - 1))
        end
        if G_mid < 0
            G_mid = 0
        end
    end

    # 6. Build the piecewise–linear boundaries over the effective region.
    boundaries = zeros(Float64, N_eff + 1)
    boundaries[1] = G_start
    if N_eff > 1
        N_mid = Int(floor(N_eff / 2))
        # First part: interpolate from G_start to G_mid over N_mid subintervals.
        for j in 1:N_mid
            boundaries[j + 1] = G_start + (j / N_mid) * (G_mid - G_start)
        end
        # Second part: interpolate from G_mid to effective_G_end over the remaining subintervals.
        #              effective_G_end is considered to be the first value of the next interval.
        M = N_eff - N_mid
        for k in 1:M
            boundaries[N_mid + 1 + k] = G_mid + (k / M) * (effective_G_end - G_mid)
        end
    else
        boundaries[2] = G_mid
    end

    # 7. Compute the effective irradiance values (midpoints) over the effective region.
    effective_irr = [(boundaries[i] + boundaries[i + 1]) / 2 for i in 1:(length(boundaries) - 1)]

    # 8. fill effective_irr with zeros at the beginning or end, depending if we are during 
    #    sunset or sunrise
    if N_eff < N_total
        if is_sunrise
            effective_irr = vcat(zeros(Float64, N_total - N_eff), effective_irr)
        elseif is_sunset
            effective_irr = vcat(effective_irr, zeros(Float64, N_total - N_eff))
        end
    end

    return t_blocks, effective_irr, effective_G_end
end

"""
 segment_profile() processes a long period (e.g. a full year) given:
 
   • begin_dt, end_dt : DateTime objects defining the overall period.
   • hourly_H         : Vector of hourly integrated irradiation values (in Wh/m²)
                        that correspond sequentially to hourly intervals starting at begin_dt.
   • dt               : a simulation timestep provided as a Millisecond period.
  
 This function iterates over each hourly interval (which may span many days).
 For each interval:
   - It computes the effective start and end (in fractional hours relative to midnight)
     taking into account partial hours at the beginning or end.
   - It obtains that day's sunrise and sunset.
   - It segments the interval and converts the fractional hours back to DateTime.
 
 In addition, the returned timestamps are shifted so that each timestamp represents
 the beginning of the upcoming dt period (i.e. the mean value is associated with the start of the interval).
"""
function segment_profile(begin_dt::DateTime,
                         end_dt::DateTime,
                         hourly_H::Vector{Float64},
                         dt::Millisecond,
                         sim_params::Dict{String,Any})
    dt_h = dt / Hour(1)  # Convert dt to hours (Float64)

    # Prepare containers for results.
    times_total = DateTime[]
    irr_total = Float64[]

    # Helper: compute fractional hour (since midnight) from a DateTime.
    function fractional_hour(dt_obj::DateTime)
        Dates.hour(dt_obj) +
        Dates.minute(dt_obj) / 60 + Dates.second(dt_obj) / 3600 + Dates.millisecond(dt_obj) / 3600000
    end

    n_hours = length(hourly_H)  # total number of hourly intervals
    G_start = 0.0  # carry-over start

    # Loop over each hourly interval.
    for i in 1:n_hours
        current_hour_dt = add_ignoring_leap_days(begin_dt, Hour(i - 1))
        current_date = Date(current_hour_dt)
        # Reset G_start at the beginning of a new day.
        if i == 1 || (i > 1 && Date(add_ignoring_leap_days(begin_dt, Hour(i - 2))) != current_date)
            G_start = 0.0
        end

        # Compute the boundaries of this hourly interval (in fractional hours).
        start_of_hour = fractional_hour(current_hour_dt)
        end_of_hour = fractional_hour(current_hour_dt + Hour(1)) == 0.0 ? 24.0 :
                      fractional_hour(current_hour_dt + Hour(1))
        interval_start = (i == 1) ? fractional_hour(begin_dt) : start_of_hour
        interval_end = (i == n_hours) ? fractional_hour(end_dt) : end_of_hour

        # Get sunrise and sunset for current_date.
        tz = +1.0  # TODO
        sunrise, sunset = get_sunset_sunrise(current_date, sim_params["latitude"], sim_params["longitude"], tz)

        # Retrieve current and next hourly integrated irradiation.
        Hn = hourly_H[i]
        if i < n_hours
            next_date = Date(add_ignoring_leap_days(begin_dt, Hour(i)))
            Hnp1 = (next_date == current_date) ? hourly_H[i + 1] : 0.0
        else
            Hnp1 = 0.0
        end

        # Segment the current hourly interval.
        t_sub, irr_sub, G_end = segment_interval(Hn, Hnp1,
                                                 interval_start, interval_end,
                                                 dt_h, G_start, sunrise, sunset)
        # Convert each fractional hour (relative to midnight) to DateTime for current_date.
        t_sub_dt = [add_ignoring_leap_days(DateTime(current_date), Millisecond(round(Int, t * 3600000)))
                    for t in t_sub]

        append!(times_total, t_sub_dt)
        append!(irr_total, irr_sub)
        G_start = G_end
    end

    return irr_total, times_total
end

# -----------------------------------------------------------------
# get_sunset_sunrise computes sunrise and sunset (in fractional hours)
# for a given Date and location.
#
# (This is one common NOAA‐based implementation.)
# -----------------------------------------------------------------
function get_sunset_sunrise(date::Date, lat::Float64, lon::Float64, tz::Float64)
    # Helper functions for degree–radian conversion (used in sunrise/sunset)
    deg2rad(deg) = deg * (π / 180)
    rad2deg(rad) = rad * (180 / π)

    n = dayofyear(date)
    lat_rad = deg2rad(lat)
    g = 2π / 365 * (n - 1)
    decl = 0.006918 - 0.399912 * cos(g) + 0.070257 * sin(g) - 0.006758 * cos(2 * g) +
           0.000907 * sin(2 * g) - 0.002697 * cos(3 * g) + 0.00148 * sin(3 * g)
    E = 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) -
                  0.014615 * cos(2 * g) - 0.040849 * sin(2 * g))
    zenith = deg2rad(90.833)  # accounts for refraction
    cos_omega = (cos(zenith) - sin(lat_rad) * sin(decl)) / (cos(lat_rad) * cos(decl))
    if cos_omega > 1
        return (NaN, NaN)   # sun never rises
    elseif cos_omega < -1
        return (NaN, NaN)   # sun never sets
    end
    omega = acos(cos_omega)
    delta_minutes = rad2deg(omega) * 4
    solar_noon = 720 - 4 * lon - E + tz * 60   # in minutes
    sunrise_min = solar_noon - delta_minutes
    sunset_min = solar_noon + delta_minutes
    return sunrise_min / 60, sunset_min / 60  # fractional hours
end

"""
    profile_linear_interpolation(values, timestamps, original_time_step, new_time_step, sim_params)

This function takes the values and corresponding timestamps and converts them from the original_time_step 
to the new_time_step using linear interpolation. It can also handle time shifts if the datestamp of
the profile does not align with the sim_params["start"] datetime.
It is used to break down profiles from a coarser timestep to a finer one (segmentation) or for a time shift only.

"""
function profile_linear_interpolation(values::Vector{Float64},
                                      timestamps::Vector{DateTime},
                                      original_time_step::Period,
                                      new_time_step::Period,
                                      sim_params::Dict{String,Any})
    ref_time = sim_params["start_date"]
    numeric_timestamps = [Dates.value(Millisecond(sub_ignoring_leap_days(dt, ref_time))) for dt in timestamps]
    interp = interpolate((numeric_timestamps,), values, Gridded(Linear()))
    new_timestamps = remove_leap_days(collect(range(sim_params["start_date"];
                                                    stop=sim_params["end_date"],
                                                    step=new_time_step)))
    new_numeric_timestamps = [Dates.value(Millisecond(sub_ignoring_leap_days(dt, ref_time))) for dt in new_timestamps]
    converted_profile = [interp(t) for t in new_numeric_timestamps]

    # add missing entries by coping the last entry
    for _ in 1:ceil(Int, (original_time_step - new_time_step) / new_time_step)
        append!(converted_profile, converted_profile[end])
    end

    return converted_profile, new_timestamps
end

"""
    profile_spread_to_segments(values, timestamps, original_time_step, new_time_step, sim_params)

This function takes the values and corresponding timestamps and converts them from the original_time_step 
to the new_time_step using a stepwise interpolation. 
It is used to break down profiles from a coarser timestep to a finer one (segmentation).

"""
function profile_spread_to_segments(values::Vector{Float64},
                                    timestamps::Vector{DateTime},
                                    original_time_step::Period,
                                    new_time_step::Period,
                                    sim_params::Dict{String,Any})
    segmentation_factor = Int(original_time_step / new_time_step)   # is always > 1 and of type {Int} as only full 
    #                                                                 dividers are allowed

    if sim_params["start_date"] in timestamps
        start_idx = findfirst(x -> x == sim_params["start_date"], timestamps)
    else
        start_idx = findfirst(x -> x >= sim_params["start_date"], timestamps) - 1
    end
    new_length = ceil(Int, sub_ignoring_leap_days(sim_params["end_date"], timestamps[start_idx]) / new_time_step) + 1

    converted_profile = zeros(new_length)

    for new_idx in 1:new_length
        old_idx = Int(start_idx + floor((new_idx - 1) / segmentation_factor))
        converted_profile[new_idx] = values[old_idx]
    end

    end_time = add_ignoring_leap_days(timestamps[start_idx], new_time_step * (new_length - 1))
    new_timestamps = remove_leap_days(collect(range(timestamps[start_idx]; stop=end_time, step=new_time_step)))

    return converted_profile, new_timestamps
end

"""
    profile_calculate_means(values, timestamps, original_time_step, new_time_step, sim_params)

This function takes the values and corresponding timestamps and converts them from the original_time_step 
to the new_time_step using the mean of several original_time_steps to get one new_time_step. 
It is used to aggregate profiles from a finer timestep to a coarser one (aggregation).

"""
function profile_calculate_means(values::Vector{Float64},
                                 timestamps::Vector{DateTime},
                                 original_time_step::Period,
                                 new_time_step::Period,
                                 sim_params::Dict{String,Any})
    aggregation_factor = Int(new_time_step / original_time_step)   # is always > 1 and of type {Int} as only full 
    #                                                                dividers are allowed
    old_length = length(timestamps)
    new_length = ceil(Int, sub_ignoring_leap_days(sim_params["end_date"], sim_params["start_date"]) / new_time_step) + 1
    start_idx = findfirst(x -> x == sim_params["start_date"], timestamps) # timeshift is already done here

    converted_profile = zeros(new_length)

    for n in 1:new_length
        old_index_start = Int(start_idx + (n - 1) * aggregation_factor)
        old_index_end = Int(min(start_idx + n * aggregation_factor - 1, old_length))
        old_number_of_steps = Int(old_index_end - old_index_start + 1)
        converted_profile[n] = sum(values[old_index_start:old_index_end]) / old_number_of_steps
    end

    end_time = add_ignoring_leap_days(sim_params["start_date"], new_time_step * (new_length - 1))
    new_timestamps = remove_leap_days(collect(range(sim_params["start_date"]; stop=end_time, step=new_time_step)))

    return converted_profile, new_timestamps
end

end

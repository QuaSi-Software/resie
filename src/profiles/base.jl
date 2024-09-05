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

The data_type is essential to make sure that the profile data is converted correcty.
It has to be set to "intensive" for values like temperatures, power, wind speed or
for profiles containting states. Only for energies, data_type has to be set to "extensive".
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
                     shift::Dates.Second=Second(0))                   # optional: timeshift for data. A positive shift 
        #                                                                         adds to the timestamp = value is earlier
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
            use_linear_segmentation = false

            file_handle = nothing

            try
                file_handle = open(abspath(file_path), "r")
                for (line_idx, line) in enumerate(readlines(file_handle))
                    line = strip(line)

                    if isempty(line) || length(line) < 2 # handle empty lines
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
                        elseif strip(splitted[1]) == "use_linear_segmentation"
                            use_linear_segmentation = parse(Bool, String(strip(splitted[2])))
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
            use_linear_segmentation = false
        end

        # remove leap days and daiylight savings from profile
        profile_timestamps_date, profile_values = remove_leap_days(profile_timestamps_date, profile_values, file_path)
        profile_timestamps_date = remove_day_light_saving(profile_timestamps_date, time_zone, file_path)
        if shift > Second(0)
            profile_timestamps_date = add_ignoring_leap_days.(profile_timestamps_date, shift)
            @info "The timestamp of the profile at $(file_path) was shifted by $shift."
        elseif shift < Second(0)
            profile_timestamps_date = sub_ignoring_leap_days.(profile_timestamps_date, -shift)
            @info "The timestamp of the profile at $(file_path) was shifted by $shift."
        end

        # add first or last entry if only one time step is missing. copy first or last value.
        if profile_timestamps_date[1] > sim_params["start_date"] &&
           profile_timestamps_date[1] - Second(profile_time_step) <= sim_params["start_date"]
            profile_timestamps_date = vcat(profile_timestamps_date[1] - Second(profile_time_step),
                                           profile_timestamps_date)
            profile_values = vcat(profile_values[1], profile_values)
            @info "The profile at $(file_path) has been extended by one timestep at the begin by doubling the fist value."
        end
        time_diff = Second(max(0, sim_params["time_step_seconds"] - profile_time_step))
        if profile_timestamps_date[end] < sim_params["end_date"] + time_diff &&
           profile_timestamps_date[end] + time_diff >= sim_params["end_date"]
            while profile_timestamps_date[end] < sim_params["end_date"] + time_diff
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
                                                            Second(Int(profile_time_step)),
                                                            Second(sim_params["time_step_seconds"]),
                                                            data_type,
                                                            file_path,
                                                            sim_params,
                                                            use_linear_segmentation)

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
    use_linear_segmentation::Bool   flag if a linear (true) or stepwise (false) interpolation should 
                                    be used for segmentation.

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
                         use_linear_segmentation::Bool)

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
        info_message *= "was converted from the profile timestep $(original_time_step) to the simulation " *
                        "timestep of $(new_time_step)."
    elseif original_time_step > new_time_step                   # segmentation
        if profile_type == :extensive
            values .*= new_time_step / original_time_step
            info_message *= "(extensive profile) "
        else
            info_message *= "(intensive profile) "
        end
        if use_linear_segmentation
            values, timestamps = profile_linear_interpolation(values,
                                                              timestamps,
                                                              original_time_step,
                                                              new_time_step,
                                                              sim_params)
            if !time_is_aligned
                info_message *= "was shifted to fit the simulation start time and "
            end
            info_message *= "was converted from the profile timestep $(original_time_step) to the simulation " *
                            "timestep of $(new_time_step) using linear interpolation."
        else # stepwise conversion
            values, timestamps = profile_spread_to_segments(values,
                                                            timestamps,
                                                            original_time_step,
                                                            new_time_step,
                                                            sim_params)
            if !time_is_aligned  # do time shift
                values, timestamps = profile_linear_interpolation(values,
                                                                  timestamps,
                                                                  new_time_step,
                                                                  new_time_step,
                                                                  sim_params)
                info_message *= "was shifted using linear interpolation to fit the simulation start time and "
            end
            info_message *= "was converted from the profile timestep $(original_time_step) to the simulation " *
                            "timestep of $(new_time_step) using stepwise segmentation."
        end
    end

    @info info_message

    return values, timestamps
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
    numeric_timestamps = [Dates.value(Second(sub_ignoring_leap_days(dt, ref_time))) for dt in timestamps]
    interp = interpolate((numeric_timestamps,), values, Gridded(Linear()))
    new_timestamps = remove_leap_days(collect(range(sim_params["start_date"];
                                                    stop=sim_params["end_date"],
                                                    step=new_time_step)))
    new_numeric_timestamps = [Dates.value(Second(sub_ignoring_leap_days(dt, ref_time))) for dt in new_timestamps]
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
        old_idx = Int(ceil((start_idx + new_idx - 1) / segmentation_factor))
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

    for n in start_idx:new_length
        old_index_start = Int((n - 1) * aggregation_factor + 1)
        old_index_end = Int(min(n * aggregation_factor, old_length))
        old_number_of_steps = Int(old_index_end - old_index_start + 1)
        converted_profile[n - start_idx + 1] = sum(values[old_index_start:old_index_end]) / old_number_of_steps
    end

    end_time = add_ignoring_leap_days(sim_params["start_date"], new_time_step * (new_length - 1))
    new_timestamps = remove_leap_days(collect(range(sim_params["start_date"]; stop=end_time, step=new_time_step)))

    return converted_profile, new_timestamps
end

end

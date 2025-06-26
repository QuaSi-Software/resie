module SolarIrradiance

using Dates
using Roots

export sun_position, irr_in_plane, get_sunrise_sunset

"""
    sun_position(dt::DateTime, time_step_seconds::Number, longitude::Number, 
                 latitude::Number, time_zone::Number=0, pressure::Number=1.0, 
                 temperature::Number=10.9)

Calculate solar position in degrees.
Time is shifted by a half time_step forward, because it all profiles in Resie are given as 
sum or mean over the next full time_step.
Based on Roberto Grena (2012), Five new algorithms for the computation of sun position
from 2010 to 2110, Solar Energy, 86(5):1323–1337, doi:10.1016/j.solener.2012.01.024.
dt must be provided in UTC, if no time_zone is provided.
time_zone is the time difference to UTC as INT in hours.
longitude and latitude should be provided in WGS84.
pressure in atm and temperature in °C are needed to apply the refraction correction which is 
relevant when the sun is low.

# Arguments
- `dt::DateTime`: Current date of the simulation as DateTime
- `time_step_seconds::UInt64`: Time_step of the simulation in seconds
- `longitude::Number`: Longitude of the location
- `latitude::Number`: Latitude of the location
- `time_zone::Number=0`: Time zone of dt
- `pressure::Number=1.0`: Average air pressure at the location
- `temperature::Number=10.9`: Average air temperature at the location

# Returns
- 'zenith::Float64': Zenith of the sun in degrees.
- 'azimuth::Float64': Azimuth of the sun in degrees.
"""
function sun_position(dt::DateTime, time_step_seconds::Number, longitude::Number, 
                      latitude::Number, time_zone::Number=0, pressure::Number=1.0, 
                      temperature::Number=10.9)

    longitude = deg2rad(longitude)
    latitude = deg2rad(latitude)

    dt_utc = dt - Hour(time_zone) + Second(round(time_step_seconds / 2))

    dt2060 = DateTime(2060, 1, 1)
    t2060 = Int64(Dates.value(dt_utc.instant - dt2060.instant)) / 86400000.0
    tt = t2060 + 1.1574e-5 * (96.4 + 0.00158 * t2060)

    rightAscension, declination, hourAngle = alg5(t2060, tt, longitude)

    sp = sin(latitude)
    cp = sqrt((1 - sp * sp))
    sd = sin(declination)
    cd = sqrt(1 - sd * sd)
    sH = sin(hourAngle)
    cH = cos(hourAngle)
    se0 = sp * sd + cp * cd * cH
    ep = asin(se0) - 6371.0 / 149597871 * sqrt(1.0 - se0 * se0)

    azimuth = atan(sH, cH * sp - sd * cp / cd)

    epr = max((0.08422 * pressure) / ((273.0 + temperature) * tan(ep + 0.003138 / (ep + 0.08919))), 0)
    # if epr > 0.0
    #     zenith = pi/2 - ep - epr
    # else
    #     zenith = pi/2 - ep
    # end
    zenith = pi/2 - ep - epr
    return rad2deg(zenith), rad2deg(azimuth)
end

function alg5(t2060::Number, tt::Number, longitude::Number)
    wtt = tt * 0.0172019715

    s1 = sin(wtt)
    c1 = cos(wtt)
    s2 = 2.0 * s1 * c1
    c2 = (c1 + s1) * (c1 - s1)
    s3 = s2 * c1 + c2 * s1
    c3 = c2 * c1 - s2 * s1

    L = 1.7527901 + 1.7202792159e-2*tt + 3.33024e-2*s1 - 2.0582e-3*c1 +
        3.512e-4*s2 - 4.07e-5*c2 + 5.2e-6*s3 - 9e-7*c3 -
        8.23e-5*s1*sin(2.92e-5*tt) + 1.27e-5*sin(1.49e-3*tt - 2.337) +
        1.21e-5*sin(4.31e-3*tt + 3.065) + 2.33e-5*sin(1.076e-2*tt - 1.533) +
        3.49e-5*sin(1.575e-2*tt - 2.358) + 2.67e-5*sin(2.152e-2*tt + 0.074) +
        1.28e-5*sin(3.152e-2*tt + 1.547) + 3.14e-5*sin(2.1277e-1*tt - 0.488)

    nu = tt * 9.282e-4 - 0.8
    dlam = 8.34e-5 * sin(nu)
    lambda = L + pi + dlam

    epsi = 4.089567e-1 - 6.19e-9 * t2060 + 4.46e-5 * cos(nu)

    sl = sin(lambda)
    cl = cos(lambda)
    se = sin(epsi)
    ce = sqrt(1 - se * se)

    rightAscension = atan(sl * ce, cl)
    if rightAscension < 0.0
        rightAscension += 2 * pi
    end

    declination = asin(sl * se)

    hourAngle = 1.7528311 + 6.300388099 * t2060 + longitude - rightAscension + 0.92 * dlam
    hourAngle = mod(hourAngle + pi, 2 * pi) - pi

    return rightAscension, declination, hourAngle
end

"""
    irr_in_plane(sim_params::Dict{String,Any}, tilt_angle::Number, azimuth_angle::Number, 
                 beam_solar_hor_irradiance::Number, diffuse_solar_hor_irradiance::Number, 
                 dni=nothing, pressure::Number=1.0, temperature::Number=20.0, 
                 ground_reflectance::Number=0.2, zenith_threshold_for_zero_dni::Number=89.5)

Calculate beam irradiance in collector plane and direct_normal_irradiance.
All angles are in degrees and irradiances in W/m².

# Arguments
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `tilt_angle::Number`: Tilt angle of the collector.
- `azimuth_angle::Number` azimuth angle of the collector.
- `beam_solar_hor_irradiance::Number`: Horizontal beam irradiance. 
                                          Either this or dni has to be given
- `diffuse_solar_hor_irradiance::Number` Horizontal diffuse irradiance.
- `dni::Floathing=nothing`: Direct normal irradiance. 
                            Either this or beam_solar_hor_irradiance has to be given
- `pressure::Number=1.0`: Average air pressure at the location.
- `temperature::Number=10.9`: Average air temperature at the location.
- `ground_reflectance::Number=0.2`: Ground reflectance of the surrounding ground.
- `zenith_threshold_for_zero_dni::Number=89.5': Threshold after which the direct normal 
                                                irradiance is set to 0. Is used to prevent 
                                                big calculation mistakes close to 90° zenith

# Returns
- `beam_solar_irradiance_in_plane::Float64`: Beam irradiance in the collector plane.
- `diffuse_solar_irradiance_in_plane::Float64`: Diffuse irradiance in the collector plane.
- `direct_normal_irradiance::Float64`: Direct_normal_irradiance.
- `aoi::Float64`: Angle of incidence on the collector.
- `aoi_l::Float64`: Lateral angle of incidence on the collector.
- `aoi_t::Float64`: Transversal angle of incidence on the collector.
- `solar_zenith::Float64`: Zenith of the sun in degrees.
- 'solar_azimuth::Float64': Azimuth of the sun in degrees.
"""
function irr_in_plane(sim_params::Dict{String,Any}, tilt_angle::Number, azimuth_angle::Number,
                      beam_solar_hor_irradiance::Number, diffuse_solar_hor_irradiance::Number, 
                      dni::Union{Float32,Nothing}=nothing, pressure::Number=1.0, temperature::Number=10.9, 
                      ground_reflectance::Number=0.2, zenith_threshold_for_zero_dni::Number=89.5)

    # get sunrise and senset to adjust time_step length from sunrise to end of time_step and 
    # start of time_step to sunset
    if "weather_data" in keys(sim_params)
        sr = sim_params["weather_data"].sunrise.data[sim_params["current_date"]]
        ss = sim_params["weather_data"].sunset.data[sim_params["current_date"]]
    else
        sr, ss = get_sunrise_sunset(sim_params["current_date"], sim_params["latitude"],
                                    sim_params["longitude"], sim_params["timezone"],
                                    pressure, temperature)
    end
    sr = convert_decimal_time_to_datetime(sim_params["current_date"], sr)
    ss = convert_decimal_time_to_datetime(sim_params["current_date"], ss)
    dt_start = sim_params["current_date"]
    dt_end = dt_start + Second(sim_params["time_step_seconds"])
    if sr >= dt_start && sr <= dt_end
        new_time_step = round(dt_end - sr, Second).value
    elseif ss >= dt_start && ss <= dt_end
        new_time_step = round(ss - dt_start, Second).value
    else
        new_time_step = sim_params["time_step_seconds"]
    end

    solar_zenith, solar_azimuth = sun_position(sim_params["current_date"], new_time_step,
                                               sim_params["longitude"], sim_params["latitude"],
                                               sim_params["timezone"], pressure, temperature)

    # calculate angle of incidence in degrees on the collector plane
    aoi = acosd(cosd(tilt_angle) * cosd(solar_zenith) +
                sind(tilt_angle) * sind(solar_zenith) *
                cosd(solar_azimuth - azimuth_angle))

    # calculate longitudinal and transversal projections of the incidence angle
    aoi_t = atand(sind(solar_zenith) * sind(azimuth_angle - solar_azimuth) / cosd(aoi))
    aoi_t_iso = atand(sind(aoi) * sind(azimuth_angle - solar_azimuth) / cosd(aoi))
    aoi_l = abs(tilt_angle - atand(tand(solar_zenith) * cosd(azimuth_angle - solar_azimuth)))
    aoi_l_iso = atand(sind(aoi) * cosd(azimuth_angle - solar_azimuth) / cosd(aoi))

    ext_terr_normal_irr = extraterrestrial_normal_irradiance(dayofyear(sim_params["current_date"]))

    # calculate direct normal irradiance from beam horrizontal irradiance
    if isnothing(dni)
        if solar_zenith >= zenith_threshold_for_zero_dni
            direct_normal_irradiance = 0
        else
            direct_normal_irradiance = max((beam_solar_hor_irradiance) / cosd(solar_zenith), 0)
            direct_normal_irradiance = min(direct_normal_irradiance, ext_terr_normal_irr)
        end
    else
        direct_normal_irradiance = dni
    end

    # calculate the beam irradiance on the collector plane
    beam_solar_irradiance_in_plane = max(direct_normal_irradiance * cosd(aoi), 0)

    # calculate the diffuse irradiance on the collector plane from the sky with the Hay 
    # model see: J. E. Hay, Calculation of monthly mean solar radiation for horizontal and  
    # inclined surfaces, Solar Energy, Jg. 23, Nr. 4, S. 301–307, 1979. 
    # doi: 10.1016/0038-092X(79)90123-3. 
    anisotropy_index = direct_normal_irradiance / ext_terr_normal_irr

    diffuse_sky_irradiance_in_plane = max(0,
                                          diffuse_solar_hor_irradiance *
                                          ((1 - anisotropy_index) * (1 + cosd(tilt_angle)) / 2 +
                                           anisotropy_index * cosd(aoi) / sind(solar_zenith)))

    # calculate the diffuse irradiance on the collector plane from the ground
    diffuse_ground_irradiance_in_plane = ground_reflectance *
                                         (beam_solar_hor_irradiance + diffuse_solar_hor_irradiance) *
                                         (1 - cosd(tilt_angle)) / 2

    diffuse_solar_irradiance_in_plane = diffuse_sky_irradiance_in_plane + diffuse_ground_irradiance_in_plane

    return beam_solar_irradiance_in_plane, diffuse_solar_irradiance_in_plane,
           direct_normal_irradiance, aoi, aoi_l, aoi_t, solar_zenith, solar_azimuth
end

"""
    extraterrestrial_normal_irradiance(day_of_year::Int)

Calculate the extraterrestrial normal irradiance in W/m² at the given day of the year. Based 
on J. W. Spencer, "Fourier series representation of the sun," Search, vol. 2, p. 172, 1971.

# Arguments
- `day_of_year::Int`: The day of the year starting with 1 for 01. January.

# Returns
- `etrn::Float64`: Extraterrestrial normal irradiance.
"""
function extraterrestrial_normal_irradiance(day_of_year::Int)
    day_angle = 2 * pi * (day_of_year - 1) / 365
    solar_constant = 1361 # new solar constant since 2015
    etrn = solar_constant * (1.000110 + 0.034221 * cos(day_angle) + 0.001280 * sin(day_angle) +
                             0.000719 * cos(2 * day_angle) + 0.000077 * sin(2 * day_angle))
    return etrn
end

"""
    get_sunrise_sunset(dt::DateTime, latitude::Number, longitude::Number, time_zone::Number, 
                       pressure::Number=1.0, temperature::Number=10.9)

Computes sunrise and sunset (in fractional hours between 0 and 24) for a given Date and 
location.

# Arguments
- `dt::DateTime`: Current date of the simulation as DateTime
- `latitude::Number`: Latitude of the location
- `longitude::Number`: Longitude of the location
- `time_zone::Number=0`: Time zone of dt
- `pressure::Number=1.0`: Average air pressure at the location
- `temperature::Number=10.9`: Average air temperature at the location

# Returns
- 'sunrise_decimal::Float64': Time of sunrise as fractional hour.
- 'sunset_decimal::Float64': Time of sunset as fractional hour.
"""
function get_sunrise_sunset(dt::DateTime, latitude::Number, longitude::Number,
                            time_zone::Number, pressure::Number=1.0, temperature::Number=10.9)

    sunrise_decimal = try
        find_zero(time_decimal -> sun_position_itt(dt, time_decimal, longitude, latitude, 
                                                   time_zone, 90, pressure, temperature),
                  (0, 12),
                  Bisection())
    catch
        0.0
    end
    sunset_decimal = try
        find_zero(time_decimal -> sun_position_itt(dt, time_decimal, longitude, latitude, 
                                                   time_zone, 90, pressure, temperature),
                  (12, 24),
                  Bisection())
    catch
        24.0
    end

    return sunrise_decimal, sunset_decimal
end

"""
    sun_position_itt(dt::DateTime, time_decimal::Number, longitude::Number, 
                     latitude::Number, time_zone::Number=0, target_zenith::Number=90, 
                     pressure::Number=1.0, temperature::Number=10.9)

Calculate solar position in degrees in an itteratrive process. 
Can be used to solve for sunrise and sunset.

# Arguments
- `dt::DateTime`: Date of the day as DateTime
- `time_decimal::Float64`: Time of the day as fractional hour between 0 and 24
- `longitude::Number`: Longitude of the location
- `latitude::Number`: Latitude of the location
- `time_zone::Number=0`: Time zone of dt
- `target_zenith::Number=90`: A target zenith which is substracted from the sun_position. 
                              Set to 90 to find sunrise and sunset.
- `pressure::Number=1.0`: Average air pressure at the location
- `temperature::Number=10.9`: Average air temperature at the location

# Returns
- 'Float64': Difference between zenith of the sun and target_zenith in degrees.
"""
function sun_position_itt(dt::DateTime, time_decimal::Number, longitude::Number, 
                          latitude::Number, time_zone::Number=0, target_zenith::Number=90, 
                          pressure::Number=1.0, temperature::Number=10.9)

    dt = convert_decimal_time_to_datetime(dt, time_decimal)
    solar_zenith, _ = sun_position(dt, 0, longitude, latitude, time_zone, pressure, temperature)
    return solar_zenith - target_zenith
end

"""
    convert_decimal_time_to_datetime(dt::DateTime, time_decimal::Number)

Quick helper function to convert time from fractional hour to datetime format.

# Arguments
- `dt::DateTime`: Date of the day as DateTime
- `time_decimal::Float64`: Time of the day as fractional hour between 0 and 24

# Returns
- 'DateTime': Combined DateTime with Date of dt and Time of time_decimal
"""
function convert_decimal_time_to_datetime(dt::DateTime, time_decimal::Number)
    h = floor(time_decimal)
    m = floor((time_decimal - h) * 60)
    s = floor(((time_decimal - h) * 60 - m) * 60)
    ns = floor((((time_decimal - h) * 60 - m) * 60 - s) * 1000000000)
    return DateTime(year(dt), month(dt), day(dt), h, m, s) + Nanosecond(ns)
end

end
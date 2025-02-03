using Dates

function sun_position(dt::DateTime, time_step_seconds, longitude::Number, latitude::Number, 
                             pressure::Number=1.0, temperature::Number=20.0)
"""
Calculate solar position in degrees.
Time is shifted by a half time_step forward, because it all profiles in Resie are given as 
sum or mean over the next full time_step.
Based on Roberto Grena (2012), Five new algorithms for the computation of sun position
from 2010 to 2110, Solar Energy, 86(5):1323–1337, doi:10.1016/j.solener.2012.01.024.
dt must be provided in UTC.
longitude and latitude should be provided in WGS84.
pressure and temperature are needed to apply the refraction correction which is relevant
when the sun is low.
"""

    longitude = deg2rad(longitude)
    latitude = deg2rad(latitude)

    dt = dt + Dates.Second(time_step_seconds / 2)

    dt2060 = DateTime(2060,1,1)
    t2060 = Int64(Dates.value(dt.instant - dt2060.instant)) / 86400000.0
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

    zenith  = pi/2 - ep
    if ep > 0.0
        zenith -= (0.08422*pressure)/((273.0+temperature)*tan(ep + 0.003138/(ep + 0.08919)))
    end

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
    ce = sqrt(1-se*se)

    rightAscension = atan(sl * ce, cl)
    if rightAscension < 0.0
        rightAscension += 2*pi
    end

    declination = asin(sl * se)

    hourAngle = 1.7528311 + 6.300388099 * t2060 + longitude - rightAscension + 0.92 * dlam;
    hourAngle = mod(hourAngle + pi, 2*pi) - pi;

    return rightAscension, declination, hourAngle
end


function beam_irr_in_plane(tilt_angle::Number, azimuth_angle::Number, solar_zenith::Number, solar_azimuth::Number, 
                           global_solar_hor_irradiance::Number, diffuse_solar_hor_irradiance::Number, dni=nothing, 
                           zenith_threshold_for_zero_dni::Number=89.5)
    """
    Calculate beam irradiance in collector plane and direct_normal_irradiance
    All angles are in degrees and irradiances in W/m²
    """
    # calculate angle of incidence in degrees on the collector plane
    aoi = acosd(
        cosd(tilt_angle) * cosd(solar_zenith) + sind(tilt_angle) * sind(solar_zenith) * 
        cosd(solar_azimuth - azimuth_angle)
        )

    # calculate longitudinal and transversal projections of the incidence angle
    aoi_t = atand(sind(solar_zenith) * sind(azimuth_angle - solar_azimuth) / cosd(aoi))
    aoi_t_iso = atand(sind(aoi) * sind(azimuth_angle - solar_azimuth) / cosd(aoi))
    aoi_l = abs(tilt_angle - atand(tand(solar_zenith) * cosd(azimuth_angle - solar_azimuth)))
    aoi_l_iso = atand(sind(aoi) * cosd(azimuth_angle - solar_azimuth) / cosd(aoi))

    # calculate direct normal irradiance from global horrizontal irradiance (GHI) and 
    # diffuse horrizontal irradiance (DHI)
    if isnothing(dni)
        if solar_zenith >= zenith_threshold_for_zero_dni
            direct_normal_irradiance = 0
        else
            direct_normal_irradiance = max(
                (global_solar_hor_irradiance - diffuse_solar_hor_irradiance) / cosd(solar_zenith), 0
                )
        end
    else
        direct_normal_irradiance = dni
    end

    # calculate the beam irradiance on the collector plane
    beam_solar_irradiance_in_plane = max(direct_normal_irradiance * cosd(aoi), 0)

    return beam_solar_irradiance_in_plane, direct_normal_irradiance, aoi, aoi_l, aoi_t
end

export sun_position, beam_irr_in_plane

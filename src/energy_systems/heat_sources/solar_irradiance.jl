using Dates

function sun_position(dt::DateTime, longitude::Float64, latitude::Float64, 
                             pressure::Float64=1.0, temperature::Float64=20.0)
"""
calculate solar position in radians
Based on Roberto Grena (2012), Five new algorithms for the computation of sun position
from 2010 to 2110, Solar Energy, 86(5):1323–1337, doi:10.1016/j.solener.2012.01.024.
"""

    longitude = deg2rad(longitude)
    latitude = deg2rad(latitude)

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

function alg5(t2060::Float64, tt::Float64, longitude::Float64)
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


function beam_irr_in_plane(tilt_angle, azimuth_angle, solar_zenith, solar_azimuth, 
                           global_solar_hor_irradiance, diffuse_solar_hor_irradiance)
    """
    All angles are in degrees and irradiances in W/m²
    """
    # calculate angle of incidence in degrees on the collector plane
    aoi = acosd(
        cosd(tilt_angle) * cosd(solar_zenith) + sind(tilt_angle) * sind(solar_zenith) * 
        cosd(solar_azimuth - azimuth_angle)
        )

    # calculate direct normal irradiance from global horrizontal irradiance (GHI) and 
    # diffuse horrizontal irradiance (DHI)
    direct_nomal_irradiance = max(
        (global_solar_hor_irradiance - diffuse_solar_hor_irradiance) / cosd(solar_zenith), 0
        )

    # calculate the beam irradiance on the collector plane
    beam_solar_irradiance_in_plane = max(direct_nomal_irradiance * cosd(aoi), 0)

    return beam_solar_irradiance_in_plane, direct_nomal_irradiance
end

export sun_position, beam_irr_in_plane

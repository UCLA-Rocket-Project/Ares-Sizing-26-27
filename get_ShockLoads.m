% get_recLoads

% Script from anthony, just made it a function with our desired outputs

% Inputs: dry mass
% Outputs: Shock Load from drogue and main deployment

function [f_drogue, f_main] = get_ShockLoads(dry_mass)
    % Parachute
    cd_drogue = 2.2; % cd from manufacturer
    cd_main = 2.2; % cd from manufacturer
    
    % Flight and Rocket Info
    mass_rocket = dry_mass/2.205; % dry mass of rocket (kg)
    weight_rocket = (mass_rocket)* 9.81; % weight of rocket (N)
    alt_drogue = 15240; % altitude in meters for drogue deployment
    alt_main = 304.8; % altitude in meters for main deployment
    elevation = 624.2; % elevation of launch site (m)
    velocity_at_drogue = 300/3.281; % Velocity at drogue deployment (m/s)
    area_drogue = 1.16745; % canopy area of drogue (m^2)
    area_main = 13.8583; % canopy area of main (m^2)

    % Force Drogue
    rho1 = density_air_func(alt_drogue + elevation);
    RmD = mass_ratio(rho1, cd_drogue, area_drogue, mass_rocket);
    [Ck1_avg, Ck1_max] = Ck_approx(RmD);
    Ck1 = [Ck1_avg, Ck1_max];
    Forces1 = F_max(rho1, velocity_at_drogue, cd_drogue, area_drogue, Ck1);% lbf
    
    if 1.75 * Forces1(1) >= 1.5 * Forces1(2)
      f_drogue = Forces1(1) * 2; % lbf
    else
       f_drogue = Forces1(2) * 2; % lbf
    end
    
    % Force Main
    rho2 = density_air_func(alt_main + elevation); % kg/m^3
    velocity_at_main = sqrt((2 * weight_rocket) / (cd_drogue * rho2 * area_drogue)); 
    RmM = mass_ratio(rho2, cd_main, area_main, mass_rocket);
    [Ck2_avg, Ck2_max] = Ck_approx(RmM);
    Ck2 = [Ck2_avg, Ck2_max];
    Forces2 = F_max(rho2, velocity_at_main, cd_main, area_main, Ck2); % lbf
    
    if 1.75 * Forces2(1) >= 1.5 * Forces2(2)
       f_main = Forces2(1) * 2; % lbf
    else
       f_main = Forces2(2) * 2; % lbf
    end
end

%% Functions

function[rho] = density_air_func(altitude)
    if altitude <= 11000
        Temp = 15.04 - 0.00649 * altitude; % Temperature in Celsius
        Pres = 101.29 * ((Temp + 273.1) / 288.08)^(5.256); % Pressure in kPa
    elseif altitude > 11000 && altitude < 25000
        Temp = -56.46; % Temperature in Celsius
        Pres = 22.65 * exp(1.73 - 0.000157 * altitude); % Pressure in kPa
    else
        Temp = -131.21 + 0.00299 * altitude; % Temperature in Celsius
        Pres = 2.488 * ((Temp + 273.1) / 216.6)^(-11.388); % Pressure in kPa
    end
    rho = Pres / (0.2869 * (Temp + 273.1)); % Air density in kg/m^3
end

function[Rm] = mass_ratio(rho, Cd, Area, mass)
    Rm = rho * (Cd * Area) ^ (3/2) / mass;
end

function[Ck_avg, Ck_max] = Ck_approx(Rm)
    Ck_max = -0.019*(log(Rm) + 5.6).^2 + 1.45;
    Ck_min = -0.01*(log(Rm) + 7.8).^2 + 0.99;
    Ck_avg = (Ck_max + Ck_min) / 2;
end

function[f_max] = F_max(rho, velocity, Cd, Area, Ck)
    f_max =  0.224809 * 0.5 * rho * ((velocity)^2) * Area * Cd * Ck; % lbf
end
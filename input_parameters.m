function p = input_parameters()

    p.Ctau_eff = 0.98;

    p.Cstar_eff = 0.90;

    p.P_amb = 14.7; % psi

    p.T_fuel = 300; % K

    p.T_ox = 94; % K

    p.P_fuel_tank = 600 * 6894; %pa

    p.P_ox_tank = 600 * 6894; %pa

    % initial pressure inputs for CoolProp density call, not used
    % Tank pressure calculated from run_press used in main
    
    p.ox_density = py.CoolProp.CoolProp.PropsSI('D','T', p.T_ox, 'P', p.P_ox_tank, 'oxygen'); % kg/m3

    p.water_density = py.CoolProp.CoolProp.PropsSI('D','T', p.T_fuel, 'P', p.P_fuel_tank, 'water'); % kg/m^3
    p.ethanol_density = py.CoolProp.CoolProp.PropsSI('D','T',p.T_fuel, 'P', p.P_fuel_tank, 'ethanol'); % kg/m^3

    %density conversions happen inside run_press, CoolProp returns SI 
    p.eth_ratio = 0.75;

    p.fuel_density = mass_fraction(p.eth_ratio, p.ethanol_density, p.water_density); % kg/m3

    p.ullage = 0.05;

    p.T_He = 298; % K initial He temp

    p.P_He_init = 4500; %psi COPV pressure

    % p.ox_feed_CdA = ; % m^2

    % p.fuel_feed_CdA = ; % m^2

    %p.channel_CdA = ; % m^2

    %p.inj_ox_CdA = ; % m^2

    %p.inj_f_CdA = ; % m^2

    % CdA conversions happen inside run_press

    %p.ac_at = ;

    %p.A_e = ; % in^2

end

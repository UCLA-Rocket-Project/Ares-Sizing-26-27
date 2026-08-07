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

    % run_press computes m^3 tank volume using kg/m^3 densities, then converts to L for get_PV_mel
    
    p.eth_ratio = 0.75;

    p.fuel_density = mass_fraction(p.eth_ratio, p.ethanol_density, p.water_density); % kg/m3

    p.ullage = 0.05;

    p.T_He = 298; % K initial He temp

    p.P_He_init = 4500; %psi COPV pressure

    p.GN2_bottle_orifice_area = 2.027E-5; % m^2 

    p.Dome_orifice_area = 2.677414E-5; % m^2 

%     p.ox_feed_CdA = 4.435E-05 ; % m^2

  %  p.fuel_feed_CdA = 6.8516E-05 ; % m^2

   % p.channel_CdA = 0.000047624 ; % m^2

 %   p.inj_ox_CdA = 0.00005837500098; % m^2

  %  p.inj_f_CdA = 0.00008144; % m^2

    % CdA conversions happen inside run_press

    p.ac_at = 4 ;

    p.A_e = 13.59; % in^2 From FRR

end

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

    p.T_He = 298; % K initial He tempterm

    p.P_He_init = 4500; %psi COPV pressure

    p.GN2_bottle_orifice_area = 1.7719053832789E-5; % m^2, from CGA V-1 2013 Page 55 (CGA 580) nipple diameter https://www.scribd.com/document/967862515/CGA-V1-2013

    p.Dome_orifice_area = 2.6804825641854E-5; % m^2, from Aqua 873-D Cv = 0.8, Diameter = 0.23in

    p.T_N2 = 293; %K

    % feed/injector CdAs
    % taken from config.py in Michael's engine balance script

    p.CdA_fuel_feed = 8.860385824464009E-5; % m^2
    p.CdA_ox_feed = 5.660629324484712E-5; % m^2
    p.CdA_fuel_injector = 0.00003687789023; % m^2
    p.CdA_ox_injector = 0.00004918255865; % m^2
 
    % Series-combined CdA (feed + injector)
    
    p.CdA_fuel = 1 / sqrt((1/p.CdA_fuel_feed)^2 + (1/p.CdA_fuel_injector)^2); % m^2
    p.CdA_ox = 1 / sqrt((1/p.CdA_ox_feed)^2 + (1/p.CdA_ox_injector)^2); % m^2

    p.ac_at = 4 ;

    p.A_e = 13.59; % in^2 From FRR

end

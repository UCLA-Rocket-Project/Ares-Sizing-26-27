% run_press

% computes required tank pressure for given chamber pressure input
% computes Helium volume needed
% Uses stiffness + measured flight system dP now, instead of using CdAs and mdot equation
% conservative feed loss overhead added, can orifice as needed 
% assume 20% stiffness (not heritage 30%, high stiffness was due to manufacturing, not planned)

% Inputs:  Prop (OF, Pc, mdot, mdot_fuel, mdot_ox), params
% Outputs: Press.tank_press [psi], Press.V_He [m^3]

function [Prop, Press] = run_press(Prop, params)

   % mdot_fuel_SI = Prop.mdot_fuel / 2.20462; % lbm/s to kg/s
   %  mdot_ox_SI = Prop.mdot_ox / 2.20462; %lbm/s to kg/s

    % Find tank pressure needed for given chamber pressure
 
    % combine fuel and ox CdA (OLD)
    % fuel_CdA = CdA_series([params.fuel_feed_CdA, params.channel_CdA, params.inj_f_CdA]); % m^2
   %  ox_CdA   = CdA_series([params.ox_feed_CdA, params.inj_ox_CdA]); % m^2

    % find dP needed using mdot equation (OLD)
    % fuel_dp = dP_from_CdA_mdot(fuel_CdA, mdot_fuel_SI, params.fuel_density); % psi
   %  ox_dp   = dP_from_CdA_mdot(ox_CdA, mdot_ox_SI, params.ox_density); % psi




   % find hydrostatic pressure head at liftoff (most conservative since this adds pressure)

   % find fuel and ox masses (no ullage) 

    prop_mass_SI = Prop.prop_mass / 2.20462; % lbm to kg
    fuel_mass_SI = prop_mass_SI / (Prop.OF + 1); % kg
    ox_mass_SI   = prop_mass_SI - fuel_mass_SI; % kg

   % tank geometry

  r_o_SI = 4 * 0.0254; % in to m 
  A_tank_SI = pi * r_o_SI ^2; % m^2

  h_fuel_liquid = fuel_mass_SI / params.fuel_density / A_tank_SI; % m, in-tank column only
  h_ox_liquid  = ox_mass_SI / params.ox_density / A_tank_SI;   % m

  % Feed line height (tank outlet to injector face), from heritage MEL

  h_fuel_feed = 28.5; 
  h_ox_feed  = 68.95;
  h_fuel_feed_SI = h_fuel_feed * 0.0254; % in to m
  h_ox_feed_SI  = h_ox_feed * 0.0254;  % in to m

  h_fuel_total = h_fuel_liquid + h_fuel_feed_SI; % m
  h_ox_total  = h_ox_liquid + h_ox_feed_SI;   % m

  g = 9.81; % m/s^2 at liftoff

  dP_head_fuel_pa = params.fuel_density * g * h_fuel_total; % Pa
  dP_head_ox_pa  = params.ox_density * g * h_ox_total;   % Pa

  dP_head_fuel = dP_head_fuel_pa * 0.000145038; % Pa to psi
  dP_head_ox  = dP_head_ox_pa * 0.000145038;  % Pa to psi

 % find dP from heritage and add conservative overhead

   injector_dP = 0.2 * Prop.Pc; % 20% stiffness for low freq instability

   fuel_feed_dP = 20 + 20; % from Pandora flight system hotfire + conservative overhead (can be orificed) 

   ox_feed_dP = 14.4 + 20; % from Pandora flight system hotifre + conservative overhead (can be orificed) 

   % find total dP across feed system

   fuel_dp = injector_dP + fuel_feed_dP;

   ox_dp = injector_dP + ox_feed_dP;

    % add dP to chamber pressure to find tank pressure
    fuel_tank_press = fuel_dp - dP_head_fuel + Prop.Pc; % psi
    ox_tank_press   = ox_dp - dP_head_ox + Prop.Pc; % psi

    Press.tank_press = max(fuel_tank_press, ox_tank_press); % psi

    % Find He volume needed 

    % find tank vollumes by using mass, density, and ullage

    fuel_tank_volume = (fuel_mass_SI / params.fuel_density) / (1 - params.ullage); % m^3
    ox_tank_volume   = (ox_mass_SI / params.ox_density) / (1 - params.ullage); % m^3
    Press.fuel_tank_volume = fuel_tank_volume;
    Press.ox_tank_volume = ox_tank_volume;

    tank_press_pa = Press.tank_press * 6894.76; % psi to Pa (CoolProp wants Pa)

    % He density for fuel side
    He_density_low = py.CoolProp.CoolProp.PropsSI('D','T', params.T_He, 'P', tank_press_pa, 'helium'); % kg/m^3

    % He density for LOx side
 
    He_density_cryo = py.CoolProp.CoolProp.PropsSI('D','T', params.T_ox, 'P', tank_press_pa, 'helium'); % kg/m^3

    % He density at COPV 

    P_He_full_pa = params.P_He_init * 6894.76; % psi to Pa
    He_density_full = py.CoolProp.CoolProp.PropsSI('D','T', params.T_He, 'P', P_He_full_pa, 'helium'); % kg/m^3

    % find mass needed by combining fuel side & ox side Helium density w/ respective volumes
    mass_He_required = (fuel_tank_volume * He_density_low) + (ox_tank_volume * He_density_cryo); % kg

    % divide by COPV side Helium density to find COPV volume

    V_He_SI = mass_He_required / He_density_full; % m^3
    Press.V_He = V_He_SI * 1000; % m^3 to L for get_PV_mel

end


%local helper functions (OLD)

 %function CdA_total = CdA_series(CdA_list)
    % combines CdAs in series (feed + channel + injector)
  %  CdA_total = 1 / sqrt(sum(1 ./ (CdA_list.^2)));
%end

%function dP = dP_from_CdA_mdot(CdA, mdot, rho)
    % uses mdot = CdA*sqrt(2*rho*dP)

    %dP_pa = (mdot^2) / (2 * rho * CdA^2); % Pa
    %dP = dP_pa * 0.000145038; % Pa to psi
%end

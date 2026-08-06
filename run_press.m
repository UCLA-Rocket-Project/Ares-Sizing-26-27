% run_press

% computes required tank pressure for given chamber pressure input
% computes Helium volume needed

% Inputs:  Prop (OF, Pc, mdot, mdot_fuel, mdot_ox), params
% Outputs: Press.tank_press [psi], Press.V_He [m^3]

function [Prop, Press] = run_press(Prop, params)

    mdot_fuel_SI = Prop.mdot_fuel / 2.20462; % lbm/s to kg/s
    mdot_ox_SI = Prop.mdot_ox / 2.20462; %lbm/s to kg/s

    % Find tank pressure needed for given chamber pressure
 
    % combine fuel and ox CdA (OLD)
    % fuel_CdA = CdA_series([params.fuel_feed_CdA, params.channel_CdA, params.inj_f_CdA]); % m^2
   %  ox_CdA   = CdA_series([params.ox_feed_CdA, params.inj_ox_CdA]); % m^2

    % find dP needed using mdot equation (OLD)
    % fuel_dp = dP_from_CdA_mdot(fuel_CdA, mdot_fuel_SI, params.fuel_density); % psi
   %  ox_dp   = dP_from_CdA_mdot(ox_CdA, mdot_ox_SI, params.ox_density); % psi

   injector_dP = 0.2 * Prop.Pc; % 20% stiffness for low freq instability

   fuel_feed_dP = 20 + 20; % from Pandora flight system hotfire + conservative overhead (can be orificed) 

   ox_feed_dP = 14.4 + 20; % from Pandora flight system hotifre + conservative overhead (can be orificed) 

   fuel_dp = injector_dP + fuel_feed_dP;

   ox_dp = injector_dP + ox_feed_dP;

    % add dP to chamber pressure to find tank pressure
    fuel_tank_press = fuel_dp + Prop.Pc; % psi
    ox_tank_press   = ox_dp + Prop.Pc; % psi

    Press.tank_press = max(fuel_tank_press, ox_tank_press); % psi

    % Find He volume needed 

    prop_mass_SI = Prop.prop_mass / 2.20462; % lbm to kg
    fuel_mass_SI = prop_mass_SI / (Prop.OF + 1); % kg
    ox_mass_SI   = prop_mass_SI - fuel_mass_SI; % kg

    fuel_tank_volume = (fuel_mass_SI / params.fuel_density) / (1 - params.ullage); % m^3
    ox_tank_volume   = (ox_mass_SI / params.ox_density) / (1 - params.ullage); % m^3

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

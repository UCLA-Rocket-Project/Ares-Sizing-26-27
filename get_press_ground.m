function [out, status] = get_press_ground(Prop, Press, params)
  % Pressurant Sizing Script for Ground / GN2 K-Bottles

% static mass balance to find bottle #
% assumes choked flow at bottle, checks if vdot_gas > vdot_prop, adds to bottle # if not
% assumes choked flow at dome, checks if vdot_gas > vdot_prop, adds to dome # if not

% added check that time for bottle from 2k psi to tank pressure > burn time
% uses worst case vdot instead of averaged now
% accounts for isentropic cooling by fixing entropy instead of temperature


  % Unpack inputs from struct
  mdot          = Prop.mdot / 2.20462;          % lbm/s to kg/s
  tank_volume         = Press.fuel_tank_volume + Press.ox_tank_volume; % m^3
  fuel_volume         = Press.fuel_tank_volume;         % m^3
  ox_volume     = Press.ox_tank_volume;     % m^3
  OF                = Prop.OF;
  tank_pressure   = Press.tank_press * 6894.76;   % Pa
  GN2_pressure      = 2000 * 6894.76;      % Pa
  Dome_orifice_area = params.Dome_orifice_area; % m^2, from Cv/SCFM
  A         = params.GN2_bottle_orifice_area;         % m^2, bottle outlet area
  eth_ratio     = params.eth_ratio;
  T_fuel              = params.T_fuel;              % K
  T_ox          = params.T_ox;          % K
  T_nitrogen            = params.T_N2;            % K initial bottle temp

  % Initialize outputs
  out.max_domes = 0;
  out.bottle_number = 0;
  out.t_blowdown = 0;
  status = 0; % 0: Success

  R_nitrogen = 296.8; % J/(kg*K)

  % find fuel density

  rho_water = py.CoolProp.CoolProp.PropsSI('D', 'T', T_fuel, 'P', tank_pressure, 'water');
  rho_ethanol = py.CoolProp.CoolProp.PropsSI('D', 'T', T_fuel, 'P', tank_pressure, 'ethanol');
  rho_fuel = 1/(((1 - eth_ratio) / rho_water) + (eth_ratio / rho_ethanol)); % kg/m^3

  rho_ox = py.CoolProp.CoolProp.PropsSI('D', 'T', T_ox, 'P', tank_pressure, 'oxygen'); % kg/m^3

  % split total mdot by OF, convert each to volumetric flow

  mdot_fuel = mdot / (1 + OF); % kg/s
  mdot_ox = mdot - mdot_fuel;    % kg/s
  vdot_fuel = mdot_fuel / rho_fuel; % m^3/s
  vdot_ox = mdot_ox / rho_ox;         % m^3/s
  vdot_tot = vdot_fuel + vdot_ox;       % m^3/

  % calc burn time
  prop_mass = (fuel_volume * rho_fuel) + (ox_volume * rho_ox); % kg
  burn_time = prop_mass / mdot; % s

  % static mass balance: finds number of bottles & mass of n2 needed to fill tank volume
  % uses worst case density after isnetropic cooling

  s_initial = py.CoolProp.CoolProp.PropsSI('S', 'T', T_nitrogen, 'P', GN2_pressure, 'nitrogen'); % J/(kg*K)

  rho_nitrogen = py.CoolProp.CoolProp.PropsSI('D', 'P', tank_pressure, 'S', s_initial, 'nitrogen'); % kg/m^3
  nitrogen_mass = rho_nitrogen * tank_volume; % kg

  raw_bottle_number = nitrogen_mass / 11; % 11 kg for K-bottle
  bottle_number = ceil(raw_bottle_number);
  fprintf('Number of GN2 bottles needed: %d\n', bottle_number);

  % Transient blowdown for ground
  % track bottle pressure by assuming choked flow at bottle

  bottle_p = zeros(1,1);
  bottle_p(1) = GN2_pressure; % Pa
  delta_t = 0.01;       % s
  m_gas_old = 11;       % kg
  t = 0;              % s
  vdot_array = [];    % m^3/s
  gamma_array = [];
  stop_pressure = tank_pressure; % Pa

  i = 1;
  while bottle_p(i) > stop_pressure

    % fixing state w/ changing pressure & constant entropy (isentropic) instead of constant temp
    rho_nitrogen_t = py.CoolProp.CoolProp.PropsSI('D', 'P', bottle_p(i), 'S', s_initial, 'nitrogen'); % kg/m^3
    Cp = py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', bottle_p(i), 'S', s_initial, 'nitrogen'); % J/(kg*K)
    Cv = py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', bottle_p(i), 'S', s_initial, 'nitrogen'); % J/(kg*K)
    gamma = Cp / Cv;
    gamma_array(i) = gamma;

    % choked flow thru bottle outlet
    mdot_gas = ((A * bottle_p(i) * sqrt(gamma)) / (sqrt(298) * sqrt(R_nitrogen))) * ((gamma + 1) / 2)^(-(gamma + 1) / (2 * (gamma - 1))); % kg/s

    vdot_gas_t = mdot_gas / rho_nitrogen_t; % m^3/s
    vdot_array(i) = vdot_gas_t;

    % subtract gas lost this step & update pressure

    mass_lost = mdot_gas * delta_t; % kg
    m_gas_new = m_gas_old - mass_lost;
    bottle_p(i+1) = bottle_p(i) * m_gas_new / m_gas_old;

    m_gas_old = m_gas_new;
    t = t + delta_t;
    i = i + 1;

  end

  fprintf('Time for pressurant bottle to reach tank pressure: %.4f seconds\n', t);

  % Duration check
  % checks if the bottle actually last whole burn

  fprintf('Burn time: %.4f s\n', burn_time);
  if t < burn_time
    fprintf('bottle reaches tank pressure before burn ends (short by %.4f s)\n', burn_time - t);
    duration_ok = false;
    status = -1; % Sentinel: Reaches tank pressure before burn ends
  else
    fprintf('bottle covers burn time (margin: %.4f s)\n', t - burn_time);
    duration_ok = true;
  end

  % Flow rate check
  % checks minimum GN2 vdot >= propellant vdot if choked at bottle
  % if not adds to bottle_number
  vdot_gas_min = bottle_number * min(vdot_array);

  if vdot_gas_min > vdot_tot
    enough_flow = true;
  else
    enough_flow = false;
  end

  while ~enough_flow

    vdot_gas_min = vdot_gas_min / bottle_number;
    bottle_number = bottle_number + 1;

    fprintf('Insufficient pressurant volumetric flow rate: increasing bottle number to %d\n', bottle_number);

    vdot_gas_min = vdot_gas_min * bottle_number;
    if vdot_gas_min > vdot_tot
      enough_flow = true;
    end
  end

  fprintf('Final bottle count: %d\n', bottle_number);
  fprintf('Duration check: %d, flow rate check : %d\n', duration_ok, enough_flow);

  % Dome sizing/check
  % sizes # of domes needed so minimum GN2 vdot >= propellant vdot if choked at dome

  number_of_domes = [];
  current_dome_number = 1;

  for j = 1:(length(bottle_p) - 1)

    rho_nitrogen_j = py.CoolProp.CoolProp.PropsSI('D', 'T', 298, 'P', bottle_p(j), 'nitrogen'); % kg/m^3

    domes_number_good = true;
    while domes_number_good

      M = py.CoolProp.CoolProp.PropsSI('M', 'T', 298, 'P', bottle_p(j), 'nitrogen'); % kg/mol
      R = 8.31446261815324 / M; % J/(kg*K)

      % choked flow through the dome's orifice

      mdot_dome = Dome_orifice_area * bottle_p(j) * sqrt(gamma_array(j) / (R * 298)) * ((gamma_array(j) + 1)/2)^(-(gamma_array(j) + 1)/(2*(gamma_array(j) - 1))); % kg/s

      vdot_dome = (mdot_dome / rho_nitrogen_j) * current_dome_number; % m^3/s

      if vdot_dome < vdot_tot
        current_dome_number = current_dome_number + 1;
      else
        number_of_domes(end+1) = current_dome_number;
        current_dome_number = 1;
        domes_number_good = false;
      end
    end
  end

  max_domes = max(number_of_domes);
  fprintf('Total number of domes needed: %d\n', max_domes);

% Graph GN2 bottle pressure over time
 t_array = (0:length(bottle_p)-1) * delta_t;

 figure;
 plot(t_array, bottle_p/6894.76, 'LineWidth', 1.5);
 yline(tank_pressure/6894.76, '--r', 'Tank pressure');
 xline(burn_time, '--k', 'Burn time');
 xline(t, ':b', 'Reaches tank pressure');
 xlabel('Time [s]'); ylabel('GN2 bottle pressure [psi]');
 title('GN2 Bottle Pressure vs Time'); grid on;

  % Pack structure output
  out.max_domes = max_domes;
  out.bottle_number = bottle_number;
  out.t_blowdown = t;
end
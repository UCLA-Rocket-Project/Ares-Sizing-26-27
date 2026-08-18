function [out, status] = get_press_ground(Prop, Press, params)
  % Pressurant Sizing Script for Ground / GN2 K-Bottles

% static mass balance to find bottle #
% assumes choked flow at bottle, checks if vdot_gas > vdot_prop, adds to bottle # if not
% assumes choked flow at dome, checks if vdot_gas > vdot_prop, adds to dome # if not

% added check that time for bottle from 2k psi to tank pressure > burn time
% uses worst case vdot instead of averaged now
% polytropic adiabatic ideal blowdown (P/rho^gamma = const), gamma recomputed each step


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
   % fixing state w/ changing pressure & constant entropy (isentropic) instead of constant temp

  s_initial = py.CoolProp.CoolProp.PropsSI('S', 'T', T_nitrogen, 'P', GN2_pressure, 'nitrogen'); % J/(kg*K)

  rho_nitrogen = py.CoolProp.CoolProp.PropsSI('D', 'P', tank_pressure, 'S', s_initial, 'nitrogen'); % kg/m^3
  nitrogen_mass = rho_nitrogen * tank_volume; % kg

  raw_bottle_number = nitrogen_mass / 11; % 11 kg for K-bottle
  bottle_number = ceil(raw_bottle_number);
  fprintf('Number of GN2 bottles needed: %d\n', bottle_number);

  % Transient blowdown for ground
  
    rho_nitrogen_at_GN2_pressure = py.CoolProp.CoolProp.PropsSI('D', 'T', T_nitrogen, 'P', GN2_pressure, 'nitrogen'); % kg/m^3
     bottle_volume = 11 / rho_nitrogen_at_GN2_pressure; % m^3
     total_bottle_volume = bottle_number * bottle_volume; % m^3, all bottles in the bank

    rho_N2_ullage = py.CoolProp.CoolProp.PropsSI('D','P', tank_pressure, 'S', s_initial, 'nitrogen'); % kg/m^3

    prop_mass_kg = Prop.prop_mass / 2.20462; % lbm to kg
    fuel_mass_t = prop_mass_kg / (1 + OF); 
    ox_mass_t = prop_mass_kg - fuel_mass_t; % kg

    V_fuel_ullage_t = fuel_volume - fuel_mass_t / rho_fuel;
    V_ox_ullage_t = ox_volume - ox_mass_t / rho_ox;

    m_bottles_initial = bottle_number * 11; % kg
    m_total = m_bottles_initial + rho_N2_ullage * V_fuel_ullage_t + rho_N2_ullage * V_ox_ullage_t;

% Loop

    rho_old = 11 / bottle_volume; % kg/m^
    T_old = T_nitrogen;
    bottle_p = GN2_pressure;

    while fuel_mass_t > 0 && ox_mass_t > 0
     fuel_mass_t = fuel_mass_t - mdot_fuel*delta_t;
    ox_mass_t = ox_mass_t - mdot_ox*delta_t;
    if fuel_mass_t <= 0 || ox_mass_t <= 0, break; end

     V_fuel_ullage_t = fuel_volume - fuel_mass_t/rho_fuel;
     V_ox_ullage_t = ox_volume - ox_mass_t/rho_ox;

     m_ullage_now = rho_N2_ullage * V_fuel_ullage_t + rho_N2_ullage * V_ox_ullage_t;
     m_bottles_now = m_total - m_ullage_now;
     if m_bottles_now <= 0, break; end

    rho_new = m_bottles_now / total_bottle_volume;
     T_new = T_old * (rho_new/rho_old)^(polytropic_n - 1);
    bottle_p = py.CoolProp.CoolProp.PropsSI('P','D',rho_new,'T',T_new,'nitrogen');

     rho_old = rho_new; T_old = T_new;

    end
    
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
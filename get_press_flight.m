function [out, status] = get_press_flight(Prop, Press, PV_mel, params)

  % Pressurant Sizing Function for Helium / Flight COPV
% Pressurant Sizing Script for Helium / Flight COPV 

% static mass balance to check if 12L holds enough Helium mass
% assumes choked flow at COPV outlet, checks if vdot_gas > vdot_prop
% assumes choked flow at dome, checks if vdot_gas > vdot_prop, adds to dome # if not

% added check that time for bottle from 4.5k to tank pressure > burn time
% uses worst case vdot instead of averaged
% accounts for isentropic cooling by fixing entropy instead of temperature

  % Unpack inputs from struct
  mdot        = Prop.mdot / 2.20462;        % lbm/s to kg/s
  tank_volume     = Press.fuel_tank_volume + Press.ox_tank_volume;     % m^3
  fuel_volume     = Press.fuel_tank_volume;     % m^3
  ox_volume     = Press.ox_tank_volume;     % m^3
  OF          = Prop.OF;
  tank_pressure   = Press.tank_press * 6894.76;   % psi to Pa
  COPV_pressure   = params.P_He_init * 6894.76;   % Pa
  COPV_volume     = PV_mel.copv_v / 1000;     % L to m^3
  Dome_orifice_area = params.Dome_orifice_area; % m^2, from Cv/SCFM
  eth_ratio     = params.eth_ratio;
  T_fuel        = params.T_fuel;        % K
  T_ox        = params.T_ox;        % K
  T_helium      = params.T_He;      % K initial COPV temp

  % Initialize outputs
  out.max_domes = 0;
  out.t_cross = 0;
  out.t_blowdown = 0;
  status = 0; % 0: Success

  R_helium = 2077; % J/(kg*K), specific gas constant for Helium

  % find fuel density
  rho_water = py.CoolProp.CoolProp.PropsSI('D', 'T', T_fuel, 'P', tank_pressure, 'water');
  rho_ethanol = py.CoolProp.CoolProp.PropsSI('D', 'T', T_fuel, 'P', tank_pressure, 'ethanol');
  rho_fuel = 1/(((1 - eth_ratio) / rho_water) + (eth_ratio / rho_ethanol)); % kg/m^3
  rho_ox = py.CoolProp.CoolProp.PropsSI('D', 'T', T_ox, 'P', tank_pressure, 'oxygen'); % kg/m^3

  % split total mdot by OF, convert each to volumetric flow

  mdot_fuel = mdot / (1 + OF);
  mdot_ox = mdot - mdot_fuel;
  vdot_fuel = mdot_fuel / rho_fuel;
  vdot_ox = mdot_ox / rho_ox;
  vdot_tot = vdot_fuel + vdot_ox; % m^3/s

  % calc burn time

  prop_mass = (fuel_volume * rho_fuel) + (ox_volume * rho_ox);
  burn_time = prop_mass / mdot; % s

  % static mass balance: checks if 12L COPV can hold enough helium mass to fill tanks (boolean)
  % uses worst case density after isnetropic cooling

  s_initial = py.CoolProp.CoolProp.PropsSI('S', 'T', T_helium, 'P', COPV_pressure, 'helium'); % J/(kg*K)
  rho_helium_tank = py.CoolProp.CoolProp.PropsSI('D', 'P', tank_pressure, 'S', s_initial, 'helium'); % kg/m^3
  helium_mass_needed = rho_helium_tank * tank_volume; % kg

  rho_helium_full = py.CoolProp.CoolProp.PropsSI('D', 'P', COPV_pressure, 'S', s_initial, 'helium'); % kg/m^3
  helium_mass_available = COPV_volume * rho_helium_full; % kg

  if helium_mass_needed > helium_mass_available
    fprintf('Insufficient helium mass: need %.4f kg, COPV holds %.4f kg\n', helium_mass_needed, helium_mass_available);
    enough_helium = false;
    status = -2; % Sentinel: Insufficient helium mass
  else
    fprintf('Sufficient helium in COPV (need %.4f kg, have %.4f kg)\n', helium_mass_needed, helium_mass_available);
    enough_helium = true;
  end

  % Transient blowdown for flight
  % track bottle pressure by assuming choked flow at COPV outlet

  bottle_p = zeros(1,1);
  bottle_p(1) = COPV_pressure; % Pa
  A = 0.00114; % m^2 COPV outlet orifice area
  delta_t = 0.0001; % s
  m_gas_old = helium_mass_available; % kg
  t = 0; % s
  vdot_array = [];
  gamma_array = [];
  temp_array = []; % K
  stop_pressure = 100*6894.76; % Pa lower than tank pressure to see full blowdown in graph

  idx_cross = []; % timestep where bottle_p first drops below tank_pressure

  i = 1;
  while bottle_p(i) > stop_pressure

     % fixing state w/ changing pressure & constant entropy (isentropic) instead of constant temp

    rho_helium_t = py.CoolProp.CoolProp.PropsSI('D', 'P', bottle_p(i), 'S', s_initial, 'helium'); % kg/m^3
    Cp = py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', bottle_p(i), 'S', s_initial, 'helium');
    Cv = py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', bottle_p(i), 'S', s_initial, 'helium');
    gamma = Cp / Cv;
    gamma_array(i) = gamma;

    T_helium_t = py.CoolProp.CoolProp.PropsSI('T', 'P', bottle_p(i), 'S', s_initial, 'helium'); % K
    temp_array(i) = T_helium_t;

    % choked flow thru COPV outlet
    mdot_gas = ((A * bottle_p(i) * sqrt(gamma)) / (sqrt(298) * sqrt(R_helium))) * ((gamma + 1) / 2)^(-(gamma + 1) / (2 * (gamma - 1))); % kg/s

    vdot_array(i) = mdot_gas / rho_helium_t;

    % subtract gas lost this step & update pressure

    mass_lost = mdot_gas * delta_t;
    m_gas_new = m_gas_old - mass_lost;
    bottle_p(i+1) = bottle_p(i) * m_gas_new / m_gas_old;

    % log timestep where COPV pressure drops below tank pressure
    if isempty(idx_cross) && bottle_p(i+1) <= tank_pressure
      idx_cross = i+1;
    end

    m_gas_old = m_gas_new;
    t = t + delta_t;
    i = i + 1;

  end

  t_cross = (idx_cross - 1) * delta_t; % s, time at which COPV drops below tank pressure
  fprintf('Time for COPV to drop below tank pressure: %.4f s\n', t_cross);
  fprintf('Time for COPV to reach stop_pressure (full blowdown): %.4f s\n', t);

  % Duration check
  % checks if 12L actually last whole burn

  fprintf('Burn time: %.4f s\n', burn_time);
  if t_cross < burn_time
    fprintf('COPV drops below tank pressure before burn ends (short by %.4f s)\n', burn_time - t_cross);
    duration_ok = false;
    if status == 0
      status = -1; % Sentinel: Reaches tank pressure before burn ends
    end
  else
    fprintf('COPV stays above tank pressure for full burn (margin: %.4f s)\n', t_cross - burn_time);
    duration_ok = true;
  end

  % Flow rate check
  % checks minimum Helium vdot >= propellant vdot if choked at COPV

  vdot_gas_min = min(vdot_array(1:idx_cross-1));

  if vdot_gas_min > vdot_tot
    enough_flow = true;
  else
    enough_flow = false;
  end

  if ~enough_flow
    fprintf('Insufficient pressurant volumetric flow rate.\n');
  end
  fprintf('Flow rate OK (worst-case, pre-crossing): %d\n', enough_flow);

  % Dome sizing/check
  % sizes # of domes needed so minimum Helium vdot >= propellant vdot if choked at dome

  number_of_domes = [];
  current_dome_number = 1;

  for j = 1:(idx_cross - 1)

    rho_helium_j = py.CoolProp.CoolProp.PropsSI('D', 'P', bottle_p(j), 'S', s_initial, 'helium'); % kg/m^3

    domes_number_good = true;
    while domes_number_good

      M = py.CoolProp.CoolProp.PropsSI('M', 'P', bottle_p(j), 'S', s_initial, 'helium'); % kg/mol
      R = 8.31446261815324 / M;

      % choked flow through the dome's orifice
      mdot_dome = Dome_orifice_area * bottle_p(j) * sqrt(gamma_array(j) / (R * temp_array(j))) * ((gamma_array(j) + 1)/2)^(-(gamma_array(j) + 1)/(2*(gamma_array(j) - 1))); % kg/s

      vdot_dome = (mdot_dome / rho_helium_j) * current_dome_number;

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

  % Graph COPV pressure & temperature
  t_array = (0:length(bottle_p)-1) * delta_t;
  t_array_temp = (0:length(temp_array)-1) * delta_t;

  figure;

  subplot(2,1,1);
  plot(t_array, bottle_p/6894.76, 'LineWidth', 1.5);
  yline(tank_pressure/6894.76, '--r', 'Tank pressure');
  xline(burn_time, '--k', 'Burn time');
  xline(t_cross, ':b', 'Crosses tank pressure');
  xlabel('Time [s]'); ylabel('COPV pressure [psi]');
  title('COPV Pressure vs Time'); grid on;

  subplot(2,1,2);
  plot(t_array_temp, temp_array, 'LineWidth', 1.5);
  xline(burn_time, '--k', 'Burn time');
  xline(t_cross, ':b', 'Crosses tank pressure');
  xlabel('Time [s]'); ylabel('COPV gas temperature [K]');
  title('COPV Temperature vs Time'); grid on;

  % Pack structure output
  out.max_domes = max_domes;
  out.t_cross = t_cross;
  out.t_blowdown = t;

end
function [out, status] = get_press_flight(Prop, Press, PV_mel, params)

  % Pressurant Sizing Function for Helium / Flight COPV
% Pressurant Sizing Script for Helium / Flight COPV 

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
  polytropic_n = params.polytropic_n; 

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

  %% Transient blowdown for flight

   % find total system helium mass

prop_mass_kg = Prop.prop_mass / 2.20462; % lbm to kg
fuel_mass_t = prop_mass_kg / (1 + OF); % kg 
ox_mass_t = prop_mass_kg - fuel_mass_t; % kg

V_fuel_ullage_t = Press.fuel_tank_volume - fuel_mass_t / rho_fuel; % m^3
V_ox_ullage_t = Press.ox_tank_volume - ox_mass_t / rho_ox; %m^3 

burn_time = prop_mass_kg / mdot; %s

rho_He_ullage = py.CoolProp.CoolProp.PropsSI('D', 'P', tank_pressure, 'S', s_initial, 'helium'); % kg/m^3

m_helium_total = helium_mass_available + (rho_He_ullage * V_fuel_ullage_t) + (rho_He_ullage * V_ox_ullage_t); 

 % initializing

   delta_t = 0.01; % s
  bottle_p = COPV_pressure; % Pa
  T_copv = T_helium; % K
  rho_helium_now = rho_helium_full; % kg/m^3

  t = 0; % s
  i = 1;

  bottle_p_array = COPV_pressure; % Pa
  temp_array = T_helium; % K
  thrust_array = []; % N
  Pc_array = []; % Pa
  mdot_array = []; % kg/s

  idx_cross = []; % timestep index where bottle_p first drops below tank_pressure
  t_cross = []; % s

  in_blowdown = false;

  % run until prop mass reaches 0 

  while fuel_mass_t > 0 && ox_mass_t > 0

    if ~in_blowdown
      % phase 1 
      % regulated at tank pressure
      % tank pressure  at dome setpoint
    
      mdot_fuel_i = mdot_fuel;
      mdot_ox_i = mdot_ox;
      thrust_array(i) = NaN; 
      Pc_array(i) = Prop.Pc * 6894.76; % Pa
      mdot_array(i) = mdot;

    else
      % phase 2
      % blowdown, tank pressure = COPV pressure
      % systemSolver called to iterate system for chamber pressure

      P_tank = bottle_p; % Pa

      [thrust_i, mdot_i, of_i, Pc_i, ~, ~] = systemSolver( ...
          rho_fuel, rho_ox, P_tank, P_tank, CdA_fuel_inj, CdA_ox_inj, ...
          params.P_amb * 6894.76, A_throat, params.Cstar_eff, params.Ctau_eff, ...
          A_exit, params);

      thrust_array(i) = thrust_i; % N
      Pc_array(i) = Pc_i; % Pa
      mdot_array(i) = mdot_i; % kg/s

      mdot_fuel_i = mdot_i / (1 + of_i); % kg/s
      mdot_ox_i = mdot_i - mdot_fuel_i;  % kg/s
    end

    % Deplete prop mass by mdot * dt

    fuel_mass_t = fuel_mass_t - mdot_fuel_i * delta_t; % kg
    ox_mass_t = ox_mass_t - mdot_ox_i * delta_t; % kg

    if fuel_mass_t <= 0 || ox_mass_t <= 0
      break
    end

    % Update ullage volumes

    V_fuel_ullage_t = Press.fuel_tank_volume - fuel_mass_t / rho_fuel; % m^3
    V_ox_ullage_t = Press.ox_tank_volume - ox_mass_t / rho_ox; % m^3

    % subtract mass required in ullages (at current tank pressure) from the fixed total, remainder is left in the COPV 

    m_fuel_tank = rho_He_ullage * V_fuel_ullage_t; % kg
    m_ox_tank = rho_He_ullage * V_ox_ullage_t; % kg
    m_copv_now = m_helium_total - (m_fuel_tank + m_ox_tank); % kg

    if m_copv_now <= 0
      break
    end

    % COPV pressure
    %  density from remaining mass in COPV & volume
    % temperature from polytropic relation w/ density ratio
    % pressure from CoolProp at (density, T)

    rho_helium_prev = rho_helium_now; % kg/m^3
    rho_helium_now = m_copv_now / COPV_volume; % kg/m^3

    T_copv = T_copv * (rho_helium_now / rho_helium_prev) ^ (polytropic_n - 1); % K

    bottle_p = py.CoolProp.CoolProp.PropsSI('P','D',rho_helium_now,'T',T_copv,'helium'); % Pa

    bottle_p_array(i+1) = bottle_p; % Pa
    temp_array(i+1) = T_copv; % K

    % Log time COPV goes into blowdown

    if ~in_blowdown && bottle_p <= tank_pressure
      in_blowdown = true;
      idx_cross = i+1;
      t_cross = (i+1) * delta_t; % s
      fprintf('Entering blowdown at t = %.4f s (COPV pressure <= tank pressure)\n', t_cross);
    end

    t = t + delta_t; % s
    i = i + 1;

  end

  fprintf('Sim complete: t = %.4f s\n', t);

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


  %% Dome sizing
  % sizes # of domes needed so minimum Helium vdot >= propellant vdot if choked at dome

    if isempty(idx_cross)
    dome_check_range = 1:length(bottle_p_array);
  else
    dome_check_range = 1:(idx_cross-1);
  end

  number_of_domes = [];
  current_dome_number = 1;

  for j = dome_check_range

    rho_helium_j = py.CoolProp.CoolProp.PropsSI('D', 'P', bottle_p_array(j), 'S', s_initial, 'helium'); % kg/m^3

    domes_number_good = true;
    while domes_number_good

          Cp_j = py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', bottle_p_array(j), 'S', s_initial, 'helium');
      Cv_j = py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', bottle_p_array(j), 'S', s_initial, 'helium');
      gamma_j = Cp_j / Cv_j;

      M = py.CoolProp.CoolProp.PropsSI('M', 'P', bottle_p_array(j), 'S', s_initial, 'helium'); % kg/mol
      R = 8.31446261815324 / M;

      % choked flow through the dome's orifice
      mdot_dome = Dome_orifice_area * bottle_p_array(j) * sqrt(gamma_j / (R * temp_array(j))) * ((gamma_j + 1)/2)^(-(gamma_j + 1)/(2*(gamma_j - 1))); % kg/s

      vdot_dome = (mdot_dome / rho_helium_j) * current_dome_number; % m^3/s

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

 %% Graph

  t_array = (0:length(bottle_p_array)-1) * delta_t;

  figure;
  subplot(3,1,1);
  plot(t_array, bottle_p_array/6894.76, 'LineWidth', 1.5);
  yline(tank_pressure/6894.76, '--r', 'Tank pressure (dome setpoint)');
  xline(burn_time, '--k', 'Nominal burn time');
  if ~isempty(t_cross)
    xline(t_cross, ':b', 'Enters blowdown');
  end
  xlabel('Time [s]'); ylabel('COPV pressure [psi]');
  title('COPV Pressure vs Time'); grid on;

  subplot(3,1,2);
  plot(t_array, temp_array, 'LineWidth', 1.5);
  xlabel('Time [s]'); ylabel('COPV temperature [K]');
  title('COPV Temperature vs Time'); grid on;

  subplot(3,1,3);
  t_solved = t_array(1:length(Pc_array));
  plot(t_solved, Pc_array/6894.76, 'LineWidth', 1.5);
  xlabel('Time [s]'); ylabel('Chamber pressure [psi]');
  title('Chamber Pressure vs Time (Phase 1 nominal, Phase 2 solved)'); grid on;

 %% Pack structure output

  out.max_domes = max_domes;
  out.t_cross = t_cross;
  out.t_blowdown = t;
  out.bottle_p = bottle_p_array; % Pa
  out.temp_array = temp_array; % K
  out.thrust_array = thrust_array; % N
  out.Pc_array = Pc_array; % Pa
  out.mdot_array = mdot_array; % kg/s
  out.duration_ok = duration_ok;

end
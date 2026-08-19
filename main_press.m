% main_press.m
%
% driver for pressurant sizing (flight COPV + ground GN2),

clear; clc;

%% Parameter inputs 
Pc = 370;         % psi
OF = 1.4;
eps = 5;
prop_mass = 100;  % lbm

%% run functions 
params = input_parameters();
params.polytropic_n_He = 1.67; 
params.polytropic_n_N2 = 1.4;

Prop = struct('OF', OF, 'Pc', Pc, 'eps', eps, 'prop_mass', prop_mass);
Prop = run_CEA(Prop, params);
[Prop, Press] = run_press(Prop, params);
PV_mel.copv_v = 12;

% geometry 
A_throat = Prop.A_t * 0.00064516; % in^2 to m^2
A_exit = params.A_e * 0.00064516; % in^2 to m^2

% CEA_obj

card_str = sprintf(['fuel C2H5OH(L) C 2 H 6 O 1\n', ...
  'h,cal=-66370.0 t(k)=298.00 wt%%=75.00\n', ...
  'fuel water H 2.0 O 1.0 wt%%=25.00\n', ...
  'h,cal=-68308. t(k)=298.00 rho,g/cc=0.9998']);
py.rocketcea.cea_obj.add_new_fuel('ETHANOL_WATER_75_25(L)', card_str);

CEA_obj = py.rocketcea.cea_obj.CEA_Obj(pyargs('oxName', 'LOX', 'fuelName', 'ETHANOL_WATER_75_25(L)'));

%% Run pressurant sizing 
[flight_out, flight_status] = get_press_flight(Prop, Press, PV_mel, params, CEA_obj, A_throat, A_exit);
[ground_out, ground_status] = get_press_ground(Prop, Press, params);

%  Summary 

fprintf('\n=== Flight (COPV) ===\n');
fprintf('status = %d\n', flight_status);
fprintf('t_cross = %.4f s, t_blowdown = %.4f s\n', flight_out.t_cross, flight_out.t_blowdown);
fprintf('max_domes = %d\n', flight_out.max_domes);

fprintf('\n=== Ground (GN2) ===\n');
fprintf('status = %d\n', ground_status);
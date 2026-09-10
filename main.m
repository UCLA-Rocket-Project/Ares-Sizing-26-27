%% 
% main

% Sweeps prop mass, OF, Pc, and area ratio (eps). 
% Pipeline order
% per iteration: run_CEA -> run_press -> get_PV_mel -> get_ShockLoads -> get_recLoads -> get_dryMass -> get_apogee

%Units 
% params are mostly imperial, some metric
% run_CEA does calcs internally in metric, converts to imperial
% outputs Prop struct all imperial
% run_press converts to SI internally for CoolProp/He calcs, but ouputs imperial 
% outputs Press struct all imperial
% get_PV_mel / get_ShockLoads / get_recLoads / get_dryMass / get_apogee are all imperial

% Plotting results: graph of vehicle length and prop mass vs apogee

% MEOP: 600 psia
% Chamber Pressure: 370 psia
% Manifold Pressure: 444
% Thrust: 1450
% OF: 1.4
% Mdot: 6.639 lbm/s
% area ratio: 5
% eth/water: 75/25
% Burn time: 15.06
% Isp: 219s
% Prop mass: 100lbs
% Ctau efficiency: 98%
% C* efficiency: 90%

%% Setup
% CSV path handling and writing

base_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(base_dir, 'out');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

data_csv_path = fullfile(out_dir, 'output.csv'); % every iter
filtered_csv_path = fullfile(out_dir, 'output_filtered.csv'); % for successful iterations (meets all filters)
log_path = fullfile(out_dir, 'log.txt');

params = input_parameters();

col_names = {'prop_mass', 'eth_ratio', 'OF', 'Pc', 'eps', 'mdot','thrust', 'Isp', 't_b','tank_press', 'V_He', 'dry_mass', ...
            'fuel_tank_length', 'ox_tank_length', 'tank_wall', 'vehicle_length', 'apogee', 'fail_code'};

%Filter summary
% fail_code = 0 means apogee computed, all filters passed
% -1 means COPV volume infeasible
% -2 means tank wall infeasible
% -3 means dry mass convergence didnt converge
% -4 means OTRS did not pass
% removed caps on mdot, tank pressure, burn time, and TWR from last year
%all variables are evaluated and checked in this script, no arbitrary threshold needed

writecell(col_names, data_csv_path);
writecell(col_names, filtered_csv_path);


%% Sweep Ranges
% mass_dist = 85:5:100; %lbm
OF_dist = 1.2:0.02:1.6;
% Pc_dist = 200:10:600; %psi
% eps_dist = 3:0.5:5;
eth_dist = 0.75:0.02:0.95;

%mass_dist = 95:5:100; %lbm
%OF_dist = 1.3:0.02:1.4;
%Pc_dist = 200:10:400; %psi
%eps_dist = 4:0.5:5;

tol = 0.1;
max_iter = 50;

%% Preload Cd vs Mach data once per prop_mass
cd_files = containers.Map('KeyType', 'double', 'ValueType', 'any');
% for m = mass_dist
%     switch m
%         case 85
%             fname = 'MvsCd_data_85.CSV';
%         case 90
%             fname = 'MvsCd_data_90.CSV';
%         case 95
%             fname = 'MvsCd_data_95.CSV';
%         otherwise
            fname = 'MvsCd_data_100.csv';
    % end
    data = readmatrix(fullfile(base_dir, 'input', fname));
    [M_u, uniq_idx] = unique(data(:,1));
    cd_files(100) = struct('M_data', M_u, 'Cd_data', data(uniq_idx, 2));
% end

%% Linear index into grid
% [MM, OO, PP, EE] = ndgrid(mass_dist, OF_dist, Pc_dist, eps_dist);
% combos = [MM(:), OO(:), PP(:), EE(:)];
[OO, EE] = ndgrid(OF_dist, eth_dist);
combos = [EE(:), OO(:)];
it_ct = size(combos, 1);

results = cell(it_ct, 1);

h = waitbar(0, 'Running');
dq = parallel.pool.DataQueue;
afterEach(dq, @(~) updateWaitbar(h, it_ct));

tic;
eth_rho = params.ethanol_density;
water_rho = params.water_density;
parfor it = 1:it_ct
    eth_ratio  = combos(it, 1);
    prop_mass = 100; % combos(it, 1);
    OF        = combos(it, 2);
    Pc        = 370; % combos(it, 3);
    eps       = 5; % combos(it, 4);
    fuel_density = mass_fraction(eth_ratio, eth_rho, water_rho); % kg/m3
    fprintf('Running iteration %d with eth_ratio = %.2f, fuel_density = %.2f\n', it, eth_ratio, fuel_density);

    fail_code = 0;
    Prop  = struct('OF', OF, 'eth_ratio', eth_ratio, 'Pc', Pc, 'eps', eps, 'prop_mass', prop_mass);
    Press = struct('tank_press', NaN, 'V_He', NaN, 'fuel_tank_volume', NaN, 'ox_tank_volume', NaN); %#ok
    dry_mass = NaN;
    apogee = NaN; %#ok

    % for ablative char depth -> mass calc
    % [injector, nozzle entrance]
    Abl = struct('T0', [NaN, NaN], 'Pr', [NaN, NaN], 'gamma', [NaN, NaN], 'Cp', [NaN, NaN], 'Pc', NaN); %#ok

    %% run_CEA
    [Prop, Abl] = run_CEA(Prop, params, eth_ratio);

    %% run_press
    [Prop, Press] = run_press(Prop, params, fuel_density);

    %% get_abl
    abl_mass = NaN;
    char_depth = NaN;
    [abl_mass, char_depth] = get_abl(Abl, Prop);

    %% get_PV_mel
    PV_mel = get_PV_mel(Press, fuel_density, Press.tank_press, Press.V_He);

    if isnumeric(PV_mel) && PV_mel == -1
        fail_code = -1;
        PV_mel = struct('fuel_l', NaN, 'ox_l', NaN, 'tank_wall', NaN);
    elseif isnumeric(PV_mel) && PV_mel == -2
        fail_code = -2;
        PV_mel = struct('fuel_l', NaN, 'ox_l', NaN, 'tank_wall', NaN);
    end

    %% get_ShockLoads & get_recLoads (iterate dry mass until unchanging)
    if fail_code == 0
        dry_mass = 120; % Initializing dry mass
        old_dryMass = 0; % Set old dry mass to zero for the first iteration
        iter = 0;

        while abs(dry_mass - old_dryMass) > tol && iter < max_iter
            old_dryMass = dry_mass; % Update old dry mass

            % Calculate new dry mass based on current shock loads
            [f_drogue, f_main] = get_ShockLoads(dry_mass);
            recLoads = get_highestLoad(f_drogue, f_main, PV_mel);
            dry_mass = abl_mass + get_dryMass(recLoads, PV_mel);
            iter = iter + 1;
        end

        if iter == max_iter
            fail_code = -3;
        end
    end

    %% get_apogee
    cd_s = cd_files(prop_mass);
    M_data  = cd_s.M_data;
    Cd_data = cd_s.Cd_data;

    if fail_code == 0
        apogee = get_apogee(Prop, params, Cd_data, M_data, dry_mass); % ft
        if apogee == -3
            fail_code = -4;
        end
    else
        apogee = fail_code;
    end

    vehicle_length = 158.5 + PV_mel.fuel_l + PV_mel.ox_l; % in
    % ^ 158.5 in is Pandora's length without tank barrels

    %% Row assembly
    results{it} = {prop_mass, eth_ratio, OF, Pc, eps, Prop.mdot, Prop.Thrust, Prop.Isp, Prop.t_b, ...
        Press.tank_press, Press.V_He, dry_mass, PV_mel.fuel_l, PV_mel.ox_l, PV_mel.tank_wall, ...
        vehicle_length, apogee, fail_code};

    send(dq, 1);
end
elapsed_time = toc;

close(h);

function updateWaitbar(h, N)
    persistent completed;
    if isempty(completed) || nargin == 0
        completed = 0;
    end
    completed = completed + 1;
    waitbar(completed / N, h, sprintf('Progress: %d / %d', completed, N));
    if completed >= N
        completed = []; % Reset for future use
    end
end


%% Serial write after the parallel loop
buffer_all = vertcat(results{:});
writecell(buffer_all, data_csv_path, 'WriteMode', 'append');

pass_mask = cellfun(@(r) r{end} == 0, results);
buffer_pass = buffer_all(pass_mask, :);
if ~isempty(buffer_pass)
    writecell(buffer_pass, filtered_csv_path, 'WriteMode', 'append');
end

fprintf('%d / %d iterations complete\n', it_ct, it_ct);
%% Optimized Output

results_filtered = readtable(filtered_csv_path);

if isempty(results_filtered)
    fprintf('No valid combinations found, check filters & ranges\n')

else
    [max_apogee, opt_idx] = max(results_filtered.apogee);

    log = sprintf([...
        '%s\n' ...
        'ARES 26-27 SIZING SWEEP COMPLETE\n' ...
        '%s\n\n' ...
        'Trial combinations: %d\n' ...
        'Valid combinations found: %d\n' ...
        'Time to complete: %.2f seconds\n\n' ...
        'Maximum apogee found: %.2f ft\n\n' ...
        'Prop Mass: %.2f lb\n' ...
        'Thrust: %.2f lbf\n' ...
        'Burn Time: %.2f s\n' ...
        'Mass Flow Rate: %.3f lbm/s\n' ...
        'Ethanol concentration: %.2f\n' ...
        'OF Ratio: %.2f\n' ...
        'Area Ratio: %.2f\n'...
        'Chamber Pressure: %.2f psi\n' ...
        'Tank Pressure: %.2f psi\n' ...
        'Isp: %.2f s\n\n' ...
        'Dry Mass: %.2f lb\n' ...
        'Fuel Barrel Length: %.2f in\n' ...
        'Ox Barrel Length: %.2f in\n' ...
        'Tank wall thickness: %.4f in\n'...
        'Vehicle Length: %.2f in\n'], ...
        repmat('=',1,80), repmat('=',1,80), ...
        it_ct, height(results_filtered), elapsed_time, ...
        max_apogee, ...
        results_filtered.prop_mass(opt_idx), results_filtered.thrust(opt_idx), ...
        results_filtered.t_b(opt_idx), results_filtered.mdot(opt_idx), ...
        results_filtered.eth_ratio(opt_idx), results_filtered.OF(opt_idx), results_filtered.eps(opt_idx), results_filtered.Pc(opt_idx), ...
        results_filtered.tank_press(opt_idx), results_filtered.Isp(opt_idx), ...
        results_filtered.dry_mass(opt_idx), results_filtered.fuel_tank_length(opt_idx), ...
        results_filtered.ox_tank_length(opt_idx), results_filtered.tank_wall(opt_idx), results_filtered.vehicle_length(opt_idx));

    fid = fopen(log_path, 'w');
    fprintf(fid, '%s', log);

    fclose(fid);

    disp(log);
end


% only ran once for optimum combo

best_prop_mass = 100; %results_filtered.prop_mass(opt_idx);
best_OF = results_filtered.OF(opt_idx);
best_Pc = 370; %results_filtered.Pc(opt_idx);
best_eps = 5; %results_filtered.eps(opt_idx);
best_eth_ratio = results_filtered.eth_ratio(opt_idx);
best_fuel_density = mass_fraction(best_eth_ratio, params.ethanol_density, params.water_density); % kg/m3

% Re-run the pipeline once for this combo to get the full Prop,Press, PV_mel structs

Prop = struct('OF', best_OF, 'Pc', best_Pc, 'eps', best_eps, 'prop_mass', best_prop_mass);
Prop = run_CEA(Prop, params, best_eth_ratio);
[Prop, Press] = run_press(Prop, params, best_fuel_density);
PV_mel = get_PV_mel(Press, best_fuel_density, Press.tank_press, Press.V_He);


%% Press sizing (WIP)

%[flight_out, flight_status] = get_press_flight(Prop, Press, PV_mel, params);
%[ground_out, ground_status] = get_press_ground(Prop, Press, params);

%fprintf('PRESSURANT SIZING RESULTS (optimal combo)\n');

%fprintf('-- Flight (COPV / Helium) --\n');
%%fprintf('Status: %d (0 = OK, -1 = COPV drops below tank pressure before burn ends, -2 = insufficient helium mass)\n', flight_status);
%fprintf('Domes needed: %d\n', flight_out.max_domes);
%%fprintf('Time COPV stays above tank pressure: %.4f s\n\n', flight_out.t_cross);
%fprintf('Full COPV time, w/ blowdown: %.4f s\n\n', flight_out.t_blowdown)

%fprintf('-- Ground (GN2 K-Bottles) --\n');
%fprintf('Status: %d (0 = OK, -1 = bottle reaches tank pressure before burn ends)\n', ground_status);
%fprintf('Bottles needed: %d\n', ground_out.bottle_number);
%%fprintf('Domes needed: %d\n', ground_out.max_domes);
%fprintf('Blowdown time: %.4f s\n', ground_out.t_blowdown);

%% Tube mass for optimum combo to back calc layers

best_dry_mass = results_filtered.dry_mass(opt_idx);
[f_drogue, f_main] = get_ShockLoads(best_dry_mass);
recLoads = get_highestLoad(f_drogue, f_main, PV_mel);
[dry_mass, tube_masses] = get_dryMass(recLoads, PV_mel);

fprintf('UBT mass: %.2f lb\n', tube_masses.ubt_m);
fprintf('LBT mass: %.2f lb\n', tube_masses.lbt_m);
fprintf('ITS mass: %.2f lb\n', tube_masses.its_m);

%% RSE File Gen

if best_prop_mass == 85
data = readmatrix(fullfile(base_dir, 'input/MvsCd_data_85.CSV'));
elseif best_prop_mass == 90
data = readmatrix(fullfile(base_dir, 'input/MvsCd_data_90.CSV'));
elseif best_prop_mass == 95
data = readmatrix(fullfile(base_dir, 'input/MvsCd_data_95.CSV'));
else
data = readmatrix(fullfile(base_dir, 'input/MvsCd_data_100.CSV'));
end
Cd_data = data(:,2);
M_data = data(:,1);
[M_data, uniq_idx] = unique(M_data);
Cd_data = Cd_data(uniq_idx);

% get_RSE_files(Prop, params, dry_mass, PV_mel, Cd_data, M_data, "Ares", "UCLA_Rocket_Project", out_dir);

%% Plotting

if ~isempty(results_filtered) 
    figure;
    plot(results_filtered.vehicle_length, results_filtered.apogee, '.');
    xlabel('Vehicle Length (in)');
    ylabel('Apogee (ft)');

    figure;
    plot(results_filtered.prop_mass, results_filtered.apogee, '.');
    xlabel('Prop mass (lb)');
    ylabel('Apogee (ft))');

end

tol = 1e-6; % float compare tolerance for matching fixed sweep values

% 1-4: fixed eps=5, prop_mass=100, Pc vs OF vs {mdot, t_b, apogee, V_He}
plot_slice(results_filtered, 'Pc','OF','mdot',  {'eps','prop_mass'}, [5,100], tol, ...
  'mdot vs Pc & OF (eps=5, prop\_mass=100)', 'Pc (psi)','OF','mdot (lbm/s)');
plot_slice(results_filtered, 'Pc','OF','t_b',   {'eps','prop_mass'}, [5,100], tol, ...
  'Burn Time vs Pc & OF (eps=5, prop\_mass=100)', 'Pc (psi)','OF','t\_b (s)');
plot_slice(results_filtered, 'Pc','OF','apogee', {'eps','prop_mass'}, [5,100], tol, ...
  'Apogee vs Pc & OF (eps=5, prop\_mass=100)', 'Pc (psi)','OF','Apogee (ft)');
plot_slice(results_filtered, 'Pc','OF','V_He',  {'eps','prop_mass'}, [5,100], tol, ...
  'V\_He vs Pc & OF (eps=5, prop\_mass=100)', 'Pc (psi)','OF','V\_He (L)');

% 5: fixed eps=5, OF=1.38, Pc vs prop_mass vs V_He
plot_slice(results_filtered, 'Pc','prop_mass','V_He', {'eps','OF'}, [5,1.38], tol, ...
  'V\_He vs Pc & Prop Mass (eps=5, OF=1.38)', 'Pc (psi)','Prop Mass (lb)','V\_He (L)');

% 6: fixed prop_mass=100, Pc vs eps vs mdot, layered over OF = 1.2, 1.3, 1.38
plot_layered_slice(results_filtered, 'Pc','eps','mdot', 'prop_mass', 100, tol, ...
  'OF', [1.2, 1.3, 1.38], ...
  'mdot vs Pc & eps, layered by OF (prop\_mass=100)', 'Pc (psi)','eps','mdot (lbm/s)');

% 7-10: fixed prop_mass=100, OF=OF_fixed, Pc vs eps vs {mdot, t_b, apogee, V_He}
OF_fixed = 1.38; 
plot_slice(results_filtered, 'Pc','eps','mdot',  {'prop_mass','OF'}, [100,OF_fixed], tol, ...
  sprintf('mdot vs Pc & eps (prop\\_mass=100, OF=%.2f)', OF_fixed), 'Pc (psi)','eps','mdot (lbm/s)');
plot_slice(results_filtered, 'Pc','eps','t_b',   {'prop_mass','OF'}, [100,OF_fixed], tol, ...
  sprintf('Burn Time vs Pc & eps (prop\\_mass=100, OF=%.2f)', OF_fixed), 'Pc (psi)','eps','t\_b (s)');
plot_slice(results_filtered, 'Pc','eps','apogee', {'prop_mass','OF'}, [100,OF_fixed], tol, ...
  sprintf('Apogee vs Pc & eps (prop\\_mass=100, OF=%.2f)', OF_fixed), 'Pc (psi)','eps','Apogee (ft)');
plot_slice(results_filtered, 'Pc','eps','V_He',  {'prop_mass','OF'}, [100,OF_fixed], tol, ...
  sprintf('V\\_He vs Pc & eps (prop\\_mass=100, OF=%.2f)', OF_fixed), 'Pc (psi)','eps','V\_He (L)');

%% Local Functions

% Returns struct with highest equivalent axial load and highest bending
% moment from recovery on each of the 3 carbon tubes
function load = get_highestLoad(f_drogue, f_main,PV_mel)
drogueLoad = get_recLoads(f_drogue,PV_mel,"drogue");
mainLoad = get_recLoads(f_main,PV_mel,"main");

% Lengths

load.ubt_l = drogueLoad.ubt_l;
load.lbt_l = drogueLoad.lbt_l;
load.its_l = drogueLoad.its_l;

% UBT Axial

load.ubt_axial = max(drogueLoad.ubt_axial, mainLoad.ubt_axial);

% UBT Bending

load.ubt_bending = max(drogueLoad.ubt_bending, mainLoad.ubt_bending);

% LBT Axial

load.lbt_axial = max(drogueLoad.lbt_axial, mainLoad.lbt_axial);

% LBT Bending

load.lbt_bending = max(drogueLoad.lbt_bending, mainLoad.lbt_bending );

% ITS Axial

load.its_axial = max(drogueLoad.its_axial, mainLoad.its_axial);

%ITS bending

load.its_bending = max(drogueLoad.its_bending, mainLoad.its_bending);


% 3d plotter

end

function plot_slice(T, xcol, ycol, zcol, fixed_cols, fixed_vals, tol, ttl, xlab, ylab, zlab)
  mask = true(height(T),1);
  for k = 1:numel(fixed_cols)
    mask = mask & abs(T.(fixed_cols{k}) - fixed_vals(k)) < tol;
  end
  Tsub = T(mask,:);

  figure;
  if isempty(Tsub)
    title([ttl ' -- no matching rows']);
    return
  end
  scatter3(Tsub.(xcol), Tsub.(ycol), Tsub.(zcol), 25, Tsub.(zcol), 'filled');
  xlabel(xlab); ylabel(ylab); zlabel(zlab);
  title(ttl);
  colorbar; grid on; view(45,25);
end


%3d plotter, layered

function plot_layered_slice(T, xcol, ycol, zcol, fixed_col, fixed_val, tol, layer_col, layer_vals, ttl, xlab, ylab, zlab)
  mask = abs(T.(fixed_col) - fixed_val) < tol;
  Tsub = T(mask,:);

  figure; hold on;
  colors = lines(numel(layer_vals));
  legend_labels = strings(numel(layer_vals),1);
  for k = 1:numel(layer_vals)
    layer_mask = abs(Tsub.(layer_col) - layer_vals(k)) < tol;
    Tl = Tsub(layer_mask,:);
    if isempty(Tl)
      continue
    end
    scatter3(Tl.(xcol), Tl.(ycol), Tl.(zcol), 25, colors(k,:), 'filled');
    legend_labels(k) = sprintf('%s = %.2f', layer_col, layer_vals(k));
  end
  xlabel(xlab); ylabel(ylab); zlabel(zlab);
  title(ttl);
  legend(legend_labels(legend_labels ~= ""), 'Location','best');
  grid on; view(45,25); hold off;
end
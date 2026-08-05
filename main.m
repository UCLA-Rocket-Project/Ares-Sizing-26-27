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

%% Setup
% CSV path handling and writing

base_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(base_dir, 'out');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

data_csv_path = fullfile(out_dir, 'output.csv'); % every iter
filtered_csv_path = fullfile(out_dir, 'output_filtered.csv'); % for succesful iterations (meets all filters)
log_path = fullfile(out_dir, 'log.txt');

params = input_parameters();

col_names = {'prop_mass', 'OF', 'Pc', 'eps', 'mdot','thrust', 'Isp', 't_b','tank_press', 'V_He', 'dry_mass', ...
            'fuel_tank_length', 'ox_tank_length', 'vehicle_length', 'apogee', 'fail_code'};

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

% CURRENTLY PLACEHOLDERS !!!!
mass_dist = 85; %lbm
OF_dist = 1.3;
Pc_dist = [250, 300]; %psi
eps_dist = 3.5; 

it_ct = numel(mass_dist) * numel(OF_dist) * numel(Pc_dist) * numel(eps_dist);

%% Sweep

tic;
it = 0;

% CSV write for a small batch of rows in memory
% cleaner than writing every iteration
% less memory use than keeping every row and writing at the end
% changed to 4 nested loops (like last year) instead of running main per combo and appending 

flush_every = 50;
buffer_all={}; 
buffer_pass ={};

tol = 0.1;
max_iter = 50;

for prop_mass = mass_dist
    for OF = OF_dist
        for Pc = Pc_dist
            for eps = eps_dist
                
                it = it + 1;
                fail_code = 0 ;

                %defaults so that row doesnt return error if fail code occurs

                Prop = struct('OF', OF, 'Pc', Pc, 'eps', eps, 'prop_mass', prop_mass);
                Press = struct('tank_press', NaN, 'V_He', NaN);
                dry_mass = NaN;
                apogee = NaN;

                %% run_CEA

                Prop = run_CEA(Prop, params);

                %% run_press

                [Prop, Press] = run_press(Prop, params);

                %% get_PV_mel

                PV_mel = get_PV_mel(prop_mass, OF, Press.tank_press, Press.V_He);

                if isnumeric(PV_mel) && PV_mel == -1
                    fail_code = -1;
                    PV_mel = struct('fuel_l', NaN, 'ox_l', NaN);

                elseif isnumeric(PV_mel) && PV_mel == -2
                    fail_code = -2;
                    PV_mel = struct('fuel_l', NaN, 'ox_l', NaN);

                end

                %% get_ShockLoads & get_recLoads
                % Iterate through loads and dry mass until unchanging

                if fail_code == 0 
                dry_mass = 120; % Initializing dry mass
                old_dryMass = 0; % Set old dry mass to zero for the first iteration
                iter = 0;

               while abs(dry_mass - old_dryMass) > tol && iter < max_iter
                   old_dryMass = dry_mass; % Update old dry mass
  
                   % Calculate new dry mass based on current shock loads
                     [f_drogue, f_main] = get_ShockLoads(dry_mass);
                 recLoads = get_highestLoad(f_drogue, f_main, PV_mel);
                  dry_mass = get_dryMass(recLoads,PV_mel); 
                  iter = iter + 1;
               end

                    if iter == max_iter

                    fail_code = -3; 

                    end

                end

                %% get_apogee 
                if prop_mass == 85
                    data = readmatrix(fullfile(base_dir, 'MvsCd_data_85.csv'));
                    Cd_data = data(:,2);
                    M_data = data(:,1);
                elseif prop_mass == 90
                    data = readmatrix(fullfile(base_dir, 'MvsCd_data_90.csv'));
                    Cd_data = data(:,2);
                    M_data = data(:,1);
                elseif prop_mass == 95
                    data = readmatrix(fullfile(base_dir, 'MvsCd_data_95.csv'));
                    Cd_data = data(:,2);
                    M_data = data(:,1);
                else
                    data = readmatrix(fullfile(base_dir,'MvsCd_data_100.csv'));
                    Cd_data = data(:,2);
                    M_data = data(:,1); 
                end

                [M_data, uniq_idx] = unique(M_data);
                Cd_data = Cd_data(uniq_idx);
            

                % Calculate apogee and vehicle length
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


                %% Row assembly and write

                row = {prop_mass, OF, Pc, eps, Prop.mdot, Prop.Thrust, Prop.Isp, Prop.t_b, Press.tank_press, Press.V_He, dry_mass, PV_mel.fuel_l ...
                    PV_mel.ox_l, vehicle_length, apogee, fail_code};

                buffer_all(end+1, :) = row;
                % every row,  gets periodically written and cleared every iteration to output.csv

                if fail_code == 0
                    buffer_pass(end+1, :) = row; 
                    % rows that pass the filters, written to filtered.csv
                    
                end

                if mod(it, flush_every) ==0 || it == it_ct
                    if ~isempty(buffer_all)
                        writecell(buffer_all, data_csv_path, 'WriteMode','append');
                        buffer_all = {};
                    end
                    
                    if ~isempty(buffer_pass)
                        writecell(buffer_pass, filtered_csv_path, 'WriteMode', 'append');
                        buffer_pass = {};
                    end
                    fprintf(' %d / %d iterations complete\n,', it, it_ct);

                end
            end
        end
    end
end
elapsed_time = toc;
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
        'OF Ratio: %.2f\n' ...
        'Chamber Pressure: %.2f psi\n' ...
        'Tank Pressure: %.2f psi\n' ...
        'Isp: %.2f s\n\n' ...
        'Dry Mass: %.2f lb\n' ...
        'Fuel Barrel Length: %.2f in\n' ...
        'Ox Barrel Length: %.2f in\n' ...
        'Vehicle Length: %.2f in\n'], ...
        repmat('=',1,80), repmat('=',1,80), ...
        it, height(results_filtered), elapsed_time, ...
        max_apogee, ...
        results_filtered.prop_mass(opt_idx), results_filtered.thrust(opt_idx), ...
        results_filtered.t_b(opt_idx), results_filtered.mdot(opt_idx), ...
        results_filtered.OF(opt_idx), results_filtered.Pc(opt_idx), ...
        results_filtered.tank_press(opt_idx), results_filtered.Isp(opt_idx), ...
        results_filtered.dry_mass(opt_idx), results_filtered.fuel_tank_length(opt_idx), ...
        results_filtered.ox_tank_length(opt_idx), results_filtered.vehicle_length(opt_idx));

    fid = fopen(log_path, 'w');
    fprintf(fid, '%s', log);

    fclose(fid);

    disp(log);
end

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

end


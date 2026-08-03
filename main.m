% main

% Important notes: 
% need Cd and Mach data to read in, currently storing apogee and length in 
% struct called "results" for graphing at end, need to combine into overall 
% function that iterates through prop parameters

% Plotting results: graph of vehicle length and prop mass vs apogee

%% Code
PV_mel = get_PV_mel(prop_mass, OF, p_operating,press_v);
if PV_mel == -1
    apogee = -1;
    results.apogee.append(apogee);
    results.length.append(158.5);
    fprintf("COPV not available for calculated pressurant volume of %.0f\n", press_v);
    % add statement about what prop parameters cause this
    return
elseif PV_mel == -2
    apogee = -2;
    results.apogee.append(apogee);
    results.length.append(158.5);
    fprintf("Tanks unable to withstand operating pressure of %.0f\n",p_operating);
    return
end

% Iterate through loads and dry mass until unchanging
dry_mass = 120; % Initializing dry mass
old_dryMass = 0; % Set old dry mass to zero for the first iteration
while dry_mass > old_dryMass
    old_dryMass = dry_mass; % Update old dry mass
    % Calculate new dry mass based on current shock loads
    [f_drogue, f_main] = get_ShockLoads(dry_mass);
    recLoads = get_highestLoad(f_drogue, f_main, PV_mel);
    dry_mass = get_dryMass(recLoads,PV_mel); 
end

% Calculate apogee and vehicle length
Cd_data = readmatrix('Cd_data.csv'); % ** need data
M_data = readmatrix('M_data.csv'); % **need data
apogee = get_apogee(Prop, params, Cd_data, M_data, dry_mass); % ft
vehicle_length = 158.5 + PV_mel.fuel_l + PV_mel.ox_l; % in
% ^ 158.5 in is Pandora's length without tank barrels
if apogee == -3
    fprintf("Not meeting OTRS.\n");
end

% Store apogee and length in an array of results for graphing at end
results.apogee.append(apogee);
results.length.append(vehicle_length);
results.prop_mass.append(prop_mass);

%% Plotting
for i = 1:length(results.apogee)
    if results.apogee(i) == -1 || results.apogee(i) == -2 || results.apogee(i) == -3
        results.apogee(i) = [];
        results.length(i = []);
        results.prop_mass(i) = [];
    end
end
plot(results.length,results.apogee);
xlabel("Vehicle Length (in)");
ylabel("Apogee (ft)");
figure;
plot(results.prop_mass, results.apogee);
xlabel("Prop Mass (lb)");
ylabel("Apogee (ft)");

%% Local Functions

% Returns struct with highest equivalent axial load and highest bending
% moment from recovery on each of the 3 carbon tubes
function load = get_highestLoad(f_drogue, f_main,PV_mel)
drogueLoad = get_recLoads(f_drogue,PV_mel,"drogue");
mainLoad = get_recLoads(f_main,PV_mel,"main");

% UBT Axial
if drogueLoad.ubt_axial > mainLoad.ubt_axial
    load.ubt_axial = drogueLoad.ubt_axial;
else
    load.ubt_axial = mainLoad.ubt_axial;
end

% UBT Bending
if drogueLoad.ubt_bending > mainLoad.ubt_bending
    load.ubt_bending = drogueLoad.ubt_bending;
else
    load.ubt_bending = mainLoad.ubt_bending;
end

% LBT Axial
if drogueLoad.lbt_axial > mainLoad.lbt_axial
    load.lbt_axial = drogueLoad.lbt_axial;
else
    load.lbt_axial = mainLoad.lbt_axial;
end

% LBT Bending
if drogueLoad.lbt_bending > mainLoad.lbt_bending
    load.lbt_bending = drogueLoad.lbt_bending;
else
    load.lbt_bending = mainLoad.lbt_bending;
end

% ITS Axial
if drogueLoad.ubt_axial > mainLoad.ubt_axial
    load.ubt_axial = drogueLoad.ubt_axial;
else
    load.ubt_axial = mainLoad.ubt_axial;
end

% ITS Bending
if drogueLoad.its_bending > mainLoad.its_bending
    load.its_bending = drogueLoad.its_bending;
else
    load.its_bending = mainLoad.its_bending;
end

end
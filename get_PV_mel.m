% get_PV_mel

% Calculates masses and lengths of pressure vessels (tanks and COPV)
% Assuming 2:1 semi ellipsoidal end caps

% Inputs: Propellant mass, OF ratio, tank MEOP, press tank volume
% Outputs: Struct of PV lengths and masses or -1 (COPV not possible), 
% -2 (tanks unable to withstand pressure)

function PV_mel = get_PV_mel(prop_mass, OF, p_operating,press_v)

% First, check if a COPV is possible (return -1 if not)
copv_volumes = [9,12]; % L
copv_masses= [4.9,6.9].* 2.205; % lbm
copv_lengths = [560,595]/25.4; % in

if press_v > copv_volumes(length(copv_volumes))
    PV_mel = -1; % COPV not possible
    return
end

for i = 1:length(copv_volumes)
    if copv_volumes(i) > press_v
        PV_mel.copv_v = copv_volumes(i);
        PV_mel.copv_m = copv_masses(i);
        PV_mel.copv_l = copv_lengths(i);
        break
    end
end

% Proceed to check tank thickness based on MEOP (return -2 if not possible)

% Constants
S_y = 35000; % yield strength of Aluminum 6061-T6, psi
r_o = 4; % in
FS_y = 1.25; % FOS Yield for pressure vessels
% FF = 1.15; % Fitting factor, currently ignoring
% E = 1; % Weld efficiency factor (can ignore)

tank_thicknesses = [0.125, 0.25, 0.375]; % in, 
% Can be found on OnlineMetals, Speedy metals, Grainger, and BuyMetal
% ^^ need to research available tubes to make array

for t = 1:length(tank_thicknesses)
    r_i = r_o - tank_thicknesses(t);
    hoop_stress = 1.5*p_operating*r_i/tank_thicknesses(t);
    MOS_y = (S_y/(FS_y*hoop_stress))-1;

    if MOS_y > 0
        PV_mel.tank_wall = tank_thicknesses(t);
        break
    end

    if t == tank_thicknesses(length(tank_thicknesses))
        PV_mel = -2; % Tanks unable to withstand pressure
        return
    end
end

r_i = r_o - PV_mel.tank_wall;

% Calculate end cap thickness (crown and knuckle)
cap_t = 0.1375; % in
a = r_o - cap_t; % in
knuckle_hoop = 1.5*p_operating*a/cap_t; % psi
MOS_cap_y = S_y/(FS_y*knuckle_hoop)-1;

while MOS_cap_y <= 0
    cap_t = cap_t + 1/32; % in
    a = r_o - cap_t; % in
    knuckle_hoop = 1.5*p_operating*a/cap_t; % psi
    MOS_cap_y = S_y/(FS_y*knuckle_hoop)-1;
end
PV_mel.cap_wall = cap_t; % in

% Calculate end cap volume and mass
a = r_o - PV_mel.cap_wall; % in
cap_vol = (1/3*pi*a^3)+(pi*r_i^2*1.5); % in^3 (crown + tangent volume)
cap_mass = 0.0975 * (pi/24*((r_o*2)^3-(a*2)^3) + pi*r_o*(r_o^2-r_i^2)); % lbm
cap_mass = cap_mass * 1.3; % mark up factor based on pandora


% Calculate barrel length and mass
fuel_mass = prop_mass/(1+OF); % lbm
ox_mass = OF*fuel_mass; % lbm
fuel_density = 0.75*0.0285 + 0.25*0.036; % lb/in^3 for 75/25 ethanol/water
ox_density = 0.041; % lb/in^3
fuel_vol = fuel_mass/fuel_density; % in^3
ox_vol = ox_mass/ox_density; % in^3
ullage = 0.05;

fuel_length = (fuel_vol + ullage*fuel_vol - 2*cap_vol)/pi/(r_i^2); % in
ox_length = (ox_vol + ullage*ox_vol - 2*cap_vol)/pi/(r_i^2); % in

% Add results to struct
PV_mel.fuel_l = fuel_length;
PV_mel.fuel_m = 0.0975*pi*fuel_length*(r_o^2-r_i^2); % lbm
PV_mel.ox_l = ox_length;
PV_mel.ox_m = 0.0975*pi*ox_length*(r_o^2-r_i^2); % lbm
PV_mel.cap_m = cap_mass;
PV_mel.pv_m = PV_mel.fuel_m + PV_mel.ox_m + 4*PV_mel.cap_m + PV_mel.copv_m;

end

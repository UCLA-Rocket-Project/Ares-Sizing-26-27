% get_dryMass

% Calculates dry mass of the vehicle based on recovery loads

% Inputs: Struct with max recovery loads, PV_mel
% Output: Dry mass of the vehicle

function [dry_mass, tube_masses] = get_dryMass(recLoads,PV_mel)

% Determine UBT layer number
n = 6;
ubt_MOS_buckling = get_Buckling(recLoads.ubt_l,n,recLoads.ubt_axial);
ubt_MOS_bending = get_Bending(n,recLoads.ubt_bending);
while ubt_MOS_buckling <= 0 || ubt_MOS_bending <= 0
    n = n + 2;
    ubt_MOS_buckling = get_Buckling(recLoads.ubt_l,n,recLoads.ubt_axial);
    ubt_MOS_bending = get_Bending(n,recLoads.ubt_bending);
end
ubt_n = n;

% Determine LBT layer number
n = 6;
lbt_MOS_buckling = get_Buckling(recLoads.lbt_l,n,recLoads.lbt_axial);
lbt_MOS_bending = get_Bending(n,recLoads.lbt_bending);
while lbt_MOS_buckling <= 0 || lbt_MOS_bending <= 0
    n = n + 2;
    lbt_MOS_buckling = get_Buckling(recLoads.lbt_l,n,recLoads.lbt_axial);
    lbt_MOS_bending = get_Bending(n,recLoads.lbt_bending);
end
lbt_n = n;

% Determine ITS layer number
n = 6;
its_MOS_buckling = get_Buckling(recLoads.its_l,n,recLoads.its_axial);
its_MOS_bending = get_Bending(n,recLoads.its_bending);
its_MOS_bearing = get_Bearing(n,recLoads.its_axial);
while its_MOS_bearing <= 0 || its_MOS_buckling <= 0 || its_MOS_bending <= 0
    n = n + 1;
    its_MOS_buckling = get_Buckling(recLoads.its_l,n,recLoads.its_axial);
    its_MOS_bending = get_Bending(n,recLoads.its_bending);
    its_MOS_bearing = get_Bearing(n,recLoads.its_axial);
end
its_n = n;

% Calculate mass of each tube
ubt_m = get_tubeMass(ubt_n,recLoads.ubt_l);
lbt_m = get_tubeMass(lbt_n,recLoads.lbt_l);
its_m = get_tubeMass(its_n,recLoads.its_l);

tube_masses.ubt_m = ubt_m;
tube_masses.lbt_m = lbt_m;
tube_masses.its_m = its_m;

% Dry mass of Pandora minus all PV things and all carbon tubes = 116 lb
dry_mass = ubt_m + lbt_m + its_m + PV_mel.pv_m + 120; % lb

end

%% Local Functions

% Buckling (input tube length, # of layers, expected axial load)
function MOS_buckling = get_Buckling(l,n,P)
r_i = 4;
r_o = 4 + 0.0083*n;
I = pi/4*(r_o^4-r_i^4);
P_crit = 4320538.26*I/(l^2);
P_exp = P;
FS_u = 1.75;
MOS_buckling = P_crit/(FS_u*P_exp)-1;
end

% Bending (input # of layers, expected bending moment)
function MOS_bending = get_Bending(n, M)
r_i = 4; % in
t = n*0.0083; % in
r_o = r_i + t; % in
I = pi/4*(r_o^4-r_i^4);
E = 1.25*10^7; % Young's Modulus (psi)
u = 0.05; % Poisson's Ratio
phi = (1/16)*sqrt(r_o/t);
gamma_bending = 1-0.713*(1-exp(-phi));
bending_stress = (gamma_bending*E)/sqrt(3*(1-u^2))*(t/r_o); % psi
M_crit = bending_stress*I/r_o; % lbf-in
M_exp = M; % lbf-in
FS_u = 1.75;
MOS_bending = M_crit/(FS_u*M_exp)-1;
end

% ITS Bearing (input # of layers, axial load)
function MOS_bearing = get_Bearing(n,P)
S_br_CF = 72000; % Carbon fiber bearing strength (psi)
t = 0.0083*n;
d_bushing = 0.625; % in
n_bolts = 8;
bearing_stress = P/(n_bolts*d_bushing*t); % psi
FS_u = 1.75;
MOS_bearing = S_br_CF/(FS_u*bearing_stress)-1;
end

% Mass (lb) of a tube (input # of layers, tube length)
function mass = get_tubeMass(n,l)
density = 1.4/27.68; % lb/in^3
r_i = 4; % in
r_o = 4 + 0.0083*n; % in
V = pi*(r_o^2-r_i^2)*l; % in^3
mass = V*density; % lb
end
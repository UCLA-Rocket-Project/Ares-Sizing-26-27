% get_abl.m
% inputs: 
    % [injector, nozzle entrance] T0 stagnation temp, Pr Prandtl number, gamma, Cp, t_b burn time, Pc chamber pressure psia
% outputs: max char depth, mass of ablative


function [abl_mass, char_depth] = get_abl(Abl, Prop)

    % Calculate char depth
    % At injector
    abl_depth_inj = get_charDepth(Abl.T0(1), Abl.Pr(1), Abl.gamma(1), Abl.Cp(1), Prop.t_b, Prop.Pc);
    % At nozzle entrance
    abl_depth_noz = get_charDepth(Abl.T0(2), Abl.Pr(2), Abl.gamma(2), Abl.Cp(2), Prop.t_b, Prop.Pc);
    char_depth = max(abl_depth_inj, abl_depth_noz);

    % Calculate mass of ablative
    abl_mass = get_ablMass(char_depth, 1.25, 0.2); % 1.25 safety factor, 0.2 margin of safety

end

% char depth from Huzel and Huang
function char_depth = get_charDepth(T0, Pr, gamma, Cp, t_b, Pc)
    c = 1.606771881; % empirical constant
    k = 9.80e-6; % thermal conductivity of char, btu/s*in*F

    R_r = 0.5; % weight fraction of pyrolyzed resin content in ablative material
    R_v = 0.41; % weight fraction of pyrolyzed resin vs total resin R_r

    C_p = Cp * 0.0002388; % J/kg-K to btu/lbm-F
    rho = 0.04046257; % lbm/in^3

    r = Pr ^ (1/3); % recovery number
    mach = 0.0; % Mach number

    T_aw_factor = (gamma - 1) / 2 * mach * mach;
    T_aw = (T0 * 1.8) * ((1 + T_aw_factor * r) / (1 + T_aw_factor)); % adiabatic wall temp, R, just converts T0 K to R?

    T_d = 1200; % pyrolysis temp, R
    L_d = 1190.7; % latent heat of pyrolysis, btu/lbm
    
    char_depth = c * ((2 * k * t_b) / (rho * C_p * R_r * R_v) * log(1 + (R_r * R_v * C_p * (T_aw - T_d)) / L_d)) ^ 0.5 * (Pc / 100) ^ 0.4; % in
end

function abl_mass = get_ablMass(char_depth, fos, mos)
    rho = 0.04046257; % lbm/in^3
    R_ch_inner = 4.76; % in, inner radius of chamber, heritage
    abl_mass = 0.5 * pi * (R_ch_inner^2 - (R_ch_inner - (2 * char_depth * (fos + mos)))^2) * char_depth * rho; % lb
end
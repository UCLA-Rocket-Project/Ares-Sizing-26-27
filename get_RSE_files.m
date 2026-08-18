% get_RSE_files

function get_RSE_files(Prop, params, dry_mass, PV_mel, Cd_data, M_data, file_name, company, out_dir)

    %% Thrust curve 
    % (copied from get_apogee

    mdot = Prop.mdot;
    cstar = Prop.C_star;
    cstar_eff = params.Cstar_eff;
    eps = Prop.eps;
    Pc = Prop.Pc;
    Ctau_vac = Prop.Ctau_vac;
    ctau_eff = params.Ctau_eff;
    prop_mass = Prop.prop_mass;

    g = 32.174; % ft/s^2
    dt = 0.01; % s
    h0 = 3000; % ft
    A = pi*16/144; % ft^2
    m = dry_mass + prop_mass; % lb

    h = h0;
    v = 0;

    t_log = [];
    thrust_N_log = [];
    t = 0;

    while m > dry_mass

        [rho, temp, p] = get_air_properties(h);
        Pamb = p * 0.145038;
        rho = rho/16.018; % lb/ft^3
        v_s = sqrt(1.4*287*(temp+273.15))*3.281;
        M = v/v_s;
        Cd = interp1(M_data,Cd_data, M,'linear','extrap');

        Ctau_t = Ctau_vac - eps * (Pamb / Pc);
        T = mdot * cstar * cstar_eff * Ctau_t * ctau_eff / g; % lbf

        t_log(end+1) = t;
        thrust_N_log(end+1) = T * 4.44822; % lbf -> N

        m = max(dry_mass, m - mdot*dt);
        drag = 0.5 * Cd * rho * (v^2) * A;
        a = (T - m - drag) * g / m;
        v = v + a*dt;
        h = h + v*dt + 0.5*a*dt^2;
        t = t + dt;

    end

    thrust_lbf_log = thrust_N_log / 4.44822;

    %% Fuel / ox mass over time
    fuel_mass0 = prop_mass / (1 + Prop.OF); % lb
    ox_mass0   = prop_mass - fuel_mass0;    % lb

    fuel_mass_lb = max(0, fuel_mass0 - Prop.mdot_fuel .* t_log);
    ox_mass_lb   = max(0, ox_mass0 - Prop.mdot_ox .* t_log);

    fuel_mass_g = fuel_mass_lb * 453.592;
    ox_mass_g   = ox_mass_lb * 453.592;

    %% Fuel / ox CG over time
    r_i = 4 - PV_mel.tank_wall;
    area_in2 = pi * r_i^2;

    fuel_density_lb_in3 = params.fuel_density * 3.61273e-5;
    ox_density_lb_in3   = params.ox_density * 3.61273e-5;

    h_fuel_in = fuel_mass_lb ./ (fuel_density_lb_in3 * area_in2);
    h_ox_in   = ox_mass_lb ./ (ox_density_lb_in3 * area_in2);

    fuel_depth_mm = (PV_mel.fuel_l - 0.5 * h_fuel_in) * 25.4;
    ox_depth_mm   = (PV_mel.ox_l - 0.5 * h_ox_in) * 25.4;

    %% Header values 
    % copied from anthony's script 

    fuel_thrust_N = 0.00001;
    ox_thrust_N = 0.00001;

    dia_mm = 25.4 * 8;
    burn_time = max(t_log);

    fuel_code = "Fuel_" + file_name;
    ox_code = "Ox_" + file_name;

    fuel_header = createEngineHeader(mean(thrust_N_log), burn_time, fuel_code, dia_mm, ...
        fuel_mass_g(1), PV_mel.fuel_l*25.4, mean(thrust_N_log), fuel_mass_g(1), company);
    ox_header = createEngineHeader(mean(thrust_N_log), burn_time, ox_code, dia_mm, ...
        ox_mass_g(1), PV_mel.ox_l*25.4, mean(thrust_N_log), ox_mass_g(1), company);

    footer = [
        '      </data>\n' ...
        '    </engine>\n' ...
        '  </engine-list>\n' ...
        '</engine-database>\n'
    ];

    thrust_code = "Thrust_" + file_name;

    thrust_dia_mm = 127; % mm

    thrust_header = createEngineHeader(mean(thrust_N_log), burn_time, thrust_code, thrust_dia_mm, ...
        0, max(PV_mel.fuel_l, PV_mel.ox_l)*25.4, max(thrust_N_log), 0, company);

    %% Write RSE files
    thrustFile = fopen(fullfile(out_dir, thrust_code + ".RSE"), 'w');
    fprintf(thrustFile, thrust_header);
    for i = 1:length(t_log)
        fprintf(thrustFile, '<eng-data cg="0.00000" f="%.5f" m="0.00000" t="%.5f"/>\n', thrust_N_log(i), t_log(i));
    end
    fprintf(thrustFile, footer);
    fclose(thrustFile);

    fuelFile = fopen(fullfile(out_dir, fuel_code + ".RSE"), 'w');
    fprintf(fuelFile, fuel_header);
    for i = 1:length(t_log)
        fprintf(fuelFile, '<eng-data cg="%.5f" f="%.5f" m="%.5f" t="%.5f"/>\n', fuel_depth_mm(i), fuel_thrust_N, fuel_mass_g(i), t_log(i));
    end
    fprintf(fuelFile, footer);
    fclose(fuelFile);

    oxFile = fopen(fullfile(out_dir, ox_code + ".RSE"), 'w');
    fprintf(oxFile, ox_header);
    for i = 1:length(t_log)
        fprintf(oxFile, '<eng-data cg="%.5f" f="%.5f" m="%.5f" t="%.5f"/>\n', ox_depth_mm(i), ox_thrust_N, ox_mass_g(i), t_log(i));
    end
    fprintf(oxFile, footer);
    fclose(oxFile);

    %% Write CSVs
    writetable(table(t_log', fuel_depth_mm', fuel_mass_g', 'VariableNames', {'t_s','cg_mm','mass_g'}), ...
        fullfile(out_dir, fuel_code + "_cg.csv"));

    writetable(table(t_log', ox_depth_mm', ox_mass_g', 'VariableNames', {'t_s','cg_mm','mass_g'}), ...
        fullfile(out_dir, ox_code + "_cg.csv"));

end

function header = createEngineHeader(avgThrust, burnTime, code, dia, initWt, len, peakThrust, propWt, mfg)
    header = sprintf([
        '<engine-database>\n' ...
        '<engine-list>\n' ...
        '<engine FDiv="10" FFix="1" FStep="-1." Isp="195.96" Itot="39.78" Type="liquid" auto-calc-cg="0" auto-calc-mass="1"\n' ...
        '    avgThrust="%.5f" burn-time="%.3f" cgDiv="10" cgFix="1" cgStep="-1." code="%s" delays="0" dia="%.1f" exitDia="0." initWt="%.2f"\n' ...
        '    len="%.1f" mDiv="10" mFix="1" mStep="-1." massFrac="36.32" mfg="%s" peakThrust="%.5f" propWt="%.2f" tDiv="10" tFix="1"\n' ...
        '    tStep="-1." throatDia="0">\n' ...
        '      <data>\n'
    ], avgThrust, burnTime, code, dia, initWt, len, mfg, peakThrust, propWt);
end

function [rho,temp,p] = get_air_properties(h)
    h = h/3.281;
    if h < 11000
        T = 15.04 - 0.00649*h;
        p = 101.29 *((T+273.15)/288.08)^(5.256);
    elseif h < 25000
        T = -56.46;
        p = 22.65*exp(1.73-0.000157*h);
    else
        T = -131.21 + 0.00299*h;
        p = 2.488 * ((T+273.15)/216.6)^(-11.388);
    end
    rho = p/(0.2869*(T+273.15));
    temp = T;
end

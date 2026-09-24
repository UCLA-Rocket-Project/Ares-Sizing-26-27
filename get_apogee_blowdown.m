function apogee = get_apogee_blowdown(thrust_array, mdot_array, dry_mass, prop_mass, Cd_data, M_data)

    g = 32.174; % ft/s^2
    dt = 0.01; % s
    h0 = 2100; % ft
    v_oftr = 100; % ft/s
    A = pi*16/144; % ft^2
    m = dry_mass + prop_mass; % lb

    h = h0; % ft
    v = 0; % ft/s
    oftr_checked = false;

    n_burn = length(thrust_array);

    for i = 1:n_burn
        T = thrust_array(i) / 4.44822; % N -> lbf
        mdot_i = mdot_array(i) * 2.20462; % kg/s -> lbm/s

       [rho, temp, ~] = get_air_properties(h);
        rho = rho/16.018/32.174;
        v_s = sqrt(1.4*287*(temp+273.15))*3.281;
        M = v/v_s;
        Cd = interp1(M_data, Cd_data, M, 'linear', 'extrap');
     
        m = max(dry_mass, m - mdot_i*dt);
        drag = 0.5 * Cd * rho * (v^2) * A;
        Force = (T - m - drag);
        a = Force * g / m;
        v = v + a*dt;
        h = h + v*dt + 0.5*a*dt^2;

        if oftr_checked == false && h >= h0 + 48
            if v < v_oftr
                apogee = -3;
                return
            end
            oftr_checked = true;
        end
    end

    while v >= 0 % coasting, same as get_apogee
        [rho, temp] = get_air_properties(h);
        rho = (rho/16.018)/32.174;
        v_s = sqrt(1.4*287*(temp+273.15))*3.281;
        M = v/v_s;
        Cd = interp1(M_data, Cd_data, M, 'linear', 'extrap');

        a = (-dry_mass - 0.5*Cd*rho*v^2*A)*g/dry_mass;
        v = v + a*dt;
        h = h + v*dt + 0.5*a*dt^2;
    end
    apogee = h - h0;
end
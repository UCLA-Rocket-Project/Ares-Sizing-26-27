% Simplified Bending Moment Script for integration with sizing script
% Cyril's Script adapted for rocket sizing

% Inputs: Recovery load (lbf), PV mel struct, chute type (drogue or main)
% Outputs: Struct with max loads on UBT, LBT, ITS for layer sizing

function recLoads = get_recLoads(P, PV_mel,chute_type)

    if chute_type == "main"
        deploy_angle = atan(8/4);
        P_radial = P*cos(deploy_angle);
        P_axial = P*sin(deploy_angle);
    else
        P_radial = P;
        P_axial = 0;
    end
    
    %% Constants
    g = 32.174;     % Gravitational Acceleration [ft/s²]
    in2ft = 12;     % Inches per Foot [in]
    dL = 0.01;      % Tiny Division size [in]
    
    %% Define Beam Length and Mass Distribution
    L_tot = 158.5 + PV_mel.ox_l + PV_mel.fuel_l; % Pandora length + tanks, may want to add extra to length of pandora for safety
    
    % Divide the Rocket Into Several Tiny Pieces of Length dL
    N_divs = round(L_tot / dL);            % Total Divisions of Rocket
    Locs = dL / 2 : dL : L_tot - dL / 2;   % Location of Each Division
    
    % Define Mass Components in lbm (nothing above recovery coupler)
    % Taken from Pandora's MEL w/ PV commented out
    masses = zeros(1,N_divs);
    
    masses = defMass(masses, 39.76, 3, dL, 2.5875); %Main Body Electronics
    masses = defMass(masses, 39.76, 4, dL, 4.725); %Recovery Electronics + Mechanism
    masses = defMass(masses, 26.61, 8.5, dL, 4.0845); %Recovery Tube
    masses = defMass(masses, 27.21, 10, dL, 6.21); %Recovery Parachutes & Shock Cords
    masses = defMass(masses, 31.11, 8, dL, 2.541); %Recovery Coupler + Fasteners
    masses = defMass(masses, 35.11, 48.16, dL, 5.9184); %Upper Body Tube
    masses = defMass(masses, 39.11, 0.65, dL, 1.449); %Lower Recovery Bulkhead + Pin + Fasteners
    masses = defMass(masses, 45.61, 22.44, dL, 1.15); %Pressurant Tank Centering Rings
    % masses = defMass(masses, 45.61, 22.44, dL, 15.96); %Pressurant Tank
    masses = defMass(masses, 68.05, 10.31, dL, 8.8725); %Pressurant Bay Plumbing
    masses = defMass(masses, 74, 75, dL, 1.15); %Fairing
    masses = defMass(masses, 74, 75, dL, 1.035); %Raceway Plumbing
    % masses = defMass(masses, 79.27, 4, dL, 2.814); %Hemi + Skirt (Upper LOX) + Nut Holders
    % masses = defMass(masses, 83.27, 22.3, dL, 7.35); %LOx Tank (No Hemi)
    % masses = defMass(masses, 83.27, 15.25, dL, 10.395); %ITS Joint + Hemis + Nut Holders
    % masses = defMass(masses, 120.82, 25.2, dL, 8.19); %Fuel Tank (No Hemi)
    % masses = defMass(masses, 146.02, 4, dL, 2.9085); %Fuel Skirt + Hemi (Singular) + Nut Holders
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l, 37, dL, 4.053); %Lower Body Tube
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+4, 23, dL, 15.12); %Engine Bay Plumbing
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+4+23+4, 0.5, dL, 0.42); %Thrust Bulkhead + Fasteners
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+4+23+4.5, 0.5, dL, 2.9925); %Injector
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+4+23+5, 18.75, dL, 22.386); %Engine
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+12, 0.3, dL, 0.525); %Engine Centering Ring
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+29, 20, dL, 12.6); %Fins
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l+37, 12, dL, 1.9845); %Boattail
    
    % Pressure Vessels
    masses = defMass(masses, 68.05-PV_mel.copv_l, PV_mel.copv_l, dL, PV_mel.copv_m); %Pressurant Tank
    masses = defMass(masses, 79.27, 4, dL, PV_mel.cap_m); % Upper Lox End Cap
    masses = defMass(masses, 83.27, PV_mel.ox_l, dL, PV_mel.ox_m); % Lox Tank
    % note that end cap mass is constant across trials
    masses = defMass(masses, 83.27+PV_mel.ox_l, 15.25, dL, 10.395); % ITS + Lower Lox and Upper Fuel End Caps
    masses = defMass(masses, 83.27+PV_mel.ox_l+15.25, PV_mel.fuel_l, dL, PV_mel.fuel_m); % Fuel Tank
    masses = defMass(masses, 83.27+15.25+PV_mel.ox_l+PV_mel.fuel_l, 4, dL, PV_mel.cap_m); % Lower Fuel End Cap
    
    %% Apply Loads
    % Define Applied Forces
    ExternalForces = zeros(1,N_divs);   % Initialize External Force Vector [lbf]
    x = 31.1; % lip of recovery coupler, stays constant
    loc_app = find(Locs > x, 1, 'first'); 
    ExternalForces(loc_app) = ExternalForces(loc_app) + P_radial;
    
    %% Calculate CG and MoI
    % Mass and MoI at Each Location
    MoIs = zeros(1,N_divs);          % Moment of inertia of each section about its center [lbm*in²]
    
    M_tot = sum(masses);
    CG_numerator = dot(Locs, masses);
    CG_loc = CG_numerator / M_tot;              % Calculate CG location
    MoIs_axis = masses .* (CG_loc - Locs) .^2;  % Calculate MoI of each section about the rocket CG
    MoI_tot = sum(MoIs_axis);                   % Sum individual MoIs to get the total rocket MoI about its CG
    
    %% Compute Accelerations
    % Calculate Sum of Forces and Moments
    F_tot = sum(ExternalForces);                        % Total external force [lbf]
    Mom_tot = sum(ExternalForces .* (CG_loc - Locs));   % Total moment [lbf*in]
    
    % Calculate Linear and Angular Acceleration
    a_lin = F_tot / (M_tot / g);
    a_ang = (Mom_tot / in2ft) / (MoI_tot / (g * in2ft ^2));
    
    %% Calculate Inertial Loading
    % acceleration at each location
    Accels = a_lin + a_ang * (CG_loc - Locs) / in2ft;   % Acceleration of each tiny section [ft/s²]
    
    InertialForces = Accels .* (masses / g);            % Inertial forces on each tiny section [lbf]
    
    %% Calculate Shear Force at Each Station
    DeltaShear = ExternalForces - InertialForces;       % Shear added at each tiny section [lbf]
    
    % Numerically Integrate DeltaShear to Get Shear(x)
    Shear = zeros(1,N_divs);
    Shear(1) = DeltaShear(1);
    for i = 2 : N_divs
        Shear(i) = Shear(i - 1) + DeltaShear(i);
    end
    
    %% Calculate Bending Moment at Each Station
    % Moment is the Integration of Shear Over the Length
    DeltaMoment = Shear * (dL / in2ft);
    
    Moment = zeros(1,N_divs);
    Moment(1) = DeltaMoment(1);
    for i = 2 : N_divs
        Moment(i) = Moment(i - 1) + DeltaMoment(i);
    end
    
    
    %% Calculate P_eq at Each Station
    r = 4;                    % Rocket Radius [in]
    P_eq = abs(2 * Moment / (r / in2ft));
    
    if abs(P_eq) >= 0.0001
        P_eq = P_eq + P_axial;
    end
    
    %% Results
    
    % Define component boundaries
    ubt_idx = round(35.11/dL):round((35.11+51.25)/dL);

    its_start_x = 83.27 + PV_mel.ox_l;
    its_end_x = its_start_x + 15.25;
    its_idx = round(its_start_x/dL):round(its_end_x/dL);

    lbt_start_x = its_end_x + PV_mel.fuel_l;
    lbt_end_x = lbt_start_x + 37;
    lbt_idx = round(lbt_start_x/dL):round(lbt_end_x/dL);

    % UBT Max Load
    recLoads.ubt_bending = max(abs(Moment(ubt_idx)*12));
    recLoads.ubt_axial = max(abs(P_eq(ubt_idx)));
    recLoads.ubt_l = 51.25;

    % ITS Max Load
    recLoads.its_bending = max(abs(Moment(its_idx)*12));
    recLoads.its_axial = max(abs(P_eq(its_idx)));
    recLoads.its_l = 15;
    
    % LBT Max Load
    recLoads.lbt_bending = max(abs(Moment(lbt_idx)*12));
    recLoads.lbt_axial = max(abs(P_eq(lbt_idx)));
    recLoads.lbt_l = 37;

end

% Calculates new vehicle mass distribution when a component is added
function Mass_new = defMass(Mass_cur, x_i, L_comp, dL, m1)
    Mass_new = Mass_cur;
    N_divs = round(L_comp / dL);
    start_index = round(x_i / dL) + 1;
    end_index = start_index + N_divs - 1;
    mass = m1;
    m_each = mass / N_divs;
    Mass_new(start_index:end_index) = Mass_new(start_index:end_index) + m_each;
end
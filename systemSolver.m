% systemSolver
% taken from Michael's engine balance script
% iterates for Pc 

% solveChamberPressure 
% engineSim


function [thrust, mdot_total, of_guess, P_chamber_guess, T_adiabatic, message] = systemSolver( ...
    rho_fuel, rho_ox, P_tank_fuel, P_tank_ox, CdA_fuel, CdA_ox, P_amb, ...
    A_throat, cstar_eff, ctau_eff, A_exit, CEA_obj, relax, chamber_guess_scale, tolerance, max_iterations)

% defining iteration 
    if nargin < 13, relax = 0.1; end
    if nargin < 14, chamber_guess_scale = 0.5; end
    if nargin < 15, tolerance = 1e-3; end
    if nargin < 16, max_iterations = 1000; end

    iteration = 0;

    % initial Pc guess

    P_chamber_guess = chamber_guess_scale * P_tank_ox;

    % mdot using Pc guess dP & CdA equation

    mdot_ox_guess = massFlowRate(CdA_ox, rho_ox, max(P_tank_ox - P_chamber_guess, 0.0));
    mdot_fuel_guess = massFlowRate(CdA_fuel, rho_fuel, max(P_tank_fuel - P_chamber_guess, 0.0));
    of_guess = mdot_ox_guess / mdot_fuel_guess;

    while iteration < max_iterations
        
        % run CEA with guessed chamber pressure
        % get total mdot w/ guessed Pc using A_t * Pc / C*

        [mdot_total, thrust, T_adiabatic] = engineSim(P_chamber_guess, P_amb, of_guess, A_throat, A_exit, cstar_eff, ctau_eff, CEA_obj);

        % solve for Pc needed to supply the mdot_total found in engineSim
        % uses CdA equation

        P_chamber_guess_2 = solveChamberPressure(mdot_total, CdA_ox, CdA_fuel, rho_ox, rho_fuel, P_tank_ox, P_tank_fuel, P_amb);

        if abs(P_chamber_guess_2 - P_chamber_guess) / P_chamber_guess < tolerance
            break
        end

        % iteration relaxation

        P_chamber_guess = (1 - relax) * P_chamber_guess + relax * P_chamber_guess_2;
        mdot_ox_guess = massFlowRate(CdA_ox, rho_ox, max(P_tank_ox - P_chamber_guess, 0.0));
        mdot_fuel_guess = massFlowRate(CdA_fuel, rho_fuel, max(P_tank_fuel - P_chamber_guess, 0.0));
        of_guess = mdot_ox_guess / mdot_fuel_guess;

        iteration = iteration + 1;

    end

    message = sprintf([ ...
        '\n    Results from systemSolver:\n\n' ...
        '    Thrust: %.2f lbf\n' ...
        '    Total Mass Flow Rate: %.4f kg/s\n' ...
        '    Oxidizer-to-Fuel Ratio (O/F): %.2f\n' ...
        '    Chamber Pressure: %.2f PSIA\n' ...
        '    Iterations: %d\n'], ...
        thrust * 0.224809, mdot_total, of_guess, P_chamber_guess * 0.000145038, iteration);

end


%% Local helper functions

function mdot = massFlowRate(CdA, rho, dP)
    % mdoot

    mdot = CdA * sqrt(2*rho*dP);
end

function P_root = solveChamberPressure(mdot_total, CdA_ox, CdA_fuel, rho_ox, rho_fuel, P_tank_ox, P_tank_fuel, P_min)
    % Returns the estimated chamber pressure (Pa) given total mdot 

    tol = 1e-3;

    % equation for difference btw mdot_total and mdot from ox & fuel side CdA equations
    flow_error = @(P) (CdA_ox * sqrt(2.0 * rho_ox * max(P_tank_ox - P, 0.0)) ...
                      + CdA_fuel * sqrt(2.0 * rho_fuel * max(P_tank_fuel - P, 0.0))) - mdot_total;

    % Pc search range

    P_low = max(P_min, 0.0);
    P_high = min(P_tank_ox, P_tank_fuel) - 1.0; % leave 1 Pa margin

    f_low = flow_error(P_low);
    f_high = flow_error(P_high);

    if f_low < 0
        error('Target mdot exceeds what the tanks can supply - Check inputs. mdot_total: %g kg/s', mdot_total);
    end
    if f_high > 0
        error('Target mdot exceeds what the tanks can supply - Check inputs. mdot_total: %g kg/s', mdot_total);
    end

    % Brent's method 
    % find the Pc thats needed for flow error = 0 (mdot from CdA equation to match mdot from Pc * At / C* equation)

    P_root = fzero(flow_error, [P_low, P_high]);

end

function [mdot, thrust, T_adiabatic] = engineSim(P_chamber, P_ambient, of_ratio, A_throat, A_exit, cstar_eff, ctau_eff, CEA_obj)
    % runs CEA to get C* 
    % finds total mdot using Pc * At / C* equation

    P_chamber_psia = P_chamber * 0.000145038; % Pa -> psia
    P_ambient_psia = P_ambient * 0.000145038; % Pa -> psia
    eps = A_exit / A_throat;

    % get_Cstar(Pc, MR) -> ft/sec
    cstar_ftps = double(CEA_obj.get_Cstar(pyargs('Pc', P_chamber_psia, 'MR', of_ratio)));
    cstar = cstar_ftps * 0.3048 * cstar_eff; % ft/s -> m/s

    % get_PambCf(Pamb, Pc, MR, eps) -> python tuple, index [1] in Python
    % (0-based) is the ambient-corrected Cf -- MATLAB py.tuple indexing
    % via {} is 1-based, so Python's [1] is MATLAB's {2}.
    cf_result = CEA_obj.get_PambCf(pyargs('Pamb', P_ambient_psia, 'Pc', P_chamber_psia, 'MR', of_ratio, 'eps', eps));
    ctau = double(cf_result{2}) * ctau_eff;

    % get_Temperatures(Pc, MR, eps) -> python tuple, index [0] in Python
    % is chamber/combustion temperature -- MATLAB py.tuple index {1}.
    temps_result = CEA_obj.get_Temperatures(pyargs('Pc', P_chamber_psia, 'MR', of_ratio, 'eps', eps));
    T_adiabatic_R = double(temps_result{1}); % Rankine
    T_adiabatic = T_adiabatic_R * (5/9); % Rankine -> Kelvin

    mdot = (A_throat * P_chamber) / cstar; % kg/s
    thrust = ctau * P_chamber * A_throat; % N

end

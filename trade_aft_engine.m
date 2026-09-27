
%  Boattail Aft and Engine Exit Diameter Optimization

% To use, change A_e in input_parameters to the engine exit diameter you
% are trading for. This script runs main.m to get the optimal prop
% parameters, and then rechecks apogee for each boattail configuration and
% returns the highest apogee found. 

% User Inputs
fname = {"input\v6\4.75_15.CSV", "input\v6\5.0_15.CSV","input\v6\5.25_15.CSV","input\v6\5.50_15.CSV","input\v6\5.75_15.CSV","input\v6\6.0_15.CSV","input\v6\6.25_15.CSV","input\v6\6.50_15.CSV","input\v6\6.75_15.CSV","input\v6\7.0_15.CSV"};   % RASAero exports
D_aft = {4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7}; % Corresponding to RASAero file configurations
data_range = 500; % Data range for AoF of 0 up to Mach 5


apogees = zeros(1,length(fname));
evalc('main');
D_e = round(2*sqrt(params.A_e/pi),2); % Engine exit diameter [in]


for i = 1:length(fname)
    if D_e + 1 <= D_aft{i}
        data = readmatrix(fname{i});
        testM_data = data(1:data_range,1);
        testCd_data = data(1:data_range,3);
        apogees(i) = get_apogee(Prop, params, testCd_data, testM_data, dry_mass);
    else
        apogees(i) = 0;
        % fprintf("Aft diameter (%4.2f in) is too small for the engine exit diameter (%4.2f in)\n", D_aft{i},D_e);
    end
end

[max_apo, max_index] = max(apogees);
fprintf("Aft diameter of %4.2f in is most optimal for a %4.2f in engine exit diameter. Apogee: %5.0f ft\n",D_aft{max_index}, D_e, max_apo);

% run_CEA

% runs NASA CEA for a given set of variable parameters

% Inputs: Pc, OF, Pamb, A_e, eps (Ae/At), Propellant Mass
% Outputs: C*, C_tau, A_t, mdot, Thrust, Isp, Burn Time

function Prop = run_CEA(Prop, params)

    Engine = CEA('problem','rocket','equilibrium','fac','ac/at',params.ac_at,'o/f',Prop.OF,'p(psi)',Pc_psi,'supersonic',Prop.eps,'reactants', ...
             'fuel','C2H5OH(L)','wt%',params.eth_ratio*100,'t(k)',params.T_fuel_inlet,'fuel','H2O(L)','wt%',(1-params.eth_ratio)*100,'t(k)', ...
            params.T_fuel_inlet,'oxid','O2(L)','wt%',100,'t(k)',params.T_ox_inlet,'output','transport','mks','end');

    % Calculate C_star

    Prop.C_star = Engine.output.eql.cstar(1);

     % Calculate C_tau

    Prop.Ctau = Engine.output.eql.cf(end);

    % Calculate Exit Pressure

    Prop.Pe = Engine.output.eql.pressure(end) * 14.5038; 

    % Calculate C_tau at Vaccuum

    Prop.Ctau_vac = Prop.Ctau + Prop.eps * (Prop.Pe / Pc);

    % Calculate throat area

    Prop.A_t = params.A_e / Prop.eps;

    % Calculate mdot

    Prop.mdot = (Prop.A_t * Prop.Pc) / (Prop.C_star * params.Cstar_eff);

    % Calculate burn time

    Prop.t_b = params.prop_mass / Prop.mdot;

    % Calculate Ctau at set ambient pressure

    Ctau_ref = Prop.Ctau_vac - Prop.eps * (params.P_amb / Pc);

    % Calculate thrust for Isp

    Thrust_ref = Prop.mdot * Prop.C_star * params.Cstar_eff * Ctau_ref * params.Ctau_eff;

    % Calculate specific impulse

    Prop.Isp = Thrust_ref / (Prop.mdot * 9.81);

end
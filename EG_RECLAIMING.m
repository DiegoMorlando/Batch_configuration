clear all
close all
clc
global p  
%% INPUT
p=[];
%  %distribution inlet [MDEA H2O NAOH HCl Formic acid AceticAcid H2SO4 EG]
p.wp_in = [0.0388 0.2508 0.0020 0.0020 0.0042 0.0040 0.0040 0.6942];
% p.wp_in = [0.4245 0.5664 0 0.0002 0.003 0.0029 0.0003+0.0029 0]; %WHEN THEY GET RID OF EG BUT CONDITIONS ARE NOT SPECIFIED
p.feed = 19194.48716; %Kg/h
p.Tin = 120.2+273.15;     %Temperature [K]
p.P_tot = 0.3*1.0135*10^2; %Pressure [KPa]
%density of the solution Kg/m^3
p.rho=998;
p.V = 1; %m^3

% [MDEA H2O NAOH HCl Formicacid AceticAcid H2SO4 EG ] Molecular weight [Kg/Kmol] 
p.PM = [119.163 18 40 36.458 46.03 60.052 98.079 62.0678 ]; %for salt a medium among sodiumformate and sodiumacetate
p.c = [286.741 74 74 29.1822 35.07215 89.92099 98.2132 180 ]; %KJ/Kmol K
%inlet mass
p.m_in = p.wp_in*p.feed;

% Kinetics parameter: Neutralization
p.k_free=1;                            %GUESS [m^3/Kmol h]
% Stoic coefficient  [MDEA H2O NAOH HCL formicacid aceticacid H2SO4 EG ] Neutralization
p.vi_free_hcl = [+1 +1 -1 -1 0 0 0 0 ];   
p.vi_free_formicacid = [+1 +1 -1 0 -1 0 0 0];
p.vi_free_aceticacid = [+1 +1 -1 0 0 -1 0 0];
p.vi_free_h2so4 = [+1 +1 -1 0 0 0 -1 0];
%Outlet flow
p.feed_mol = p.feed*p.wp_in./p.PM;
p.feed_mol = sum(p.feed_mol);
%defining mole and evaluating inlet molar fraction
p.n_in = p.m_in./p.PM;
p.ml_in = p.n_in/p.feed_mol;

%Heat of reaction and heat of vaporization (i use the value as first approximation from the master
%thesis i found need to check on deltah_r)  [MDEA H2O NAOH HCL FORMICACID ACETICACID SOLFORICACID EG]
p.lambda = [73*1000 44000 0 0 0 0 0 65.6*1000]; % KJ/Kmol 
p.deltah_r  = 0; % don't take that into account
% WASTE
p.V_OVER_FEED = 0.8737; 
p.waste = p.feed_mol*(1-p.V_OVER_FEED);

%% Dynamic simulation 
%giving the time span and initial condition
t_span = linspace(0,5,10^3); %[h]
%[MDEA H2O NAOH HCL FORMICACID ACETICACID H2SO4 EG] 
x0 = [p.wp_in]; 
%solving the system of dae
option_dyn =odeset('RelTol',10^-10,'AbsTol',10^-10*ones(1,length(x0)));
[t, sol_dyn] = ode15s(@(t,x) ds_reclaiming_eg_mass(t,x,p.P_tot,p.waste,p.V),t_span,x0,option_dyn);
%getting other result
for i = 1:length(t_span)
[der,Tss,Heat_duty(i),y_mass,Vapor_flow_rate] = ds_reclaiming_eg_mass(t,sol_dyn(i,:)',p.P_tot,p.waste,p.V);
end
wf_ss=sol_dyn(end,:);
DMEA_recovery = Vapor_flow_rate*y_mass(1)/(p.m_in(1));
disp('The equilibrium temperature is = ' + string(Tss-273.15) + '°C')
disp('The necessary heat duty is equal to  = ' + string(Heat_duty(end)/(3600*1000)) + ' MW')
disp('DMEA recovery is = ' + string(DMEA_recovery*100) + '%')

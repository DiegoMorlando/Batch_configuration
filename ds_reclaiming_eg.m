%function for dynamic simulation
function [f_dyn,Heat_duty,wf] = ds_reclaiming_eg(t,x,P_tot,waste,V)
global p     
f_dyn=zeros(length(p.PM)+1,1); %EACH COMPONENTS PLUS TEMPERATURE
% Defining variable in terms of molar fraction for each components [MDEA H2O NAOH ACIDS EG]
mf_mdea=x(1);
mf_h2o=x(2);
mf_naoh=x(3);
mf_acid = x(4);
mf_eg = x(5);
T = x(6); %K 
Mtot = V*p.rho;
mf = [mf_mdea mf_h2o mf_naoh mf_acid mf_eg];
%Evaluation of molar content in the system 
N = Mtot*1/sum(p.PM.*mf);                                                       
%Vapour-Liquid Equilibrium
vpmdea = vp_mdea(T);        %KPa
vph2o =vp_h2o(T);          %Kpa
vpeg = vp_eg(T);
psat_mdea=vpmdea.*(mf(1));                     %saturation pressure of MEA 
psat_water=vph2o.*(mf(2));                     %saturation pressure of H2O 
psat_eg = vpeg.*mf(5);                         %Saturation pressure of EG
%molar fraction in the gas phase
y_mdea = psat_mdea/P_tot;
y_h2o = psat_water/P_tot;
y_eg = psat_eg/P_tot;

Vapor_flow_rate = p.feed_mol - waste;
%Evaluating the vapor rate [MDEA H2O NAOH ACIDS EG]
y = [y_mdea y_h2o 0 0 y_eg];

%concentration conversion for kinetics [Kmol/m^3 h]

c = (N*mf)/V; 
c_naoh = c(3);
c_acids = c(4);
r_free=p.k_free*c_naoh*c_acids;                              %Kmol/m^3 h     

%I USE PTOT = PMEA+PH2O+PCO2 
f_dyn(end) = P_tot - psat_mdea -psat_water -psat_eg;
%defining material balance

for i = 1:length(mf)
     f_dyn(i) = (p.n_in(i)-Vapor_flow_rate*y(i)+ r_free*V*p.vi_free(i)-p.waste*mf(i))/N;
end
%Evaluation of the Heat Duty
Q_in = sum(p.n_in.*p.c.*(p.Tin)));
Q_waste = sum(waste.*p.c.*mf);
Q_ev = sum(Vapor_flow_rate*p.lambda.*y);
Heat_duty = -Q_in + Q_ev + Q_waste;
%Massic fraction
M_i = N*mf./p.PM; %mass of each components
wf = M_i/sum(M_i); %Mass fraction
end                                                                                                                                             






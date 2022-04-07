function [f,T,Vapor_flow_rate,y_mass_mea,y_mass_h2o,y_mass_co2,viscosity,Mea_recovery_rate,r_deg,feed_water] = batch_configuration_function(t,x,P_tot,Qin,beta)
global p  
f = zeros(length(x),1);
% Defining variable in terms of molar fraction for each components
wf_mea=x(1);
wf_h2o=x(2);
wf_hef=x(3);
wf_naoh=x(4);
wf_hss=x(5);
wf_hepo=x(6);
wf_co2=x(7);
wf_nasalt =x(8);
V = x(9);
Mtot = V*p.rho;
wf = [wf_mea wf_h2o wf_hef wf_naoh wf_hss wf_hepo wf_co2 wf_nasalt];
N_c = Mtot*wf./(p.PM);
N = sum(N_c);
mf = N_c/N;
alfa_loading = mf(7)/(mf(1));                                         %nCO2/nMEA 
%concentration conversion for kinetics [Kmol/m^3 h]
c = (N_c)/V; 
c_mea =c(1);
c_co2 = c(7);
c_hef=c(3);
c_naoh=c(4); 
c_salt = sum(c)-c_mea-c_co2-c(2);
% PTOT = PMEA+PH2O+PCO2 
opts = optimset("Display","Off");
T = fsolve(@(T) P_tot-vp_mea(T)*mf(1)-vp_h2o(T)*mf(2)-psat_co2(T,alfa_loading),370,opts);
%molar fraction in the gas phase
y_mea = vp_mea(T)*mf(1)/P_tot;
y_h2o = vp_h2o(T)*mf(2)/P_tot;
y_co2 = psat_co2(T,alfa_loading)/P_tot;

%Evaluating the vapor rate [MEA H2O HEF NAOH DEGP HEPO CO2 Nasalt]
y = [y_mea y_h2o 0 0 0 0 y_co2 0];
y_mass = y.*p.PM/(sum(p.PM.*y));
y_mass_mea=y_mass(1);
y_mass_h2o = y_mass(2);
y_mass_co2 = y_mass(7);
%viscosity 
viscosity = viscosity_Weiland_Nielsen(wf_mea,30+273.15,alfa_loading,c_salt,beta);
toll=5/100;
% Take the time when the viscosity reach the value that i don't want to
% reach
if viscosity < p.viscosity_upper_limit*(1+toll) && viscosity > p.viscosity_upper_limit*(1-toll)
    p.t_start = [p.t_start,t];
elseif isempty(p.t_start)
    p.t_start = 100;
end
%Eliminate the first element which is equal to 100 when i find a second one 
if length(p.t_start) > 1 && p.t_start(1) == 100
    p.t_start(1) = [];
end
% when we reach the viscosity add water for t_flux_water
if t < p.t_start(1)+p.t_flux_water && t > p.t_start(1)
   feed_water = p.feed_water;
else
   feed_water = 0;
end
%ENERGY BALANCE
Q_water_in = feed_water*4.186*(p.T_in_water);
Vapor_flow_rate_mol = (Qin+Q_water_in)/(sum(p.lambda.*y));
Vapor_flow_rate = sum(Vapor_flow_rate_mol*y.*p.PM);

%Declaring GEN/CONS term
r_free=p.k_free*c_naoh*c_hef;                              %Kmol/m^3 h     
r_deg1=p.k_ref*exp(-p.Ea/p.R*(1/T-1/p.T_ref))*c_mea*c_co2;
r_deg2=p.k_ref2*c_mea*c_co2;
r_deg = r_deg1+r_deg2;
%VOLUME dVdt = - V/rho
f(end) = (feed_water-Vapor_flow_rate)/p.rho;
%defining material balance
for i = 1:length(mf)
     f(i) = (feed_water*p.wf_feed_water(i)-Vapor_flow_rate*y_mass(i)+ r_free*V*p.vi_free(i)*p.PM(i)+p.vi_deg1(i)*r_deg1*V*p.PM(i)+r_deg2*p.vi_deg2(i)*V*p.PM(i)-p.rho*wf(i)*f(end))/(V*p.rho);
end
%Evaluation of MEA Recovery
Mea_recovery_rate= Vapor_flow_rate_mol*y_mea/(p.V*p.rho*p.wp_in(1)/p.PM(1));
                  
end






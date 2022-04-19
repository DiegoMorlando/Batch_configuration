%function for dynamic simulation
function [f_dyn,T,Q_used,y_mass,Vapor_flow_rate] = ds_reclaiming_eg_mass(t,x,P_tot,waste,V)
global p     
f_dyn=zeros(length(p.PM),1); %EACH COMPONENTS PLUS TEMPERATURE
% Defining variable in terms of molar fraction for each components
wf_mdea=x(1);
wf_h2o=x(2);
wf_naoh=x(3);
wf_hcl = x(4);
wf_formicacid = x(5);
wf_aceticacid = x(6);
wf_h2so4 = x(7);
wf_eg = x(8);
Mtot = V*p.rho;
wf = [wf_mdea wf_h2o wf_naoh wf_hcl wf_formicacid wf_aceticacid wf_h2so4 wf_eg ];
N_c = Mtot*wf./(p.PM);
N_tot = sum(N_c);
mf = N_c/N_tot;
waste_mass = sum(waste*mf.*p.PM);
Vapor_flow_rate = p.feed-waste_mass;

%concentration conversion for kinetics [Kmol/m^3 h]

c = (N_c)/V; 
c_naoh = c(3);
c_hcl = c(4);
c_formicacid =c(5);
c_aceticacid = c(6);
c_h2so4 = c(7);
r_free_hcl=p.k_free*c_naoh*c_hcl;    %Declaring GEN/CONS term
r_free_formicacid = p.k_free*c_naoh*c_formicacid;
r_free_aceticacid = p.k_free*c_naoh*c_aceticacid;
r_free_h2so4 = p.k_free*c_naoh*c_h2so4;
%I USE PTOT = PMDEA+PH2O+PEG 
opts = optimoptions('fsolve','Display','off');
T = fsolve(@(T) P_tot-vp_mdea(T)*mf(1)-vp_h2o(T)*mf(2)-vp_eg(T)*mf(8),400,opts);
y_mdea = vp_mdea(T)*mf(1)/p.P_tot;
y_h2o  = vp_h2o(T)*mf(2)/p.P_tot;
y_eg   = vp_eg(T)*mf(8)/p.P_tot;
%Evaluating the vapor rate [MDEA H2O NAOH HCL FORMICACID ACETICACID H2SO4 EG ]
y = [y_mdea y_h2o 0 0 0 0 0 y_eg ];
y_mass = y.*p.PM/(sum(p.PM.*y));
%defining material balance

for i = 1:length(mf)
     f_dyn(i) = (p.m_in(i)-Vapor_flow_rate*y_mass(i)-waste_mass*wf(i)+ r_free_hcl*V*p.vi_free_hcl(i)*p.PM(i)+r_free_formicacid*V*p.vi_free_formicacid(i)*p.PM(i)+r_free_aceticacid*V*p.vi_free_aceticacid(i)*p.PM(i)+r_free_h2so4*V*p.vi_free_h2so4(i)*p.PM(i))/Mtot;
end

%Evaluation of the Heat Duty
Q_in = sum((p.n_in.*p.c.*(p.Tin)));
Q_waste = sum(waste_mass./p.PM.*p.c.*mf);
Q_ev = sum(Vapor_flow_rate*p.lambda.*y./p.PM);
Q_used = Q_ev  + Q_waste -Q_in;
end                                                                                                                                             






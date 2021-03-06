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

%k sensitivity
k_ss = linspace(0.1,100,10^3);

%% Solving
%defining the first try
x0 = p.wp_in;
%defining vector to store solution
wf_ss =[];
T_ss = [];
Vapor_flow_rate_ss = [];
y_mass_ss =[];
MDEA_recovery_ss =[];
%options for fsolve
opts = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
for i = 1:length(k_ss)
    if i == 1
        x_try = x0;
    else
        x_try = sol_ss;
    end
    sol_ss = fsolve(@(x) ss_sensitivity(x,p.P_tot,p.waste,p.V,k_ss(i)),x_try,opts);
    %getting additional result
    [f_ss,T,Q_used,y_mass,Vapor_flow_rate,MDEA_recovery] = ss_sensitivity(sol_ss,p.P_tot,p.waste,p.V,k_ss(i));
    wf_ss = [wf_ss;sol_ss];
    T_ss = [T_ss;T];
    Vapor_flow_rate_ss = [Vapor_flow_rate_ss;Vapor_flow_rate];
    y_mass_ss = [y_mass_ss;y_mass];
    MDEA_recovery_ss = [MDEA_recovery_ss,MDEA_recovery];
end
%% Plotting and evaluating Error
wf_industrial = [0.1068 0.0151 0.0131 0.0117 0.032 0.0238 0.0146+0.0008 0.7781];
y_industrial = [0.02576 0.2958 0 0 0 0 0 0.6776];
T_industrial = 142.5+273.15;
Vapor_flow_rate_industrial = 16123;
MDEA_recovery_industrial = Vapor_flow_rate_industrial*y_industrial(1)/(p.m_in(1));

%Defining absolute error
%defining error on liquid mass fraction
AE_liquidfraction = abs((wf_ss-wf_industrial)*100);
%defining error on vapor mass fraction
AE_vaporfraction = abs(y_mass_ss-y_industrial)*100;
%Defining error temperature
AE_temperature = abs(T_ss-T_industrial)';
%Defining error recovery
AE_recovery = abs(MDEA_recovery_ss-MDEA_recovery_industrial)*100;

%plotting result Liquid fraction
subplot(4,2,1)
plot(k_ss,AE_liquidfraction(:,1),'linewidth',1.2) %MDEA
title('Absolute error [\%] MDEA','interpreter','latex')
subplot(4,2,2)
plot(k_ss,AE_liquidfraction(:,2),'linewidth',1.2) %H2O
title('Absolute error [\%] H2O','interpreter','latex')
subplot(4,2,3)
plot(k_ss,AE_liquidfraction(:,3),'linewidth',1.2) %NaOH
title('Absolute error [\%] NaOH','interpreter','latex')
subplot(4,2,4)
plot(k_ss,AE_liquidfraction(:,4),'linewidth',1.2) %HCl
title('Absolute error [\%] HCl','interpreter','latex','Fontsize',12)
subplot(4,2,5)
plot(k_ss,AE_liquidfraction(:,5),'linewidth',1.2) %Formicacid
title('Absolute error [\%] Formic Acid','interpreter','latex','Fontsize',12)
subplot(4,2,6)
plot(k_ss,AE_liquidfraction(:,6),'linewidth',1.2) %Aceticacid
title('Absolute error [\%] Acetic Acid','interpreter','latex','Fontsize',12)
subplot(4,2,7)
plot(k_ss,AE_liquidfraction(:,7),'linewidth',1.2) %solforicacid
title('Absolute error [\%] Solforic Acid ','interpreter','latex','Fontsize',12)
subplot(4,2,8)
plot(k_ss,AE_liquidfraction(:,8),'linewidth',1.2) %EG
title('Absolute error [\%] EG','interpreter','latex','Fontsize',12)
%Plotting Vapor fraction
figure
subplot(3,1,1)
plot(k_ss,AE_vaporfraction(:,1),'linewidth',1.2) %MDEA
title('Absolute error [\%] MDEA Vapor','interpreter','latex','Fontsize',12)
subplot(3,1,2)
plot(k_ss,AE_vaporfraction(:,2),'linewidth',1.2) %H2O
title('Absolute error [\%] H2O Vapor','interpreter','latex','Fontsize',12)
subplot(3,1,3)
plot(k_ss,AE_vaporfraction(:,8),'linewidth',1.2) %EG
title('Absolute error [\%] EG Vapor','interpreter','latex','Fontsize',12)
figure
%Plotting error temperature and recovery
subplot(2,1,1)
plot(k_ss,AE_temperature,'-r','linewidth',1.2) %Temperature
title('Absolute error Temperature','interpreter','latex','Fontsize',12)
subplot(2,1,2)
plot(k_ss,AE_recovery,'-r','linewidth',1.2) %Recovery
title('Absolute error Recovery [\%]','Interpreter','latex','Fontsize',12)














    
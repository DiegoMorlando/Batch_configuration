% SEMI FED BATCH SIMULATION : INPUT: HEAT DUTY AND PRESSURE)
clear all
close all
clc
global p 
%% INPUT

p=[];
%  %initial distribution [MEA H2O HEF NAOH DEGP HEPO CO2 Nasalt]
p.wp_in = [0.485 0.454 0.013 0.01 0.001 0.002 0.035 0];
%specific heat [MEA H2O HEF NaOH DEGP HEPO CO2] KJ/Kmol K
p.P_tot = 1*1.0135*10^2; %Pressure [KPa]
%density of the solution Kg/m^3
p.rho=998;
p.V = 1; %m^3
p.Mtot = p.V*p.rho;
p.feed_water = 200; %Kg/h
p.T_in_water = 40+273.15;
p.wf_feed_water = [0 1 0 0 0 0 0 0];

% [MEA H2O HEF NAOH DEGP HEPO CO2 ] Molecular weight [Kg/Kmol] 
p.PM = [61.08 18 89.08 40 68 144.17 44 63];
p.c = [155.8 74 100 90 83 100 0 100]; %KJ/Kmol K

% Kinetics parameter: Neutralization
p.k_free=10^7;                            %GUESS [m^3/Kmol h] Instantenous
% Stoic coefficient  [MEA H2O HEF NAOH DEGP HEPO CO2 Nasalt] Neutralization
p.vi_free = [1 0 -1 -1 0 0 0 +1];   
p.vi_deg1  = [-1 +0.5 0 0 +0.5 0 0 0]; %2MEA --> HEEDA + H2O
p.vi_deg2  = [-1 +0.5 0 0 +0.5 0 0 0]; %2MEA --> BHEU + H2O


% Antoine coefficient NIST Database
p.A = [4.29252 3.55959]; %MEA and H2O
p.B = [1408.873 643.748];
p.C = [-116.093 -198.043];

%Heat of reaction and heat of vaporization (i use the value as first approximation from the master
%thesis i found need to check on deltah_r) 
p.lambda = [50438 44000 0 0 0 0 16.7*1000 0]; % KJ/Kmol 
p.deltah_r  = -10^4; %

%Degradation properties MEA to HEEDA
p.k_ref = 1.599*10^-11*1000*3600; %m^3/Kmol h
p.Ea = 151.1*10^3; %KJ/Kmol
p.R = 8.31; %KJ/Kmol K 
%Degradation properties MEA to BHEU
p.k_ref2 = 1.281*10^-12*1000*3600; %m^3/Kmol h
p.T_ref = 400; %K

%Heat duty
p.Q_in = 10^5; %KJ/h

%How much do we want to continue the flux of water
p.t_flux_water = 1; %[h]
%Viscosity Parameter
p.beta = 1;
p.viscosity_upper_limit = 500;
p.t_start =[];
%% Integration
%giving the time span and initial condition
t_span = linspace(0,15,1000); %[h]
% [MEA H2O HEF NAOH DEGP HEPO CO2 Nasalt Volume]
x0 = [p.wp_in, p.V]; 
options = odeset('RelTol',10^-10,'AbsTol',10^-10*length(x0));
[t, sol_dyn] = ode15s(@(t,x) batch_configuration_function(t,x,p.P_tot,p.Q_in,p.beta),t_span,x0,options);
%getting other results
for i = 1:length(t)
    [~,T(i),Vapor_flow_rate(i),y_mea(i),y_h2o(i),y_co2(i),viscosity(i),Mea_recovery_rate(i),r_deg(i),feed_water(i)] = batch_configuration_function(t(i), sol_dyn(i,:)',p.P_tot,p.Q_in,p.beta);
end

%% Plotting
%Integration Mea recovery
dMea_recovery_result=[];
for i = 1 :length(t)
    if i == length(t)
    dMea_recovery=Mea_recovery_rate(i)*(t(i)-t(i-1));
    else
    dMea_recovery=Mea_recovery_rate(i)*(t(i+1)-t(i));
    end
    dMea_recovery_result = [dMea_recovery_result,dMea_recovery];
    Mea_recovery = sum(dMea_recovery_result)*100;
    subplot(4,2,1)
    hold on
    plot(t(i),Mea_recovery,'o-b','MarkerSize',3)
    title('Mea Recovery [\%]','FontSize',12,'Interpreter','latex','FontWeight','bold')
end
disp(Mea_recovery(end));
%Integration Mea Degradated

dMea_degradated_result=[];
for i = 1 :length(t)
    if i == length(t)
    dMea_degradated=r_deg(i)*sol_dyn(i,end)*(t(i)-t(i-1))*p.PM(1); %Kg
    else
    dMea_degradated=r_deg(i)*sol_dyn(i,end)*(t(i+1)-t(i))*p.PM(1);  %Kg
    end
    dMea_degradated_result = [dMea_degradated_result,dMea_degradated];
    Mea_degradated = sum(dMea_degradated_result); %g
    subplot(4,2,2)
    hold on
    plot(t(i),Mea_degradated/(p.V*p.wp_in(1)*p.rho)*100,'o-r','MarkerSize',3)
    xlabel('time [h]','FontSize',12,'Interpreter','latex','FontWeight','bold')
    title('Fraction of MEA degradated [\%]','FontSize',12,'Interpreter','latex','FontWeight','bold')
end
% Mea_degradated(end)/(p.V*p.wp_in(1)*p.rho)*100
% T(end)-273.15
% viscosity(end)
% p.feed_water*(t(end)-p.t_start(1))
% p.t_start(1)

subplot(4,2,3)
plot(t,sol_dyn(:,1),'-r','linewidth',1.2)
hold on
plot(t,sol_dyn(:,2),'-b','linewidth',1.2)
plot(t,sol_dyn(:,7),'-g','linewidth',1.2)
title('Mass liquid fraction','FontSize',12,'Interpreter','latex','FontWeight','bold')
legend('MEA','H2O','CO2','Location','best','Orientation','vertical')
hold off
subplot(4,2,4)
plot(t,y_mea,'-r','linewidth',1.2)
hold on
plot(t,y_h2o,'-b','linewidth',1.2)
plot(t,y_co2,'-g','linewidth',1.2)
title('Mass vapor fraction','FontSize',12,'Interpreter','latex','FontWeight','bold')
legend('MEA','H2O','CO2','Location','best','Orientation','vertical')
    subplot(4,2,5)
plot(t,viscosity,'-k','LineWidth',1.2)
title_viscosity = 'Viscosity $[mPa \cdot s]$ at 30  $^{\circ}$ C $\beta$  = ' + string(p.beta);
title(title_viscosity,'FontSize',12,'Interpreter','latex','FontWeight','bold')
    subplot(4,2,6)
plot(t,sol_dyn(:,end),'-b','linewidth',1.2)
title('Volume $[m^3]$','FontSize',12,'Interpreter','latex','FontWeight','bold')
    subplot(4,2,7)
plot(t,Vapor_flow_rate,'-r','linewidth',1.2)
title('Mass Vapor flow rate $[\frac{Kg}{h}]$','FontSize',12,'Interpreter','latex','FontWeight','bold')
xlabel('time [h]','FontSize',12,'Interpreter','latex','FontWeight','bold')
    subplot(4,2,8)
plot(t,T-273.15,'-r','linewidth',1.2)
title('Temperature $^{\circ}$ C','FontSize',12,'Interpreter','latex','FontWeight','bold')
xlabel('time [h]','FontSize',12,'Interpreter','latex')
figure
plot(t,sol_dyn(:,3),'-m','linewidth',1.2)
hold on
plot(t,sol_dyn(:,4),'-r','linewidth',1.2)
plot(t,sol_dyn(:,5)+sol_dyn(:,6),'-b','linewidth',1.2)
xlabel('time [h]','FontSize',12,'Interpreter','latex')
ylabel('Mass liquid fraction','FontSize',12,'Interpreter','latex')
legend('HEF','NaOH','DEGCOMPOUNDS')
figure
plot(t,feed_water,'-b','linewidth',1.5)


%% Sensitivity on beta
sensitivity_beta()








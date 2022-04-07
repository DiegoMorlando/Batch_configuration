% Sensitivity on beta parameter 
function [] = sensitivity_beta()
global p beta viscosity_result
%giving the time span and initial condition
t_span = linspace(0,15,10^3); %[h]
% [MEA H2O HEF NAOH DEGP HEPO CO2 Nasalt Volume]
x0 = [p.wp_in, p.V]; 
options = odeset('RelTol',10^-10,'AbsTol',10^-10*length(x0));
beta = linspace(0,2,20);
%storing the final viscosity result
viscosity_result = [];
for i = 1:length(beta)
[t, sol_dyn] = ode15s(@(t,x) batch_configuration_function(t,x,p.P_tot,p.Q_in,beta(i)),t_span,x0,options);
for j= 1:length(t)
[~,~,~,~,~,~,viscosity(j),~,~,~] = batch_configuration_function(t(j), sol_dyn(j,:)',p.P_tot,p.Q_in,beta(i));
end
viscosity_result = [viscosity_result,viscosity(end)];
end
plot(beta,viscosity_result,'*b','Marker','*','MarkerSize',9)
xlabel('$\beta$ ','Interpreter','latex','FontSize',12)
ylabel('$\mu$ Dynamic Viscosity $[mPa \cdot s]$','Interpreter','latex','FontSize',12)
msg = "T = 30 $^{\circ}$ C after " + string(t(end)) + " $[h]$";
title(msg,'interpreter','latex','Fontsize',13)
end














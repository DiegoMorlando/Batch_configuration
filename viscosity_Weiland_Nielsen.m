function viscosity = viscosity_Weiland_Nielsen(wf_mea,T,alfa,c_salt,beta)
% T [=] K c = [M]
a = 0;
b = 0;
c = 21.186;
d = 2373;
e = 0.01015;
f = 0.0093;
g = -2.2589;
omega = wf_mea*100;
ln_mi_over_mi_h2o = ((a*omega+b)*T+(c*omega+d))*(alfa*(e*omega+f*T+g)+1)*omega/T^2;
A = 1.856*10^-11;
B = 4209;
C= 0.04527;
D = -3.376*10^-5;
mi_water = A*exp(B/T+C*T+D*T^2);
mi_amine_clean = mi_water*exp(ln_mi_over_mi_h2o);
%degradation Nielsen
viscosity = exp(beta*c_salt)*mi_amine_clean; 

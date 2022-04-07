function psat_co2_vle= psat_co2(T,alfa)
A_vle = 1.8;
B_vle = 10;
k1 = -9155.95*(1/T)+28.03;
k2 = exp(-6146.18/T+15);
k3 = 7527.04*(1/T)-16.94;
% dovrei usare l
if isnan(alfa) == 1 
    psat_co2_vle = 0;
else
    psat_co2_vle =exp(A_vle*log(alfa)+k1+B_vle/(1+k2*exp(-k3*log(alfa))));
end


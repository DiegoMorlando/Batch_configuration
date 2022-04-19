function vp_h2o = vp_h2o(T)
A = 3.55959;
B = 643.748;
C = -198.043;
vp_h2o = 10^(A-B/(T+C))*1.013*10^2;

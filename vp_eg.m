function vp_eg = vp_eg(T)
A = 4.97012;
B = 1914.951;
C = -84.996;
vp_eg = 10^(A-B/(T+C))*1.013*10^5/10^3;
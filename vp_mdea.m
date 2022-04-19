function vp_mdea = vp_mdea(T)
c1 = 23.3224;
c2 = -5126.49;
c3 = -81.8475;
vp_mdea= exp(c1+c2/(c3+T))/10^3;
function vp_h2o = vp_h2o(T)
global p 
vp_h2o=10.^(p.A(2)-p.B(2)./(T+p.C(2)))*1.013*100;          %Kpa
end
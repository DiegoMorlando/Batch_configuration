function vp_mea = vp_mea(T)
global p 
vp_mea=10.^(p.A(1)-p.B(1)./(T+p.C(1)))*1.013*100;  %Kpa
end
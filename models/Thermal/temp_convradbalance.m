function [Tw] = temp_convradbalance (T_h, rho_h, v, vehicle, const_sb)

k = 1.7415e-4;  % Earth
% k = 1.9027e-4;  % Mars
q_conv = k*(rho_h/vehicle.Rnose)^0.5*v^3;
Tw = (T_h^4+q_conv/(const_sb*vehicle.Epsilon))^0.25;
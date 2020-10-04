function w_p = w_p_gyro(B,m,q)
% W_P_GYRO - angular (charged) particle gyro frequency (non-relativistic)
%   
% Calling: 
%  w_p = w_p_gyro(B,m,q)
% Input:
%  B magnetic field [T] (can also be a vector with differernt values of B),
%  m mass of particle [kg]
%  q charge of particle [C]

w_p = q.*B/m;

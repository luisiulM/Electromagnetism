function v = v_of_E(E)
% Transform from electron energy (eV) to velocity (m/s)
%
% Calling:
%  v = v_of_E(E)
% Input:
%  E - electron energy (eV)  double [ N x M ]
% Output:
%  v - electron velocity (m/s) double [ N x M ]

m_e     = 9.10939e-31;
q_e     = 1.6021773e-19;

v = (2*q_e*E/m_e).^(1/2);

function drdtdvdt = ode_ExB(t,rv,B,E,m,nq)
% ODE_EXB - momentum eq for particles in E and B fields
%  ode_ExB is the equation of motion for charged particles in
%  static and homogenous E and B-fields, suitable for integration
%  with matlab's ode-suite (ode23, ode45, ode113, ode15s, ode23s,
%  ode23t, ode23tb)
%  
%  The momentum equation 
%   d²r/dt² = q/m*( E + v x B) 
%  is here rewritten as 6 first-order ODEs as is required for the 
%  odeXX-functions
%
% Calling:
%  drdtdvdt = ode_ExB(t,rv,B,E,m,nq)
% Input:
%  t  - time, double scalar (s)
%  rv - position and velocity-vector, double array [6 x 1] (m, and m/s)
%  B  - magnetic field vector, double array [3 x 1] (T)
%  E  - Electric field vector, doublde array [3 x 1] (V/m)
%  m  - (optional) particle mass, double scalar, (kg), defaults to
%       electron mass
%  nq - (optional) charge state, scalar integer, defaults to 1
% Output:
%  drdtdvdt - velocity and acceleration, double array [6 x 1] (m/s,
%             m/s^2)


m_e     = 9.10938291e-31;         % electron rest mass [kg]
q_e     = 1.602176565e-19;        % elementary charge [C]


if nargin < 6 || isempty(nq)
  nq = 1;
end
if nargin < 5 || isempty(m)
  m = m_e;
end

v = rv(4:end);

drdtdvdt = rv(:);
drdtdvdt(1:3) = v;

a = nq*q_e/m*(E + cross(v,B));
drdtdvdt(4:6) = a;




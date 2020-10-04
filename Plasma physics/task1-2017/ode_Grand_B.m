function drdtdvdt = ode_Grand_B(t,rv,B0,Grand_B,E,q,m,nq)
% ODE_Grand_B - momentum eq for particles in E and B fields
%  ode_Grand_B is the equation of motion for charged particles in
%  static and nonuniform B-fields, suitable for integration
%  with matlab's ode-suite (ode23, ode45, ode113, ode15s, ode23s,
%  ode23t, ode23tb)
%  
%  The momentum equation 
%   d²r/dt² = q/m*( E + v x (B0 + y Grand_B)) 
%  is here rewritten as 6 first-order ODEs as is required for the 
%  odeXX-functions
%
% Calling:
%  drdtdvdt = ode_Grand_B(t,rv,B0,Grand_B,E,q,m,nq)
% Input:
%  t  - time, double scalar (s)
%  rv - position and velocity-vector, double array [6 x 1] (m, and m/s)
%  B0 - Start magnetic field vector, double array [3 x 1] (T)
%  Grand_B - gradient of the magnetic field, scaler (T)
%  E  - Electric field vector, doublde array [3 x 1] (V/m)
%  q  - (optional) electric charge of the particle, scalar (C), defaults to
%       electron charge
%  m  - (optional) particle mass, double scalar (kg), defaults to
%       electron mass
%  nq - (optional) charge state, scalar integer, defaults to 1
% Output:
%  drdtdvdt - velocity and acceleration, double array [6 x 1] (m/s,
%             m/s^2)


m_e     = 9.10938291e-31;          % electron rest mass [kg]
q_e     = -1.602176565e-19;        % electron charge [C]


if nargin < 8 || isempty(nq)
  nq = 1;
end
if nargin < 7 || isempty(m)
  m = m_e;
end
if nargin < 6 || isempty(q)
  q = q_e;
end

v = rv(4:end);
y = rv(2);

drdtdvdt = rv(:);
drdtdvdt(1:3) = v;

a = nq*q/m*(E + cross(v,(B0 + y*Grand_B*sign(B0))));
drdtdvdt(4:6) = a;


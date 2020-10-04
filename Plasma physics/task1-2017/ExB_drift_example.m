%% Example for calculation of electron ExB trajectory
%  This script calculates the ExB drift of an 1 Ev electron in a
%  50000 nT magnetic field with an electric field of 1 V/m
%  Here we will use one of the ODE-integrating functions in matlab
%  ODE23T. To integrate the equation of motion we need to write that
%  function as a set of first-order ODEs, see ode_ExB.
clear all

%% B and E-fields

B = 5e-5*[0 0 1];        % Magnetic field || z
E = 1 *[0 1 0];          % Electric field || y
m_e = 9.10939e-31;       % electron rest mass [kg]
q_e = 1.6021773e-19;     % elementary charge [C]

%% Time-span of simulation
w_gyro = w_p_gyro(norm(B),m_e,q_e); % Angular electron gyro frequency
T_gyro = (2*pi)/w_gyro;
n_gyrations = 30;
n_per_orbit = 100;
T_span = linspace(0,n_gyrations*T_gyro,n_gyrations*n_per_orbit+1);

%% Initial conditions
v0 = v_of_E(1);        % speed of 1 eV electron
r0v0 = [0;0;0;v0;0;0]; % Starting position and velocity

%% Integration of the Momentum Equation:
[t23,rv23] = ode23t(@(t,rv)  ode_ExB(t,rv,B,E),T_span,r0v0);


%% Inspection of results 1, trajectory
figT = figure;
plot(rv23(:,1),rv23(:,2),'linewidth',2)
xlabel('X-Distance (m)')
ylabel('Y-Distance (m)')

%% Inspection of results 1, energy conservation

% For kinetic energy we also need particle mass
m_e     = 9.10938291e-31;         % electron rest mass [kg]
q_e     = 1.602176565e-19;        % elementary charge [C]

figE = figure;
ph1 = plot(t23,m_e/2*sum(rv23(:,4:end).^2,2),'linewidth',2);
hold on
ph2 = plot(t23,-q_e*E(2)*rv23(:,2),'linewidth',2);
ph3 = plot(t23,m_e/2*sum(rv23(:,4:end).^2,2)-q_e*E(2)*rv23(:,2),'linewidth',2);
legend([ph1,ph2,ph3],'E_{kin}','E_E','E_{tot}')
%% Looks reasonable


%% Check that the integrating function makes the particle end up at
%  starting point when integrating ME backwards:

[t23B,rv23B] = ode23t(@(t,rv)  ode_ExB(t,rv,B,E),T_span(end:-1:1),rv23(end,:)');

figure(figT)
phT1 = plot(rv23(:,1),rv23(:,2),'linewidth',2);
hold on
phT2 = plot(rv23B(:,1),rv23B(:,2),'linewidth',2);

xlabel('X-Distance (m)')
ylabel('Y-Distance (m)')
legend([phT1,phT2],'Forward trajectory','Reversed trajectory')

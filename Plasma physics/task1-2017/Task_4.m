%% Task 4 from Fall 2017 Home Exam. candidate number: 8
% Reference: 
% This code was heavelly inpered by ExB_drift_example.m. Found on the 
% fronter rom ''FYS-2009 - Innf. i plasmafys. - 1. Sem (H2017)'' under 
% ''Arkiv'', inside ''Programming task 1''.

% This script calculates the gradiant-B drift of an 1 eV electron in a
% 50000 nT magnetic field with 6 different cases of magnetic field strength
% gradiants in the y-direction. Here we use both the gyro-centre 
% approximation and the ODE-integrating functions in matlab
% ODE23T, to calculate the motion of the electron. To integrate the 
% equation of motion we need to write the function as a set of first-order 
% ODEs, see ode_Grand_B.

% Assumtion: 
% Since there's no mention of electric field, we assume it to be zero. 

%% Known and defined variable
EV = 1;                           % Number of Electron volts [eV]
E = [0 0 0];                    % Electric field [V/m]

m_e = 9.10939e-31;                % electron mass [kg]
q_e = -1.6021773e-19;              % electron charge [C]

B0 = 5e-5*[0 0 1];                     % Start magnetic field [T] || z
Grand_B = [5;50;350;370;380;500]*1e-6; % Magnetic field strength gradient for 6 cases

%% Initial conditions
v0 = v_of_E(EV);                % Initial speed
r0v0 = [0;0;0;0;v0;0];          % Starting position and velocity
 
%% Time-span of simulation for all 6 cases
y_mid = (v0*m_e)./(abs(q_e)*norm(B0)); % approximately halfway between the max and min value 
B_mid = (B0(3) + y_mid.*Grand_B);

w_gyro = -w_p_gyro(B_mid,m_e,q_e);     % Angular electron gyro frequency for all 6 cases
T_gyro = (2*pi)./abs(w_gyro);          % Gyro-period for all 6 cases
n_gyro = 5;                            % Number of gyro-periods
n_per_orbit = 100;                     % Number of points per orbit 

T_span1 = linspace(0,n_gyro*T_gyro(1),n_gyro*n_per_orbit);
T_span2 = linspace(0,n_gyro*T_gyro(2),n_gyro*n_per_orbit);
T_span3 = linspace(0,n_gyro*T_gyro(3),n_gyro*n_per_orbit);
T_span4 = linspace(0,n_gyro*T_gyro(4),n_gyro*n_per_orbit);
T_span5 = linspace(0,n_gyro*T_gyro(5),n_gyro*n_per_orbit);
T_span6 = linspace(0,n_gyro*T_gyro(6),n_gyro*n_per_orbit);

%% Integration of the Momentum Equation for all 6 cases:

[t1,rv1] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(1),E,q_e),T_span1,r0v0);
[t2,rv2] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(2),E,q_e),T_span2,r0v0);
[t3,rv3] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(3),E,q_e),T_span3,r0v0);
[t4,rv4] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(4),E,q_e),T_span4,r0v0);
[t5,rv5] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(5),E,q_e),T_span5,r0v0);
[t6,rv6] = ode23t(@(t,rv)  ode_Grand_B(t,rv,B0,Grand_B(6),E,q_e),T_span6,r0v0);

t = [t1 t2 t3 t4 t5 t6];
%% Gyro-centre approximation for all 6 cases

v_perp = sqrt( r0v0(4)^2 + r0v0(5)^2 );             % velocity perpendicular to the B-field
r_L = v_perp./w_gyro;                               % Larmor radius
v_D = @(c) (v0.*r_L(c).*Grand_B(c))./(2.*norm(B0)); % Gradiant_B Drift velocity for each case (s)

N = length(t(:,1));
x1 = zeros(1,N);   y1 = zeros(1,N);
x2 = zeros(1,N);   y2 = zeros(1,N);
x3 = zeros(1,N);   y3 = zeros(1,N);
x4 = zeros(1,N);   y4 = zeros(1,N);
x5 = zeros(1,N);   y5 = zeros(1,N);
x6 = zeros(1,N);   y6 = zeros(1,N);

for i = 1:N
    % In order to start at x(0) = 0, we need to subtract the x-position
    % equation by the Larmor radius. Due to the fact that cos(0) = 1, 
    % giving us a starting position equal to the radius.
    
    x1(i) = -r_L(1) + r_L(1).*cos(w_gyro(1)*t(i,1)) + v_D(1).*t(i,1);
    y1(i) = r0v0(2) + r_L(1).*sin(w_gyro(1)*t(i,1));
    
    x2(i) = -r_L(2) + r_L(2).*cos(w_gyro(2)*t(i,2)) + v_D(2).*t(i,2);
    y2(i) = r0v0(2) + r_L(2).*sin(w_gyro(2)*t(i,2));
    
    x3(i) = -r_L(3) + r_L(3).*cos(w_gyro(3)*t(i,3)) + v_D(3).*t(i,3);
    y3(i) = r0v0(2) + r_L(3).*sin(w_gyro(3)*t(i,3));
    
    x4(i) = -r_L(4) + r_L(4).*cos(w_gyro(4)*t(i,4)) + v_D(4).*t(i,4);
    y4(i) = r0v0(2) + r_L(4).*sin(w_gyro(4)*t(i,4));
    
    x5(i) = -r_L(5) + r_L(5).*cos(w_gyro(5)*t(i,5)) + v_D(5).*t(i,5);
    y5(i) = r0v0(2) + r_L(5).*sin(w_gyro(5)*t(i,5));
    
    x6(i) = -r_L(6) + r_L(6).*cos(w_gyro(6)*t(i,6)) + v_D(6).*t(i,6);
    y6(i) = r0v0(2) + r_L(6).*sin(w_gyro(6)*t(i,6));
end

%% Ploting of the numerical and gyro-centre approximation of the x- and y-direction
figXY = figure;

subplot(2,3,1)
plot(rv1(:,1),rv1(:,2))
hold on
plot(x1,y1)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')
legend('Numerical integration','Gyro-centre approx')

subplot(2,3,2)
plot(rv2(:,1),rv2(:,2))
hold on
plot(x2,y2)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')

subplot(2,3,3)
plot(rv3(:,1),rv3(:,2))
hold on
plot(x3,y3)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')

subplot(2,3,4)
plot(rv4(:,1),rv4(:,2))
hold on
plot(x4,y4)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')

subplot(2,3,5)
plot(rv5(:,1),rv5(:,2))
hold on
plot(x5,y5)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')

subplot(2,3,6)
plot(rv6(:,1),rv6(:,2))
hold on
plot(x6,y6)
xlabel('X-Distance [m]')
ylabel('y-Distance [m]')

%% Ploting of the numerical and gyro-centre approximation of the x-direction as a function of time
figXt = figure;

subplot(2,3,1)
plot(t1,rv1(:,1))
hold on
plot(t1,x1)
xlabel('Time [s]')
ylabel('X-Distance [m]')
legend('Numerical integration','Gyro-centre approx')

subplot(2,3,2)
plot(t2,rv2(:,1))
hold on
plot(t2,x2)
xlabel('Time [s]')
ylabel('X-Distance [m]')

subplot(2,3,3)
plot(t3,rv3(:,1))
hold on
plot(t3,x3)
xlabel('Time [s]')
ylabel('X-Distance [m]')

subplot(2,3,4)
plot(t4,rv4(:,1))
hold on
plot(t4,x4)
xlabel('Time [s]')
ylabel('X-Distance [m]')

subplot(2,3,5)
plot(t5,rv5(:,1))
hold on
plot(t5,x5)
xlabel('Time [s]')
ylabel('X-Distance [m]')

subplot(2,3,6)
plot(t6,rv6(:,1))
hold on
plot(t6,x6)
xlabel('Time [s]')
ylabel('X-Distance [m]')

%% Ploting of the numerical and gyro-centre approximation of the y-direction as a function of time
figYt = figure;

subplot(2,3,1)
plot(t1,rv1(:,2))
hold on
plot(t1,y1)
xlabel('Time [s]')
ylabel('y-Distance [m]')
legend('Numerical integration','Gyro-centre approx')

subplot(2,3,2)
plot(t2,rv2(:,2))
hold on
plot(t2,y2)
xlabel('Time [s]')
ylabel('y-Distance [m]')

subplot(2,3,3)
plot(t3,rv3(:,2))
hold on
plot(t3,y3)
xlabel('Time [s]')
ylabel('y-Distance [m]')

subplot(2,3,4)
plot(t4,rv4(:,2))
hold on
plot(t4,y4)
xlabel('Time [s]')
ylabel('y-Distance [m]')

subplot(2,3,5)
plot(t5,rv5(:,2))
hold on
plot(t5,y5)
xlabel('Time [s]')
ylabel('y-Distance [m]')

subplot(2,3,6)
plot(t6,rv6(:,2))
hold on
plot(t6,y6)
xlabel('Time [s]')
ylabel('y-Distance [m]')
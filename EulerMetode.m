% Startverdi til differensial ligning x(t0) = x0
% tids interval [t0,t1] 

function x = EulerMetode(t0,x0,t1)
% Definerer fysiske størrelser
w = [5000 10000 20000]; % verdier til omega
L = 5*10^(-3);          % Spolen          
R = 50;                 % Motstand
V0 = 1;                 % Spennings styrke
tau = L/R;

% Tids perioder til hver frekvens
dt1 = 1/(4*2*pi*w(1));
time1 = t0:dt1:t1;

dt2 = 1/(4*2*pi*w(2));
time2 = t0:dt2:t1;

dt3 = 1/(4*2*pi*w(3));
time3 = t0:dt3:t1;

% Differensial ligning til hver frekvens
dx1 = @(t,x) 1/tau*( V0*sign( sin(w(1)*t) ) - x );
dx2 = @(t,x) 1/tau*( V0*sign( sin(w(2)*t) ) - x );
dx3 = @(t,x) 1/tau*( V0*sign( sin(w(3)*t) ) - x );

% inisial verdier
x1(1) = x0;
x2(1) = x0;
x3(1) = x0;

% Euler metode for alle tre frekvenser
for i=1:length(time1)-1
  x1(i+1)= x1(i) + dt1*dx1(time1(i),x1(i));
end
for i=1:length(time2)-1
  x2(i+1)= x2(i) + dt2*dx2(time2(i),x2(i));
end
for i=1:length(time3)-1
  x3(i+1)= x3(i) + dt3*dx3(time3(i),x3(i));
end

plot(time1,x1,'b',time2,x2,'r',time3,x3,'g');
legend('w = 5000 rad/s', 'w = 10000 rad/s', 'w = 20000 rad/s')
xlabel('Tid [s]')
ylabel('V_R [V]')

end

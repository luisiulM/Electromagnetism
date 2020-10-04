clear all
% MATLAB skript som visualiserer elektrisk felt og potensiale fra en enkeltladning

% definerer fysiske stoerrelser
   
    q1 = 1E-9;           % stoerrelse paa til ladning 1
    q2 = -1E-9;          % stoerrelse paa til ladning 2
    
    r0 = 0.05;           % radius for ladning
    
    e0 = 8.854188E-12;  % permittivitet i vakuum
    Ke = 1/(4*pi*e0);   % elektrostatis konstant i Coulombs lov

% definerer ytre grenser for omraadet som skal vises i xz-planet

    xmin = -1;
    xmax = 1;    
    zmin = -1;
    zmax = 1;
   
% definerer avstand mellom punkter (samples) i x- og z-retningen

    dx = 0.1;
    dz = 0.1;
   
% definerer vektorer for samplingspunktene langs x- og yz-aksen    
    
    xv = xmin:dx:xmax;
    zv = zmin:dz:zmax;

% bruker meshgrid funksjon til aa lage 2D matiser fra x- og z-samples
% dette gjoer det lett aa beregne felter og potensialer vha matriseoperajoner. 

    [X,Z] = meshgrid(xv,zv);
    
    Xp1 = 0;     % x posisjon til partikkel 1
    Zp1 = 0.2;   % z posisjon til partikkel 1
    
    Xp2 = 0;     % x posisjon til partikkel 2
    Zp2 = -0.2;  % z posisjon til partikkel 2
    
    r1 = sqrt((X-Xp1).*(X-Xp1)+(Z-Zp1).*(Z-Zp1));% beregner avstand fra sentrum av partikkel 1   
    r2 = sqrt((X-Xp2).*(X-Xp2)+(Z-Zp2).*(Z-Zp2));% beregner avstand fra sentrum av partikkel 2
    
    ind=find(r1<r0);     % finner indekser (posisjoner) til punkter der r1<r0  
    r1(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r2<r0);     % finner indekser (posisjoner) til punkter der r2<r0
    r2(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    V1 = Ke*q1./r1;       % beregner potensiale til partikkel 1
    V2 = Ke*q2./r2;       % beregner potensiale til partikkel 2
    Total_V = V1 + V2;          % beregner den totale potensiale
    
    Ex1 = Ke*q1*(X-Xp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på x-aksen
    Ez1 = Ke*q1*(Z-Zp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på z-aksen
    
    Ex2 = Ke*q2*(X-Xp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på x-aksen
    Ez2 = Ke*q2*(Z-Zp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på z-aksen
    Ex = Ex1 + Ex2;                 % beregner totale elektrisk felt på x-aksen
    Ez = Ez1 + Ez2;                 % beregner totale elektrisk felt på z-aksen
    
% plotter felter i xz-planet

    figure(1)
    subplot(1,1,1)
    
    contour(X,Z,Total_V);      % viser potensialet som et contourplot
    hold on
    quiver(X,Z,Ex,Ez)    % viser det elektriske feltet som et vektorplot 
    xlabel('x-aksen')
    ylabel('z-aksen')
    hold off
    axis equal
    
%% oppgave 1 c)

    plot(Z,Total_V)            % Plotter Coulombfeltet for V som funksjon av z langs z-aksen
    xlabel('z-aksen')
    ylabel('Elektrisk potensiale V')

%% oppgave 1 d) I

clear all
% MATLAB skript som visualiserer elektrisk felt og potensiale fra en enkeltladning

% definerer fysiske stoerrelser
   
    q1 = 1E-9;           % stoerrelse paa til ladning 1
    q2 = -1E-9;          % stoerrelse paa til ladning 2
    q3 = 1E-9;           % stoerrelse paa til ladning 3
    
    r0 = 0.05;           % radius for ladning
    
    e0 = 8.854188E-12;  % permittivitet i vakuum
    Ke = 1/(4*pi*e0);   % elektrostatis konstant i Coulombs lov

% definerer ytre grenser for omraadet som skal vises i xz-planet

    xmin = -0.4;
    xmax = 0.4;    
    zmin = -0.4;
    zmax = 0.4;
   
% definerer avstand mellom punkter (samples) i x- og z-retningen

    dx = 0.02;
    dz = 0.02;
   
% definerer vektorer for samplingspunktene langs x- og yz-aksen    
    
    xv = xmin:dx:xmax;
    zv = zmin:dz:zmax;

% bruker meshgrid funksjon til aa lage 2D matiser fra x- og z-samples
% dette gjoer det lett aa beregne felter og potensialer vha matriseoperajoner. 

    [X,Z] = meshgrid(xv,zv);
    
    Xp1 = 0;     % x posisjon til partikkel 1
    Zp1 = 0.2;   % z posisjon til partikkel 1
    
    Xp2 = 0;     % x posisjon til partikkel 2
    Zp2 = 0;  % z posisjon til partikkel 2
    
    Xp3 = 0;     % x posisjon til partikkel 3
    Zp3 = -0.2;  % z posisjon til partikkel 3
    
    r1 = sqrt((X-Xp1).*(X-Xp1)+(Z-Zp1).*(Z-Zp1));% beregner avstand fra sentrum av partikkel 1   
    r2 = sqrt((X-Xp2).*(X-Xp2)+(Z-Zp2).*(Z-Zp2));% beregner avstand fra sentrum av partikkel 2
    r3 = sqrt((X-Xp3).*(X-Xp3)+(Z-Zp3).*(Z-Zp3));% beregner avstand fra sentrum av partikkel 3
    
    ind=find(r1<r0);     % finner indekser (posisjoner) til punkter der r1<r0  
    r1(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r2<r0);     % finner indekser (posisjoner) til punkter der r2<r0
    r2(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r3<r0);     % finner indekser (posisjoner) til punkter der r3<r0
    r3(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    V1 = Ke*q1./r1;         % beregner potensiale til partikkel 1
    V2 = Ke*q2./r2;         % beregner potensiale til partikkel 2
    V3 = Ke*q3./r3;         % beregner potensiale til partikkel 3
    Total_V = V1 + V2 + V3; % beregner den totale potensiale
    
    Ex1 = Ke*q1*(X-Xp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på x-aksen
    Ez1 = Ke*q1*(Z-Zp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på z-aksen
    
    Ex2 = Ke*q2*(X-Xp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på x-aksen
    Ez2 = Ke*q2*(Z-Zp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på z-aksen
    
    Ex3 = Ke*q3*(X-Xp3)./(r3.^3);   % beregner elektrisk felt til partikkel 3 på x-aksen
    Ez3 = Ke*q3*(Z-Zp3)./(r3.^3);   % beregner elektrisk felt til partikkel 3 på z-aksen
    
    Total_Ex = Ex1 + Ex2 + Ex3;    % beregner totale elektrisk felt på x-aksen
    Total_Ez = Ez1 + Ez2 + Ez3;    % beregner totale elektrisk felt på z-aksen
    
% plotter felter i xz-planet

    figure(1)
    subplot(1,1,1)
    
    contour(X,Z,Total_V);      % viser potensialet som et contourplot
    hold on
    quiver(X,Z,Total_Ex,Total_Ez)    % viser det elektriske feltet som et vektorplot 
    xlabel('x-aksen')
    ylabel('z-aksen')
    hold off
    axis equal
    
%% oppgave 1 e I
    
    %subplot(1,2,1)
    plot(Z,Total_V)            % Plotter Coulombfeltet for V som funksjon av z langs z-aksen
    xlabel('z-aksen')
    ylabel('Elektrisk potensiale V')
    hold on

%% oppgave 1 d) II
    
 clear all
% MATLAB skript som visualiserer elektrisk felt og potensiale fra en enkeltladning

% definerer fysiske stoerrelser
   
    q1 = 1E-9;           % stoerrelse paa til ladning 1
    q2 = -1E-9;          % stoerrelse paa til ladning 2
    q3 = 1E-9;           % stoerrelse paa til ladning 3
    q4 = -1E-9;          % stoerrelse paa til ladning 4
    
    r0 = 0.05;           % radius for ladning
    
    e0 = 8.854188E-12;  % permittivitet i vakuum
    Ke = 1/(4*pi*e0);   % elektrostatis konstant i Coulombs lov

% definerer ytre grenser for omraadet som skal vises i xz-planet

    xmin = -0.4;
    xmax = 0.4;    
    zmin = -0.4;
    zmax = 0.4;
   
% definerer avstand mellom punkter (samples) i x- og z-retningen

    dx = 0.02;
    dz = 0.02;
   
% definerer vektorer for samplingspunktene langs x- og yz-aksen    
    
    xv = xmin:dx:xmax;
    zv = zmin:dz:zmax;

% bruker meshgrid funksjon til aa lage 2D matiser fra x- og z-samples
% dette gjoer det lett aa beregne felter og potensialer vha matriseoperajoner. 

    [X,Z] = meshgrid(xv,zv);
    
    Xp1 = 0;     % x posisjon til partikkel 1
    Zp1 = 0.2;   % z posisjon til partikkel 1
    
    Xp2 = -0.2;  % x posisjon til partikkel 2
    Zp2 = 0;     % z posisjon til partikkel 2
    
    Xp3 = 0;     % x posisjon til partikkel 3
    Zp3 = -0.2;  % z posisjon til partikkel 3
    
    Xp4 = 0.2;   % x posisjon til partikkel 4
    Zp4 = 0;     % z posisjon til partikkel 4
    
    r1 = sqrt((X-Xp1).*(X-Xp1)+(Z-Zp1).*(Z-Zp1));% beregner avstand fra sentrum av partikkel 1   
    r2 = sqrt((X-Xp2).*(X-Xp2)+(Z-Zp2).*(Z-Zp2));% beregner avstand fra sentrum av partikkel 2
    r3 = sqrt((X-Xp3).*(X-Xp3)+(Z-Zp3).*(Z-Zp3));% beregner avstand fra sentrum av partikkel 3
    r4 = sqrt((X-Xp4).*(X-Xp4)+(Z-Zp4).*(Z-Zp4));% beregner avstand fra sentrum av partikkel 4
    
    ind=find(r1<r0);     % finner indekser (posisjoner) til punkter der r1<r0  
    r1(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r2<r0);     % finner indekser (posisjoner) til punkter der r2<r0
    r2(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r3<r0);     % finner indekser (posisjoner) til punkter der r3<r0
    r3(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    ind=find(r4<r0);     % finner indekser (posisjoner) til punkter der r4<r0
    r4(ind) = r0;        % og setter disse til r0 for aa ungaa divisjon med 0
    
    V1 = Ke*q1./r1;                 % beregner potensiale til partikkel 1
    V2 = Ke*q2./r2;                 % beregner potensiale til partikkel 2
    V3 = Ke*q3./r3;                 % beregner potensiale til partikkel 3
    V4 = Ke*q4./r4;                 % beregner potensiale til partikkel 4
    Total_V = V1 + V2 + V3 + V4;    % beregner den totale potensiale
    
    Ex1 = Ke*q1*(X-Xp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på x-aksen
    Ez1 = Ke*q1*(Z-Zp1)./(r1.^3);   % beregner elektrisk felt til partikkel 1 på z-aksen
    
    Ex2 = Ke*q2*(X-Xp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på x-aksen
    Ez2 = Ke*q2*(Z-Zp2)./(r2.^3);   % beregner elektrisk felt til partikkel 2 på z-aksen
    
    Ex3 = Ke*q3*(X-Xp3)./(r3.^3);   % beregner elektrisk felt til partikkel 3 på x-aksen
    Ez3 = Ke*q3*(Z-Zp3)./(r3.^3);   % beregner elektrisk felt til partikkel 3 på z-aksen
    
    Ex4 = Ke*q4*(X-Xp4)./(r4.^3);   % beregner elektrisk felt til partikkel 4 på x-aksen
    Ez4 = Ke*q4*(Z-Zp4)./(r4.^3);   % beregner elektrisk felt til partikkel 4 på z-aksen
    
    Total_Ex = Ex1 + Ex2 + Ex3 + Ex4;    % beregner totale elektrisk felt på x-aksen
    Total_Ez = Ez1 + Ez2 + Ez3 + Ez4;    % beregner totale elektrisk felt på z-aksen
    
% plotter felter i xz-planet

    figure(1)
    subplot(1,1,1)
    
    contour(X,Z,Total_V);      % viser potensialet som et contourplot
    hold on
    quiver(X,Z,Total_Ex,Total_Ez)    % viser det elektriske feltet som et vektorplot 
    xlabel('x-aksen')
    ylabel('z-aksen')
    hold off
    axis equal
    
 %% oppgave 1 e II
 
    plot(Z,Total_V)            % Plotter Coulombfeltet for V som funksjon av z langs z-aksen
    xlabel('z-aksen')
    ylabel('Elektrisk potensiale V')
    
 %%
 H = 1./abs(Z).^3;
 
 plot(Z,H)
 xlabel('z-aksen')
 ylabel('1/|z|^3')
 
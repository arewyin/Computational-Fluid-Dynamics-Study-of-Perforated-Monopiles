function [DifZLoca, F_Time_Series] = MFCentroid(H,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% inputs
h = 35/80;                                   % Water depth
rho = 998;                                % density (kg/m^3)
g = 9.80665;                               % gravity acceleration
dtheta = 5;                                   % d angle (degree)
a = 5/80;                                   % Radius
ka = 2*pi()/L*a;                           % ka
Nb = 100;                                   % Number of bessel function calculation 

theta = (0:dtheta:360-dtheta)/360*(2*pi);  % angle in radian
dz = h/500;
zloca = (H/2:-dz:-h);                     % submergence depth

kw = ka/a;                                 % Wave number k
w = sqrt(g*kw.*tanh(kw*h));                % Wave Frequency
T = 2*pi./w;                               % Wave period T
div = 100;                                  % Time division
dt = (2*pi./w)/div;                        % dt

t = (0:dt:2*T);
t = transpose(t);

Nz = size(zloca,2);
Nt = size(t,1);
Ntheta = size(theta,2);

%% Functions
% Definition of beta_m
bm = zeros(1,Nb);
for m = 0:Nb
    if m == 0
        bm(m+1) = 1;
    else 
        bm(m+1) = 2*1i^m;
    end
end

% Bessel function and Hankel function and its derivative
% [1] https://en.wikipedia.org/wiki/Bessel_function
% [2] https://math.stackexchange.com/questions/2204475/derivative-of-bessel
% -function-of-second-kind-zero-order
% [3] https://sharif.edu/~aborji/25735/files/bessel%20functions.pdf

for i = 1:Nb+1
    Jm(i,:) = besselj(i-1,ka);  % Bessel function first kind  => ~ e^-ikx ~ cos(kx)
    Ym(i,:) = bessely(i-1,ka);  % Bessel function second kind => ~ +-i*e^-ikx ~ sin(kx)
    Hm(i,:) = besselh(i-1,ka);  % Hankel function first kind  => ~ Jm + i*Ym ~ cos(kx)+sin(kx)
end
% DJm = derivatie of Bessel function first kind Jm'(x) 
% DYm = derivatie of Bessel function second kind Ym'(x)
% DHm = derivate of Hankel function of first kind Hm'(x)
for i = 1:Nb
    if i == 1 % zeroth
        DJm(i,:) = -Jm(i+1,:);
        DYm(i,:) = -Ym(i+1,:);
    else      % higher
        DJm(i,:) = 0.5*(Jm(i-1,:)-Jm(i+1,:));
        DYm(i,:) = 0.5*(Ym(i-1,:)-Ym(i+1,:));
    end
    DHm(i,:) = DJm(i,:)+1i*DYm(i,:);

end

%% Incident wave profile 
% calculation incident wave profile ita0 around the angle. 
itai0 = zeros(Ntheta,1);

for i = 1:Ntheta
    temp = 0;
    for m = 1:Nb
        temp = temp + H/2*(bm(m)*Jm(m,:)*cos((m-1)*theta(i)));
    end
    itai0(i,:) = temp; % Time independent incident wave elevation
end

itai = zeros(Nt,Ntheta);
for j = 1:Ntheta
    for k = 1:Nt
        itai(k,j) = itai0(j)*exp(-1i*w(1)*t(k));  % Time history of incident wave elevation
    end
end

%% Force Calc
alpha = atan(DJm(2)/DYm(2));

F = zeros(Nt,Nz);
for i = 1:Nt
    for ii = 1:Nz
        if zloca(ii)<=real(itai(i,1))
            temp = 2*rho*g*H/kw*cosh(kw*(zloca(ii)+h))/cosh(kw*h)/...
                (sqrt(DJm(2)^2+DYm(2)^2))*cos(w*t(i)-alpha);
            F(i,ii) = temp*dz;
        end
    end
end

for i = 1:Nt
    F_Time_Series(i) = sum(F(i,:),"all");
end

F_Max = max(F_Time_Series,[],'all');
%% Find Center of force for wave hitting monopile at 0 deg
for i = 1:Nt
    for k = 1:Nz
        before = sum(F(i,1:k));
        m = k+1;
        after = sum(F(i,m:end));
        difference(i,k) = before - after;
    end
end
[MinDif, IndexMinDif] = min(abs(difference),[],2);
DifZLoca = zloca(IndexMinDif);

% figure()
% plot(t,DifZLoca)
% title('Center of Force vs. Time')
% xlabel('Time(s)')
% ylabel('Center of Force from MWL(m)')
end
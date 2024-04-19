function [L0,L,C0,C,umax,wmax,umaxaccel] = ...
    WaveSolver(Period, Depth, Height)
%WaveSolver
%   Input Variables: Period(s), Height(m)
%   Output Variables: L0 (Array deep water wavelength), L (Wavelength), C0
%   (Array of deep water celerity), C (Array of celerity), n0 (Array of
%   deep water group numbers), n (Array of group numbers), Cg0 (Array of
%   deep water celerity), Cg (Array of water celerity), E0 (Array of deep
%   wave energy), E (Array of wave energy), EF (Array of engergy flux),
%   umax (Array of maximum horizontal velocity), wmax (Array of maximum 
%   vertical velocity). 

%MSL
AvgDepth = Depth;



%Loops through each wave to determine L0, L, C0, C, n0, n, Cg0, Cg, E0, E,
%   EF, umax, wmax
for i = 1:length(Period)
    %Determinine wavelength for current wave
    L0(i) = 9.81*(Period(i))^2/(2*pi());
    L(i) = L0(i)*tanh(2*pi()*AvgDepth/L0(i));
    Lt(1) = L0(i);
    j=2;
    while abs(Lt-L(i)) > 0.001
        Lt(j) = L(i);
        L(i) = L0(i)*tanh(2*pi()*AvgDepth/Lt(j));
        j=j+1;
    end

    
    %Determining wave celerity
     C0(i) = sqrt(9.81*L(i)/(2*pi()));
     C(i) = sqrt(9.81*L(i)/(2*pi())*tanh(2*pi()*AvgDepth/L(i)));
% 
%     %Determining group number
%     n(i) = 1/2*(1+4*pi()*AvgDepth/(L(i))/sinh(4*pi()*AvgDepth/L(i)));
%     n0(i) = 1/2*(1+4*pi()*AvgDepth/(L0(i))/sinh(4*pi()*AvgDepth/L0(i)));
% 
%     %Determining celerity
%     Cg(i) = C(i)*n(i);
%     Cg0(i) = C0(i)*n0(i);
% 
%     %Determining wave energy
%     E(i) = 1/8*(1026)*9.81*(Height(i));
%     E0(i) = E(i)*Cg(i)/Cg0(i);
% 
%     %Energy Flux
%     EF(i) = E(i)*Cg(i);
% 
    %Velocities
    Sigma = 2*pi()/(Period(i));
    k = 2*pi()/L(i);
    a = Height(i)/2;
    umax(i,1) = a*Sigma*cosh(k*(AvgDepth+Height(i,1)))/sinh(k*AvgDepth)*...
        cos(k*C(i)*4/(Period(i,1))-Sigma*(Period(i,1))/4);
    wmax(i,1) = a*Sigma*sinh(k*(AvgDepth))/sinh(k*AvgDepth)*...
        sin(k*C(i)*2/(Period(i,1))-Sigma*(Period(i,1))/2);
    umaxaccel(i,1) = a*Sigma*cosh(k*(AvgDepth+Height(i,1)))/...
        sinh(k*AvgDepth)*k*sin(k*C(i)*4/(Period(i,1))-...
        Sigma*(Period(i,1))/4);
end



function [Waves, Height, Period, AvgH, AvgT, Hrms, Hs, Ts, psd1,...
    Tm01, Tm02, Hm0, f_hobo, Tp, AvgH_Er, AvgT_Er, Ts_Er] = ...
    ZeroUp(DepthD,Time)
%ZeroUp Analysis
%   Input Variables: Trimmed detrended depth, Time
%   Output Variables: Waves (Index of wave locations), Height, Period, 
%       Length, AvgH (Average Wave Height), AvgT (Average Period), Hrms 
%       (Root mean sum of wave heights), Hs (Average of top 1/3 of
%       heights), Ts (Average of top 1/3 of wave periods), psd (Power
%       Spectral Density), Tm01 (Zero-crossing period), 
%       Tm02 (Peak Period), Hm0 (average wave height from
%       0th spectral moment), f_hobo (Array of frequency), Tp (Peak period),
%       AvgH_Er (Error between average heights), AvgT_Er (Error between
%       average periods), Ts_Er (Error between averages of top 1/3 periods)


%Deterimining all values above 0
surzero = find(DepthD > 0);

%% Check for Downcrossing
n=0;
Waves = [];
%Determine limit of for loop
if surzero(length(surzero)) == length(DepthD)
    limit = length(surzero)-1;
else
    limit = length(surzero);
end

%For loop
for l = 1:limit
    if DepthD(surzero(l)+1) <= 0 %Check if next value is below zero
        n=n+1;
        Waves(n,1) = surzero(l) + 1;
    end
end

%% If height is empty search for waves by upcrossings
if isempty(Waves) == 1 | length(Waves) ==1
    %Deterimining all values below 0
    subzero = find(DepthD < 0);

    %Check for Upcrossing
    n=0;
    %Determine limit of for loop
    if length(subzero) == length(DepthD)
        limit = length(subzero)-1;
    else
        limit = length(subzero);
    end 

    for l = 1:limit
        if DepthD(subzero(l)+1) >= 0 %Check if next value is above zero
            n=n+1;
            Waves(n,1) = subzero(l) + 1;
        end
    end
end

%% Finding Wave Heights for Experimental Data
for l = 1:length(Waves)-1
    Hobo_Sec = DepthD(Waves(l):Waves(l+1));
    Sub = find(Hobo_Sec < 0);
    Over = find(Hobo_Sec >= 0);
    Height(l,1) = max(Hobo_Sec(Over)) - min(Hobo_Sec(Sub));
    Period(l,1) = Time(Waves(l+1)) - Time(Waves(l));
end


%% Determining wave variables
AvgH = sum(Height)/length(Height);
AvgT = time2num(sum(Period)/length(Period));
Hrms = 0; %(sumsqr(Height)/length(Height))^(1/2);
SortedH = sort(Height,"Descend");
Hs = sum(SortedH(1:round(length(SortedH)/3)))/(length(SortedH)/3);
SortedT = sort(Period,"Descend");
Ts = time2num(sum(SortedT(1:round(length(SortedT)/3)))/(length(SortedT)/3));


%% Generating Frequency Vector
fs_hobo = 1;
N_real = length(DepthD(:,1));           %Length of the time vector
df_hobo = fs_hobo/N_real;                       %Difference of consecutive samples
f_hobo = [0:df_hobo:fs_hobo/2];                 %Frequency vector
f_hobo = f_hobo.';                              %Transposing

%% Fast Fourier Transform
HoboFFT = fft(DepthD);
HoboFFT = HoboFFT(1:length(f_hobo(:,1)));

psd2 = (1/(N_real*fs_hobo))*(abs(HoboFFT)).^2;  %Two-sided
psd1 = psd2;                                %One-sided
psd1(2:end-1,1)=2.*psd1(2:end-1,1);

%% Calculating spectral moments
m0 = sum(psd1.*f_hobo.^0.*df_hobo); %Zero-moment
m1 = sum(psd1.*f_hobo.^1.*df_hobo); %First-moment
m2 = sum(psd1.*f_hobo.^2.*df_hobo); %Second-moment

%% Wave Height
Hm0 = 4*sqrt(m0);

%% Calculating wave periods
Tm01 = m0/m1;       %mean period
Tm02 = (m0/m2)^0.5; %zero-crossing period

%% Calculating peak period
[psdMax MaxIndex] = max(psd1(:,1));         %Locating spectral peak
Fp = f_hobo(MaxIndex,1);                    %Peak frequency
Tp = 1/Fp;                                  %Peak period

%% Errors
AvgH_Er = abs(AvgH-Hm0)/AvgH*100;
AvgT_Er = abs(AvgT-Tm01)/AvgT*100;
Ts_Er = abs(Ts-Tp)/Ts*100;

end
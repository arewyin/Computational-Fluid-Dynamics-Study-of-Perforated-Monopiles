function [Waves] = StrainZeroCrossings(Strain)
%% StrainZeroCrossings Analysis
%   Checks for downcrossing then upcrossings
%   Input Variables: Zeroed Strain
%   Output Variables: Waves (Index of wave locations)


%Deterimining all values above 0
surzero = find(Strain > 0);

%% Check for Downcrossing
n=0;
Waves = [];
%Determine limit of for loop
if surzero(length(surzero)) == length(Strain)
    limit = length(surzero)-1;
else
    limit = length(surzero);
end

%For loop
for l = 1:limit
    if Strain(surzero(l)+1) <= 0 %Check if next value is below zero
        n=n+1;
        Waves(n,1) = surzero(l) + 1;
    end
end

%% If height is empty search for waves by upcrossings
if isempty(Waves) == 1 | length(Waves) ==1
    %Deterimining all values below 0
    subzero = find(Strain < 0);

    %Check for Upcrossing
    n=0;
    %Determine limit of for loop
    if length(subzero) == length(Strain)
        limit = length(subzero)-1;
    else
        limit = length(subzero);
    end 

    for l = 1:limit
        if Strain(subzero(l)+1) >= 0 %Check if next value is above zero
            n=n+1;
            Waves(n,1) = subzero(l) + 1;
        end
    end
end
end
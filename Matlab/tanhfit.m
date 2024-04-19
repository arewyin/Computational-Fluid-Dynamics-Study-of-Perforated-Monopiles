function [s] = tanhfit(x,y)
% Sine fit applies  a sinusoidal funtion to a set of x and y coordinates
% and returns the parameters of the sinusoidal function to be used in
% another script.
%% INPUTS:
%   x = Array of x data points
%   y = Array of y data points
%% OUTPUTS:
%   s = Array of parameters

%% Applying a hyperbolic tangent function
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
fit = @(b,x)  b(1).*tanh(b(2).*x);    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr;  per;  -1;  ym]);                      % Minimise Least-Squares
xp = linspace(min(x),max(x));

%% Plot
% figure(1)
% plot(x,y, xp, fit(s,xp),'r')
% grid
% close(1)
end
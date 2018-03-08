function [t, y] = eulerode(fh, t, y0)
% eulerode: uses Euler's method to integrate an ODE
% Input:

% Auxiliar function to compute y_dot
f = str2func(fh);

lt = length(t);
ly = length(y0);

y = zeros(ly, lt);

% Insert itinial condition
y(:, 1) = y0;  

% Iterate over time vector
for i=1:lt-1
    % Time step
    h = t(i+1) - t(i);
    
    % Eurler Formula
    y(:, i+1) = y(:, i) + f(t(i),y(:, i)) .* h;
end

y = y';

end
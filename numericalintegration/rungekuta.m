function [t, y] = rungekuta(fh, t, y0)

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
    
    % Runge-Kutta Formulae
    
    k1 = f(t(i),y(:, i));
    
    k2 = f(t(i) + h/2, y(:, i) + k1/2 .* h) ;
    
    k3 = f(t(i) + h/2, y(:, i) + k2/2 .* h);
    
    k4 = f(t(i) + h/2, y(:, i) + k3 .* h);
    
    y(:, i+1) = y(:, i) + (k1 + 2*k2 + 2*k3 + k4) .* h/6;
end

y = y';

end
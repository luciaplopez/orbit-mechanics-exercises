function yp = yprime(t, y)
% Input: position and velocity (r, v)
% Output: derivative for both vectors (rdot, dot)

ae = 6378 * 1e3;
u = 3.986004418e+14;

% Create the vector
yp = zeros(6, 1);

% Derivative of the position
yp(1:3, 1) = y(4:6, 1);

r = y(1:3, 1);
R = norm(r, 2);

% Derivative of the velocity
yp(4:6, 1) = -(u/R^3) * r;

end


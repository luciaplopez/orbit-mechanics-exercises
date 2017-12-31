function yp = yprime(t, y, u)
% Input: position and velocity (r, v)
% Output: derivative for both vectors (rdot, dot)

yp = zeros(6, 1);

% Derivative of the position
yp(1:3, 1) = y(4:6, 1);

R = norm(y(1:3, 1), 2);

% Derivative of the velocity
yp(4:6, 1) = -(u/R^3) * y(1:3, 1);

end


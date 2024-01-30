function [u_xt] = heat_sol_brute(t, x, dx, dt, D)
% The solution to the heat equation in point t and x
% u(t,x) = int_0^t int_R g(s, y, t, x) dW(s,y)
% Approximated by "brute force" with definition of the "Walsh-esque" 
% multivariate stochastic integral w.r.t Brownian sheet W(s,y)

% Input:
% t and x point to be evaluated
% dx and dt step size
% D diameter of space, "should be" = inf since integral over R for space
mu_A = dx * dt;

g = @(s, y) exp(-((x - y).^2) ./ (4*(t - s))) ./ ((4*pi*(t - s)).^(1/2));

s_nbr_points = t/(dt) + 1;
s_points = linspace(0, t, s_nbr_points);
% + 1 makear sense, rita upp grid så får du se
% kom ihåg vi vill evaluera i övre högra punkten
y_nbr_points = 2*D/(dx) + 1;
y_points = linspace(-D,D,y_nbr_points);

u_xt = 0;
% t = 0 är u = 0 så börja från 2 och gå inte till slutet för då t - s = 0
for i = 2:(s_nbr_points - 1)
    % j = 2 för vi tar i övre högra hörnet
    for j = 2:y_nbr_points
        s_point = s_points(i);
        y_point = y_points(j);
        u_xt = u_xt + g(s_point, y_point) * sqrt(mu_A) * normrnd(0, 1);
    end
end
end
function [u_t, t_points] = u_time_simulation(x_point, dy, dt, D, T, theta)
mu_A = dt * dy;
G = @(t, y, s) exp(-((x_point - y).^2) ./ (4*theta*(t - s))) ./(2*sqrt(theta) * ((pi*(t - s)).^(1/2)) );
t_nbr_points = T/dt;
t_points = linspace(0, T, t_nbr_points);
y_nbr_points = 2*D/(dy);
y_points = linspace(-D,D,y_nbr_points);
white_noise = sqrt(mu_A)*normrnd(zeros(t_nbr_points, y_nbr_points), 1);
u_t = zeros(1, t_nbr_points);

for j = 2:t_nbr_points
    sum = 0;
    j
    for k = 1:y_nbr_points
        sum = sum + ...
        G(t_points(j), y_points(k), t_points(j-1)) * ...
        white_noise(j,k);
    end
    u_t(j) = u_t(j-1) + sum;
end
end

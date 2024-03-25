%% Estimation of the drift parameter SHE

%% Path of time for a given x
dy = 1/1000;
dt = 1/5000; 
D = 20; %diameter of integration
x_point = 1; 
T = 1;
theta = 0.5;
t_nbr_points = T/dt;

[u_time, t_points] = u_time_simulation(x_point, dy, dt, D, T, theta);
%% Plots 
close all
plot(t_points, u_time)

%% Drift parameter estimation
alpha = 2;
c1alpha = (1/(2*pi*(alpha - 1)))*gamma(1/alpha);
c2alpha = sqrt(c1alpha * 2^(1-(1/alpha))); % Kvadreras eller ej?? Mycket typos och slarv i artikeln
mu_4 = 3; % tror jag
factor = c2alpha^4 * 2 * mu_4 * T;% antar alpha = 2 ekv 20 https://arxiv.org/pdf/1912.07917.pdf
sum = 0;
for j = 1:t_nbr_points - 1
    sum = sum + (u_time(j + 1) - u_time(j))^4;
end
theta = 1/(sum/factor)

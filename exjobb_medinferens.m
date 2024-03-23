

n = 10;
T = 1;
D = 8;
t = linspace(0,T,1000);
x = linspace(-D,D,1000);
%sizebs = size(t,2)*size(x,2);
bs = zeros(size(t,2), size(x,2));
%% Brownian Sheet pointwise
% Skapar pointwise lakan m.h.a KL-expansionWN = zeros(size(t,2), size(x,2));
for i = 1:size(t,2)
    i
    for j = 1:size(x,2)
        bs(i,j) = wn(t(i), x(j), n, T, D);
    end
end
%% White Noise approximation (term-wise differentiation of BS-KL)
%sizebs = size(t,2)*size(x,2);
WN = zeros(size(t,2), size(x,2));
for i = 1:size(t,2)
    i
    for j = 1:size(x,2)
        WN(i,j) = w_hat(t(i), x(j), n, T, D);
    end
end

%%

h = surf(t,x , bs)
set(h,'LineStyle','none')
%a = load('WN_100x100_100_1_1.mat');

%% Stokastisk integral brute force (riktig seg och dålig)
dx = 1/20;
dt = 1/20;
D = 100;
T = 1;
x_points = linspace(-D, D, 100);
t_points = linspace(0, T, 100);

u = zeros(size(t_points,2), size(x_points,2));

for i = 2:size(t_points,2)
    i
    for j = 1:size(x_points,2)
        a = heat_sol_brute(t_points(i), x_points(j), dx, dt, D);
        u(i,j) = a;
    end
end

%% plotting
figure
surf(u)

%% Euler-Maruyama
dx = 1/1000;
dt = 1/1000; % enl Walsh
mu_A = dt * dx;
D = 8; %diameter of integration
X = 1; % x limits
T = 1;
G = @(x, t, y, s) exp(-((x - y).^2) ./ (2*(t - s))) ./ ((2*pi*(t - s)).^(1/2));
t_nbr_points = round(T/(dt));
t_points = linspace(0, T, t_nbr_points);
y_nbr_points = 2*D/(dx);
y_points = linspace(-D,D,y_nbr_points);
x_nbr_points = 2*X/(dx);
x_points = linspace(-X,X,x_nbr_points);
%%
white_noise = sqrt(mu_A)*normrnd(zeros(t_nbr_points, y_nbr_points), 1);
u = zeros(x_nbr_points, t_nbr_points);
%%

for i = 1:x_nbr_points
    i
    for j = 2:t_nbr_points
        sum = 0;
        for k = 1:y_nbr_points
            sum = sum + ...
            G(x_points(i), t_points(j), y_points(k), t_points(j-1)) * ...
            white_noise(j,k);
        end
        u(i,j) = u(i, j-1) + sum;
    end
end
%% only one x point simulation
i = round(x_nbr_points/2);
i =2;
for j = 2:t_nbr_points
    j
    sum = 0;
    for k = 1:y_nbr_points
        sum = sum + ...
        G(x_points(i), t_points(j), y_points(k), t_points(j-1)) * ...
        white_noise(j,k);
    end
    u(i,j) = u(i, j-1) + sum;
end
%%
alpha = 2;
c1alpha = (1/(2*pi*(alpha - 1)))*gamma(1/alpha);
c2alpha = sqrt( c1alpha * 2^(1-(1/alpha)));
mu_4 = 3; % tror jag
factor = c2alpha^4 * 2 * mu_4 * T;% antar alpha = 2 ekv 20 https://arxiv.org/pdf/1912.07917.pdf
sum = 0;
x_point_index = i;
u_t = u(x_point_index, :);
for j = 2:t_nbr_points - 1
    sum = sum + (u_t(j + 1) - u(j))^4;
end
theta = 1 / (sum / factor)

%% Plotting
close all
figure;
set(gcf, 'Color', 'w');  

%v = VideoWriter('temperature_over_time.avi');
%open(v);
for t = 1:size(u, 1)
    p = plot(x_points, u(t, :), 'LineWidth', 2);
    p.Color = [1, 0.592, 0];  
    
    title(['Time = ', num2str(t_points(t))], 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('Space', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Temperature', 'Interpreter', 'latex', 'FontSize', 14);
    ylim([-0.04 0.04])
    set(gca, 'FontSize', 12); 
    set(gca, 'LineWidth', 1.5);
    box on;  
    drawnow; 
    pause(0.1);
%    frame = getframe(gcf);  % Capture the plot as a frame
%    writeVideo(v, frame);
end
%close(v);  % Close the video file
%% 
close all
figure
h = surf(t_points, x_points, u)
set(h,'LineStyle','none')
title('$u(t,x) = \int_0^t \int_{R}\frac{e^{-\frac{|x-y|^2}{2(t-s)}}}{(2\pi|t-s|)^{1/2}}W(dyds)$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('t', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('x', 'Interpreter', 'latex', 'FontSize', 14);

%% Plotting
close all
figure;
set(gcf, 'Color', 'w');  

%v = VideoWriter('temperature_over_time.avi');
%open(v);
for x = 1:size(u, 2)
    p = plot(t_points, u(:, x), 'LineWidth', 2);
    p.Color = [1, 0.592, 0];  
    
    title(['Space = ', num2str(x_points(x))], 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('Time', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Temperature', 'Interpreter', 'latex', 'FontSize', 14);
    ylim([-4 4])
    set(gca, 'FontSize', 12); 
    set(gca, 'LineWidth', 1.5);
    box on;  
    drawnow; 
    pause(0.1);
%    frame = getframe(gcf);  % Capture the plot as a frame
%    writeVideo(v, frame);
end
%close(v);  % Close the video file


%% covariance

CovM = @(t,s,x,y) (1/sqrt(2*pi))* (sqrt(t + s)*exp((abs(x-y)^2)/(2*(t+s))) ...
- sqrt(abs(t - s))*exp(-(abs(x-y)^2)/(2*abs(t-s)))) ...
+ (x-y)*(cdf('Normal', (x-y)/sqrt(s+t), 0, 1) - cdf('Normal', (x-y)/sqrt(abs(t-s)), 0, 1));

Cov = @(t,s,tau) (1/sqrt(2*pi)) .* (sqrt(t + s).*exp(-(abs(tau).^2)/(2.*(t+s))) ...
- sqrt(abs(t - s)).*exp(-(abs(tau).^2)/(2.*abs(tau)))) ...
+ (tau).*(cdf('Normal', (tau)./sqrt(s+t), 0, 1) - cdf('Normal', (tau)./sqrt(abs(t-s)), 0, 1));

t_p = 0;
s_p = 1;
D = 10;
tau_points = linspace(-D, D, 10000);
figure
plot(tau_points, Cov(t_p,s_p, tau_points))

%% white-coloured noise
dx = 1/10;
dt = 1/10;
mu_A = dt * dx;
D = 5; %diamter av integration
X = 1; % x gränser för utvärdering
T = 1;
G = @(t, s, x, y) exp(-((x - y).^2) ./ (4*(t - s))) ./ ((4*pi*(t - s)).^(1/2));
alpha = 3;
d = 1;
f = @(x,y) 2.^(d-alpha).*pi.^(d/2).*(gamma((d-alpha)/2))/(gamma(alpha/2)) .* abs(x-y).^(-d+alpha);
t_nbr_points = T/(dt);
t_points = linspace(0, T, t_nbr_points);
% + 1 kanske inte behövs
y_nbr_points = 2*D/(dx) + 1;
y_points = linspace(-D,D,y_nbr_points);
x_nbr_points = 2*X/(dx) + 1;
x_points = linspace(-X,X,x_nbr_points);

%% calculating variances of coloured noise
coloured_variances = zeros(t_nbr_points, y_nbr_points);
for i = 2:t_nbr_points
    for j = 1:(y_nbr_points -1)
        coloured_variances(i,j) = t_points(i)*integral2(f, y_points(j), y_points(j + 1), y_points(j), y_points(j + 1));
    end
end

%% Theta estimation test (time)

%load("u_dx150_dt200_D8_X1_T10.mat")


%% Theta estimation test (space)

%load("u_dx150_dt200_D8_X1_T10.mat")
alpha = 2;
c1alpha = (1/(2*pi*(alpha - 1)))*gamma(1/alpha);
c2alpha = sqrt( c1alpha * 2^(1-(1/alpha)));
mu_4 = 3; % tror jag
factor = c2alpha^4 * 2 * mu_4 * T;% antar alpha = 2 ekv 20 https://arxiv.org/pdf/1912.07917.pdf
sum = 0;
x_point_index = 1000;
u_t = u(x_point_index, :);
for j = 2:t_nbr_points - 1
    sum = sum + (u_t(j + 1) - u(j))^4;
end
theta = 1 / (sum / factor)

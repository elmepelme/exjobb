

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
dx = 1/150;
dt = 1/200;
mu_A = dt * dx;
D = 8; %diamter av integration
X = 1; % x gränser för utvärdering
T = 10;
G = @(t, s, x, y) exp(-((x - y).^2) ./ (4*(t - s))) ./ ((4*pi*(t - s)).^(1/2));
t_nbr_points = T/(dt);
t_points = linspace(0, T, t_nbr_points);
% + 1 kanske inte behövs
y_nbr_points = 2*D/(dx) + 1;
y_points = linspace(-D,D,y_nbr_points);
x_nbr_points = 2*X/(dx) + 1;
x_points = linspace(-X,X,x_nbr_points);
%%
white_noise_A = normrnd(zeros(t_nbr_points, y_nbr_points), 1);
%%
u = zeros(t_nbr_points, x_nbr_points);
% vi evaluerar funktionen i nedre högra punkten pga singulariteten
for i = 1:x_nbr_points
    i
    for j = 2:t_nbr_points
        sum = 0;
        for k = 1:y_nbr_points
            sum = sum + ... 
                G(t_points(j), t_points(j-1), x_points(i), y_points(k)) * ...
                sqrt(dx*dt)*white_noise_A(j,k);
        end
        u(j,i) = u(j-1, i) + sum;
    end
end


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
h = surf(x_points, t_points, u)
set(h,'LineStyle','none')
title('$u(t,x) = \int_0^t \int_{R}\frac{e^{-\frac{|x-y|^2}{4(t-s)}}}{(4\pi|t-s|)^{n/2}}W(dyds)$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('x', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('t', 'Interpreter', 'latex', 'FontSize', 14);

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

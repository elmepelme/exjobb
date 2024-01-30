%% Brownian Sheet pointwise
% Skapar pointwise lakan m.h.a KL-expansion

n = 10;
T = 1;
D = 1;
t = linspace(0,T,1000);
x = linspace(0,D,1000);
%sizebs = size(t,2)*size(x,2);
bs = zeros(size(t,2), size(x,2));
for i = 1:size(t,2)
    i
    for j = 1:size(x,2)
        bs(i,j) = wn(t(i), x(j), n, T, D);
    end
end


%%
figure
surf(bs)
a = load('bs_100x100_100_1_1.mat');

%% Stokastisk integral brute force
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


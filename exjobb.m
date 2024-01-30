%% Simulera lösning

n = 100;
T = 1;
D = 1;
t = linspace(0,T,100);
x = linspace(0,D,100);
%sizebs = size(t,2)*size(x,2);
bs = zeros(size(t,2), size(x,2));
for i = 1:size(t,2)
    for j = 1:size(x,2)
        i/size(t,2)
        bs(i,j) = wn(t(i), x(j), n, T, D);
    end
end


%%
figure
surf(bs)
a = load('bs_100x100_100_1_1.mat')
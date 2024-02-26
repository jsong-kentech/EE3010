%% Random walk 
% 2024.02.25
% Ju Song
% KENTECH

% Contents
% (1) Config 
% (2) 2D Random walk (to visualize)
% (3) 1D Random walk (to analyze)

clear; clc; close all
%% Config
dx = 1;
dt = 1;

N = 10^3; % number of random walkers
t_end = 10000; % time period to simulate
t_vec = 0:dt:t_end;
M = length(t_vec);


%% 1D Random walk

% initialize
x_mat = zeros(N,M);

% random walks
% time step loop
for m = 1:M
    
    % generate random walk for N particles
    dx_now = 2*unidrnd(2,N,1)-3;
    
    % calculate position
    if m ==1
       x_mat(:,m) = zeros(N,1);
    else
        x_mat(:,m) = x_mat(:,m-1) + dx_now;
    end
end


% analyze: root mean distance
figure(5)
plot(t_vec,x_mat)
xlabel('time')
ylabel('x-position')

d_mat = sqrt(x_mat.^2);
d_avg = mean(d_mat,1);
figure(6) 
plot(t_vec,d_avg)
xlabel('time')
ylabel('avg distance')

% analyze: flux
    % define concentration
edges = linspace(-2*d_avg(end),2*d_avg(end),21);
positions = (edges(2:end) + edges(1:end-1))/2;
c_mat = zeros(length(edges)-1,M);

for m = 1:M

    c_mat(:,m) = histcounts(x_mat(:,m),edges);

end


figure(7)
for m = 1:1000:M
    
    plot(positions',c_mat(:,m)); hold on

end


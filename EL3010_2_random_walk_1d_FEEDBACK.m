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
    dx_now = dx*(2*unidrnd(2,N,1)-3);
    
    % calculate position
    if m ==1
       x_mat(:,m) = zeros(N,1);
    else
        x_mat(:,m) = x_mat(:,m-1) + dx_now;
    end
end


% analyze: average distance
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

figure(7) 
plot(t_vec.^0.5,d_avg)
xlabel('time^{0.5}')
ylabel('avg distance')


% 


% analyze: flux - 
    % define concentration
x_vec = -2*round(d_avg(end)):10:2*round(d_avg(end));
c_mat = zeros(length(x_vec),M);
dc_mat = zeros(length(x_vec)-1,M);

for m = 1:M

    c_mat(:,m) = histcounts(x_mat(:,m),[x_vec,x_vec(end)+dx]);
    
    if ismember(m,[1000, 5000, 10000])
        figure(7)
        bar(x_vec,c_mat(:,m),'facealpha',0.6); hold on
    end

    dc_mat(:,m) = diff(c_mat(:,m)); 

end

f_mat = zeros(length(x_vec)-1,M);
for m = 2:M

    for l = 1:length(x_vec)-1
            

        if l == length(x_vec)-1
            f_mat(l,m) = nnz(...
                    x_vec(l) <= x_mat(:,m-1) && x_mat(:,m-1) < x_vec(l+1) ...
                    && x_vec(l+1) <= x_mat(:,m) && x_mat(:,m) < x_vec(l+2)) ...
                     - nnz(...
                     x_vec(l+1) <= x_mat(:,m-1) && x_mat(:,m-1) < x_vec(l+1) ...
                    && x_vec(l+1) <= x_mat(:,m) && x_mat(:,m) < x_vec(l));
        
   

        else
            f_mat(l,m) = nnz(...
                    (x_vec(l) <= x_mat(:,m-1)) .* (x_mat(:,m-1) < x_vec(l+1)) ...
                    && (x_vec(l+1) <= x_mat(:,m)) .* (x_mat(:,m) < x_vec(l+2))) ...
                     - nnz(...
                     x_vec(l+1) <= x_mat(:,m-1) && x_mat(:,m-1) < x_vec(l+1) ...
                    && x_vec(l+1) <= x_mat(:,m) && x_mat(:,m) < x_vec(l));
        
        
        end
    
    end

end




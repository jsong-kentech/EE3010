
clear; clc; close all
% Parameters

dx = 1;
dt = 1;

t_simulate = 10^4;
t_vec = 0:dt:t_simulate;

N = 4; 
M = length(t_vec);


x_mat = zeros(N,M);
y_mat = zeros(N,M);

for m = 2:M
    
    dx_now = dx*2*unidrnd(2,4,1)-3;
    dy_now = dy*2*unidrnd(2,4,1)-3;

    
    x_mat(:,m) = x_mat(:,m-1)+dx_now;
    y_mat(:,m) = y_mat(:,m-1)+dy_now;

end 


figure(1)
for n = 1:N
plot(x_mat(n,:),y_mat(n,:)); hold on
end

c_mat = lines(N);
figure(2)
for m = 1:M

    for n = 1:N
    plot(x_mat(n,1:m),y_mat(n,1:m),'-','color',c_mat(n,:)); hold on
    plot(x_mat(n,m),y_mat(n,m),'o','MarkerFaceColor',c_mat(n,:))
    end
    hold off
    pause(0.1)
end
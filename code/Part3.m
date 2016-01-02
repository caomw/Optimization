% 4 Function minimization in nD

clear all
close all
clc

%% Groundtruth data creation

a = 5;
b = 3;
xc = 1;
yc = 1;
phi = pi/4;
n = 500;
global data;
data = dtEllipse_noisy(a,b,xc,yc,phi,n);
figure;
scatter(data(1,:),data(2,:));
grid on;
title('Ellipse data points');

%% Levenberg Marquardt data fitting

clc

% x = [a; b; xc; yc; phi] -> defines the ellipse

% Initial guess
x0 = [1; 1; 0; 0; 0];
u = 1;
ku = 2^0.25;

% Calculation parameters
precision = 0.001;
counter = 0;
iteration_limit = 1000;
I = eye(5);

% Calculation variables
x = x0;
prev = x0 + 1;

% Approximate centre for faster convergence
x(3) = mean(data(1,:));
x(4) = mean(data(2,:));

while norm(prev - x)>precision && counter<iteration_limit
    F = dtF(x);
    J = dtGrad(@dtF,x,5)';
    H = J'*J;
    h = (H+u*I)\(-J'*F); %%%
    
    if dtF(x+h)>= F
        u = u*ku;
    else
        u = u/ku;
        prev = x;
        x = x+h;
    end
    
    counter = counter+1;
end

x_min = x
y_min = dtF(x)
num_iterations = counter

% Plot results
figure;
scatter(data(1,:),data(2,:));
grid on;
hold on;
result = dtEllipse(x_min);
plot(result(1,:),result(2,:),'r','LineWidth',2);
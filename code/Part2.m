% 3. Function Minimization in 2D

% 3.1 Minimization of a quadratic form

clear all
close all
clc

%% 1 - Arbitrary Line Search

clc

% Directions
d1 = [1 0];
d1 = d1/norm(d1); % Just in case of change in direction
d2 = [0 1];
d2 = d2/norm(d2); % Just in case of change in direction

% Starting point
% x0 = [1000 1000000]; % Absurd point to test stability
x0 = [1.5 1.5];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;

% Calculation variables
x = x0;
prev_1 = x0 + d1;
prev_2 = x0 + d2;
coords = x;

while norm(prev_1 - x)>precision && norm(prev_2 - x)>precision && counter<iteration_limit
    % Computes 2 steps at a time
    prev_2 = x;
    prev_1 = dtBrent(@f2,d1,x);
    coords = [coords; prev_1];
    x = dtBrent(@f2,d2,prev_1);
    coords = [coords; x];
    counter = counter + 1;
end

xp = -2:0.1:2;
yp = -2:0.1:2;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2(x)
num_iterations = counter*2 % Because two computations are done per loop

%% 2 - Steepest Descent Method

clc

% Starting point
% x0 = [1000 1000000]; % Absurd point to test stability
x0 = [1.5 1.5];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;
dx = 0.001;
dx1 = [dx 0];
dx2 = [0 dx];

% Calculation variables
x = x0;
grad = [1 1];
prev = x0 + grad;
coords = x;

while norm(prev - x)>precision && counter<iteration_limit
    prev = x;
    deriv_x1 = (f2(x+dx1)-f2(x))/dx;
    deriv_x2 = (f2(x+dx2)-f2(x))/dx;
    grad = [deriv_x1 deriv_x2];
    grad = grad/norm(grad);
    x = dtBrent(@f2,grad,x);
    coords = [coords; x];
    counter = counter + 1;
end

xp = -2:0.1:2;
yp = -2:0.1:2;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2(x)
num_iterations = counter


%% 3 - Powell's Method

clc

% Starting point
% x0 = [1000 1000000]; % Absurd point to test stability
x0 = [1.5 1.5];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;

% Calculation variables
P0 = x0;
P1 = x0;
P2 = x0;
u1 = [1 0];
u2 = [0 1];
prev = x0 + u1;
coords = P0;

while norm(prev - P0)>precision && counter<iteration_limit
    prev = P0;
    
    P1 = dtBrent(@f2,u1,P0);
    coords = [coords; P1];
    P2 = dtBrent(@f2,u2,P1);
    coords = [coords; P2];
    u1 = u2;
    u2 = P2-P0;
    u2 = u2/norm(u2);
    P0 = dtBrent(@f2,u2,P2);
    coords = [coords; P0];
    
    counter = counter + 1;
end

xp = -2:0.1:2;
yp = -2:0.1:2;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2(x)
num_iterations = counter


%% 4 - Conjugate Gradients
% Shewchuk - B4. Nonlinear Conjugate Gradients with Newton-Raphson and
% Fletcher-Reeves

clc

% Starting point
% x0 = [1000; 1000000]; % Absurd point to test stability
x0 = [1.5; 1.5];

% Calculation parameters
precision = 0.001;
counter = 0;
iteration_limit = 100;

% Calculation variables
x = x0;
prev = x+1;

coords = x;

% Begin
k = 0;
r = -dtGrad(@f2,x,2);
d = r;
delta = norm(r)^2;

while norm(prev - x)>precision && counter<iteration_limit
    counter = counter + 1;
    prev = x;
    
    j = 0;
    delta_d = norm(d)^2;
    alpha = 1;
    while j<100 && (alpha^2)*delta_d > precision^2
        grad = dtGrad(@f2,x,2);
        H = dtHess2d(@f2,x);
        alpha = -(grad'*d)/(d'*H*d);
        x = x + alpha*d;
        j = j+1;
    end
    r = -dtGrad(@f2,x,2);
    delta_old = delta;
    delta = norm(r)^2;
    beta = delta/delta_old;
    d = r + beta*d;
    k = k+1;
    if k == 2 || r'*d<=0
        d = r;
        k = 0;
    end
    
    coords = [coords x];
end

coords = coords';
xp = -2:0.1:2;
yp = -2:0.1:2;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2(x)
num_iterations = counter

%% 3.2 Minimization of a non quadratic form

clear all
close all
clc

%% 1 - Arbitrary line search

clc

% Directions
d1 = [1 0];
d1 = d1/norm(d1); % Just in case of change in direction
d2 = [0 1];
d2 = d2/norm(d2); % Just in case of change in direction

% Starting point
x0 = [-2 -2];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;

% Calculation variables
x = x0;
prev_1 = x0 + d1;
prev_2 = x0 + d2;
coords = x;

while (norm(prev_1 - x)>precision || norm(prev_2 - x)>precision) && counter<iteration_limit
    % Computes 2 steps at a time
    prev_2 = x;
    prev_1 = dtBrent(@f2_2,d1,x);
    coords = [coords; prev_1];
    x = dtBrent(@f2_2,d2,prev_1);
    coords = [coords; x];
    counter = counter + 1;
end

xp = -2:0.1:0.2;
yp = -2:0.1:-0.8;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2_2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2_2(x)
num_iterations = counter*2 % Because two computations are done per loop


%% 2 - Steepest descent method

clc

% Starting point
x0 = [-1.5 -1.5];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;
grad_multiplier = 5;
dx = 0.001;
dx1 = [dx 0];
dx2 = [0 dx];

% Calculation variables
x = x0;
grad = [1 1];
prev = x0 + grad;
coords = x;

while norm(prev - x)>precision && counter<iteration_limit
    prev = x;
    deriv_x1 = (f2_2(x+dx1)-f2_2(x))/dx;
    deriv_x2 = (f2_2(x+dx2)-f2_2(x))/dx;
    grad = [deriv_x1 deriv_x2];
    grad = grad/norm(grad);
    x = dtBrent(@f2_2,grad,x);
    coords = [coords; x];
    counter = counter + 1;
end

xp = -2:0.1:0.2;
yp = -2:0.1:0;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2_2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2_2(x)
num_iterations = counter


%% Quasi-Newton Method

clc

% Starting point
x0 = [0; 0];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;
grad_multiplier = 5;
dx = 0.001;
dx1 = [dx 0];
dx2 = [0 dx];

% Calculation variables
x = x0;
prev = x0 + 1;
% H = eye(2);
H = dtHess2d(@f2_2,x);
coords = x;

while norm(prev - x)>precision && counter<iteration_limit
    prev = x;
    
    h = -H*dtGrad(@f2_2,x,2);
    x = dtBrent(@f2_2,h,x);
    coords = [coords x];
    
    % BFGS
    f = dtGrad(@f2_2,x,2);
    u = x/(x'*f)-(H*f)/(f'*H*f);
    H_temp = H + (x*(x'))/(x'*f) - ((H*f)*((H*f)'))/(f'*H*f);
    H = H_temp + (f'*H*f)*(u*(u'));
    
%     H = H+((x-H*f)/(x'*H*f))*(x'*H);
    
    counter = counter + 1;
end

coords = coords';
xp = -2:0.1:0.5;
yp = -2:0.1:0.5;
[xx,yy] = meshgrid(xp,yp);
zz = 0*xx;
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        zz(i,j) = f2_2([xx(i,j) yy(i,j)]);
    end
end
figure;
contour(xx,yy,zz);
grid on;
hold on;
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

x_min = x
y_min = f2_2(x)
num_iterations = counter
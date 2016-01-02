% 2. Function Minimization in 1D

clear all
close all
clc

% To run code properly, please execute cells in order of appearence 

%% 1 - Brute Force Approach

clc

x = -10:0.0001:10; % Stoopid array of points
y = f1(x); % Stoopid evaluation of function on billion points
s = size(x);
s = s(2); % Number of points
index = 1; % Initialize minimum index
min = y(1); % Initialize minimum value

% Find smalles value or values
for i = 2:s
    if y(i)<min
        min = y(i);
        index = i;
    elseif y(i) == min
        % If more than one points are minima, return array
        index = [index i];
    end
end

% Print answers on terminal
x_minimum = x(index)
minimum = min

figure;
plot(x,y);
hold on;
grid on;
stem(x(index),y(index),'r');

%% 2 - Golden Search

% Note: a < b < c

clc

% Initial Points
a = -10+randn();
c = 10+randn();
b = -0.4;

% Constant
w = 0.38197; 

% Distances
d1 = b - a;
d2 = c - b;
d = c - b;

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;

coords = [b f1(b)];

% First, check if there is a minimum inside interval
if f1(b)<f1(a) && f1(b)<f1(c)
    
    % Apply Golden Search
        % Loop stops when precision is obtained or when iterations reach an
        % iteration number limit
    while d>precision && counter<iteration_limit
        if d1>d2
            x = b - w*(b - a);
            if f1(x)>f1(b)
                a = x;
            else
                c = b;
                b = x;
            end
        else
            x = b + w*(c - b);
            if f1(x)>f1(b)
                c = x;
            else
                a = b;
                b = x;
            end
        end
        
        coords = [coords; b f1(b)];
        
        % Update distances and counter to check if loop should continue
        d1 = b - a;
        d2 = c - b;
        d = max([d1 d2]);
        counter = counter + 1;
    end
else
    'No guaranteed minimum in defined interval'
end

xp = -0.5:0.001:1.5;
yp = f1(xp);
figure;
plot(xp,yp);
hold on;
grid on;
% stem(coords(:,1),coords(:,2),'r');
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

% Print answers and counter on terminal
x_minimum = b
minimum = f1(b)
number_of_iterations = counter

%% 3 - Brent's Method

% Note: a < c < b

clc

% Initialize values
a = -7;
b = 6.2;
c = (a+b)/2;

coords = [c f1(c)];

% Calculation parameters
precision = 0.00001;
counter = 0;
iteration_limit = 100;
d = precision + 1;

% Check if there's a minimum inside interval
if f1(c)<f1(a) && f1(c)<f1(b)
    
    % Apply Brent's Method
    while d>precision && counter<iteration_limit
        
        % Calculate new point
        num = ((b-a)^2)*(f1(b)-f1(c)) - ((b-c)^2)*(f1(b)-f1(a));
        den = 2*((b-a)*(f1(b)-f1(c)) - (b-c)*(f1(b)-f1(a)));
        x = b - num/den;

        % Update points' positions
        if x > c
            if f1(x) > f1(c)
                b = x;
            else
                a = c;
                c = x;
            end
            
        elseif x < c
            if f1(x) > f1(c)
                a = x;
            else
                b = c;
                c = x;
            end

        end
        
        coords = [coords; c f1(c)];

        % Update distance and counter
        d = b - a;
        counter = counter + 1;
    end
else
    'No guaranteed minimum in defined interval'
end

xp = -0.5:0.001:1.5;
yp = f1(xp);
figure;
plot(xp,yp);
hold on;
grid on;
% stem(coords(:,1),coords(:,2),'r');
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

% Print results and counter
x_minimum = c
minimum = f1(c)
number_of_iterations = counter


%% 4 - Gauss-Newton Method

clc

% Initial value
x = -0.3;

% Calculation parameters
dx = 0.001;
precision = 0.00001;
counter = 0;
iteration_limit = 30;
d = precision + 1;

coords = [x f1(x)];

% Apply Newton's Method
while d>precision && counter<iteration_limit
    
    % Calculate derivatives to both sides
    deriv_1 = (f1(x)-f1(x-dx))/dx;
    deriv_2 = (f1(x+dx)-f1(x))/dx;
    
    % Average derivatives to minimize numerical error
    first_deriv = (deriv_2+deriv_1)/2;
    
    % Calculate second derivative
    sec_deriv = (deriv_2-deriv_1)/dx;
    
    % Store old point for distance calculation
    x_old = x;
    
    % Calculate new point
    x = x - (first_deriv/sec_deriv);
    
    coords = [coords; x f1(x)];
    
    % Update distance and counter
    d = abs(x - x_old);
    counter = counter + 1;
end

xp = -0.5:0.001:0.8;
yp = f1(xp);
figure;
plot(xp,yp);
hold on;
grid on;
% stem(coords(:,1),coords(:,2),'r');
plot(coords(:,1),coords(:,2),'--ro','LineWidth',2);

% Print solutions and counter
x_extremum = x
extremum = f1(x)
number_of_iterations = counter


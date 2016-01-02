function H = dtHess2d(f,x)
    H = zeros(2);
    dx = 0.0001;
    u1 = [dx; 0];
    u2 = [0; dx];
    H(1,1) = (feval(f,x+u1)-2*feval(f,x)+feval(f,x-u1))/(dx^2);
    H(1,2) = (((feval(f,x+u1+u2)-feval(f,x-u1+u2))/(2*dx))-((feval(f,x+u1-u2)-feval(f,x-u1-u2))/(2*dx)))/(2*dx);
    H(2,1) = (((feval(f,x+u2+u1)-feval(f,x+u1-u2))/(2*dx))-((feval(f,x-u1+u2)-feval(f,x-u1-u2))/(2*dx)))/(2*dx);
    H(2,2) = (feval(f,x+u2)-2*feval(f,x)+feval(f,x-u2))/(dx^2);
end
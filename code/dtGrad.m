function grad = dtGrad(f,x0,dim)
    grad = zeros(dim,1);
    
    dx = 0.0001;
    for i = 1:dim
        u = zeros(dim,1);
        u(i) = dx;
        grad(i) = -feval(f,x0+2*u)+8*feval(f,x0+u)-8*feval(f,x0-u);
        grad(i) = (grad(i)+feval(f,x0-2*u))/(12*dx);
    end
end
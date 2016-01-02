function cost = dtF(parameters)
    a = parameters(1);
    b = parameters(2);
    xc = parameters(3);
    yc = parameters(4);
    phi = parameters(5);
    global data;
    n = size(data,2);
    cost = 0;
    for i = 1:n
        cost = cost + 0.5*dtError(a,b,xc,yc,phi,data(:,i))^2;
    end
end
function points = dtEllipse(parameters)
    a = parameters(1);
    b = parameters(2);
    xc = parameters(3);
    yc = parameters(4);
    phi = parameters(5);
    % Random angle
    theta = linspace(0,2*pi,100);
    % Ellipse point
    x = xc + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi);
    y = yc + a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi);
    points = [x; y];
end
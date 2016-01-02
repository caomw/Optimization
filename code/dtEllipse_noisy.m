function points = dtEllipse_noisy(a,b,xc,yc,phi,n)
    points = zeros(2,n);
    for i = 1:n
        % Random angle
        theta = 2*pi*rand();
        % Ellipse point
        x = xc + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi);
        y = yc + a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi);
        points(:,i) = [x; y];
    end
    % Noise
    d = (a+b)/100; % Controls noise variance
    noise = d*randn(2,n);
    points = points + noise;
end
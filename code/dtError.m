function dist_sq = dtError(a,b,xc,yc,phi,point)
    % Calculates the square of the distance between point and the a point
    % in the ellipse with the fiven parameters and the same angle from the
    % centre
    diff = point - [xc; yc];
    theta = atan2(diff(2),diff(1))-phi;
    x = xc + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi);
    y = yc + a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi);
    estimated = [x; y];
    dist_sq = norm(point - estimated)^2;
end
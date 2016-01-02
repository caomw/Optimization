function out = f2_2(in)
    x = in(1);
    y = in(2);
    out = (1+((x+y+1)^2)*(19-14*x+3*x*x-14*y+6*x*y+3*y*y));
    out = out*(30+((2*x-3*y)^2)*(18-32*x+12*x*x+48*y-36*x*y+27*y*y));
end
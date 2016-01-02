function [f_of_x] = f1(x)
    f_of_x = exp(-(x+1).^2)-exp(-x.^2);
end
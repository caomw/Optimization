function J = dtJacob(f,n,data)
    s = size(data,2);
    J = zeros(s,n);
    for i = 1:s
        J(i,:) = dtGrad(f,data(
    end
end
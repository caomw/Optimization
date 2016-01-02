function x_min = dtBrent(f,d,x0)
    %%%%% N-dimensional Brent's Method minimization %%%%%
    % f -> fucntion handler (@f2 for example)
    % d -> direction of minimization
    % x0 -> starting point

    a = 0;
    b = 0;
    c = 0;
    
    d = d/norm(d); % Just in case
    
    % Initialize values
    counter = 0;
    found  = 0;
    while counter<10000 && found==0
        a = -exp((counter-2000)/100)*tan(rand()*pi/2);
        b = exp((counter-2000)/100)*tan(rand()*pi/2);
        
        if feval(f,x0+c*d)<feval(f,x0+a*d) && feval(f,x0+c*d)<feval(f,x0+b*d)
            found = 1;
        end
        counter = counter+1;
    end
    
    % Calculation parameters
    precision = 0.00001;
    counter = 0;
    iteration_limit = 1000;
    dist = precision + 1;

    % Check if there's a minimum inside interval
    if feval(f,x0+c*d)<feval(f,x0+a*d) && feval(f,x0+c*d)<feval(f,x0+b*d)

        % Apply Brent's Method
        while dist>precision && counter<iteration_limit

            % Calculate new point
            num = ((b-a)^2)*(feval(f,x0+b*d)-feval(f,x0+c*d)) - ((b-c)^2)*(feval(f,x0+b*d)-feval(f,x0+a*d));
            den = 2*((b-a)*(feval(f,x0+b*d)-feval(f,x0+c*d)) - (b-c)*(feval(f,x0+b*d)-feval(f,x0+a*d)));
            x = b - num/den;

            % Update points' positions
            if x > c
                if feval(f,x0+x*d) > feval(f,x0+c*d)
                    b = x;
                else
                    a = c;
                    c = x;
                end

            elseif x < c
                if feval(f,x0+x*d) > feval(f,x0+c*d)
                    a = x;
                else
                    b = c;
                    c = x;
                end

            end

            % Update distance and counter
            dist = b - a;
            counter = counter + 1;
        end
    else
        % 'No guaranteed minimum in defined interval'
    end
    
    x_min = c*d+x0;
    
end
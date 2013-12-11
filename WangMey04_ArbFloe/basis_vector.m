function N_d = basis_vector(d, x, y)
% Initializing Vector of Basis Functions
    N_d = zeros(1,12);
    
    global a;   % Panel Size
    
    x_temp(1) = -1;
    x_temp(2) = 1;
    y_temp(1) = -1;
    y_temp(2) = 1;
    
    [ x_1 , y_1 , x_2 , y_2 , mp_x , mp_y ] = integration_points( d );

    % Points within Panel d
     x = ( x - mp_x ) / a;
     y = ( y - mp_y ) / a;
   
    count = 1;
    % Running through Panels: bottom left -> top left -> bottom right -> top right
    for i = 1 : 2
        for j = 1 : 2
            N_d(count) = 1/8 * (1 + x_temp(i)*x) * (1 + y_temp(j)*y) * ( 2 + x_temp(i)*x + y_temp(j)*y - x^2 - y^2 );
            N_d(count+1) = a/16 * (1 + x_temp(i)*x) * (y_temp(j) + y) * (y^2 - 1);
            N_d(count+2) = -a/16 * (x_temp(i) + x) * (x^2 - 1) * (1 + y_temp(j)*y);
            count = count + 3;
        end
    end
end
function g = green_matrix( d , e)
% Calculation of Green matrix [g]_de for every panel combination
    
    global a;
    u = a^2 * [1; 1; 1; 1]; % Vector of Integration Weights for 4 Integration Points
 	v = a^2 * [25/81; 40/81; 25/81; 40/81; 64/81; 40/81; 25/81; 40/81; 25/81]; % Vector of Integration Weights for 9 Integration Points
    
    % number of integration points per panel
	N = 4;     % = length(u);
	M = 9;     % = length(v);
     
	N_1 = zeros(12,N);
	N_2 = zeros(M,12);
	
    % x-y-coordinates of integration points & panel centre
	[ x_1_d , y_1_d , x_2_d , y_2_d , mp_x_d , mp_y_d ] = integration_points( d );
	[ x_1_e , y_1_e , x_2_e , y_2_e , mp_x_e , mp_y_e ] = integration_points( e );
    
     
	for i = 1 : N
        N_1(1:12,i) = u(i) * transpose(basis_vector(d, x_1_d(i), y_1_d(i)));
%         N_1(1:12,i) = u(i) * transpose(basis(d, x_1_d(i), y_1_d(i)));
    end
    
	for i = 1 : M
        N_2(i,1:12) = v(i) * basis_vector(d, x_2_d(i), y_2_d(i));
%         N_2(i,1:12) = v(i) * basis(d, x_2_d(i), y_2_d(i));
    end
     
	if d == e
        green = zeros(N,M);
        for i = 1 : N
            for j = 1 : M
                % green(i,j) = green_finite_meylan(x_1_d(i), y_1_d(i), x_2_d(j), y_2_d(j)); % d=e
                green(i,j) = green_finite(x_1_d(i), y_1_d(i), x_2_d(j), y_2_d(j)); % d=e
                % green(i,j) = green_infinite(x_1_d(i), y_1_d(i), x_2_d(j), y_2_d(j)); % d=e
            end
        end
        g = N_1*green*N_2;
    end
    
    if d ~= e
        green = zeros(N,N);
        for i = 1 : N
            for j = 1 : N
                % green(i,j) = green_finite_meylan(x_1_d(i), y_1_d(i), x_1_e(j), y_1_e(j)); % d != e
                green(i,j) = green_finite(x_1_d(i), y_1_d(i), x_1_e(j), y_1_e(j)); % d != e
                % green(i,j) = green_infinite(x_1_d(i), y_1_d(i), x_1_e(j), y_1_e(j)); % d != e
            end
        end
        g = N_1*green*transpose(N_1);
    end     
end
function g = green_finite_meylan(x, y, xi, eta)
% Calculation of Green Function for finite Water Depth (like Meylan)

% Parameters
    global H;
    global k;
    global alpha;
    
    N = 20; %100 % number of summands for infinite sum

	disp_roots = dispersion_free_surface(N+1,H); % Note: disp_roots(1) is purely imaginary
    
    g = 0; % Initializing of Green Function
    pkt_abst = sqrt((x-xi)^2 + (y-eta)^2); % Euclidean Distance between Observation Point and Source Point
    
    hlp = alpha^2 - k^2; % = w^4/g^2 - k^2
    g = -.5*1i*(hlp / ( hlp*H - alpha))*(cosh(k*H))^2 * besselh(0,1,k*pkt_abst);    
    
    C = 0;    
    for j = 1 : N
        k_j = disp_roots(j+1);
        hlp = k_j^2 + alpha^2;
        C = C + (hlp / ( hlp*H - alpha))*(cos(k_j*H))^2 * besselk(0,k_j*pkt_abst);
        % besselk(nu,Z) computes the modified Bessel function of the second kind, K?(z), 
        % for each element of the array Z. The order nu need not be an integer, but must be real. 
        % The argument Z can be complex. The result is real where Z is positive.
    end
    g = g - 1 / pi * C;
end
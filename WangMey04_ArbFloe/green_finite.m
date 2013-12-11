function g = green_finite(x, y, xi, eta)
% Calculation of Green Function for finite Water Depth

% Parameters
    global H;
    N = 10; %100 % number of summands for infinite sum

	disp_roots = dispersion_free_surface(N+1,H); % Note: disp_roots(1) is purely imaginary
    
    g = 0; % Initializing of Green Function
    pkt_abst = sqrt((x-xi)^2 + (y-eta)^2); % Euclidean Distance between Observation Point and Source Point
    
    for j = 1 : N
        C = H/2 * ( 1 + ( sin(2*disp_roots(j+1)*H) ) / ( 2*disp_roots(j+1)*H ) );        
        bessel = besselk(0,disp_roots(j+1)*pkt_abst);
        % besselk(nu,Z) computes the modified Bessel function of the second kind, K?(z), 
        % for each element of the array Z. The order nu need not be an integer, but must be real. 
        % The argument Z can be complex. The result is real where Z is positive.
        
        g = g - ((bessel/(2*pi*C))*(cos(disp_roots(j+1)*H))^2);
    end
    
    
%   global alpha;
%   global k;
% 	g = -1i/2 * ( alpha^2 - k^2 ) / ( (alpha^2 - k^2)*H - alpha^2 ) * (cosh(k*H))^2 * besselh(0, k * pkt_abst);
% 	for j = 1 : N
%         k_j = disp_roots(j+1);
%         C = 1/pi * ( (k_j)^2 + alpha^2 ) / ( ((k_j)^2 + alpha^2)*H - alpha^2 );
%         g = g - C * (cos(k_j * H))^2 * besselk(0, k_j * pkt_abst);
%     end
end
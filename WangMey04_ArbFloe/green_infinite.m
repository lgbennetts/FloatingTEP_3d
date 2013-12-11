function g = green_infinite(x, y, xi, eta)
% Calculation of Green Function for infinite Water Depth

    global alpha;   
    
    pkt_abst = sqrt((x-xi)^2 + (y-eta)^2); % Euclidean Distance between Observation Point and Source Point
    
    % Calculation of Green Function
    g = -1/(4*pi) * (2/pkt_abst - alpha*pi*(struve(0,alpha*pkt_abst) + bessely(0,alpha*pkt_abst) - 2*1i*besselj(0,alpha*pkt_abst)));

%     g = -1/(4*pi) * (2/pkt_abst - alpha*pi*(struve(0,alpha*pkt_abst) +
%     bessely(0,alpha*pkt_abst) - 2*1i*pi*besselj(0,alpha*pkt_abst))); % with  pi
end
function pot_in_dis = potential_in()
% Calculation of incident potential

    global A_in;    % Amplitude of incident waves
    global k;       % wave number
    global theta;   % angle of incident waves
    global node_no; % number of nodes
    A_in = 1.4*A_in;
    % Function Handles for Potential und x- & y-Derivative
    pot_in = @(x,y) A_in * exp(1i*k*(x*cos(theta) + y*sin(theta)));
    pot_in_x = @(x,y) A_in * 1i * k * cos(theta) * exp(1i*k*(x*cos(theta) + y*sin(theta)));
    pot_in_y = @(x,y) A_in * 1i * k * sin(theta) * exp(1i*k*(x*cos(theta) + y*sin(theta)));

    cnt = 1;
    [ coor_x , coor_y ] = coordinates(); % Coordinates of Nodes
    
    pot_in_dis = zeros(3*node_no,1);
    for i = 1 : node_no
        pot_in_dis(cnt) = pot_in(coor_x(i), coor_y(i));
        pot_in_dis(cnt+1) = pot_in_x(coor_x(i), coor_y(i));
        pot_in_dis(cnt+2) = pot_in_y(coor_x(i), coor_y(i));
        cnt = cnt + 3;
    end
end
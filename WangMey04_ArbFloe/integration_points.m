function [ x_1 , y_1 , x_2 , y_2 , mp_x , mp_y] = integration_points( d)
% Calculation of x-y-coordinates of Integration Points for Panel d

    global a;
    global panels_per_row;
    
    % x-y-coordinates of centre for every panel
    mp_x = ( -2 + a ) + ( ceil(d / panels_per_row) - 1 ) * (2 * a);
    mp_y = ( -2 + a ) +  mod(d-1,panels_per_row) * (2 * a);

    
% %   4x4 Ice Floe (by hand)
%     switch d
%         case 1
%             mp_x = -1.5;
%             mp_y = -1.5;
%         case 2
%             mp_x = -1.5;
%             mp_y = -0.5;
%     	case 3
%             mp_x = -1.5;
%             mp_y = 0.5;
%         case 4
%             mp_x = -1.5;
%             mp_y = 1.5;
%         case 5
%             mp_x = -0.5;
%             mp_y = -1.5;
%     	case 6
%             mp_x = -0.5;
%             mp_y = -0.5;
%         case 7
%             mp_x = -0.5;
%             mp_y = 0.5;
%         case 8
%             mp_x = -0.5;
%             mp_y = 1.5;
%     	case 9
%             mp_x = 0.5;
%             mp_y = -1.5;
%         case 10
%             mp_x = 0.5;
%             mp_y = -0.5;
%         case 11
%             mp_x = 0.5;
%             mp_y = 0.5;
%         case 12
%             mp_x = 0.5;
%             mp_y = 1.5;
%         case 13
%             mp_x = 1.5;
%             mp_y = -1.5;
%         case 14
%             mp_x = 1.5;
%             mp_y = -0.5;
%         case 15
%             mp_x = 1.5;
%             mp_y = 0.5;
%         case 16
%             mp_x = 1.5;
%             mp_y = 1.5;
%     end    


    % Integration Weights from "Prof. Suvranu De - Numerical Integration in
    % 2D"
    % Weights are 1!

    x_1(1) = mp_x - a*(1/sqrt(3));
	x_1(2) = mp_x - a*(1/sqrt(3));
	x_1(3) = mp_x + a*(1/sqrt(3));
	x_1(4) = mp_x + a*(1/sqrt(3));
        
	y_1(1) = mp_y - a*(1/sqrt(3));
	y_1(2) = mp_y + a*(1/sqrt(3));
	y_1(3) = mp_y - a*(1/sqrt(3));
	y_1(4) = mp_y + a*(1/sqrt(3));
    
    
    % numbered from bottom left -> top left -- from left -> right
    x_2(1) = mp_x - a*sqrt(3/5);
	x_2(2) = mp_x - a*sqrt(3/5);
	x_2(3) = mp_x - a*sqrt(3/5);
    x_2(4) = mp_x;
	x_2(5) = mp_x;
	x_2(6) = mp_x;
    x_2(7) = mp_x + a*sqrt(3/5);
	x_2(8) = mp_x + a*sqrt(3/5);
	x_2(9) = mp_x + a*sqrt(3/5);
    
    y_2(1) = mp_y - a*sqrt(3/5);
    y_2(2) = mp_y;
    y_2(3) = mp_y + a*sqrt(3/5);
    y_2(4) = mp_y - a*sqrt(3/5);
    y_2(5) = mp_y;
    y_2(6) = mp_y + a*sqrt(3/5);
    y_2(7) = mp_y - a*sqrt(3/5);
    y_2(8) = mp_y;
    y_2(9) = mp_y + a*sqrt(3/5);
end
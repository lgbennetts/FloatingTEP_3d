function [ coor_x , coor_y ] = coordinates()
% Coordinates for every Node -> Distance 2 for Panel Size 2a with a=1
    
    global a;
    global panels_per_row;
    
    cnt = 1;
    v_help = [-2 : 2*a : 2]';
    for i = -2 : 2*a : 2
        coor_x(cnt : cnt+panels_per_row,1) = i*ones(panels_per_row+1,1);
        coor_y(cnt : cnt+panels_per_row,1) = v_help;
        cnt = cnt + panels_per_row+1;
    end   


% 10x10 ice floe (by hand)
%     coor_x = [ -1 ; -1 ; -1 ; -1 ; -1 ; -0.5 ; -0.5 ; -0.5 ; -0.5 ; -0.5 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0.5 ; 0.5 ; 0.5 ; 0.5 ; 0.5 ; 1 ; 1 ; 1 ; 1 ; 1 ];
%     coor_y = [ -1 ; -0.5 ; 0 ; 0.5 ; 1 ; -1 ; -0.5 ; 0 ; 0.5 ; 1 ; -1 ; -0.5 ; 0 ; 0.5 ; 1 ; -1 ; -0.5 ; 0 ; 0.5 ; 1 ; -1 ; -0.5 ; 0 ; 0.5 ; 1 ];
end
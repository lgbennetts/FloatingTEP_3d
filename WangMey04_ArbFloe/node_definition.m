function Q = node_definition()
% Specification of ice floe geometry
    global no_panels;
	global panels_per_row;
    q_1 = zeros(1,no_panels);
    q_2 = zeros(1,no_panels);
    q_3 = zeros(1,no_panels);
    q_4 = zeros(1,no_panels);
    
    
    for i = 0 : no_panels-1
        i_h = i+1;
        q_1(i_h) = i_h + floor(i/panels_per_row);
        q_2(i_h) = q_1(i_h) + 1;
        q_3(i_h) = q_1(i_h) + ( panels_per_row + 1 );
        q_4(i_h) = q_3(i_h) + 1;
    end

%   10x10 ice floe (by hand)
%     q_1(1) = 1;
%     q_2(1) = q_1(1) + 1;
%     q_3(1) = 6;
%     q_4(1) = q_3(1) + 1;
% 
%     q_1(2) = 2;
%     q_2(2) = q_1(2) + 1;
%     q_3(2) = 7;
%     q_4(2) = q_3(2) + 1;
% 
%     q_1(3) = 3;
%     q_2(3) = q_1(3) + 1;
%     q_3(3) = 8;
%     q_4(3) = q_3(3) + 1;
% 
%     q_1(4) = 4;
%     q_2(4) = q_1(4) + 1;
%     q_3(4) = 9;
%     q_4(4) = q_3(4) + 1;
% 
%     q_1(5) = 6;
%     q_2(5) = q_1(5) + 1;
%     q_3(5) = 11;
%     q_4(5) = q_3(5) + 1;
% 
%     q_1(6) = 7;
%     q_2(6) = q_1(6) + 1;
%     q_3(6) = 12;
%     q_4(6) = q_3(6) + 1;
% 
%     q_1(7) = 8;
%     q_2(7) = q_1(7) + 1;
%     q_3(7) = 13;
%     q_4(7) = q_3(7) + 1;
% 
%     q_1(8) = 9;
%     q_2(8) = q_1(8) + 1;
%     q_3(8) = 14;
%     q_4(8) = q_3(8) + 1;
% 
%     q_1(9) = 11;
%     q_2(9) = q_1(9) + 1;
%     q_3(9) = 16;
%     q_4(9) = q_3(9) + 1;
%     
%     q_1(10) = 12;
%     q_2(10) = q_1(10) + 1;
%     q_3(10) = 17;
%     q_4(10) = q_3(10) + 1;
%     
%     q_1(11) = 13;
%     q_2(11) = q_1(11) + 1;
%     q_3(11) = 18;
%     q_4(11) = q_3(11) + 1;
%     
%     q_1(12) = 14;
%     q_2(12) = q_1(12) + 1;
%     q_3(12) = 19;
%     q_4(12) = q_3(12) + 1;
%     
%     q_1(13) = 16;
%     q_2(13) = q_1(13) + 1;
%     q_3(13) = 21;
%     q_4(13) = q_3(13) + 1;
%     
%     q_1(14) = 17;
%     q_2(14) = q_1(14) + 1;
%     q_3(14) = 22;
%     q_4(14) = q_3(14) + 1;
%     
%     q_1(15) = 18;
%     q_2(15) = q_1(15) + 1;
%     q_3(15) = 23;
%     q_4(15) = q_3(15) + 1;
%     
%     q_1(16) = 19;
%     q_2(16) = q_1(16) + 1;
%     q_3(16) = 24;
%     q_4(16) = q_3(16) + 1;
    
    Q = [q_1',q_2',q_3',q_4'];
end
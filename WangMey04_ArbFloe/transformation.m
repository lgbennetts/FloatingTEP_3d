function O = transformation(d,Q)
% Calculation of transformation matrices [o]_d for panel d
    global node_no;
    
    q_1 = Q(:,1);
    q_2 = Q(:,2);
    q_3 = Q(:,3);
    q_4 = Q(:,4);

    O = zeros(12,3*node_no);
	O(1,3*q_1(d)-2) = 1;
	O(2,3*q_1(d)-1) = 1;
	O(3,3*q_1(d)) = 1;
	O(4,3*q_2(d)-2) = 1;
	O(5,3*q_2(d)-1) = 1;
	O(6,3*q_2(d)) = 1;
	O(7,3*q_3(d)-2) = 1;
	O(8,3*q_3(d)-1) = 1;
	O(9,3*q_3(d)) = 1;
	O(10,3*q_4(d)-2) = 1;
	O(11,3*q_4(d)-1) = 1;
	O(12,3*q_4(d)) = 1;
end
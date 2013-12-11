function [mass_mtrx , steif_mtrx] = FEM_matrices()

    global a;   % Edge Length of Panels
	global nu;  % Poisson number
    b = a;

    % Initializing of Mass Matrix [m]_d

    % improved compared to Wang
    mass_mtrx = zeros(12,12);

    mass_mtrx_11 = [3454,       922/2*b,    -922/2*a;
                    922/2*b,    320/4*b^2,  -252/4*a*b;
                    -922/2*a,   -252/4*a*b, 320/4*a^2];

    mass_mtrx_22 = [3454,       -922/2*b,   -922/2*a;
                    -922/2*b,   320/4*b^2,  252/4*a*b;
                    -922/2*a,   252/4*a*b,  320/4*a^2];

    mass_mtrx_33 = [3454,       922/2*b,    922/2*a;
                    922/2*b,    320/4*b^2,  252/4*a*b;
                    922/2*a,    252/4*a*b,  320/4*a^2];

    mass_mtrx_44 = [3454,       -922/2*b,   922/2*a;
                    -922/2*b,   320/4*b^2,  -252/4*a*b;
                    922/2*a,    -252/4*a*b, 320/4*a^2];

    mass_mtrx_21 = [1226,       548/2*b,    -398/2*a;
                    -548/2*b,   -240/4*b^2, 168/4*a*b;
                    -398/2*a,   -168/4*a*b, 160/4*a^2];

    mass_mtrx_31 = [1226,       398/2*b,    -548/2*a;
                    398/2*b,    160/4*b^2,  -168/4*a*b;
                    548/2*a,    168/4*a*b,  -240/4*a^2];

    mass_mtrx_41 = [394,        232/2*b,    -232/2*a;
                    -232/2*b,   -120/4*b^2, 112/4*a*b;
                    232/2*a,    112/4*a*b,  -120/4*a^2];

    mass_mtrx_32 = [394,        -232/2*b,   -232/2*a;
                    232/2*b,    -120/4*b^2, -112/4*a*b;
                    232/2*a,    -112/4*a*b, -120/4*a^2];

    mass_mtrx_42 = [1226,       -398/2*b,   -548/2*a;
                    -398/2*b,   160/4*b^2,  168/4*a*b;
                    548/2*a,    -168/4*a*b, -240/4*a^2];

    mass_mtrx_43 = [1226,       548/2*b,    398/2*a;
                    -548/2*b,   -240/4*b^2, -168/4*a*b;
                    398/2*a,    168/4*a*b,  160/4*a^2];

    mass_mtrx_12 = transpose(mass_mtrx_21);
    mass_mtrx_13 = transpose(mass_mtrx_31);
    mass_mtrx_14 = transpose(mass_mtrx_41);
    mass_mtrx_23 = transpose(mass_mtrx_32);
    mass_mtrx_24 = transpose(mass_mtrx_42);
    mass_mtrx_34 = transpose(mass_mtrx_43);

    mass_mtrx(1:3,1:3) = mass_mtrx_11;
    mass_mtrx(1:3,4:6) = mass_mtrx_12;
    mass_mtrx(1:3,7:9) = mass_mtrx_13;
    mass_mtrx(1:3,10:12) = mass_mtrx_14;

    mass_mtrx(4:6,1:3) = mass_mtrx_21;
    mass_mtrx(4:6,4:6) = mass_mtrx_22;
    mass_mtrx(4:6,7:9) = mass_mtrx_23;
    mass_mtrx(4:6,10:12) = mass_mtrx_24;

    mass_mtrx(7:9,1:3) = mass_mtrx_31;
    mass_mtrx(7:9,4:6) = mass_mtrx_32;
    mass_mtrx(7:9,7:9) = mass_mtrx_33;
    mass_mtrx(7:9,10:12) = mass_mtrx_34;

    mass_mtrx(10:12,1:3) = mass_mtrx_41;
    mass_mtrx(10:12,4:6) = mass_mtrx_42;
    mass_mtrx(10:12,7:9) = mass_mtrx_43;
    mass_mtrx(10:12,10:12) = mass_mtrx_44;

    mass_mtrx = a*b/6300 * mass_mtrx;


%     mass_mtrx_11 = [3454    , 922*b     , -922*a    , 1226      , -548*b    , -398*a    ;
%                     922*b   , 320*b^2   , -252*a*b  , 548*b     , -240*b^2  , 168*a*b   ;
%                     -922*a  , -252*a*b  , 320*a^2   , -398*a    , 168*a*b   , 160*a^2   ;
%                     1226	, 548*b     , -398*a    , 3454      , -922*b    , -922*a    ;
%                     -548*b  , -240*b^2  , 168*a*b   , -922*b    , 320*b^2   , 252*a*b   ;
%                     -398*a  , 168*a*b   , 160*a^2   , -922*a    , 252*a*b   , 320*a^2  ];
%                 
%     
%     mass_mtrx_21 = [1226    , 398*b     , -548*a    , 394       , -232*b    , -232*a    ;
%                     398*b   , 160*b^2   , -168*a*b  , 232*b     , -120*b^2  , -112*a*b  ;
%                     548*a   , 168*a*b   , -240*a^2  , 232*a     , -112*a*b  , -120*a^2  ;
%                     394     , 232*b     , -232*a    , 1226      , -398*b    , -548*a    ;
%                     -232*b  , -120*b^2  , 112*a*b   , -398*b    , 160*b^2   , 168*a*b   ;
%                     232*a   , 112*a*b   , -120*a^2  , 548*a     , -168*a*b  , -240*a^2 ];
%                 
%                 
% 	mass_mtrx_22 = [3454    , 922*b     , 922*a     , 1226      , -548*b    , 398*a     ;
%                     922*b   , 320*b^2   , 252*a*b   , 548*b     , -240*b^2  , 168*a*b   ;
%                     922*a   , 252*a*b   , 320*a^2   , 398*a     , -168*a*b  , 160*a^2   ;
%                     1226    , 548*b     , 398*a     , 3453      , -922*b    , 922*a     ;
%                     -548*b  , -240*b^2  , -168*a*b  , -922*b    , 320*b^2   , -252*a*b  ;
%                     398*a   , 168*a*b   , 160*a^2   , 922*a     , -252*a*b  , 320*a^2  ];
%                 
%     mass_mtrx(1:6,1:6) = mass_mtrx_11;
%     mass_mtrx(7:12,1:6) = mass_mtrx_21;
%     mass_mtrx(1:6,7:12) = transpose(mass_mtrx_21);
%     mass_mtrx(7:12,7:12) = mass_mtrx_22;            
%                 
%     mass_mtrx = a*b/6300 * mass_mtrx;


%     % Initialisierung der Steifigkeitsmatrix [k]_d
% 
%     steif_mtrx = zeros(12,12);
% 
%     steif_mtrx_11 = [4*((b^2)/(a^2) + (a^2)/(b^2)) + 0.4*(7-2*nu),  (2*((a^2)/(b^2)) + 0.2*(1+4*nu))*b,         -(2*((b^2)/(a^2)) + 0.2*(1+4*nu))*a;
%                     (2*((a^2)/(b^2)) + 0.2*(1+4*nu))*b,             (4/3*((a^2)/(b^2)) + 4/15*(1-nu))*(b^2),    -nu*a*b;
%                     -(2*((b^2)/(a^2)) + 0.2*(1+4*nu))*a,            -nu*a*b,                                    (4/3*((b^2)/(a^2)) + 4/15*(1-nu))*(a^2)];
% 
%     eye_2 = eye(3);
%     eye_2(3,3) = -1;
% 
%     eye_3 = eye(3);
%     eye_3(2,2) = -1;
% 
%     eye_4 = eye(3);
%     eye_4(1,1) = -1;
% 
%     steif_mtrx_22 = eye_2*steif_mtrx_11*eye_2;
% 
%     steif_mtrx_33 = eye_3*steif_mtrx_11*eye_3;
% 
%     steif_mtrx_44 = eye_4*steif_mtrx_11*eye_4;
% 
%     steif_mtrx_21 = [-(2*(2*(b^2/a^2)-(a^2/b^2)) + 0.4*(7-2*nu)),   -(2*((a^2)/(b^2)) + 0.2*(1-nu)) * b,        -(((b^2)/(a^2)) - 0.2*(1+4*nu))*a;
%                     (((a^2)/(b^2)) - 0.2*(1+4*nu))*b,               (2/3*((a^2)/(b^2)) - 4/15*(1-nu))*(b^2),    0;
%                     -(2*((b^2)/(a^2)) + 0.2*(1-nu))*a,              0,                                          (2/3*((b^2)/(a^2)) - 1/15*(1-nu))*(a^2)];
% 
%     steif_mtrx_31 = [(2*(((b^2)/(a^2)) - 2*((a^2)/(b^2))) - 0.4*(7-2*nu)),  -(2*((a^2)/(b^2)) + 0.2*(1-nu))*b,          -(((b^2)/(a^2)) - 0.2*(1+4*nu))*a;
%                     (((a^2)/(b^2))-0.2*(1+4*nu))*b,                         (2/3*((a^2)/(b^2)) - 4/15*(1-nu))*(b^2),    0;
%                     -(2*((b^2)/(a^2)) + 0.2*(1-nu))*a,                      0,                                          (2/3*((b^2)/(a^2)) - 1/15*(1-nu))*(a^2)];
% 
%     steif_mtrx_41 = [-(2*(((b^2)/(a^2)) + ((a^2)/(b^2))) - 0.4*(7-2*nu)),   (((a^2)/(b^2)) - 0.2*(1-nu))*b,             -(((b^2)/(a^2)) - 0.2*(1-nu))*a;
%                     (-((a^2)/(b^2)) + 0.2*(1-nu))*b,                        (1/3*((a^2)/(b^2)) + 1/15*(1-nu))*(b^2),    0;
%                     (((b^2)/(a^2)) - 0.2*(1-nu))*a,                         0,                                          (1/3*((b^2)/(a^2)) + 1/15*(1-nu))*(a^2)];
% 
%     steif_mtrx_32 = eye_2*steif_mtrx_41*eye_2;
% 
%     steif_mtrx_42 = eye_3*steif_mtrx_31*eye_3;
% 
%     steif_mtrx_43 = eye_2*steif_mtrx_21*eye_2;
% 
%     steif_mtrx(1:3,1:3) = steif_mtrx_11;
%     steif_mtrx(1:3,4:6) = steif_mtrx_21;
%     steif_mtrx(1:3,7:9) = steif_mtrx_31;
%     steif_mtrx(1:3,10:12) = steif_mtrx_41;
% 
%     steif_mtrx(4:6,1:3) = steif_mtrx_21;
%     steif_mtrx(4:6,4:6) = steif_mtrx_22;
%     steif_mtrx(4:6,7:9) = steif_mtrx_32; % Dreher in Vorlage!!
%     steif_mtrx(4:6,10:12) = steif_mtrx_42; % Dreher in Vorlage!!
% 
%     steif_mtrx(7:9,1:3) = steif_mtrx_31;
%     steif_mtrx(7:9,4:6) = steif_mtrx_32;
%     steif_mtrx(7:9,7:9) = steif_mtrx_33;
%     steif_mtrx(7:9,10:12) = steif_mtrx_43;  % Dreher in Vorlage!!
% 
%     steif_mtrx(10:12,1:3) = steif_mtrx_41;
%     steif_mtrx(10:12,4:6) = steif_mtrx_42;
%     steif_mtrx(10:12,7:9) = steif_mtrx_43;
%     steif_mtrx(10:12,10:12) = steif_mtrx_44;
% 
%     steif_mtrx = steif_mtrx/(4*(1-(nu^2))*a*b);

% Initialisierung der Steifigkeitsmatrix [k]_d

    steif_mtrx = zeros(12,12);

    a_h = a / b;
    b_h = b / a;
    
    steif_mtrx_11 = [4*(b_h^2 + a_h^2) + 0.4*(7-2*nu)   , 2*(2*a_h^2 + 0.2*(1+4*nu))*b      , 2*(-2*b_h^2 - 0.2*(1-4*nu))*a ; 
                     2*(2*a_h^2 + 0.2*(1+4*nu))*b       , 4*(4/3*a_h^2 + 4/15*(1-nu))*b^2   , -4*nu*a*b                     ;
                     2*(-2*b_h^2 - 0.2*(1-4*nu))*a      , -4*nu*a*b                         , 4*(4/3*b_h^2 + 1/15*(1-nu))*a^2 ];
        
    steif_mtrx_31 = [-(2*(2*b_h^2 - a_h^2) + 0.4*(7-2*nu))  , 2*(a_h^2 - 0.2*(1+4*nu))*b        , 2*(2*b_h^2 + 0.2*(1-nu))*a ;
                     2*(a_h^2 - 0.2*(1+4*nu))*b             , 4*(2/3*a_h^2 - 4/15*(1-nu))*b^2   , 0                          ;
                     -2*(2*b_h^2 + 0.2*(1-nu))*a            , 0                                 , 4*(2/3*b_h^2 - 1/15*(1-nu))*a^2];
                 
    steif_mtrx_41 = [-(2*(b_h^2 + a_h^2) - 0.4*(7-2*nu))    , 2*(-a_h^2 + 0.2*(1-nu))*b         , 2*(b_h^2 - 0.2*(1-nu))*a ;
                     2*(a_h^2 - 0.2*(1-nu))*b               , 4*(1/3*a_h^2 + 1/15*(1-nu))*b^2   , 0                          ;
                     2*(-b_h^2 + 0.2*(1-nu))*a              , 0                                 , 4*(1/3*b_h^2 + 1/15*(1-nu))*a^2];
        
    steif_mtrx_21 = [2*(b_h^2 - 2*a_h^2) - 0.4*(7-2*nu)     , 2*(-2*a_h^2 - 0.2*(1-nu))*b       , 2*(-b_h^2 + 0.2*(1+4*nu))*a ;
                     2*(2*a_h^2 + 0.2*(1-nu))*b             , 4*(2/3*a_h^2 - 1/15*(1-nu))*b^2   , 0                          ;
                     2*(-b_h^2 + 0.2*(1+4*nu))*a            , 0                                 , 4*(2/3*b_h^2 - 4/15*(1-nu))*a^2];

                 
    eye_3 = eye(3);
    eye_3(3,3) = -1;

    eye_2 = eye(3);
    eye_2(2,2) = -1;

    eye_1 = eye(3);
    eye_1(1,1) = -1;

    steif_mtrx_33 = eye_3*steif_mtrx_11*eye_3;
    
    steif_mtrx_43 = eye_3*steif_mtrx_21*eye_3;
    steif_mtrx_44 = eye_1*steif_mtrx_11*eye_1;
    
    steif_mtrx_32 = eye_3*steif_mtrx_41*eye_3;
    steif_mtrx_42 = eye_1*steif_mtrx_31*eye_1;
    steif_mtrx_22 = eye_2*steif_mtrx_11*eye_2;
    
%     steif_mtrx_12 = steif_mtrx_21;
%     steif_mtrx_13 = steif_mtrx_31;
%     steif_mtrx_14 = steif_mtrx_41;
%     
%     steif_mtrx_23 = steif_mtrx_32;
%     steif_mtrx_24 = steif_mtrx_42;
%     
%     steif_mtrx_34 = steif_mtrx_43; 

    steif_mtrx(1:3,1:3) = steif_mtrx_11;
    steif_mtrx(1:3,4:6) = steif_mtrx_21;
    steif_mtrx(1:3,7:9) = steif_mtrx_31;
    steif_mtrx(1:3,10:12) = steif_mtrx_41;

    steif_mtrx(4:6,1:3) = steif_mtrx_21;
    steif_mtrx(4:6,4:6) = steif_mtrx_22;
    steif_mtrx(4:6,7:9) = steif_mtrx_32;
    steif_mtrx(4:6,10:12) = steif_mtrx_42;

    steif_mtrx(7:9,1:3) = steif_mtrx_31;
    steif_mtrx(7:9,4:6) = steif_mtrx_32;
    steif_mtrx(7:9,7:9) = steif_mtrx_33;
    steif_mtrx(7:9,10:12) = steif_mtrx_43;

    steif_mtrx(10:12,1:3) = steif_mtrx_41;
    steif_mtrx(10:12,4:6) = steif_mtrx_42;
    steif_mtrx(10:12,7:9) = steif_mtrx_43;
    steif_mtrx(10:12,10:12) = steif_mtrx_44;

    steif_mtrx = steif_mtrx/(4*(1-(nu^2))*a*b);
end
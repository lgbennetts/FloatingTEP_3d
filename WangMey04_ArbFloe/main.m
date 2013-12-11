% MAIN Function for solution of ice floe problem with arbitrary geometry
clear all;
% Parameters
beta = 0.0001; % 0.0001; % 0.01; % 0.02; % 0.02; % 0.01; % stiffness of ice floe
gamma = 0;% 0; % mass of ice floe
lambda = 4;    % wave length

global H; H = 2; % 100;                % non-dimensional water depth
% -> for infinite water depth: adapt green_matrix()
global alpha; alpha = pi/4; % alpha = sqrt(w)/g -> value out of 2004 Wang&Meylan

global A_in; A_in = 1;          % Amplitude of incident waves  -->  look at potential_in
global k; k = 2*pi/lambda;      % wave number
global theta; theta = pi/6; % 0; % pi/2; % pi/4;    % angle of incident waves
global nu; nu = 0.31;    % 1/3;                     % Poisson number (between 0 and 0.5)
% in 2009 Peter & Meylan p. 1583 nu is 1/3

%% Geometry of Ice Floes
global panels_per_row; panels_per_row = 10;
global no_panels; no_panels = panels_per_row^2;	% number of panels
global node_no; node_no = (panels_per_row+1)^2;	% number of nodes
global a; a = 2/panels_per_row;                 % edge length of panels

% str = sprintf('Panels = %d , H = %.1f , Amplitude = %.2f , lambda = %d , alpha = %.1f , beta = %.3f , gamma = %.1f , theta = %.1f' , panels_anz, H, A_in, lambda, alpha, beta, gamma, theta);
str = sprintf('Panels = %d , H = %.1f , Amplitude = %.1f , lambda = %d , beta = %.4f , gamma = %.3f , angle = %.2f' , panels_per_row, H, A_in, lambda, beta, gamma, theta );

% -> coordinates() % coordinates of nodes
% -> integration_points( d) % calculation of x-y-coordinates of integration points for panel d
% -> node_definition() % geometry of ice floe

% Calculation of incident potential
pot_in = potential_in();

% Calculation of Green matrix, Mass matrix and stiffness matrix
[G,M,K] = matrices();

% Solution of LS for displacement
mtr_help = inv(M - alpha*G); % Hilfsmatrix
LHS = beta * K + (1 - alpha*gamma) * M + alpha * M * mtr_help * G;
RHS = 1i * sqrt(alpha) * M * mtr_help * M * pot_in;
w = LHS\RHS;
w_plot = w(1:3:end);  % read every third element of solution vector (displacement)


figure;

% SET bdry = 1 to cut the boundary panels
% otherwise plot with boundary
bdry = 1;
if ( bdry == 1)
    w_plot(node_no-panels_per_row : node_no)=[];
    w_plot(1 : panels_per_row+1)=[];
    w_plot(1 : panels_per_row+1 : end)=[];
    w_plot(panels_per_row : panels_per_row : end)=[];
    a_h = 2 / (panels_per_row - 2);
    node_no_h = (panels_per_row-1)^2;
else
	a_h = 2 / panels_per_row;
    node_no_h = (panels_per_row+1)^2;
end


[X,Y] = meshgrid(-2:2*a_h:2);

Z = reshape(w_plot,sqrt(node_no_h),sqrt(node_no_h));
mesh(X,Y,real(Z));
xlabel('x','FontSize',24);
ylabel('y','FontSize',24);
zlabel('Re(w)','FontSize',24);

set(gca,'XTick',-2 : 1 : 2 , 'FontSize',16);
set(gca,'YTick',-2 : 1 : 2 , 'FontSize',16);
% set(gca,'ZTick',-2 : 1 : 2 , 'FontSize',16);

% set(gca,'XTick',-1.3 : 0.5 : 1.3 , 'FontSize',16)
% set(gca,'XTick',-1 : 0.5 : 1 , 'FontSize',16)
% set(gca,'YTick',-1 : 0.5 : 1 , 'FontSize',16)
% % zlim([-A_in A_in])
% set(gca,'ZTick',-A_in : A_in/2 : A_in,'FontSize',16)
% set(gca,'FontSize',20)
% str = sprintf('Panels = %d , H = %.1f , Amplitude = %.1f , lambda = %d , alpha = %.1f , beta = %.3f , gamma = %.1f , theta = %.1f' , no_panels, H, A_in, lambda, alpha, beta, gamma, theta);
title(str);
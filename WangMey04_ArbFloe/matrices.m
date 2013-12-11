function [ G , M , K ] = matrices()
% Calculation of Matrices (Green Matrix G, Mass Matrix M, Stiffness Matrix
% K)

% Reference: Petyt - Introduction to FE Vibration Analysis (S.196)

    global no_panels;
    global node_no;
    
    O_d = zeros(12,3*node_no,no_panels); % Definition of Transformation Matrix
    Q = node_definition(); % Geometry of Ice Floe

    for d = 1 : no_panels
        O_d(:,:,d) = transformation(d,Q);   % Initializing of Transformation Matrix for Panel d
    end

    % Calculation of Green Matrix G
    G = zeros(3*node_no,3*node_no);
    for d = 1 : no_panels
        A = zeros(12,3*node_no);
        for e = 1 : no_panels
            % [m,n] = size(green_matrix(d,e)*O_d(:,:,e))
            A = A + green_matrix(d,e)*O_d(:,:,e);
        end
        G = G + transpose(O_d(:,:,d))*A;
    end
    
    [mass_mtrx , steif_mtrx] = FEM_matrices(); % -> Matrices [m]_d & [k]_d

    M = zeros(3*node_no,3*node_no);
    K = zeros(3*node_no,3*node_no);
    % Calculation of M and K
    for d = 1 : no_panels
        M = M + transpose(O_d(:,:,d)) * mass_mtrx * O_d(:,:,d);
        K = K + transpose(O_d(:,:,d)) * steif_mtrx * O_d(:,:,d);
    end    
end
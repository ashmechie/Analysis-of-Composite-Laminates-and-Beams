function [A, B, D, ABD] = ABDmatrix(Qbar, t)
% ===============================================================
% Function to compute A, B, D, and ABD matrices for a laminate
% ===============================================================
%
% Inputs:
%   Qbar_matrices - 3x3xN array of Q̅ matrices for each ply
%   t              - thickness of each ply (same for all)
%
% Outputs:
%   A, B, D        - 3x3 laminate stiffness matrices
%   ABD            - 6x6 assembled stiffness matrix

n = size(Qbar, 3); % Number of plies

% z-coordinates (mid-plane at 0)
z = zeros(n+1, 1);
total_thickness = n * t;
z(1) = -total_thickness / 2;
for k = 2:n+1
    z(k) = z(k-1) + t;
end

% Initialize
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

% Compute matrices
for k = 1:n
    Qb = Qbar(:,:,k);
    A = A + Qb * (z(k+1) - z(k));
    B = B + 0.5 * Qb * (z(k+1)^2 - z(k)^2);
    D = D + (1/3) * Qb * (z(k+1)^3 - z(k)^3);
end

% Assemble ABD matrix
ABD = [A B; B D];

end

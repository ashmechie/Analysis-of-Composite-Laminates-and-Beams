function loads = inputLoads()
% ===============================================================
% Function to input in-plane and bending loads at the laminate level
% ===============================================================
%
% Outputs:
%   loads - (6x1) column matrix containing:
%           [Nxx; Nyy; Nxy; Mxx; Myy; Mxy]

disp('Enter the laminate loads (lb/in and lb-in/in):');

Nxx = input('Enter Nxx (Normal force per unit width in x-direction): ');
Nyy = input('Enter Nyy (Normal force per unit width in y-direction): ');
Nxy = input('Enter Nxy (Shear force per unit width): ');
Mxx = input('Enter Mxx (Bending moment per unit width about x-axis): ');
Myy = input('Enter Myy (Bending moment per unit width about y-axis): ');
Mxy = input('Enter Mxy (Twisting moment per unit width): ');

% Combine into a 6x1 column vector
loads = [Nxx; Nyy; Nxy; Mxx; Myy; Mxy];

end

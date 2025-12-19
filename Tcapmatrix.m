function Tcap = Tcapmatrix(theta_deg)
% TCAPMATRIX  Returns the strain transformation matrix [T̄]
% Input:  theta_deg - ply orientation angle in degrees
% Output: Tcap - 3x3 transformation matrix for strain

theta = deg2rad(theta_deg); % Convert to radians
c = cos(theta);
s = sin(theta);

Tcap = [ c^2      s^2      s*c;
         s^2      c^2     -s*c;
        -2*s*c   2*s*c   (c^2 - s^2) ];
end

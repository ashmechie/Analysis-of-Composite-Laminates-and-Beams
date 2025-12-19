function [R_initial, R_fail_at_ply, failure_order] = computeInitialAndProgressiveFailure( ...
            n, theta_deg, z_mid, z_interfaces, Q, Qbar_all, ...
            F1, F2, F11, F22, F66, F12, loads, eps0, kappa)
% Computes initial Tsai-Wu R-values and performs progressive failure (zero-stiffness)
% Inputs:
%   n, theta_deg, z_mid, z_interfaces, Q, Qbar_all,
%   F1,F2,F11,F22,F66,F12, loads, eps0, kappa
% Outputs:
%   R_initial (n x 1)
%   R_fail_at_ply (n x 1) - stored R when ply removed
%   failure_order - m x 3 [iteration, ply_number, R_fail]

tol = 1e-9;

%% ---------- Initial R-values ----------
R_initial = NaN(n,1);

for i = 1:n
    strain_global = eps0 + z_mid(i) * kappa;     % 3x1
    Tcap = Tcapmatrix(theta_deg(i));             % 3x3
    strain_local = Tcap * strain_global;         % 3x1
    stress_local = Q * strain_local;             % 3x1 (psi)

    s11 = stress_local(1)/1000;  % convert to ksi for Tsai-Wu inputs if strengths in ksi
    s22 = stress_local(2)/1000;
    t12 = stress_local(3)/1000;

    Acoef = F11*s11^2 + 2*F12*s11*s22 + F22*s22^2 + F66*t12^2;
    Bcoef = F1*s11 + F2*s22;
    Ccoef = -1;
    disc = Bcoef^2 - 4*Acoef*Ccoef;
    if disc >= 0 && abs(Acoef) > eps
        R1 = (-Bcoef + sqrt(disc)) / (2*Acoef);
        R2 = (-Bcoef - sqrt(disc)) / (2*Acoef);
        R_initial(i) = max([R1,R2,0]);
    else
        R_initial(i) = NaN;
    end
end

%% ---------- Progressive Failure Loop ----------
active = true(n,1);
R_fail_at_ply = NaN(n,1);
failure_order = [];
iteration = 1;

fprintf('\nStarting progressive failure loop (retain positions, zero stiffness for failed plies)...\n');

while true
    % Build modified Qbar stack (zero for failed plies)
    Qbar_mod = zeros(3,3,n);
    for i = 1:n
        if active(i)
            Qbar_mod(:,:,i) = Qbar_all(:,:,i);
        else
            Qbar_mod(:,:,i) = zeros(3);
        end
    end

    % Recompute ABD from Qbar_mod
    A = zeros(3); B = zeros(3); D = zeros(3);
    for i = 1:n
        z1 = z_interfaces(i);
        z2 = z_interfaces(i+1);
        Qb = Qbar_mod(:,:,i);
        A = A + Qb * (z2 - z1);
        B = B + 0.5 * Qb * (z2^2 - z1^2);
        D = D + (1/3) * Qb * (z2^3 - z1^3);
    end
    ABD = [A B; B D];

    if rcond(ABD) < 1e-12
        warning('ABD singular or nearly singular. Stopping progressive failure loop.');
        break;
    end

    % Recompute global strains under same loads
    strain_curvature = ABD \ loads;
    eps0_iter = strain_curvature(1:3);
    kappa_iter = strain_curvature(4:6);

    % Compute R per ply for this iteration
    R_iter = NaN(n,1);
    for i = 1:n
        strain_global = eps0_iter + z_mid(i) * kappa_iter;
        Tcap = Tcapmatrix(theta_deg(i));
        strain_local = Tcap * strain_global;
        stress_local = Q * strain_local;

        s11 = stress_local(1)/1000;
        s22 = stress_local(2)/1000;
        t12 = stress_local(3)/1000;

        Acoef = F11*s11^2 + 2*F12*s11*s22 + F22*s22^2 + F66*t12^2;
        Bcoef = F1*s11 + F2*s22;
        Ccoef = -1;
        disc = Bcoef^2 - 4*Acoef*Ccoef;
        if disc >= 0 && abs(Acoef) > eps
            R1 = (-Bcoef + sqrt(disc)) / (2*Acoef);
            R2 = (-Bcoef - sqrt(disc)) / (2*Acoef);
            R_iter(i) = max([R1,R2,0]);
        else
            R_iter(i) = NaN;
        end
    end

    active_indices = find(active);
    active_Rs = R_iter(active_indices);
    valid_mask = ~isnan(active_Rs);

    if ~any(valid_mask)
        fprintf('No valid R among active plies. Stopping.\n');
        break;
    end

    minR = min(active_Rs(valid_mask));

    to_fail_mask = (R_iter <= minR + tol) & active;
    to_fail_indices = find(to_fail_mask);

    if isempty(to_fail_indices)
        fprintf('No plies identified to fail this iteration (numerical). Stopping.\n');
        break;
    end

    for idx = to_fail_indices'
        R_fail_at_ply(idx) = minR;
        failure_order = [failure_order; iteration, idx, minR];
    end

    IterTable = table((1:n)', theta_deg', z_mid, R_iter, active, ...
        'VariableNames', {'Ply_Number','Angle_deg','z_mid_in','R_value','Active'});
    fprintf('\nIteration %d results:\n', iteration);
    disp(IterTable);
    fprintf('Plies failing this iteration: %s  (stored R_fail = %g)\n', mat2str(to_fail_indices'), minR);

    active(to_fail_indices) = false;

    if ~any(active)
        fprintf('All plies have failed. End of progressive failure.\n');
        break;
    end

    iteration = iteration + 1;
end

fprintf('\n=============================================\n');
if ~isempty(failure_order)
    fprintf('Ply Failure Order Summary:\n');
    FailureSummary = array2table(failure_order, ...
        'VariableNames', {'Iteration','Ply_Number','R_fail'});
    disp(FailureSummary);
else
    fprintf('No ply failures recorded.\n');
end
fprintf('=============================================\n');

end

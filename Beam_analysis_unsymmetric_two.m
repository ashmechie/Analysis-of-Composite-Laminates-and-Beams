
clear; clc; close all;

fprintf('\n=== COMPOSITE I-BEAM SECTION: Modulus-weighted properties (UNSYMMETRIC-ONLY) ===\n\n');

%% -------------------------
% 1) Geometry inputs (I-section only)
% -------------------------
b_top  = input('Enter TOP flange width b_top (in) (left edge at Z=0): ');
b_bot  = input('Enter BOTTOM flange width b_bot (in) (left edge at Z=0): ');
h_enter = input('Enter TOTAL web vertical extent h_enter (in). NOTE: include half top & bottom flange thicknesses as you normally do: ');
z_web_center = input('Enter horizontal location of web centerline from leftmost point of flange (Z, in) [user system, + right]: ');

%% -------------------------
% 2) Material selection (user functions must exist)
% -------------------------
fprintf('\n--- Material selection ---\n');
[E11, E22, G12, v12, v21] = getMaterialProperties();
[vf, SLt, SLc, STt, STc, SLTs] = getMaterialStrengths();
fprintf('Material properties loaded.\n\n');

%% -------------------------
% 3) Layups for top, web, bottom
% -------------------------
sameLayup = upper(input('Use SAME layup for top flange, web, bottom flange? (Y/N): ','s'));
if isempty(sameLayup), sameLayup = 'Y'; end

if strcmp(sameLayup,'Y')
    fprintf('\nEnter layup for ALL regions:\n');
    L_all = getLayupFromUser();
    L_top = L_all; L_web = L_all; L_bot = L_all;
else
    fprintf('\nEnter layup for TOP flange:\n');   L_top = getLayupFromUser();
    fprintf('\nEnter layup for WEB:\n');          L_web = getLayupFromUser();
    fprintf('\nEnter layup for BOTTOM flange:\n'); L_bot = getLayupFromUser();
end

% Extract laminate thicknesses (out-of-plane thickness)
t_top = L_top.totalThickness;
t_web = L_web.totalThickness; % this is the web thickness in Z-direction (out-of-plane)
t_bot = L_bot.totalThickness;

fprintf('\nLaminate thicknesses (in): top=%.4f  web(thick)=%.4f  bottom=%.4f\n', t_top, t_web, t_bot);

%% -------------------------
% 4) Compute effective web vertical height (h_web_eff)
% -------------------------
h_web_eff = h_enter - (t_top/2) - (t_bot/2);
if h_web_eff <= 0
    error('Computed effective web height <= 0. Check inputs and layup thicknesses.');
end

%% -------------------------
% 5) Place region geometry in (Y,Z) using the USER coordinate system:
%    (Origin bottom-left, Z positive RIGHT)
% -------------------------
% Bottom flange: 
z_bot_min = 0; z_bot_max = b_bot;
y_bot_min = 0; y_bot_max = t_bot;
A_bot = b_bot * t_bot;
z_bot_cent = (z_bot_min + z_bot_max)/2;   % right-positive coordinate
y_bot_cent = (y_bot_min + y_bot_max)/2;

% Web: thickness in Z direction 
z_web_min = z_web_center - t_web/2;
z_web_max = z_web_center + t_web/2;
y_web_min = t_bot;
y_web_max = t_bot + h_web_eff;
A_web = t_web * h_web_eff;
z_web_cent = (z_web_min + z_web_max)/2;
y_web_cent = (y_web_min + y_web_max)/2;

% Top flange: 
z_top_min = 0; z_top_max = b_top;
y_top_min = t_bot + h_web_eff;
y_top_max = y_top_min + t_top;
A_top = b_top * t_top;
z_top_cent = (z_top_min + z_top_max)/2;
y_top_cent = (y_top_min + y_top_max)/2;

%% -------------------------
% 6) Compute region Exx and Gxy from layups (use calcABDForLayup)
% -------------------------
calc_Exx_fromA = @(A_mat,h) (A_mat(1,1)/h) * (1 - (A_mat(1,2)^2)/(A_mat(1,1)*A_mat(2,2)));
calc_Gxy_fromA = @(A_mat,h) A_mat(3,3)/h;

% Top flange
[A_top_mat, ~, ~, ~] = calcABDForLayup(L_top, E11, E22, G12, v12, v21);
if any(isnan(A_top_mat),'all')
    Exx_top = NaN; Gxy_top = NaN;
else
    Exx_top = calc_Exx_fromA(A_top_mat, t_top);
    Gxy_top = calc_Gxy_fromA(A_top_mat, t_top);
end

% Web
[A_web_mat, ~, ~, ~] = calcABDForLayup(L_web, E11, E22, G12, v12, v21);
if any(isnan(A_web_mat),'all')
    Exx_web = NaN; Gxy_web = NaN;
else
    Exx_web = calc_Exx_fromA(A_web_mat, t_web);
    Gxy_web = calc_Gxy_fromA(A_web_mat, t_web);
end

% Bottom flange
[A_bot_mat, ~, ~, ~] = calcABDForLayup(L_bot, E11, E22, G12, v12, v21);
if any(isnan(A_bot_mat),'all')
    Exx_bot = NaN; Gxy_bot = NaN;
else
    Exx_bot = calc_Exx_fromA(A_bot_mat, t_bot);
    Gxy_bot = calc_Gxy_fromA(A_bot_mat, t_bot);
end

%% -------------------------
% 7) Build Region list (geometry stored in USER coordinate system: right-positive)
% -------------------------
RegionList = {};

RegionList{end+1} = struct('name','Bottom flange','A',A_bot,'Exx',Exx_bot,'Gxy',Gxy_bot, ...
    'z_min',z_bot_min,'z_max',z_bot_max,'y_min',y_bot_min,'y_max',y_bot_max, ...
    'z_cent',z_bot_cent,'y_cent',y_bot_cent,'width',b_bot,'height',t_bot,'thickness',t_bot,'Lay',L_bot);
RegionList{end+1} = struct('name','Web','A',A_web,'Exx',Exx_web,'Gxy',Gxy_web, ...
    'z_min',z_web_min,'z_max',z_web_max,'y_min',y_web_min,'y_max',y_web_max, ...
    'z_cent',z_web_cent,'y_cent',y_web_cent,'width',t_web,'height',h_web_eff,'thickness',t_web,'Lay',L_web);
RegionList{end+1} = struct('name','Top flange','A',A_top,'Exx',Exx_top,'Gxy',Gxy_top, ...
    'z_min',z_top_min,'z_max',z_top_max,'y_min',y_top_min,'y_max',y_top_max, ...
    'z_cent',z_top_cent,'y_cent',y_top_cent,'width',b_top,'height',t_top,'thickness',t_top,'Lay',L_top);

%% -------------------------
% 8) Modulus-weighted centroid (Y_bar, Z_bar) using Exx*Area weighting
% -------------------------
numY = 0; numZ = 0; den = 0;
for k = 1:length(RegionList)
    R = RegionList{k};
    if isnan(R.Exx)
        warning('Region %s has NaN Exx; it will be excluded from weighted centroid sums.', R.name);
        continue;
    end
    numY = numY + R.Exx * R.A * R.y_cent;
    numZ = numZ + R.Exx * R.A * R.z_cent;
    den = den + R.Exx * R.A;
end
if den == 0
    error('Denominator for modulus-weighted centroid is zero. Check Exx values.');
end
Y_bar = numY / den;    % Y_bar in user system (right-positive Y is same direction as internal Y)
Z_bar = numZ / den;    % Z_bar in user system (right-positive)

fprintf('\nModulus-weighted centroid (expressed in USER coords): Y_bar = %.6f in, Z_bar = %.6f in (user origin bottom-left, Z + right)\n', Y_bar, Z_bar);

%% -------------------------
% 9) Compute modulus-weighted A* and inertia tensor components (Iyy, Izz, Iyz)
% -------------------------
Er = 10e6; % reference modulus (psi)

A_star = 0; Iyy_star = 0; Izz_star = 0; Iyz_star = 0;
for k = 1:length(RegionList)
    R = RegionList{k};
    if isnan(R.Exx), continue; end
    Ei = R.Exx; Ai = R.A;
    dy = R.y_cent - Y_bar;                
    dz_internal = (Z_bar - R.z_cent);     

    b_region = R.width;   
    h_region = R.height;  

    Iyy_cent = h_region * (b_region^3) / 12;  
    Izz_cent = b_region * (h_region^3) / 12;  

    Iyy_shift = Iyy_cent + Ai * (dz_internal^2);
    Izz_shift = Izz_cent + Ai * (dy^2);
    Iyz_shift = Ai * (dy * dz_internal);  

    A_star = A_star + (Ei/Er) * Ai;
    Iyy_star = Iyy_star + (Ei/Er) * Iyy_shift;
    Izz_star = Izz_star + (Ei/Er) * Izz_shift;
    Iyz_star = Iyz_star + (Ei/Er) * Iyz_shift;
end

den_inertia = Iyy_star * Izz_star - Iyz_star^2;
if abs(den_inertia) < 1e-15
    warning('Denominator for unsymmetric bending/shear coupling is very small (near singular). Results may be ill-conditioned.');
end

fprintf('Modulus-weighted A* = %.6e in^2 (Er = %.3e psi)\n', A_star, Er);
fprintf('Modulus-weighted Iyy = %.6e in^4\n', Iyy_star);
fprintf('Modulus-weighted Izz = %.6e in^4\n', Izz_star);
fprintf('Modulus-weighted Iyz = %.6e in^4\n', Iyz_star);

%% -------------------------
% NEW: Print region Exx and Gxy and total geometric area (user-helpful)
% -------------------------
fprintf('\nRegion elastic parameters (from layups):\n');
fprintf('  %-15s  %-12s  %-12s  %-12s\n','Region','Area (in^2)','Exx (psi)','Gxy (psi)');
for k = 1:length(RegionList)
    R = RegionList{k};
    fprintf('  %-15s  %10.6f   %10.4e   %10.4e\n', R.name, R.A, R.Exx, R.Gxy);
end
A_total = A_top + A_web + A_bot;
fprintf('\nGeometric total area (A_total) = %.6e in^2 (previous P/A calculation used this) \n', A_total);
fprintf('Modulus-weighted area (A_star) = %.6e in^2 (used for new axial distribution)\n', A_star);

%% -------------------------
% 10) Plot cross-section and centroid (plot in USER coordinates)
% -------------------------
figure; hold on;
title('I-section with modulus-weighted centroid (USER coords: origin bottom-left, Z + right)');
xlabel('Z (in) - horizontal (positive to RIGHT for plot)'); ylabel('Y (in) - vertical (positive UP)');
set(gca,'YDir','normal'); axis equal; grid on; box on;
% draw rectangles (Position: [x y w h] where x=Zmin (USER coords), y=Ymin)
rectangle('Position',[z_bot_min, y_bot_min, b_bot, t_bot],'FaceColor',[0.8 0.85 1]);
rectangle('Position',[z_web_min, y_web_min, t_web, h_web_eff],'FaceColor',[0.7 0.9 0.7]);
rectangle('Position',[z_top_min, y_top_min, b_top, t_top],'FaceColor',[0.8 0.85 1]);
% plot centroid in USER coords:
plot(Z_bar, Y_bar, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
text(Z_bar, Y_bar, sprintf('  (Y_{bar}=%.4f, Z_{bar}=%.4f) (USER coords)', Y_bar, Z_bar), 'Color','r');
xlim([min([z_bot_min,z_web_min,z_top_min]) - 0.1, max([z_bot_max,z_web_max,z_top_max]) + 0.1]);
ylim([-0.1, y_top_max + 0.1]);
hold off;

%% -------------------------
% UNSYMMETRIC PATH (always used)
% -------------------------
fprintf('\nUsing UNSYMMETRIC analysis path (symmetric branch removed).\n');

Myy = input('Enter bending moment Myy (in-lb) about Y axis (vertical): ');
Mzz = input('Enter bending moment Mzz (in-lb) about Z axis (horizontal): ');
Vy  = input('Enter shear Vy (lb) along Y direction (vertical shear component): ');
Vz  = input('Enter shear Vz (lb) along Z direction (horizontal shear component): ');
P_axial = input('Enter axial load P (lb). Positive = tension. Enter 0 if none: ');
if isempty(P_axial), P_axial = 0; end

% New axial handling:
if A_star == 0
    error('A_star is zero: cannot distribute axial load with modulus-weighted area.');
end
% global axial scalar (to be multiplied by E_region/Er for region stress)
sigma_axial_global = P_axial / A_star;    % units: lb / in^2 (psi)

Den = Iyy_star * Izz_star - Iyz_star^2;
if abs(Den) < 1e-15
    warning('Denominator for unsymmetric shear coupling is very small (near singular). Results may be ill-conditioned.');
end

% Prepare storage for region loads
computedRegionLoads = cell(length(RegionList),1); % store {coord_list, Nxy_local_array, Qstructs}

% Loop over regions and request ONLY Z (for flanges) or ONLY Y (for web).
for r = 1:length(RegionList)
    R = RegionList{r};
    if isnan(R.Exx)
        computedRegionLoads{r} = {[], NaN, struct('Qy',[], 'Qz',[], 'cutPoints',[])};
        fprintf('%s: Exx NaN -> skipped\n', R.name);
        continue;
    end

    wi = R.Exx / Er;    % modulus weight for area contributions
    ti = R.thickness;   % thickness used to convert tau->Nxy

    % Prepare storage for multiple cut evaluations for this region
    coord_list = [];    % will store rows [cutCoordZ, cutCoordY, Qy, Qz, Nxy]
    Qy_arr = []; Qz_arr = []; Nxy_arr = [];
    cutPts = [];

    if strcmpi(R.name,'Top flange') || strcmpi(R.name,'Bottom flange')
        % FLANGE: ask only for Zcut. flange is split LEFT/RIGHT at Zcut.
        fprintf('\nRegion "%s" (FLANGE). Enter Zcut values (one per line). Enter ''000'' to finish.\n', R.name);
        while true
            s = input('  Enter Zcut (or 000 to finish) >> ','s'); s = strtrim(s);
            if isempty(s)
                Zcut = R.z_cent;
                fprintf('   empty => using centroid Z=%.4f\n', Zcut);
            elseif strcmp(s,'000')
                break;
            else
                Zcut = str2double(s);
                if isnan(Zcut)
                    fprintf('  invalid numeric; try again.\n'); continue;
                end
            end

            % clamp cut to region extents
            Zcut = min(max(Zcut, R.z_min), R.z_max);

            % LEFT piece: [z_min -> Zcut]
            width_left = max(0, Zcut - R.z_min);
            if width_left > 0
                A_left = R.thickness * width_left;
                zcent_left = R.z_min + width_left/2;
                ycent_left = R.y_cent;
                Ai_star = wi * A_left;
                z_int = Z_bar - zcent_left;   % internal left-positive sign convention
                y_int = ycent_left - Y_bar;   % up-positive
                Qy_star = Ai_star * z_int;
                Qz_star = Ai_star * y_int;
                Nxy_val = ( ((-Vy * Iyy_star +Vz * Iyz_star) / Den) * Qz_star - ((Vz*Izz_star - Vy * Iyz_star) / Den) * Qy_star );
                coord_list = [coord_list; Zcut, NaN, Qy_star, Qz_star, Nxy_val]; %#ok<AGROW>
                Qy_arr = [Qy_arr; Qy_star]; Qz_arr = [Qz_arr; Qz_star]; Nxy_arr = [Nxy_arr; Nxy_val];
                cutPts = [cutPts; Zcut, NaN];
                fprintf(' %s LEFT @ Z=%.4f: Qy*=%.6e Qz*=%.6e Nxy=%.6e\n', R.name, Zcut, Qy_star, Qz_star, Nxy_val);
            end

            % RIGHT piece: [Zcut -> z_max]
            width_right = max(0, R.z_max - Zcut);
            if width_right > 0
                A_right = R.thickness * width_right;
                zcent_right = Zcut + width_right/2;
                ycent_right = R.y_cent;
                Ai_star = wi * A_right;
                z_int = Z_bar - zcent_right;
                y_int = ycent_right - Y_bar;
                Qy_star = Ai_star * z_int;
                Qz_star = Ai_star * y_int;
                Nxy_val = ( ((-Vy * Iyy_star +Vz * Iyz_star) / Den) * Qz_star - ((Vz*Izz_star - Vy * Iyz_star) / Den) * Qy_star );
                coord_list = [coord_list; Zcut, NaN, Qy_star, Qz_star, Nxy_val]; %#ok<AGROW>
                Qy_arr = [Qy_arr; Qy_star]; Qz_arr = [Qz_arr; Qz_star]; Nxy_arr = [Nxy_arr; Nxy_val];
                cutPts = [cutPts; Zcut, NaN];
                fprintf(' %s RIGHT @ Z=%.4f: Qy*=%.6e Qz*=%.6e Nxy=%.6e\n', R.name, Zcut, Qy_star, Qz_star, Nxy_val);
            end

            if strcmp(s,'000'), break; end
        end

    else
        % WEB: ask only for Ycut. compute ABOVE (web above + top flange) and BELOW (web below + bottom flange)
        fprintf('\nRegion "%s" (WEB). Enter Ycut values (one per line). Enter ''000'' to finish.\n', R.name);
        while true
            s = input('  Enter Ycut (or 000 to finish) >> ','s'); s = strtrim(s);
            if isempty(s)
                Ycut = R.y_cent;
                fprintf('   empty => using centroid Y=%.4f\n', Ycut);
            elseif strcmp(s,'000')
                break;
            else
                Ycut = str2double(s);
                if isnan(Ycut)
                    fprintf('  invalid numeric; try again.\n'); continue;
                end
            end

            % clamp Ycut to web extents
            Ycut = min(max(Ycut, R.y_min), R.y_max);

            % ----- ABOVE: include TOP flange entire if present + web slice above Ycut -----
            Qy_above = 0; Qz_above = 0; Nxy_above = 0;
            idxTop = find(cellfun(@(r) strcmp(r.name,'Top flange'), RegionList),1);
            if ~isempty(idxTop)
                Rtop = RegionList{idxTop};
                if ~isnan(Rtop.Exx)
                    A_top = Rtop.width * Rtop.thickness;
                    Ai_star = (Rtop.Exx/Er) * A_top;
                    z_int = Z_bar - Rtop.z_cent;
                    y_int = Rtop.y_cent - Y_bar;
                    Qy_above = Qy_above + Ai_star * z_int;
                    Qz_above = Qz_above + Ai_star * y_int;
                end
            end
            % web area above cut
            h_above = max(0, R.y_max - Ycut);
            if h_above > 0
                A_web_above = R.width * h_above;
                ycent_web_above = (R.y_max + Ycut)/2;
                zcent_web_above = R.z_cent;
                Ai_star = wi * A_web_above;
                z_int = Z_bar - zcent_web_above;
                y_int = ycent_web_above - Y_bar;
                Qy_above = Qy_above + Ai_star * z_int;
                Qz_above = Qz_above + Ai_star * y_int;
            end
            if (Qy_above~=0) || (Qz_above~=0)
                Nxy_above = ( ((-Vy * Iyy_star +Vz * Iyz_star) / Den) * Qz_above - ((Vz*Izz_star - Vy * Iyz_star) / Den) * Qy_above );
                coord_list = [coord_list; NaN, Ycut, Qy_above, Qz_above, Nxy_above]; %#ok<AGROW>
                Qy_arr = [Qy_arr; Qy_above]; Qz_arr = [Qz_arr; Qz_above]; Nxy_arr = [Nxy_arr; Nxy_above];
                cutPts = [cutPts; NaN, Ycut];
                fprintf(' %s ABOVE @ Y=%.4f: Qy*=%.6e Qz*=%.6e Nxy=%.6e\n', R.name, Ycut, Qy_above, Qz_above, Nxy_above);
            end

            % ----- BELOW: include BOTTOM flange entire if present + web slice below Ycut -----
            Qy_below = 0; Qz_below = 0; Nxy_below = 0;
            idxBot = find(cellfun(@(r) strcmp(r.name,'Bottom flange'), RegionList),1);
            if ~isempty(idxBot)
                Rbot = RegionList{idxBot};
                if ~isnan(Rbot.Exx)
                    A_bot = Rbot.width * Rbot.thickness;
                    Ai_star = (Rbot.Exx/Er) * A_bot;
                    z_int = Z_bar - Rbot.z_cent;
                    y_int = Rbot.y_cent - Y_bar;
                    Qy_below = Qy_below + Ai_star * z_int;
                    Qz_below = Qz_below + Ai_star * y_int;
                end
            end
            % web area below cut
            h_below = max(0, Ycut - R.y_min);
            if h_below > 0
                A_web_below = R.width * h_below;
                ycent_web_below = (R.y_min + Ycut)/2;
                zcent_web_below = R.z_cent;
                Ai_star = wi * A_web_below;
                z_int = Z_bar - zcent_web_below;
                y_int = ycent_web_below - Y_bar;
                Qy_below = Qy_below + Ai_star * z_int;
                Qz_below = Qz_below + Ai_star * y_int;
            end
            if (Qy_below~=0) || (Qz_below~=0)
                Nxy_below = (((-Vy * Iyy_star +Vz * Iyz_star) / Den) * Qz_below - ((Vz*Izz_star - Vy * Iyz_star)/ Den) * Qy_below );
                coord_list = [coord_list; NaN, Ycut, Qy_below, Qz_below, Nxy_below]; %#ok<AGROW>
                Qy_arr = [Qy_arr; Qy_below]; Qz_arr = [Qz_arr; Qz_below]; Nxy_arr = [Nxy_arr; Nxy_below];
                cutPts = [cutPts; NaN, Ycut];
                fprintf(' %s BELOW @ Y=%.4f: Qy*=%.6e Qz*=%.6e Nxy=%.6e\n', R.name, Ycut, Qy_below, Qz_below, Nxy_below);
            end

            if strcmp(s,'000'), break; end
        end
    end

    % store computed RegionLoads: coord_list Nx5 [Zcut, Ycut, Qy, Qz, Nxy]
    % IMPORTANT: store region-specific axial stress (sigma_axial_reg) in qstruct so progressive failure uses it.
    sigma_axial_reg = (R.Exx / Er) * sigma_axial_global;  % region-specific axial stress (psi)
    computedRegionLoads{r} = {coord_list, Nxy_arr, struct('Qy',Qy_arr, 'Qz',Qz_arr, 'cutPoints',cutPts, 'sigma_axial', sigma_axial_reg, 'A_total', A_total)};

    if isempty(coord_list)
        fprintf('%s: no Q/Nxy points entered for this region.\n', R.name);
    else
        fprintf('%s: stored %d Q/Nxy evaluation(s).\n', R.name, size(coord_list,1));
    end
end % region loop (unsymmetric Q/Nxy)

%% -------------------------
% Bending sigma_xx evaluation (user loop) - unsymmetric
% -------------------------
term1 = (Myy * Izz_star + Mzz * Iyz_star) / Den;
term2 = (Mzz * Iyy_star + Myy * Iyz_star) / Den;

for r = 1:length(RegionList)
    R = RegionList{r};
    if isnan(R.Exx)
        continue;
    end

    fprintf('\nRegion "%s": enter coordinates where sigma_xx should be evaluated. Use format "Y Z" (in inches, USER system: Z + right). Enter ''000'' to finish for this region.\n', R.name);
    coord_list_bend = [];
    while true
        s = input(sprintf(' Enter Y Z (or 000 to stop) >> '),'s');
        s = strtrim(s);
        if isempty(s)
            Ys = R.y_cent;
            Zs = R.z_cent;
            fprintf(' Using region centroid (Y=%.4f, Z=%.4f)\n', Ys, Zs);
        elseif strcmp(s,'000')
            break;
        else
            toks = strsplit(s);
            if length(toks) < 2
                fprintf(' Invalid entry; please enter two numbers "Y Z" or 000.\n');
                continue;
            end
            Ys = str2double(toks{1});
            Zs = str2double(toks{2});
            if isnan(Ys) || isnan(Zs)
                fprintf('Invalid numeric values. Try again.\n');
                continue;
            end
        end

        % Convert Zs (USER right-positive) to INTERNAL left-positive relative coordinate for bending:
        Z_local = Z_bar - Zs;
        Y_local = Ys - Y_bar;  % Y direction same sign in both systems

        sigma_xx_bend = (R.Exx/Er) * ( term1 * Z_local - term2 * Y_local );

        % region-specific axial stress
        sigma_axial_region = (R.Exx / Er) * sigma_axial_global;

        % add axial stress
        sigma_xx_total = sigma_xx_bend + sigma_axial_region;
        Nxx_reg_point = sigma_xx_total * R.thickness;
        fprintf('  At (USER Y=%.4f, Z=%.4f): sigma_xx (bending) = %.4e psi ; sigma_axial_region = %.4e psi ; sigma_xx (total) = %.4e psi ; Nxx = %.4e lb/in\n', ...
            Ys, Zs, sigma_xx_bend, sigma_axial_region, sigma_xx_total, Nxx_reg_point);

        coord_list_bend = [coord_list_bend; Ys, Zs, sigma_xx_bend, sigma_axial_region, sigma_xx_total, Nxx_reg_point]; %#ok<AGROW>
    end

    % Merge or store bend coords (we'll keep separate from shear coords)
    existing = computedRegionLoads{r};
    if isempty(existing)
        % create a placeholder with bend list in first cell
        computedRegionLoads{r} = {coord_list_bend, NaN, struct('Qy',[], 'Qz',[], 'sigma_axial', (R.Exx/Er)*sigma_axial_global)};
    else
        shearCoords = existing{1};
        Nxyvals = existing{2};
        qstruct = existing{3};
        % store bend coords as additional field in qstruct
        qstruct.bendCoords = coord_list_bend;
        % ensure stored sigma_axial is region-specific (already stored previously)
        qstruct.sigma_axial = (R.Exx/Er)*sigma_axial_global;
        computedRegionLoads{r} = {shearCoords, Nxyvals, qstruct};
    end

end % bending user input loop

%% ------------------------------------------------------------------
% 12) Run progressive failure for each region using USER-SELECTED coordinates
% ------------------------------------------------------------------
fprintf('\n\n=== Running Point-Specific Progressive Failure ===\n');
[F1,F2,F11,F22,F66,F12] = getTsaiWuCoefficients(SLt, SLc, STt, STc, SLTs);
Q11 = E11 / (1 - v12 * v21);
Q12 = v12 * E22 / (1 - v12 * v21);
Q22 = E22 / (1 - v12 * v21);
Q66 = G12;
Q_mat = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66 ];

for r = 1:length(RegionList)
    reg = RegionList{r};
    fprintf('\n--- Region: %s ---\n', reg.name);
    
    % Request specific coordinate for analysis
    fprintf('Enter evaluation point (USER system) for this region.\n');
    Y_pt = input(sprintf('  Enter Y coord (range %.3f to %.3f) >> ', reg.y_min, reg.y_max));
    Z_pt = input(sprintf('  Enter Z coord (range %.3f to %.3f) >> ', reg.z_min, reg.z_max));
    
    % Calculate exact Nxx at this point
    Z_local = Z_bar - Z_pt;
    Y_local = Y_pt - Y_bar;
    sigma_xx_bend = (reg.Exx/Er) * ( ((Myy*Izz_star + Mzz*Iyz_star)/Den)*Z_local - ((Mzz*Iyy_star + Myy*Iyz_star)/Den)*Y_local );
    sigma_xx_total = sigma_xx_bend + (reg.Exx/Er)*sigma_axial_global;
    Nxx_final = sigma_xx_total * reg.thickness;
    
    % Calculate exact Nxy at this point (Shear Flow)
    % Using the same logic as the "cut" section but automated for the point
    wi = reg.Exx / Er;
    if strcmpi(reg.name,'Top flange') || strcmpi(reg.name,'Bottom flange')
        % Integration from nearest free edge
        dist_from_left = Z_pt - reg.z_min;
        dist_from_right = reg.z_max - Z_pt;
        if dist_from_left < dist_from_right
            A_eff = reg.thickness * dist_from_left;
            zcent_eff = reg.z_min + dist_from_left/2;
        else
            A_eff = reg.thickness * dist_from_right;
            zcent_eff = Z_pt + dist_from_right/2;
        end
        Qy_p = wi * A_eff * (Z_bar - zcent_eff);
        Qz_p = wi * A_eff * (reg.y_cent - Y_bar);
    else
        % Web logic: must include flange contributions
        h_above = reg.y_max - Y_pt;
        A_w_above = reg.width * h_above;
        ycent_w = (reg.y_max + Y_pt)/2;
        % Add top flange Q*
        idxT = find(cellfun(@(x) strcmp(x.name,'Top flange'), RegionList),1);
        Qt_y = (RegionList{idxT}.Exx/Er)*RegionList{idxT}.A*(Z_bar - RegionList{idxT}.z_cent);
        Qt_z = (RegionList{idxT}.Exx/Er)*RegionList{idxT}.A*(RegionList{idxT}.y_cent - Y_bar);
        Qy_p = Qt_y + wi*A_w_above*(Z_bar - reg.z_cent);
        Qz_p = Qt_z + wi*A_w_above*(ycent_w - Y_bar);
    end
    Nxy_final = (((-Vy*Iyy_star + Vz*Iyz_star)/Den)*Qz_p - ((Vz*Izz_star - Vy*Iyz_star)/Den)*Qy_p);

    % Execute Progressive Failure
    loads_p = [Nxx_final; 0; Nxy_final; 0; 0; 0];
    Lay = reg.Lay;
    n_p = Lay.totalPlies;
    z_int_p = linspace(-Lay.totalThickness/2, Lay.totalThickness/2, n_p+1);
    z_mid_p = 0.5*(z_int_p(1:end-1)+z_int_p(2:end));
    Qbar_p = zeros(3,3,n_p);
    for ii=1:n_p, Qbar_p(:,:,ii) = calcQbar(Q_mat, Lay.angles(ii)); end
    [~,~,~,ABD_p] = calcABDForLayup(Lay, E11, E22, G12, v12, v21);
    sc_p = ABD_p \ loads_p;
    
    [R_init, R_fail, ~] = computeInitialAndProgressiveFailure(n_p, Lay.angles, z_mid_p, z_int_p, Q_mat, Qbar_p, F1, F2, F11, F22, F66, F12, loads_p, sc_p(1:3), sc_p(4:6), false);

    % --- Plotting Ply Number vs R and Stresses ---
    figure('Name', ['Failure Analysis: ' reg.name]);
% Removed subplot(2,1,1) to let this plot fill the figure
hold on; grid on;

stem(1:n_p, R_init, 'b', 'LineWidth', 1.5, 'DisplayName', 'Initial R');
stem(1:n_p, R_fail, 'r--', 'DisplayName', 'Final R (Prog)');

% Add a horizontal line at R=1 for failure reference
yline(1, 'k-', 'Failure Threshold', 'LineWidth', 2); 

xlabel('Ply Number'); 
ylabel('Strength Ratio R'); 
title(['Tsai-Wu R-Value per Ply: ' reg.name]);
legend('Location', 'northeast'); % Moves legend to a clear spot
end

%% ------------------------------------------------------------------
% 14) Shear Center Calculation
% ------------------------------------------------------------------
fprintf('\n--- Calculating Shear Center (Unsymmetric) --\n');
% Moment about web centerline (z_web_center) due to unit Vy=1
% Force in Top Flange = Integral of shear flow
Qz_top_total = (Exx_top/Er)*(b_top*t_top)*(y_top_cent - Y_bar);
F_top_flange = (1/Izz_star) * (Exx_top/Er) * (t_top * b_top^3 / 12) * Vy; % simplified component
% e = (Moment of shear flows) / V
e_eccentricity = (Qz_top_total * h_web_eff * b_top) / (2 * Izz_star); 
Z_sc = z_web_center - e_eccentricity;
fprintf('Shear Center Z-coordinate: %.4f in (USER system)\n', Z_sc);

%% ------------------------------------------------------------------
% 15) Plot 3D Heatmaps and Specific Point Markers
% ------------------------------------------------------------------
z_vec = linspace(0, max([b_top, b_bot]), 120);
y_vec = linspace(0, y_top_max, 120);
[Z_grid, Y_grid] = meshgrid(z_vec, y_vec);
Nxx_grid = NaN(size(Z_grid)); Nxy_grid = NaN(size(Z_grid));

% Core Calculation Loop
for i = 1:numel(Z_grid)
    zy = [Z_grid(i), Y_grid(i)];
    for k = 1:3
        R = RegionList{k};
        if zy(2)>=R.y_min && zy(2)<=R.y_max && zy(1)>=R.z_min && zy(1)<=R.z_max
            zl = Z_bar - zy(1); yl = zy(2) - Y_bar;
            Nxx_grid(i) = ((R.Exx/Er)*(((Myy*Izz_star+Mzz*Iyz_star)/Den)*zl - ((Mzz*Iyy_star+Myy*Iyz_star)/Den)*yl) + (R.Exx/Er)*sigma_axial_global)*R.thickness;
            wi = R.Exx/Er;
            if strcmpi(R.name,'Web')
                h_ab = R.y_max - zy(2);
                idxT = find(cellfun(@(x) strcmp(x.name,'Top flange'), RegionList),1); RT = RegionList{idxT};
                Qz_val = (RT.Exx/Er)*RT.A*(RT.y_cent - Y_bar) + wi*(R.width*h_ab)*((R.y_max+zy(2))/2 - Y_bar);
                Qy_val = (RT.Exx/Er)*RT.A*(Z_bar - RT.z_cent) + wi*(R.width*h_ab)*(Z_bar - R.z_cent);
            else
                dist = min(zy(1)-R.z_min, R.z_max-zy(1));
                zc = (zy(1)<z_web_center)*(R.z_min+dist/2) + (zy(1)>=z_web_center)*(R.z_max-dist/2);
                Qz_val = wi*(R.thickness*dist)*(R.y_cent-Y_bar); Qy_val = wi*(R.thickness*dist)*(Z_bar-zc);
            end
            Nxy_grid(i) = (((-Vy*Iyy_star + Vz*Iyz_star)/Den)*Qz_val - ((Vz*Izz_star - Vy*Iyz_star)/Den)*Qy_val);
        end
    end
end

% Helper to get values at specific points
getVal = @(z,y) [interp2(Z_grid, Y_grid, Nxx_grid, z, y), interp2(Z_grid, Y_grid, Nxy_grid, z, y)];

% Define Interest Points [Z, Y, Label]
pts = [0, y_top_cent, "TF Left"; b_top, y_top_cent, "TF Right"; ...
       0, y_bot_cent, "BF Left"; b_bot, y_bot_cent, "BF Right"; ...
       z_web_center, y_bot_max, "Web-Bot"; z_web_center, y_top_min, "Web-Top"; ...
       z_web_center, y_web_cent, "Web-Mid"];

figure('Color','w','Position',[100 100 1200 500]);
titles = ["Nxx Distribution", "Nxy Distribution"]; grids = {Nxx_grid, Nxy_grid};
for m = 1:2
    subplot(1,2,m); p = pcolor(Z_grid, Y_grid, grids{m}); shading interp; set(p,'EdgeColor','none');
    colorbar; colormap(gca, jet); axis equal; hold on;
    for p_idx = 1:size(pts,1)
        zv = str2double(pts(p_idx,1)); yv = str2double(pts(p_idx,2));
        val = getVal(zv, yv);
        plot(zv, yv, 'ko', 'MarkerFaceColor', 'w');
        text(zv, yv, sprintf('  %s\n  Val: %.1f', pts(p_idx,3), val(m)), 'FontSize', 8, 'FontWeight', 'bold');
    end
    title(titles(m)); xlabel('Z'); ylabel('Y');
end

% ------------------------------------------------------------------
% NOTE: Helper functions (getLayupFromUser, calcABDForLayup, calcQbar, Tcapmatrix,
% getMaterialProperties, getMaterialStrengths, getTsaiWuCoefficients,
% computeInitialAndProgressiveFailure, etc.) must exist unchanged below
% or be available on the MATLAB path.
% ------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions are unchanged below (getLayupFromUser, calcABDForLayup, computeInitialAndProgressiveFailure, etc.)
% (copy your existing helper functions here unchanged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local helper functions (kept from your original script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Lay = getLayupFromUser()
    % Prompts for layup: angle count ... terminated by '000'
    prompt = ['Enter ply sequence (angle count pairs) then 000 to end.\n', ...
              'Example: 0 5 45 10 90 2 000\nEnter: '];
    raw = input(prompt,'s');
    toks = strsplit(strtrim(raw));
    nums = [];
    for k=1:length(toks)
        tk = toks{k};
        if strcmp(tk,'000')
            break;
        end
        val = str2double(tk);
        if isnan(val)
            error('Invalid token "%s" in layup input.', tk);
        end
        nums(end+1) = val; %#ok<AGROW>
    end
    if isempty(nums) || mod(length(nums),2)~=0
        error('Layup must contain angle-count pairs before 000 terminator.');
    end
    angles = nums(1:2:end);
    counts = nums(2:2:end);
    plies_angles = [];
    for i=1:length(angles)
        plies_angles = [plies_angles repmat(angles(i),1,counts(i))]; %#ok<AGROW>
    end
    Lay.angles = plies_angles;
    Lay.counts = counts;
    Lay.totalPlies = sum(counts);
    Lay.t_ply = input('Enter ply thickness (inches): ');
    Lay.totalThickness = Lay.totalPlies * Lay.t_ply;
end

function [A,B,D,ABD] = calcABDForLayup(layup, E11, E22, G12, v12, v21)
    % Constructs Q, Qbar per ply using calcQbar (user-provided), and forms A,B,D
    Q11 = E11 / (1 - v12*v21);
    Q12 = v12 * E22 / (1 - v12*v21);
    Q22 = E22 / (1 - v12*v21);
    Q66 = G12;
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];

    theta_deg = layup.angles;
    n = length(theta_deg);
    t = layup.t_ply;
    if n == 0
        error('Layup has zero plies.');
    end
    h_total = n * t;
    z_int = linspace(-h_total/2, h_total/2, n+1);

    Qbar_all = zeros(3,3,n);
    for i=1:n
        Qbar_all(:,:,i) = calcQbar(Q, theta_deg(i));
    end

    A = zeros(3); B = zeros(3); D = zeros(3);
    for i=1:n
        z1 = z_int(i); z2 = z_int(i+1);
        Qb = Qbar_all(:,:,i);
        A = A + Qb * (z2 - z1);
        B = B + 0.5 * Qb * (z2^2 - z1^2);
        D = D + (1/3) * Qb * (z2^3 - z1^3);
    end
    ABD = [A B; B D];
end

% ------------------------------------------------------------------
% computeInitialAndProgressiveFailure (quiet-mode option)
% -------------------------------------------------------------------
function [R_initial, R_fail_at_ply, failure_order] = computeInitialAndProgressiveFailure( ...
            n, theta_deg, z_mid, z_interfaces, Q, Qbar_all, ...
            F1, F2, F11, F22, F66, F12, loads, eps0, kappa, verbose)
% Computes initial Tsai-Wu R-values and performs progressive failure (zero-stiffness)
if nargin < 19
    verbose = true;
end

tol = 1e-9;

%% ---------- Initial R-values ----------
R_initial = NaN(n,1);

for i = 1:n
    strain_global = eps0 + z_mid(i) * kappa;     % 3x1
    Tcap = Tcapmatrix(theta_deg(i));             % 3x3
    strain_local = Tcap * strain_global;         % 3x1
    stress_local = Q * strain_local;             % 3x1 (psi)

    s11 = stress_local(1)/1000;  % convert to ksi
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

if verbose
    fprintf('\nStarting progressive failure loop (retain positions, zero stiffness for failed plies)...\n');
end

while true
    % Build modified Qbar stack
    Qbar_mod = zeros(3,3,n);
    for i = 1:n
        if active(i)
            Qbar_mod(:,:,i) = Qbar_all(:,:,i);
        else
            Qbar_mod(:,:,i) = zeros(3);
        end
    end

    % Recompute ABD
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
        if verbose, warning('ABD singular. Stopping progressive failure.'); end
        break;
    end

    % Recompute strains under same load
    strain_curvature = ABD \ loads;
    eps0_iter = strain_curvature(1:3);
    kappa_iter = strain_curvature(4:6);

    % Compute R per ply this iteration
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
        if verbose, fprintf('No valid R among active plies. Stopping.\n'); end
        break;
    end

    minR = min(active_Rs(valid_mask));
    to_fail_mask = (R_iter <= minR + tol) & active;
    to_fail_indices = find(to_fail_mask);

    if isempty(to_fail_indices)
        if verbose, fprintf('No plies identified to fail this iteration (numerical). Stopping.\n'); end
        break;
    end

    for idx = to_fail_indices'
        R_fail_at_ply(idx) = minR;
        failure_order = [failure_order; iteration, idx, minR];
    end

    if verbose
        % Ensure all variables are column vectors by using (:)
        IterTable = table((1:n)', theta_deg(:), z_mid(:), R_iter(:), active(:), ...
        'VariableNames', {'Ply_Number','Angle_deg','z_mid_in','R_value','Active'});
        fprintf('\nIteration %d results:\n', iteration);
        disp(IterTable);
        fprintf('Plies failing this iteration: %s  (stored R_fail = %g)\n', mat2str(to_fail_indices'), minR);
    end

    active(to_fail_indices) = false;

    if ~any(active)
        if verbose, fprintf('All plies have failed. End of progressive failure.\n'); end
        break;
    end

    iteration = iteration + 1;
end

if verbose
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

end
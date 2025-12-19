function [E11, E22, G12, v12, v21] = getMaterialProperties()

    % Function to select composite material properties
    % Returns E11, E22, G12, v12, v21

    disp('Select a material for the composite laminate:');
    disp('1 - Gr/Ep (T300/5208)');
    disp('2 - Boron/Ep');
    disp('3 - Gr/Ep (AS/3501)');
    disp('4 - Scotchply 1002');
    disp('5 - Kevlar49/Epoxy');

    material_choice = input('Enter the material number (1-5): ');

    switch material_choice
        case 1
            E11 = 26.25e6;   
            E22 = 1.49e6;
            G12 = 1.04e6;
            v12 = 0.28;
        case 2
            E11 = 29.59e6;
            E22 = 2.68e6;
            G12 = 0.81e6;
            v12 = 0.23;
        case 3
            E11 = 20.01e6;
            E22 = 1.3e6;
            G12 = 1.03e6;
            v12 = 0.30;
        case 4
            E11 = 5.6e6;
            E22 = 1.2e6;
            G12 = 0.6e6;
            v12 = 0.26;
        case 5
            E11 = 11.02e6;
            E22 = 0.8e6;
            G12 = 0.33e6;
            v12 = 0.34;
        otherwise
            error('Invalid selection. Please enter a number between 1 and 5.');
    end

    % Calculate v21
    v21 = v12 * E22 / E11;

    % Display assigned properties
    fprintf('\nSelected Material Properties:\n');
    fprintf('E11 = %.2e Psi\n', E11);
    fprintf('E22 = %.2e Psi\n', E22);
    fprintf('G12 = %.2e Psi\n', G12);
    fprintf('v12 = %.3f\n', v12);
    fprintf('v21 = %.3f\n', v21);
end

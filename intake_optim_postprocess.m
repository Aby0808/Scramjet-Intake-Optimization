% Post-Processing Script for Optimized Scramjet Intake Designs
% close all
% clc

%% Load neccessary data
% load bl shape paramenter tables
global BL_SHAPE_PARAM FRSTM_TH_PARAM

BL_SHAPE_PARAM = table2array(readtable("bl_shap_param.dat"));

% Load the optimized results
load('optimization_results.mat'); % Contains x (designs) and fval (objective values)

% Load the freestream and some fixed throat requirements

FRSTM_TH_PARAM = FUN_set_freestream_throat_params();

alpha = FRSTM_TH_PARAM(1);
M_oo = FRSTM_TH_PARAM(2);
P_oo = FRSTM_TH_PARAM(3);
T_oo = FRSTM_TH_PARAM(4);
% M_th = 2.1;
m_dot = FRSTM_TH_PARAM(5);
h_th = FRSTM_TH_PARAM(6);
T_th = FRSTM_TH_PARAM(7);

%% Extract parameters for comparison

fprintf('\nExtracting design parameter data...');

% Extracting objectives
Drag = fval(:,1);
PressureRecovery = -fval(:,2);

% Extracting PR, Tth, Mth, length, L/D, SI 

num_des = size(fval(:,1),1);
des_param = zeros(num_des,6);
for i = 1:num_des
    [~,~,~,param] = FUN_objective_function(alpha, M_oo, P_oo, T_oo, m_dot, h_th, T_th, ...
        x(i,:), 'n', 10);
    des_param(i,:) = param;
    if mod(i,20)==0
        fprintf('.');
    end
end

fprintf('\tdone!\n');

% % 1. Pareto Front: Drag vs Intake Length
% figure;
% scatter(Drag, IntakeLength, 'filled');
% xlabel('Drag'); ylabel('Intake Length');
% title('Pareto Front: Drag vs Intake Length'); grid on;
% 
% % 2. Pareto Front: Drag vs Intake Width
% figure;
% scatter(Drag, IntakeWidth, 'filled');
% xlabel('Drag'); ylabel('Intake Width');
% title('Pareto Front: Drag vs Intake Width'); grid on;
% 
% % 3. Pareto Front: 3D View (Drag, Length, Width, Color = Pressure Recovery)
% figure;
% scatter3(Drag, IntakeLength, IntakeWidth, 50, PressureRecovery, 'filled');
% xlabel('Drag'); ylabel('Intake Length'); zlabel('Intake Width');
% title('3D Pareto Front (Color = Pressure Recovery)');
% colorbar; grid on;

% Filter out penalty cases
% Normalize each column for better visualization
%normalized_fval = normalize(filtered_fval, 'range'); % Normalize to [0, 1] range

% Define labels for each column
column_labels = {'Drag', 'TPR', 'PR', 'Tth', 'Mth', 'lenght', 'L/D', 'SI'};

% Plot the parallel coordinates
figure;
parallelplot([fval';des_param']', 'CoordinateTickLabels', column_labels, 'LineWidth', 1.5);
title('Parallel Coordinate Plot for Optimized Designs');
%grid on;


% % 5. Extracting Best Trade-Off Designs
% best_index = find(Drag < 0.5 * max(Drag) & IntakeLength < 0.5 * max(IntakeLength) & IntakeWidth < 0.5 * max(IntakeWidth));
% best_designs = x(best_index, :);
% best_tradeoff_values = fval(best_index, :);
% 
% disp('Best Trade-Off Designs (Ramp Angles and Corresponding Trade-Off Values):');
% for i = 1:length(best_index)
%     disp(['Design ', num2str(i), ':']);
%     disp('Ramp Angles:');
%     disp(best_designs(i, :));
%     disp('Trade-Off Values (Drag, L/D, Pressure Ratio, Length, Width, Pressure Recovery):');
%     disp(best_tradeoff_values(i, :));
%     disp('-----------------------------------');
% end
% 
% % 6. Extracting the Absolute Best Design for Minimum Drag
% [~, min_drag_index] = min(Drag);
% best_drag_design = x(min_drag_index, :);
% best_drag_values = fval(min_drag_index, :);
% 
% disp('Best Design for Minimum Drag:');
% disp('Ramp Angles:');
% disp(best_drag_values);

% Post-Processing Script for Optimized Scramjet Intake Designs
close all
clc
% Load the optimized results
load('optimization_results1.mat'); % Contains x (designs) and fval (objective values)
% Define penalty threshold (adjust based on your penalty logic)
penalty_threshold = 1e5; 
valid_indices = all(fval(:,1) < 12000, 2); % Exclude rows where any column exceeds threshold
fval = fval(valid_indices, :); % Only valid designs
%fval(:,7) = round(fval(:,7));
x = x(valid_indices, :);

% Extracting objectives
Drag = fval(:,1);
Lift_to_Drag = -fval(:,2);
PressureRatio = -fval(:,3);
IntakeLength = fval(:,4);
IntakeWidth = fval(:,5);
PressureRecovery = -fval(:,6);

% 1. Pareto Front: Drag vs Intake Length
figure;
scatter(Drag, IntakeLength, 'filled');
xlabel('Drag'); ylabel('Intake Length');
title('Pareto Front: Drag vs Intake Length'); grid on;

% 2. Pareto Front: Drag vs Intake Width
figure;
scatter(Drag, IntakeWidth, 'filled');
xlabel('Drag'); ylabel('Intake Width');
title('Pareto Front: Drag vs Intake Width'); grid on;

% 3. Pareto Front: 3D View (Drag, Length, Width, Color = Pressure Recovery)
figure;
scatter3(Drag, IntakeLength, IntakeWidth, 50, PressureRecovery, 'filled');
xlabel('Drag'); ylabel('Intake Length'); zlabel('Intake Width');
title('3D Pareto Front (Color = Pressure Recovery)');
colorbar; grid on;

% Filter out penalty cases
% Normalize each column for better visualization
%normalized_fval = normalize(filtered_fval, 'range'); % Normalize to [0, 1] range

% Define labels for each column
column_labels = {'Drag', 'L/D', 'PR', 'Length', 'Width', '1st Ramp Length', 'TPR', 'side_fence_flag'};

% Plot the parallel coordinates
figure;
parallelplot([fval, x(:,5)], 'CoordinateTickLabels', column_labels, 'LineWidth', 1.5);
title('Parallel Coordinate Plot for Optimized Designs');
%grid on;


% 5. Extracting Best Trade-Off Designs
best_index = find(Drag < 0.5 * max(Drag) & IntakeLength < 0.5 * max(IntakeLength) & IntakeWidth < 0.5 * max(IntakeWidth));
best_designs = x(best_index, :);
best_tradeoff_values = fval(best_index, :);

disp('Best Trade-Off Designs (Ramp Angles and Corresponding Trade-Off Values):');
for i = 1:length(best_index)
    disp(['Design ', num2str(i), ':']);
    disp('Ramp Angles:');
    disp(best_designs(i, :));
    disp('Trade-Off Values (Drag, L/D, Pressure Ratio, Length, Width, Pressure Recovery):');
    disp(best_tradeoff_values(i, :));
    disp('-----------------------------------');
end

% 6. Extracting the Absolute Best Design for Minimum Drag
[~, min_drag_index] = min(Drag);
best_drag_design = x(min_drag_index, :);
best_drag_values = fval(min_drag_index, :);

disp('Best Design for Minimum Drag:');
disp('Ramp Angles:');
disp(best_drag_values);

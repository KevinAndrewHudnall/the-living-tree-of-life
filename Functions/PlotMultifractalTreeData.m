% PlotMultifractalTreeData

% This script generates the plots used in the manuscript. To generate all
% plots, "BuildMultifractalTree.m" or BuildMultifractalTreeForVisual.m"
% must first be ran and the workspace variables retained.
% "CalculateFractalDims.m" must also first be ran and the workspace
% variables retained.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Leaf_Scale_Matrix S, the paths %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on 
for i = 1:size(S, 2)
    plot(1:size(S, 1), S(1:end, i), 'LineWidth', 1);
end
set(gca, 'yscale', 'log');
title('Scale random variable \beta');
xlabel('Iteration');
ylabel('scale \beta');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Scale of leaves in final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
hold on
plot(1:size(S, 2), S(end, :));
title('Scale of leaves in final iterate');
xlabel('Leaf i');
ylabel('Scale \beta');
set(gca, 'yscale', 'log');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot H, the path entropies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
hold on 
for i = 1:size(H, 2)
    plot(1:size(H, 1), H(1:end, i));
end
title('Pathwise entropy of random variable \beta');
xlabel('Iteration');
ylabel('Path entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Entropies of leaves in final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
hold on
plot(1:size(H, 2), H(end, :));
title('Entropies of leaves in final iterate');
xlabel('Leaf i');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation dimension of every path at the final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
plot(1:size(D_F, 2), D_F(end, :));
title('correlation dimension of every path');
xlabel('Leaf i');
ylabel('D_F');
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical plots for correlation dimension %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the regression for the randomly chosen path from script
% "CalculatFractalDims.m"
figure(6)
hold on
scatter(log(ep_Rand_Member), log(C_ep_Rand_Member), 20, 'filled');
% Get coefficients of a linear fit to the data.
Coeffs = polyfit(log(ep_Rand_Member), log(C_ep_Rand_Member), 1);
% Create x-axis
xFit = linspace(min(log(ep_Rand_Member)), max(log(ep_Rand_Member)), 1000);
% Get the estimated yFit value.
yFit = polyval(Coeffs , xFit);
plot(xFit, yFit, '--', 'LineWidth', 2); % Plot fitted line.
title('Regression to determine D_F for a randomly chosen path')
xlabel('log(\epsilon)');
ylabel('log(C(\epsilon))');

% Plot a histogram of the coefficients of determination at the final
% iterate
figure(7)
histogram(R_Squared(end, :));
xlabel('R^2')

% Plot the average coefficient of determination for at each iterate
figure(8)
errorbar(1:size(Avg_R2_Vector(3:end), 2), Avg_R2_Vector(3:end),...
    Std_Devs(3:end), Std_Devs(3:end))  
title('Average coefficient of determination at each iterate');
xlabel('Number of points N');
ylabel('Average R^2');

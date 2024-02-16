% PlotMultifractalTreeData

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
xlabel('Leaf number');
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
xlabel('Leaf number');
ylabel('Entropy [nats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Information dimension of every path at the final iterate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
plot(1:size(D_F, 2), D_F(end, :));
title('Information dimension of every path');
xlabel('Leaf number');
ylabel('D_F');
axis tight



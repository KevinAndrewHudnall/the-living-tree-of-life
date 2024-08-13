function [S, P, H, Leaf_Coords] = ...
    BuildMultifractalTreeVisualFn(Max_Offspring, Max_Gens, ITERATIONS, ...
    ScaleIntervals)

% This function builds a multifractal tree by implementing the iterated
% function system. The main objects returned are S, a matrix of scales
% representing all paths in the tree of life; P, a matrix indicating
% the progeny had by each member at each iteration; and H, a matrix of
% entropy values. S, P, and H will all have as many rows as there are
% iterates of the system. P will have as many columns as there are paths in
% the second to last iterate since the final iterate does not birth any
% children. S and H will be the same size. The function takes as inputs the
% maximum number of offspring possible and the maximum number of 
% generations possible for an iteration of the Galton-Watson branching
% process, as well as the number of iterations. The generated tree is
% mapped to the unit interval and plotted.

%%%%%%%%%%%%%%%%%%%%% LEVEL 1: generate the first tree %%%%%%%%%%%%%%%%%%%%

% Determine the scale of the first tree
if (~isempty(ScaleIntervals))
    Scale = ScaleIntervals(1) * rand(1);
else
    Scale = rand(1, 1);
end

% Build the first tree
Tree = MakeRandomTreeForVisual(Scale, Max_Offspring, Max_Gens);

% Inititalize cell arrays to hold tree data
Leaves = cell(1, ITERATIONS); % Holds the scales of all leaves of all iterates.

Progenies = cell(1, ITERATIONS); % Holds the number of progeny of all
                                 % leaves at each iterate

% Store the data for the first tree (iterate 1)
Leaves{1} = ones(1, size(Tree.Leaf_Coords, 1)) * Scale; 

Leaf_Coords{1, 1} = Tree.Leaf_Coords; % Holds the unit square coordinates
                                      % of terminal nodes.
                                      
Node_Coords{1, 1} = Tree.Node_Coords; % Holds the unit square coordinates
                                      % of all nodes (terminal and internal)
                                      
Parent{1, 1} = Tree.Parent; % Holds the linear index of the parent members
                            % of the cells of Node_Coords

Progenies{1} = size(Tree.Leaf_Coords, 1); 

% Continue building trees until the final iterate is reached
n = 1;
while (n < ITERATIONS)
    for i = 1:size(Leaves{n}, 2) % loop over the leaves in the nth iterate
        
        % Randomly select a scale for each leaf in the nth iterate
        if(~isempty(ScaleIntervals))
            Scale = ScaleIntervals(n) * rand(1);
        else
            Scale = Leaves{n}(i) * rand(1, 1);
        end
        
        % Generate a random tree for each one of these scales
        [Tree] =...
            MakeRandomTreeForVisual(Scale, Max_Offspring, Max_Gens);
        
        % Store the leaves of the newly generated trees as iterate (n + 1)        
        Leaves{n + 1} =...
            [Leaves{n + 1}, ones(1, size(Tree.Leaf_Coords, 1)) * Scale];
        
        % Update the number of progeny along each path as iterate (n + 1)
        Progenies{n + 1} = [Progenies{n + 1}, size(Tree.Leaf_Coords, 1)];
        
        % Store additional tree info of new trees as iterate (n + 1)
        Leaf_Coords{i, n + 1} = Tree.Leaf_Coords;
        Node_Coords{i, n + 1} = Tree.Node_Coords;
        Parent{i, n + 1} = Tree.Parent;

    end
    n = n + 1;
end

% Move the downscale trees so they take the positions of upscale leaves
for i = 1:(size(Leaf_Coords, 2) - 1)
    for j = 1:sum(~cellfun(@isempty, Node_Coords(:, i)), 1)
        
        if (j > 1)
            Marker = Marker + size(Leaf_Coords{j - 1, i}, 1);
        elseif (j == 1)
            Marker = 0;
        end

        for k = 1:size(Leaf_Coords{j, i}, 1)
            
            % Get the distance from the upscale leaf to the downscale root
            Distance = Leaf_Coords{j, i}(k, :) - ...
                Node_Coords{k + Marker, i + 1}(1, :);
            
            % Move the downsacale tree by changing the coordinates
            Node_Coords{k + Marker, i + 1}(:, 1) = ...
                Node_Coords{k + Marker, i + 1}(:, 1) + Distance(1, 1);
            Node_Coords{k + Marker, i + 1}(:, 2) = ...
                Node_Coords{k + Marker, i + 1}(:, 2) + Distance(1, 2);
            Leaf_Coords{k + Marker, i + 1}(:, 1) = ...
                Leaf_Coords{k + Marker, i + 1}(:, 1) + Distance(1, 1);
            Leaf_Coords{k + Marker, i + 1}(:, 2) = ...
                Leaf_Coords{k + Marker, i + 1}(:, 2) + Distance(1, 2);
        end
    end
end

% At this point, sufficient data has been generated to plot the structure
% on the unit square. For data analysis however, it is useful to obtain a
% scale matrix S and a progeny matrix P

% Organize Scales and Progenies into matrices S and P
tally = cell(size(Progenies, 2), 1);
for i = 1:size(Progenies, 2) % Loop over all rows in Progenies
    for j = 1:size(Progenies{i}, 2) % Loop over all columns in Progenies
        
        % Count where we start column-wise in our repititions
        count_1 = sum(Progenies{i}(1:(j - 1))) + 1;
        
        % Count where we stop column-wise in our repititions relative to
        % our start location
        count_2 = sum(Progenies{i}(1:j));
        
        % Count the total number of repitions we need; initialize to 1
        count_3 = 1;
        
        clear value{i}
        
        % Loop over all cells of Progenies in iterate (i + 1)
        for k = (i + 1):(size(Progenies, 2))

            if (j == 1) % If we're on the first column
                
                % The total number of repeats equals the stop location
                count_3 = count_2;
                
                % Determine the next start-relative stop location
                count_2 = sum(Progenies{k}(count_1:count_2));
                
                % Start location is 1
                count_1 = 1;
                
            else
                % Number of repeats equals the start-relative stop location
                count_3 = count_2;
                
                temp = count_1; % Set a temporary variable 
                
                if (k < size(Progenies, 2)) % If we're not at the last row
                    
                    % Update start location
                    count_1 = value{i}(1, k + 1) + 1;
                    
                    % Update start-relative stop location
                    count_2 = count_1 + sum(Progenies{k}(temp:count_2))...
                        - 1;
                end
            end
            % Store the repeats required in a cell array
            value{i}(1, k) = count_3;
        end
        % Store the repeats required in a cell array
        tally{i}(1, j) = count_3;
    end
end

% Now use tally to build P
P = zeros(size(Progenies, 2), size(Progenies{end}, 2)); % Initialize P

% Generate the first row of P
P(1, 1:size(Progenies{end}, 2)) = ones(1, size(Progenies{end}, 2)) *...
    Progenies{1}(1, 1);

% Generate rows 2 through (last - 1) of P
for i = 2:(size(Progenies, 2) - 1)
    Continue = [];
    Begin = ones(1, tally{i}(1)) * Progenies{1, i}(1, 1);
    for j = 2:size(tally{i}, 2)
        
        Continue = [Continue, ones(1, (tally{i}(j) - tally{i}(j - 1))) *...
           Progenies{1, i}(1, j)];
        
    end
    P(i, 1:size(Progenies{end}, 2)) = [Begin, Continue]; 
end

% Generate the last row of P
P(end, 1:size(Progenies{end}, 2)) = Progenies{end}(1, :);

% Use Progenies to build S
S = cell(size(Leaves, 2), 1);
Temp_1 = Leaves;
for p = 0:(size(Leaves, 2) - 2)
    Temp_2 = cell(1, size(Leaves, 2) - p);
    for i = (size(Leaves, 2) - p):-1:2
        for j = 1:size(Progenies{1, i + p}, 2)

            Temp_2{i - 1} =...
                [Temp_2{i - 1}, ones(1, Progenies{1, i + p}(j)) *...
                Temp_1{1, i - 1}(j)];
        end
    end
    Temp_1 = Temp_2;
    S{(size(Leaves, 2) - (p + 1)), 1} =...
        Temp_2{1, (size(Leaves, 2) - (p + 1))};
    clear Temp_2
end
S{size(Leaves, 2), 1} = Leaves{1, end};
S = cell2mat(S);

% Every leaf in a final subtree is at the same scale, so lump them together
% to eliminate duplicates
S = unique(S', 'rows', 'stable')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Entropy at each level for each path %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = zeros(size(S));
for i = 1:size(S, 2)
    for j = 2:size(S, 1)
        H(j, i) = log(S(j - 1, i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the multifractal tree %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
sz = 50;    % Pick a plot size for the leaves
for i = 1:size(Leaf_Coords, 1)
    if(~isempty(Leaf_Coords{i, ITERATIONS}))
        
        %Plot the leaves
        scatter(Leaf_Coords{i, ITERATIONS}(:, 1),...
            Leaf_Coords{i, ITERATIONS}(:, 2), sz, 'filled', 'k');
    end
    
    for j = 1:size(Leaf_Coords, 2)
        if(~isempty(Leaf_Coords{i, j}))

            %Plot the edges
            for w = 2:size(Parent{i, j}, 2)

                plot([Node_Coords{i, j}(Parent{i, j}(w), 1),...
                    Node_Coords{i, j}(w, 1)],...
                    [Node_Coords{i, j}(Parent{i, j}(w), 2),...
                    Node_Coords{i, j}(w, 2)], 'k');
                set(findall(gca, 'Type', 'Line'), 'LineWidth', 2);

            end
        end
    end
end
axis equal
end
function Tree =...
    MakeRandomTreeForVisualization(Scale, Max_Offspring, Max_Gens)
% PRECONDITION:
% 
%       Scale: The scale of the entire generated tree
%       Max_Offspring: Max number of offspring possible at each generation
%       Max_Gens: Max number of generations possible
% 
% POSTCONDITION: 
% 
%       Tree: a structure containing fields:
% 
%                   Node_Coords: Coordinates of internal nodes.
%                   Root_Coord: Coordinates of the root node.
%                   Leaf_Coords: Coordinates of the leaves
%                   Parent: Linear indices of parents

Done = 0;
while(Done == 0)
    
    % Generate offspring probabilies
    Probs = rand(1, Max_Offspring + 1);
    
    Cum_Probs = [cumsum(Probs) 1]; %Cumulative probabilities

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the branching %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    Parent = 0;
    Gens = randi(Max_Gens); % Determine number of generations
    Child_Count = 1;
    Next_Gen = cell(1, Gens);
    Gens_To_Extinction = 0;

% Galton-Watson Branching.
% Modified from Ingemar Kaj Raimundas Gaigalas(2024)
% https://www.mathworks.com/matlabcentral/fileexchange/2516-random-trees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:Gens % Branching at each generation
        
        % Determine the linear indices of parents in previous generation.
        Parent_Index = length(Parent) - Child_Count + 1:length(Parent);
        
        % Update Parent_Count to continue getting linear indices of parents
        Parent_Count = Child_Count;
        
        % Generate the offspring distribution
        Child_Dist = rand(1, Parent_Count);

        Child_Count = 0;
        for j = 1:Max_Offspring
            
            % Set Index to contain all values of Parent_Index that are
            % greater than the jth value of Cum_Probs but less than the
            % (j+1)th value of Cum_Probs.
            Index = Parent_Index((Child_Dist > Cum_Probs(j)) &...
                (Child_Dist <= Cum_Probs(j + 1)));

            if (~isempty(Index)) % If at least one parent birthed j kids
                
                % Get the number of offspring for each parent by making a
                % vector of Parent concatenated with the vector Index
                % repeated j times.
                Parent = [Parent repmat(Index, 1, j)];
                
                % Update Child_Count with additional children
                Child_Count = Child_Count + length(Index) * j;
            end
        end
        Next_Gen{i} = Parent;
        
        % If the tree has gone extinct, record the generations to
        % extinction and break from the loop.
        if (Child_Count == 0)
            Gens_To_Extinction = i;
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Process tree results %
    %%%%%%%%%%%%%%%%%%%%%%%%

    % If the tree didn't go extinct in the first generation
    if (Gens_To_Extinction ~= 1)
        
        [x, y, ~, ~] = treelayout(Parent); % Plot the tree
        Node_Coords = cat(2,x',y'); % Get the node coordinates
        
        Tree = struct; % Initialize Tree as a structure
        
        Node_Coords = Scale * Node_Coords; % Scale the tree

        Tree.Node_Coords = Node_Coords; % Record the node coordinates
        
        Tree.Root_Coord = Node_Coords(1, :); % Record the root coordinates
        
        % Determine leaves
        Parent = rem(Parent + length(Parent), length(Parent) + 1) + 1; % change all 0s to n+1s
        Tree.Parent = Parent;
        Isaleaf = ones(1, length(Parent) + 1);
        Isaleaf(Parent) = zeros(length(Parent), 1);
        
        % Get the leaf coordinates
        Leaf_Coords = zeros(size(Node_Coords));
        for j = 1:length(Isaleaf)
            if (Isaleaf(j) == 1)
                Leaf_Coords(j, 1) = Node_Coords(j, 1);
                Leaf_Coords(j, 2) = Node_Coords(j, 2);
            end
        end
        
        Leaf_Count = 1;
        while (size(Leaf_Coords, 1) > Leaf_Count)
            if (Leaf_Coords(Leaf_Count, 1) == 0 && Leaf_Coords(Leaf_Count, 2) == 0)
                Leaf_Coords(Leaf_Count, :) = [];
            else
                Leaf_Count = Leaf_Count + 1;
            end
        end

        
        Tree.Leaf_Coords = Leaf_Coords; % Record the leaf coordinates
        
        Done = 1;
    end
    
end
end
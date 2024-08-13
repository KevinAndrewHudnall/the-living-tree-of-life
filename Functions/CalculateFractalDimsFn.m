function [D_F, R_Squared, Avg_R2_Vector, Std_Devs, Random_Member, C_ep_Rand_Member,...
    ep_Rand_Member, D_F_Rand_Member, R_Squared_Rand_Member] = CalculateFractalDimsFn(S)

% This script calculates the fractal dimension (i.e., the correlation
% dimension) for each path through the tree of life. The script calculates the
% correlation dimension of each path at each iterate, though values only
% become reliable once the number of iterates becomes large.The outputs are
% as follows:
%
%            D_F: a matrix the same size as S holding all the fractal
%            dimensions of each path at each iterate.
%
%            R_Squared: a matrix the same size as S holding all the
%            coefficients of determination of each path at each iterate
%
%            Avg_R2_Vector: a vector with as many columns as there are
%            iterations, giving the coefficient of determination averaged
%            over all paths at each iterate 
%
%            Random_Member: a randomly selected path
%
%            C_ep_Rand_Member: Correlation integral values for the randomly
%            selected path.
%
%            ep_Rand_Member: vector of epsilon values (spacings between
%            points) for the randomly selected member.
%
%            D_F_Rand_Member: Fractal dimension for the randomly selected
%            member.
%
%            R_Squared_Rand_Member: Coefficient of determination for the
%            randomly seleceted member.

% Compute the minimum spacing between leaves for each column
temp = min(abs(diff(S)));
% Compute the overall minimum spacing
Min_Spacing = min(temp);
% Compute the adjusted minimum spacing
Min_Spacing = abs(floor(log10(Min_Spacing))) + 2;
 
% Generate a vector of the logarithmically spaced epsilon values
ep_Orig = [logspace(10, 1, 100) / 1.0e+10, logspace(10, 1, 100) /...
    1.0e+20, logspace(10, 1, 100) / 10^Min_Spacing];

D_F = zeros(size(S));
ep_Orig = ep_Orig(2:end);
R_Squared = zeros(size(S));
for k = 1:size(S, 2)
    for p = 2:size(S, 1)
        
        % For each epsilon value determine the number of pairs of points
        % whose separation distance is less than epsilon
        for i = 1:size(ep_Orig, 2)
            All_Combos = nchoosek(S(1:p, k), 2);
            Differences = All_Combos(:, 1) - All_Combos(:, 2);
            
            % Store the number in C_ep
            C_ep(i) = sum(Differences < ep_Orig(i));
        end
        C_ep = C_ep(C_ep ~= 0); % Eliminate entries equal to zero.
        C_ep = C_ep / p^2; % Divide by the number of points squared
        ep_Altered = ep_Orig(C_ep ~= 0); % Get a truncated version of
                                         % ep_orig to match the truncated
                                         % version of C_ep
                                         
        Lin_Model = fitlm(log(ep_Altered), log(C_ep)); % Fit a linear model
        D_F(p, k) = Lin_Model.Coefficients{2, 1}; % The fractal dimension
                                                  % is equal to the slope.
        % Get the coefficient of determination for each path.                                          
        R_Squared(p, k) = Lin_Model.Rsquared.Ordinary; 
    end
end

% Get the average R-squared values for the paths at all iterates, as well
% as the standard deviations
Avg_R2_Vector = zeros(1, size(S, 1));
Std_Devs = zeros(1, size(S, 1)); 
for i = 1:size(R_Squared, 1)
    temp = unique(R_Squared(i, :), 'stable');
    Std_Devs(i) = std(temp);
    Avg_R2_Vector(i) = sum(temp) / size(temp, 2);
end

% Look at the fit for a randomly chosen path
Random_Member = randi(size(S, 2)); % Choose a random path
for i = 1:size(ep_Orig, 2)
    All_Combos = nchoosek(S(:, Random_Member), 2);
    Differences = All_Combos(:, 1) - All_Combos(:, 2);
    C_ep_Rand_Member(i) = sum(Differences < ep_Orig(i));
end
C_ep_Rand_Member = C_ep_Rand_Member(C_ep_Rand_Member ~= 0); % Eliminate entries equal to zero.
C_ep_Rand_Member = C_ep_Rand_Member / (size(S, 1)^2); % Divide by the number of points squared
ep_Rand_Member = ep_Orig(C_ep_Rand_Member ~= 0); % truncate
Lin_Model = fitlm(log(ep_Rand_Member), log(C_ep_Rand_Member)); % Regression
D_F_Rand_Member = Lin_Model.Coefficients{2, 1}; % Get fractal dimension
R_Squared_Rand_Member = Lin_Model.Rsquared.Ordinary; % Get the coefficient
                                                     % of determination
end
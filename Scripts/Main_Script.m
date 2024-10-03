% Main_Script

% This script generates a multifractal tree of life and calculates the
% entropy of all its paths as well as the fractal dimensions of all its
% paths. The tree is not mapped to the unit square, meaning visualization
% of the structure on the unit square is not carried out.

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%

% The only parameters are the maximum number of offspring and the maximum
% number of generations for the Galton-Watson branching process, as well as
% total number of system iterations. The structure blows up rapidly,
% setting practical limitations on the parameters.

MaxOffspring = 3; % Max number of offspring for every generated random tree
MaxGens = 2; % Max number of generations for every random tree
ITERATIONS = 20; % Total number of system iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the multifractal tree %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S,P, H] = BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the fractal dimensions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D_F, R_Squared, Avg_R2_Vector, Std_Devs, Random_Member, ...
    C_ep_Rand_Member, ep_Rand_Member, D_F_Rand_Member, ...
    R_Squared_Rand_Member] = CalculateFractalDimsFn(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the the generalized Hurst indices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the q-values to calculate the generalized Hurst exponent
q_Values = [-5, -3, -1, 0, 1, 2, 3, 5]; 

% Define different box sizes to be analyzed
Box_Sizes = [3, 4, 5, 6, 7];

% Call function MfDfaFn
Generalized_Hurst_values = MfDfaFn(S, q_Values, Box_Sizes);


% Main_Script_Visualization

% This script builds a multifractal tree that can be visualized on the unit
% square. There are two options. One can execute the script and build a
% multifractal tree where the scales of the iterates are decreasing but
% otherwise completely random. Or one can specifiy scale intervals at which
% the iterations are generated as random decreasing multiples of those
% scale intervals. The first choice corresponds to plotting a tree as
% generated in script BuildMultifractalTree. The second case allows for a
% clean visualization since the iterates are widely spaced enough so that
% branches of one generation do not cross the branches of another
% generation (that is, setting widely spaced scale intervals allows the
% structure to be generated so that downscale trees appear completely
% contained within upscale leaves). Setting variable ScaleIntervals to
% empty will execute the first option (i.e., entirely random). Settng 
% variable ScaleIntervals to the scales shown on line 28 will execute the
% widely spaced option.

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%

Max_Offspring = 3; % Max number of offspring for each generated random tree
Max_Gens = 2; % Max number of generations for every random tree
ITERATIONS = 10; % Total number of system iterations

% Set spacing of iterations if desired; set to empty otherwise
ScaleIntervals =...
    [10000, 1000, 100, 10, 1, 0.1, 0.01, 0.0001, 0.00001, 0.000001];
% ScaleIntervals = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the multifractal tree and plot it %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S, P, H, Leaf_Coords] = BuildMultifractalTreeVisualFn(Max_Offspring, ...
    Max_Gens, ITERATIONS, ScaleIntervals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate fractal dimensions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom through the multifractal tree %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine path along which to zoom. 
% The final leaf that is the target of the zoom must be selected. Here we
% select one randomly from Leaf_Coords in the final
% iterate, but the selection does not have to be random.

% Randomly select a final leaf from Leaf_Coords
Random_1 = randi(size(Leaf_Coords, 1));
Random_2 = randi(size(Leaf_Coords{Random_1, size(Leaf_Coords, 2)}, 1));
Point_0 = Leaf_Coords{Random_1, size(Leaf_Coords, 2)}(Random_2, :);

% A contraction factor (i.e, the zoom step size) must be specified.
Contraction_Factor = 0.95;

% Execute the zoom
FractalZoomFn(Point_0, Contraction_Factor);

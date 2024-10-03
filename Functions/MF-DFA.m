% Multifractal DeTrended Fluctuation Analysis (MF-DFA) 

% Define the q-values to calculate the generalized Hurst exponent
q_Values = [-5, -3, -1, 0, 1, 2, 3, 5]; 

% Define different box sizes to be analyzed
Box_Sizes = [3, 4, 5, 6, 7];

% Get the size of the matrix S
[N, Num_Series] = size(S);

% Initialize a matrix to store the generalized fluctuation function for
% each time series, q, and box size
MF_Fq_Values = zeros(Num_Series, length(Box_Sizes), length(q_Values));

% Loop over each column in S
for j = 1:Num_Series
    % Extract the j-th column
    X = S(:, j);
    
    % Remove the mean from the time series to make it zero-mean
    X = X - mean(X);
    
    % Calculate the cumulative sum
    Y = cumsum(X);
    
    % Loop over each q value
    for q_Idx = 1:length(q_Values)
        q = q_Values(q_Idx);  % Current q value
        
        % Loop over each box size
        for Idx = 1:length(Box_Sizes)
            n = Box_Sizes(Idx);
            Num_Boxes = floor(N / n);  % Number of boxes
            
            % Initialize the fluctuation function for the current box size
            F_nq = zeros(Num_Boxes, 1);
            
            % Loop over each box
            for i = 1:Num_Boxes
                % Define the indices for the current box
                Start_Idx = (i-1) * n + 1;
                End_Idx = i * n;
                
                % Extract the data for the current box
                Box_Data = Y(Start_Idx:End_Idx);
                
                % Fit a linear Trend to the data in the current box
                t = (1:n)';
                p = polyfit(t, Box_Data, 1);  % Linear fit
                Trend = polyval(p, t);
                
                % Calculate the detrended fluctuation for the current box
                F_n = sqrt(mean((Box_Data - Trend).^2));  % RMS of residuals
                
                % Calculate generalized fluctuation function for each q
                if q == 0
                    % For q = 0, use log averaging to avoid singularities
                    F_nq(i) = exp(0.5 * mean(log(F_n.^2)));
                else
                    F_nq(i) = (mean(F_n.^q))^(1/q);
                end
            end
            
            % Avg the generalized fluctuation function over all boxes of size n
            MF_Fq_Values(j, Idx, q_Idx) = mean(F_nq);
        end
    end
end

% Perform log-log fitting to estimate the generalized Hurst exponent 
% for each column of S and q-value
Generalized_Hurst_values = zeros(Num_Series, length(q_Values));

for j = 1:Num_Series
    for q_Idx = 1:length(q_Values)
        log_Box_Sizes = log(Box_Sizes);

        % Extract Fq values for this series and q
        log_Fq = log(squeeze(MF_Fq_Values(j, :, q_Idx)));  
        
        % Fit a linear model to the log-log plot
        p = polyfit(log_Box_Sizes, log_Fq, 1);
        
        % The slope of the log-log plot is the generalized Hurst exponent
        % for the j-th series and q-th moment
        Generalized_Hurst_values(j, q_Idx) = p(1);
    end
end

% Extract AR coefficient matrices from the cell array
A1 = EstMdl.AR{1}; % Coefficients for lag 1
A2 = EstMdl.AR{2}; % Coefficients for lag 2
A3 = EstMdl.AR{3}; % Coefficients for lag 3

% A1 = RestrictedMdl.AR{1}; % Coefficients for lag 1
% A2 = RestrictedMdl.AR{2}; % Coefficients for lag 2
% A3 = RestrictedMdl.AR{3}; % Coefficients for lag 3

k = size(A1, 1); % Number of variables (8 in this case)
p = 3;           % Number of lags (3 in this case)

% Initialize the companion matrix
C = zeros(k*p, k*p);

% Place AR coefficient matrices in the first block row
C(1:k, 1:k*p) = [A1, A2, A3];

% Fill the sub-diagonal with identity matrices
if p > 1
    C(k+1:end, 1:end-k) = eye(k*(p-1));
end

% Compute eigenvalues of the companion matrix
eigenvalues = eig(C);

% Check if all eigenvalues are within the unit circle
if all(abs(eigenvalues) < 1)
    disp('The VAR model is stable.');
else
    disp('The VAR model is unstable.');
end

% Display eigenvalues
disp('Eigenvalues of the companion matrix:');
disp(eigenvalues);
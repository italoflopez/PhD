Mdl = varm(size([out.fac(:,[1 5:8]),out.fac(:,2:4)],2),3);%,diff(out.fac(:,2:4))
EstMdl = estimate(Mdl,[out.fac(:,[1 5:8]),out.fac(:,2:4)]);%,diff(out.fac(:,2:4))

% Copy the AR coefficients from the estimated model
RestrictedAR = EstMdl.AR;

% Restrict specific coefficients to zero
% For example, restrict the second variable's effect on the first variable at lag 1
RestrictedAR{1}(1, [6 7 8]) = 0; % Lag 1, row 1, column 2
RestrictedAR{2}(2, [6 7 8]) = 0; % Lag 2, row 2, column 3
RestrictedAR{3}(2, [6 7 8]) = 0;

% Assign the restricted coefficients back to the model
EstMdl.AR = RestrictedAR;

% Create a new VAR model with restricted parameters
RestrictedMdl = varm(size([out.fac(:,[1 5:8]), out.fac(:,2:4)], 2), 3);
RestrictedMdl.AR = RestrictedAR;

% Assign other required properties from the estimated model
RestrictedMdl.Constant = EstMdl.Constant;    % Constant term
RestrictedMdl.Covariance = EstMdl.Covariance; % Covariance matrix
% Number of variables in the VAR model
numVars = size(out.fac, 2);

% Define presample responses (zeros for simplicity)
Y0 = zeros(3, numVars); % 3

% Specify the sample size
sampleSize = size([out.fac(:,[1 5:8]), out.fac(:,2:4)], 1);
% Compute the IRFs
[IRF, time] = irf(RestrictedMdl, 'NumObs', 180, 'Y0', Y0, 'SampleSize', sampleSize, 'NumPaths', 1);

% Plot the IRFs
figure;
plot(IRF(:, 1, 1)); % IRF for the first variable
%legend('Response to shock 1', 'Response to shock 2');
title('Impulse Response Function');
xlabel('Time');
ylabel('Response');
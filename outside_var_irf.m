Mdl = varm(size([out.fac(:,[1 5:8]),out.fac(:,2:4)],2),3);%,diff(out.fac(:,2:4))
EstMdl = estimate(Mdl,[out.fac(:,[1 5:8]),out.fac(:,2:4)]);%,diff(out.fac(:,2:4))

% Compute IRFs
[IRF, time] = irf(EstMdl, 'NumObs', 60); % Specify the number of periods for the IRF

% Plot the IRFs
figure;
plot(IRF(:, 1, 1)); % IRF for the first variable
%legend('Response to shock 1', 'Response to shock 2');
title('Impulse Response Function');
xlabel('Time');
ylabel('Response');


function data = replaceOutliersWithMovingAverage(data)
    % Detect outliers using the 3-sigma rule
    meanData = mean(data, 'omitnan');
    stdData = std(data, 'omitnan');
    outliers = abs(data - meanData) > 3 * stdData;

    % Calculate the moving average (e.g., window size = 5)
    movingAvg = movmean(data, 30, 'omitnan');

    % Replace outliers with the moving average
    data(outliers) = movingAvg(outliers);
end
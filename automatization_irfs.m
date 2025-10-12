% Define save path
savePath = 'C:\Users\Italo\OneDrive\Documents\PhD\PhD\First\resultados_muestra_total_betas_en_nivel_transformaciones_interanuales_3_rezagos.docx'; % Change this to your desired path

% Create a Word document at the specified path
import mlreportgen.dom.*;
doc = Document(savePath, 'docx');

% Define variables and their corresponding indices
variables = {'Inflation','IP',  'Capacity utilization, total','Employment',...
    'Industrial and Commercial Loans','Non-revolving consumer credit','Unemployment Rate',...
    'Average Hourly Earnings','Depreciation of the US Dollar against the UK Sterling Pound',...
    'M1','M2','S&P 500','Commodities Price Index','Finance Rate 24-month Loans',...
    '30-year mortgage rate minus 10-year Treasury yield'}; % Example for standard IRFs
varIndices = [2, 3, 15,19,76,78,67,63,86,88,89,98,100,94,97]; % Corresponding indices

restrictedVars = {'Inflation','IP',  'Capacity utilization, total','Employment',...
    'Industrial and Commercial Loans','Non-revolving consumer credit','Unemployment Rate',...
    'Average Hourly Earnings','Depreciation of the US Dollar against the UK Sterling Pound',...
    'M1','M2','S&P 500','Commodities Price Index','Finance Rate 24-month Loans',...
    '30-year mortgage rate minus 10-year Treasury yield'}; % Example for restricted IRFs
restrictedIndices = [2, 3, 15,19,76,78,67,63,86,88,89,98,100,94,97]; % Corresponding indices

factorVars = {'Yield curve Level (Nelson-Siegel Beta 1)','Yield Curve Slope (Negative of Nelson-Siegel Beta 2)','Yield curve Curvature (Nelson-Siegel Beta 3)'}; % Example for factor IRFs
factorIndices = [6,7,8]; % Corresponding indices

%% Standard IRFs
for i = 1:length(variables)
    varNum = varIndices(i);
    
    % Create figure
    figure;
    hold on;
    
    % Plot Point Estimate
    pointEstimate = plot(0.25 * squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(varNum,3, :)), '-b');
    
    % Confidence Bands
    upper1 = plot(0.25 * squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(varNum,3, :)) + ...
                  0.25 * squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
    
    lower1 = plot(0.25 * squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(varNum,3, :)) - ...
                  0.25 * squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
    
    upper2 = plot(0.25 * squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(varNum,3, :)) + ...
                  0.25 * squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
    
    lower2 = plot(0.25 * squeeze(irf_vdecomp_out.imp_y_fac_mat_scl(varNum,3, :)) - ...
                  0.25 * squeeze(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');

    % Titles and labels
    title(['Impulse Response Function: ', variables{i}]);
    xlabel('Time period');
    ylabel('Effect');
    legend([pointEstimate, upper1, upper2], ...
           'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
           'Location', 'Best');
    yline(0,'--black','HandleVisibility','off');

    % Save figure as image
    imgFilename = fullfile(tempdir, [variables{i}, '.png']);
    saveas(gcf, imgFilename);
    close(gcf);

    % Add plot to Word document
    append(doc, Paragraph(['Impulse Response Function: ', variables{i}]));
    img = Image(imgFilename);
    img.Width = '6in';  % Set width to 4 inches
    img.Height = '6in'; % Set height to 3 inches
    append(doc, img);

    append(doc, PageBreak());
end

%% Restricted IRFs
for i = 1:length(restrictedVars)
    varNum = restrictedIndices(i);
    
    % Create figure
    figure;
    hold on;
    
    % Plot Point Estimate
    pointEstimate = plot(0.25 * squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(varNum,3, :)), '-b');
    
    % Confidence Bands
    upper1 = plot(0.25 * squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(varNum,3, :)) + ...
                  0.25 * squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
    
    lower1 = plot(0.25 * squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(varNum,3, :)) - ...
                  0.25 * squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
    
    upper2 = plot(0.25 * squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(varNum,3, :)) + ...
                  0.25 * squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
    
    lower2 = plot(0.25 * squeeze(irf_vdecomp_out_restrictions.imp_y_fac_mat_scl(varNum,3, :)) - ...
                  0.25 * squeeze(se_irf_vdecomp_out_restrictions.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');

    % Titles and labels
    title(['Impulse Response Function: ', restrictedVars{i}]);
    xlabel('Time period');
    ylabel('Effect');
    legend([pointEstimate, upper1, upper2], ...
           'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
           'Location', 'Best');
    yline(0,'--black','HandleVisibility','off');

    % Save figure as image
    imgFilename = fullfile(tempdir, [restrictedVars{i}, '.png']);
    saveas(gcf, imgFilename);
    close(gcf);

    % Add plot to Word document
    append(doc, Paragraph(['Impulse Response Function: ', restrictedVars{i}]));
    img = Image(imgFilename);
    img.Width = '6in';  % Set width to 4 inches
    img.Height = '6in'; % Set height to 3 inches
    append(doc, img);

    append(doc, PageBreak());
end

%% Factor IRFs
for i = 1:length(factorVars)
    varNum = factorIndices(i);
    
    if varNum == 7
        % Create figure
        figure;
        hold on;
        
        % Plot Point Estimate
        pointEstimate = plot(-0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)), '-b');
        
        % Confidence Bands
        upper1 = plot(-0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) + ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
        
        lower1 = plot(-0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) - ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
        
        upper2 = plot(-0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) + ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
        
        lower2 = plot(-0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) - ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
    
        % Titles and labels
        title(['Impulse Response Function: ', factorVars{i}]);
        xlabel('Time period');
        ylabel('Effect');
        legend([pointEstimate, upper1, upper2], ...
               'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
               'Location', 'Best');
        yline(0,'--black','HandleVisibility','off');
    
        % Save figure as image
        imgFilename = fullfile(tempdir, [factorVars{i}, '.png']);
        saveas(gcf, imgFilename);
        close(gcf);
    
        % Add plot to Word document
        append(doc, Paragraph(['Impulse Response Function: ', factorVars{i}]));
        img = Image(imgFilename);
        img.Width = '6in';  % Set width to 4 inches
        img.Height = '6in'; % Set height to 3 inches
        append(doc, img);
        append(doc, PageBreak());
    else  
        % Create figure
        figure;
        hold on;
        
        % Plot Point Estimate
        pointEstimate = plot(0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)), '-b');
        
        % Confidence Bands
        upper1 = plot(0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) + ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
        
        lower1 = plot(0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) - ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 1.5, '-y');
        
        upper2 = plot(0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) + ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
        
        lower2 = plot(0.25 * squeeze(irf_factors.imp_y_fac_mat_scl(varNum,3, :)) - ...
                      0.25 * squeeze(se_irf_factors.se_imp_y_fac_mat_scl(varNum,3, :)) * 2, '-r');
    
        % Titles and labels
        title(['Impulse Response Function: ', factorVars{i}]);
        xlabel('Time period');
        ylabel('Effect');
        legend([pointEstimate, upper1, upper2], ...
               'Point Estimate', 'Confidence Bands: 1.5 standard errors', 'Confidence Bands: 2 standard errors', ...
               'Location', 'Best');
        yline(0,'--black','HandleVisibility','off');
    
        % Save figure as image
        imgFilename = fullfile(tempdir, [factorVars{i}, '.png']);
        saveas(gcf, imgFilename);
        close(gcf);
    
        % Add plot to Word document
        append(doc, Paragraph(['Impulse Response Function: ', factorVars{i}]));
        img = Image(imgFilename);
        img.Width = '6in';  % Set width to 4 inches
        img.Height = '6in'; % Set height to 3 inches
        append(doc, img);
        append(doc, PageBreak());        
end

% Save and close the document at the specified path
close(doc);
disp(['Word document saved at: ', savePath]);

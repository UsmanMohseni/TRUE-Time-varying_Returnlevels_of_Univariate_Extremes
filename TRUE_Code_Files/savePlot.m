
% Function to save a plot in multiple formats
function savePlot(plotNumber, fileName, outputFolder)
    % Full file path without extension
    fullPath = fullfile(outputFolder, sprintf('%s_%d', fileName, plotNumber));
    
    % Save in PNG format (high resolution)
    print(fullPath, '-dpng', '-r300'); 

    % Save in TIFF format
    print(fullPath, '-dtiff', '-r300'); 

    % Save in FIG format (MATLAB figure)
    savefig(fullPath); 

    % Save in EPS format (for vector graphics)
    hgexport(gcf, [fullPath, '.eps']); 

    % Save in EMF format (for compatibility with MS Office)
    print(fullPath, '-dmeta');  
end

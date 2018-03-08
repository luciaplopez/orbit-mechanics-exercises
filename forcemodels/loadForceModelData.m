function [data] = loadForceModelData(fileName)

folderName = 'integrator_results';

f = fullfile(folderName, fileName);

data = dlmread(f);

end


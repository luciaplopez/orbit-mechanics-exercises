function [pos] = positionForceModel(fileName)

data = loadForceModelData(fileName);

pos = data(:, 1:4);

end


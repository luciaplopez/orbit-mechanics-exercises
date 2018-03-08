function [] = plotDensityForceModel(fileName1, cases)

data1 = loadForceModelData(fileName1);

t = data1(:, 1) / 3600;

density = data1(:, 11);

fig = figure('visible','off');
hold on

plot(t, density(:, 1), 'LineWidth', 1.2)

title(['Air density for case ', cases])

xlabel('time (h)')
ylabel('air density (kg/m^3)')

grid on
grid minor

folderName = 'figs';
figName = [cases, '.png'];

f = fullfile(folderName, figName);

saveas(fig, f)

end


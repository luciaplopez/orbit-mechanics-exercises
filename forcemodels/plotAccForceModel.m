function [] = plotAccForceModel(fileName1, fileName2, cases)

data1 = loadForceModelData(fileName1);

t = data1(:, 1) / 3600;

if fileName2 ~= ''
    data2 = loadForceModelData(fileName2);
    difference = data2(:, 8:10) - data1(:, 8:10);
else
    difference = data1(:, 8:10);
end

fig = figure('visible','off');
hold on

m = sqrt(sum(difference(:, :).^2, 2));

plot(t, difference(:, 1), 'LineWidth', 1.2)
plot(t, difference(:, 2), 'LineWidth', 1.2)
plot(t, difference(:, 3), 'LineWidth', 1.2)
plot(t, m, 'LineWidth', 1.2)

if fileName2 ~= ''
    title(['Difference of accelerations for ', cases, ' cases'])
else
    title(['Accelerations for ', cases])
end


xlabel('time (h)')
ylabel('acceleration (m/s^2)')

legend('A_x','A_y','A_z','Magnitude');

grid on
grid minor


folderName = 'figs';
figName = [cases, '.png'];

f = fullfile(folderName, figName);

saveas(fig, f)

end


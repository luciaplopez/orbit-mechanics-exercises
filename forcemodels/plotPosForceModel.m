function [] = plotPosForceModel(fileName1, fileName2, cases)

data1 = loadForceModelData(fileName1);
data2 = loadForceModelData(fileName2);

if strcmp(fileName2, 'U_GRACE.txt')
   data2(:, 2:4) = data2(:, 2:4) + 1;
elseif strcmp(fileName2, 'V_GRACE.txt')
   data2(:, 5:7) = data2(:, 5:7) + 0.001;
end    

t = data1(:, 1) / 3600;

difference = data2(:, 1:7) - data1(:, 1:7);

fig = figure('visible','off');

subplot(2,1,1);
hold on

m = sqrt(sum(difference(:, 2:4).^2, 2));

plot(t, difference(:, 2), 'LineWidth', 1.2)
plot(t, difference(:, 3), 'LineWidth', 1.2)
plot(t, difference(:, 4), 'LineWidth', 1.2)
plot(t, m, 'LineWidth', 1.2)

title(['Difference in position ', cases, ' cases'])

xlabel('time (h)')
ylabel('position (m)')

legend('R_x','R_y','R_z','Magnitude');

grid on
grid minor

subplot(2,1,2);
hold on

m = sqrt(sum(difference(:, 5:7).^2, 2));

plot(t, difference(:, 5), 'LineWidth', 1.2)
plot(t, difference(:, 6), 'LineWidth', 1.2)
plot(t, difference(:, 7), 'LineWidth', 1.2)
plot(t, m, 'LineWidth', 1.2)

title(['Difference in velocity ', cases, ' cases'])

xlabel('time (h)')
ylabel('velocity (m/s)')

legend('V_x','V_y','V_z','Magnitude');

grid on
grid minor

folderName = 'figs';
figName = [cases, '.png'];

f = fullfile(folderName, figName);

saveas(fig, f)

end


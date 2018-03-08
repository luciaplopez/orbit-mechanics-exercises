function [] = getPerturbationPosVel(fileName1, fileName2, cases)

data1 = loadForceModelData(fileName1);
data2 = loadForceModelData(fileName2);
difference = data2(:, 1:7) - data1(:, 1:7);

m = sqrt(sum(difference(:, 2:4).^2, 2));

max_m = max(m);
min_m = min(m);
rms_m = rms(m);

disp(['Case considered: ', cases])
disp(['Max of the diff position magnitude: ', num2str(max_m), ' m'])
disp(['Min of the diff position magnitude: ', num2str(min_m), ' m'])
disp(['RSM of the diff position magnitude: ', num2str(rms_m), ' m'])

m = sqrt(sum(difference(:, 5:7).^2, 2));

max_m = max(m);
min_m = min(m);
rms_m = rms(m);

disp([' '])
disp(['Max of the diff velocity magnitude: ', num2str(max_m), ' m/s'])
disp(['Min of the diff velocity magnitude: ', num2str(min_m), ' m/s'])
disp(['RSM of the diff velocity magnitude: ', num2str(rms_m), ' m/s'])

end


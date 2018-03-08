function [] = getPerturbationAcc(fileName1, fileName2, cases)

data1 = loadForceModelData(fileName1);
data2 = loadForceModelData(fileName2);
difference = data2(:, 8:10) - data1(:, 8:10);

m = sqrt(sum(difference(:, :).^2, 2));

max_m = max(m);
min_m = min(m);
rms_m = rms(m);

disp(['Case considered: ', cases])
disp(['Max of the diff acceleration magnitude: ', num2str(max_m), ' m/s^2'])
disp(['Min of the diff acceleration magnitude: ', num2str(min_m), ' m/s^2'])
disp(['RSM of the diff acceleration magnitude: ', num2str(rms_m), ' m/s^2'])

end


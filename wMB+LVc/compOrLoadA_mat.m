function A_mat = compOrLoadA_mat(A_matPath, c, n, image_width, t, radius, angle_sensor, n_angles)

%%% check for model matrix if saved, otherwise build & save it
if exist(A_matPath, 'file')
    disp(['Model matrix found on disk (' A_matPath ') .. Loading into RAM..']);
    A_mat = load(A_matPath);
    A_mat = A_mat.A_mat;
else
    disp(['Model matrix not found on disk (' A_matPath ') .. Calculating & saving to disk..']);
    A_mat = Calculate_MatrixMB_Luis(c, n, image_width, t, radius, angle_sensor, n_angles);
    save(A_matPath, 'A_mat', '-v7.3');
end
disp_last_line("done.");
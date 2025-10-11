close all
clear;

two_cubes

params.label='on';
params.text_dc=1.05;
params.arrow_size=0.1;

params.number=Patch_1_number;
figure(1,'position',[20,100,920,900]);
plot_ctrl_points_3d(Patch_1,params);
view(15,45);

params.number=Patch_2_number;
figure(2,'position',[950,100,920,900]);
plot_ctrl_points_3d(Patch_2,params);
view(15,45);


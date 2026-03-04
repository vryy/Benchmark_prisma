% fid = fopen('time.txt', 'r');
% [step, timedata] = fscanf(fid, '%d %f', [10, Inf]);
% fclose(fid);

clear all;

format long
t = 0.0
hx = 1.0/72
kappa = 2.0e6
rho = 997.5
delta_time = 0.1*hx/sqrt(kappa/rho)
total_time = 5.0

step = 0;
sample_output = 100;
while t < total_time
    t = t + delta_time;

    if mod(step, sample_output) == 0
        filename = ['mesh_72x72_q4_' num2str(t*1.0e3, 12) ".h5"]
        load(filename)

        particles_figure = figure('visible', 'off');
        plot(particles_x, particles_y, 'o', 'MarkerSize', 4);
        %plot(particles_x, particles_y, '.');
        title('particles location');
        grid on;
        %set(gcf, 'Position', [100, 100, 500, 500])
        set(gcf, 'Position', [100, 100, 1500, 1500])
        axis([0, 6, 0, 6]);
        set(gca,'xtick',[0:0.5:6.0]);
        set(gca,'ytick',[0:0.5:6.0]);
        axis equal;
        print(particles_figure, ['particles_' num2str(t*1.0e3, 11) '.jpg']);
        close(particles_figure);
    end

    step = step + 1;
end

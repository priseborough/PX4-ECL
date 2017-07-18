function PlotTrajectory(output,gps_data,folder)

if ~exist(folder,'dir')
    mkdir(folder);
end
plotDimensions = [0 0 210*3 297*3];
figure('Units','Pixels','Position',plotDimensions,'PaperOrientation','portrait');
h=gcf;
set(h,'PaperOrientation','portrait');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

plot(output.position_NED(:,2),output.position_NED(:,1));
grid on;
axis equal;

if isfield(gps_data,'pos_ned')    
    hold on;
    plot(gps_data.pos_ned(:,2),gps_data.pos_ned(:,1),'r');
    hold off;
    legend('estimate','GPS');
end

xlabel('East Position (m)');
ylabel('North Position (m)');

fileName='trajectory.png';
fullFileName = fullfile(folder, fileName);
saveas(h,fullFileName);

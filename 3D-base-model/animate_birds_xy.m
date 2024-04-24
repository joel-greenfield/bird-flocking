%% Animate bird flock on XY plane with obstacles and save to .mp4
function animate_birds_xy(x_h, l_h, center, radius, dt)
    n = 100;
    obstacle = nsidedpoly(1000, 'Center', center, 'Radius', radius);
    figure
    fig = gcf;
    
    video = VideoWriter('xy_obstacle', 'MPEG-4');
    video.FrameRate = 10; 
    open(video)
    
    for i = 1:round(size(x_h,3)/n)
        clf(fig)
        x_l = x_h;
        x_f = x_h;
        idx_l = l_h(:,:,i*n) == 1;
        idx_f = l_h(:,:,i*n) == 0;

        x_l(idx_f,:,i*n) = nan;
        x_f(idx_l,:,i*n) = nan;
        
        scatter(x_l(:,1,i*n),x_l(:,2,i*n),"filled", 'MarkerFaceColor', 'red')
        hold on
        scatter(x_f(:,1,i*n),x_f(:,2,i*n),"filled", "MarkerFaceColor", 'black')
        plot(obstacle, 'FaceColor', 'blue')
        plot(center (1), center(2), "Marker",".", "MarkerSize",30, "Color","blue")
        

        plot3(x_h(:,1,i*100),x_h(:,2,i*100),x_h(:,3,i*100),".",'MarkerSize',30);

        xlim([min(min(x_h(:,1,:))) max(max(x_h(:,1,:)))]); % Set x-axis limits
        ylim([min(min(x_h(:,2,:))) max(max(x_h(:,2,:)))]); % Set y-axis limits
        %xlim([-10e50 10e50])
        %ylim([-10e50 10e50])
        
        xlabel('X')
        ylabel('Y')
        title(['t = ' num2str((i-1)*dt*n)]);
        
        pause(0.0001);
        frame = getframe(gcf); %get frame
        writeVideo(video, frame);

        hold off
    end
    close(video)
end
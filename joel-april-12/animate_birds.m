function animate_birds(x_h, l_h, dt)
    n = 100;
    figure
    fig = gcf;

    video = VideoWriter('flock_video', 'MPEG-4');
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
        
        plot3(x_l(:,1,i*n),x_l(:,2,i*n),x_l(:,3,i*n),".",'MarkerSize',24, 'Color', 'red')
        hold on
        plot3(x_f(:,1,i*n),x_f(:,2,i*n),x_f(:,3,i*n),".",'MarkerSize',24, 'Color', 'black')
        %plot3(x_h(:,1,i*100),x_h(:,2,i*100),x_h(:,3,i*100),".",'MarkerSize',30);

        xlim([min(min(x_h(:,1,:))) max(max(x_h(:,1,:)))]); % Set x-axis limits
        ylim([min(min(x_h(:,2,:))) max(max(x_h(:,2,:)))]); % Set y-axis limits
        zlim([min(min(x_h(:,3,:))) max(max(x_h(:,3,:)))]); % Set z-axis limits

        % l_b = [min(x_h(:,1,:)) min(x_h(:,2,:)) min(x_h(:,3,:))];
        % u_b = [max(x_h(:,1,:)) max(x_h(:,2,:)) max(x_h(:,3,:))];
        % xlim([min(min(l_b)) max(max(u_b))]); % Set x-axis limits
        % ylim([min(min(l_b)) max(max(u_b))]); % Set y-axis limits
        % zlim([min(min(l_b)) max(max(u_b))]); % Set z-axis limits
        
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title(['t = ' num2str((i-1)*dt*n)]);
        
        pause(0.0001);
        frame = getframe(gcf); %get frame
        writeVideo(video, frame);

        hold off
    end
    close( video)
end
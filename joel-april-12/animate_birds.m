function animate_birds(x_h,dt)
    figure
    fig = gcf;
    
    for i = 1:round(size(x_h,3)/100)
        clf(fig)
        
        plot3(x_h(:,1,i*100),x_h(:,2,i*100),x_h(:,3,i*100),".",'MarkerSize',30);
        hold on

        xlim([min(min(x_h(:,1,:))) max(max(x_h(:,1,:)))]); % Set x-axis limits
        ylim([min(min(x_h(:,2,:))) max(max(x_h(:,2,:)))]); % Set y-axis limits
        zlim([min(min(x_h(:,3,:))) max(max(x_h(:,3,:)))]); % Set z-axis limits
        
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title(['t = ' num2str((i-1)*dt*100)]);
        
        pause(0.0001);

        hold off
    end
end
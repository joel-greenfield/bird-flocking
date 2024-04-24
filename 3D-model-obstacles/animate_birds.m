function animate_birds(x_h,dt,filename)
    vid = VideoWriter(filename,'MPEG-4');
    open(vid)

    figure
    fig = gcf;
    
    circle1 = circle([0 1000],200);
    circle2 = circle([1000 0],200);
    circle3 = circle([0 -1000],200);
    circle4 = circle([-1000 0],200);    

    for i = 1:round(size(x_h,3)/100)
        clf(fig)
        
        plot3(x_h(:,1,i*100),x_h(:,2,i*100),x_h(:,3,i*100),'.','MarkerSize',30);
        view(2)
        hold on

        plot(circle1(:,1),circle1(:,2),'r.','MarkerSize',10)
        plot(circle2(:,1),circle2(:,2),'r.','MarkerSize',10)
        plot(circle3(:,1),circle3(:,2),'r.','MarkerSize',10)
        plot(circle4(:,1),circle4(:,2),'r.','MarkerSize',10)

        xlim([min(-1200,min(min(x_h(:,1,:)))) max(1200,max(max(x_h(:,1,:))))]); % Set x-axis limits
        ylim([min(-1200,min(min(x_h(:,2,:)))) max(1200,max(max(x_h(:,2,:))))]); % Set y-axis limits
        zlim([min(min(x_h(:,3,:))) max(max(x_h(:,3,:)))]); % Set z-axis limits
        
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title(['t = ' num2str((i-1)*dt*100)]);
        
        pause(0.1);

        hold off

        writeVideo(vid,getframe(gcf))
    end

    close(vid)
end

function circle_pts = circle(center,radius)
    theta = 0:0.01:2*pi;
    circle_pts = radius*[cos(theta);sin(theta)]'+center;
end
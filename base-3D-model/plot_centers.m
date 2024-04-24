%% Plot trajectories of the center of the flock
function plot_centers(x_h)
    x_avg = squeeze((mean(x_h)));

    r = 100;
    base = [300 300 0];
    [X,Y,Z] = cylinder(r);
    h = 2000;
    Z = Z*h;
    X=X+base(1); Y=Y+base(2); Z=Z+base(3);


    figure
    
    scatter3(x_avg(1,1),x_avg(2,1),x_avg(3,1),'filled','MarkerFaceColor','#77AC30','HandleVisibility','off')
    hold on
    scatter3(x_avg(1,1),x_avg(2,1),0,'filled','MarkerFaceColor','#77AC30','HandleVisibility','off')
    scatter3(x_avg(1,1),0,x_avg(3,1),'filled','MarkerFaceColor','#77AC30','HandleVisibility','off')
    scatter3(0,x_avg(2,1),x_avg(3,1),'filled','MarkerFaceColor','#77AC30','HandleVisibility','off')
    
    scatter3(x_avg(1,end),x_avg(2,end),x_avg(3,end-1),'filled','MarkerFaceColor','#A2142F','HandleVisibility','off')
    scatter3(x_avg(1,end),x_avg(2,end),0,'filled','MarkerFaceColor','#A2142F','HandleVisibility','off')
    scatter3(x_avg(1,end),0,x_avg(3,end),'filled','MarkerFaceColor','#A2142F','HandleVisibility','off')
    scatter3(0,x_avg(2,end),x_avg(3,end),'filled','MarkerFaceColor','#A2142F','HandleVisibility','off')

    plot3(x_avg(1,:),x_avg(2,:),x_avg(3,:),'Color','#0072BD','LineWidth',3,'DisplayName','xyz')
    plot3(x_avg(1,:),x_avg(2,:),zeros(size(x_avg(3,:))),'Color','#D95319','LineWidth',1,'DisplayName','xy plane')
    plot3(x_avg(1,:),zeros(size(x_avg(2,:))),x_avg(3,:),'Color','#7E2F8E','LineWidth',1,'DisplayName','xz plane')
    plot3(zeros(size(x_avg(1,:))),x_avg(2,:),x_avg(3,:),'Color','#4DBEEE','LineWidth',1,'DisplayName','yz plane')

    % s = surf(X,Y,Z,'FaceColor','#80B3FF', 'EdgeColor','none');
    % alpha(s,.2)

    xlabel("x")
    ylabel("y")
    zlabel("z")
    legend
end
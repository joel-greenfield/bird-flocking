function plot_velocity(v_h, dt)
    T = size(v_h, 3);
    tspan = linspace(1, T*dt, T);
    n = size(v_h, 1);
    v_mag = zeros(n, T);
    v_mag(:,:) = sqrt(v_h(:,1,:).^2 + v_h(:,2,:).^2);
    v_avg = mean(v_mag, 1);
    v_std = std(v_mag, 1);

    v_upper = v_avg + v_std;
    v_lower = v_avg - v_std;

    figure
    plot(tspan, v_avg)
    hold on
    plot(tspan, v_lower, "LineStyle", "--", "Color", "#0072BD")
    plot(tspan, v_upper, "LineStyle", "--", "Color", "#0072BD")
    xlim([0 T*dt])
    ylabel("Average velocity")
    xlabel("time")
    legend("$\mu$", "$\mu \pm \sigma$", 'interpreter', 'latex')
    hold off

end
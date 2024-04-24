
centers_file = 'myCenters.txt';
M = csvread(centers_file);
ending_idx = find(M(:,1)==0, 1)-1;
x_centers = M(1:ending_idx,1);
y_centers = M(1:ending_idx,2);

spread_file = 'mySpread.txt';
sig = csvread(spread_file);
sig = sig(1:ending_idx);

figure
hold on

c = linspace(1,10,length(x_centers));
scatter(x_centers, y_centers,sig,c, 'filled')
colorbar
colormap cool

%yline(120, '--' );
%xline(200, '--');

%xlim([min(x_centers)-5,max(x_centers)+5])
%ylim([min(y_centers)-5,max(y_centers)+5])

xlabel('x')
ylabel('y')
title('Flock Center')

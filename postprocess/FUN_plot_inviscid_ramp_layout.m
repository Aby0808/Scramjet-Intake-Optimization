function FUN_plot_inviscid_ramp_layout(x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7)
% Plot the 2D inviscid ramp and shock layout for post-processing.

figure
plot([x1,x2,x3,x4,x6],[y1,y2,y3,y4,y6],Color='k')
hold on
plot([0,x5,x7],[0,y5,y7], Color='k')
plot([x1,0,x2],[y1,0,y2],Color='b')
plot([0,x3],[0,y3],Color='b')
plot([x4,x5],[y4,y5],Color='b')
plot([x6,x7],[y6,y7],Color='b')
xlabel('x (m)')
ylabel('y (m)')
grid on
axis equal
end

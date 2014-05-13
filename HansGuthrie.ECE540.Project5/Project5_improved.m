%% Defined some slick custom colors
orange = [0.8, 0.8, 0];
red = [0.8, 0, 0.2];
blue = [0.2, 0, 0.6];
%% Kp: 1, Ki: 0,2,3

close all

%defined files as variables to make copying for other plots easier
file1 = load( 'FS_output_Kp_1_Ki_0.csv' );
file2 = load( 'FS_output_Kp_1_Ki_2.csv' );
file3 = load( 'FS_output_Kp_1_Ki_3.csv' );

figure(1);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5);
hold on
title( 'Kp = 1, Varying Ki (0, 2, 3)', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Ki = 0', 'Ki = 2', 'Ki = 3' );
print -dpng Kp_1.png

%% Kp: 2, Ki: 0,4,6
close all

file1 = load( 'FS_output_Kp_2_Ki_0.csv' );
file2 = load( 'FS_output_Kp_2_Ki_4.csv' );
file3 = load( 'FS_output_Kp_2_Ki_6.csv' );

figure(2);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on
title( 'Kp = 2, Varying Ki (0, 4, 6)', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Ki = 0', 'Ki = 4', 'Ki = 6' );
print -dpng Kp_2.png

%% Kp: 4, Ki: 0,8,12
close all

file1 = load( 'FS_output_Kp_4_Ki_0.csv' );
file2 = load( 'FS_output_Kp_4_Ki_8.csv' );
file3 = load( 'FS_output_Kp_4_Ki_12.csv' );

figure(3);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on
title( 'Kp = 4, Varying Ki (0, 8, 12)', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Ki = 0', 'Ki = 8', 'Ki = 12' );
print -dpng Kp_4.png

%% Kp: 1,2,4 Ki: 0*Kp
close all

file1 = load( 'FS_output_Kp_1_Ki_0.csv' );
file2 = load( 'FS_output_Kp_2_Ki_0.csv' );
file3 = load( 'FS_output_Kp_4_Ki_0.csv' );

figure(4);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on
title( 'Varying Kp (1, 2, 4), Ki = 0*Kp', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Kp = 1', 'Kp = 2', 'Kp = 4' );
print -dpng Ki_0xKp.png
%% Kp: 1,2,4, Ki: 2*Kp
close all

file1 = load( 'FS_output_Kp_1_Ki_2.csv' );
file2 = load( 'FS_output_Kp_2_Ki_4.csv' );
file3 = load( 'FS_output_Kp_4_Ki_8.csv' );

figure(5);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on
title( 'Varying Kp (1, 2, 4), Ki = 2*Kp', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Kp = 1', 'Kp = 2', 'Kp = 4' );
print -dpng Ki_2xKp.png

%% Kp: 1,2,4, Ki: 3*Kp
close all

file1 = load( 'FS_output_Kp_1_Ki_3.csv' );
file2 = load( 'FS_output_Kp_2_Ki_6.csv' );
file3 = load( 'FS_output_Kp_4_Ki_12.csv' );

figure(6);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on
title( 'Varying Kp (1, 2, 4), Ki = 3*Kp', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Kp = 1', 'Kp = 2', 'Kp = 4' );
print -dpng Ki_3xKp.png

%% Realistic Sawtooth Wave, Kp: 1,2,4, Ki: 2*kp, Period = 4
close all

file1 = load( 'FS_RST_output_Kp_1_Ki_2_T_4.csv' );
file2 = load( 'FS_RST_output_Kp_2_Ki_4_T_4.csv' );
file3 = load( 'FS_RST_output_Kp_4_Ki_8_T_4.csv' );

figure(7);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on

title( 'Sawtooth Wave: Varying Kp (1, 2, 4), Ki = 2*Kp', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Kp = 1', 'Kp = 2', 'Kp = 4' );
print -dpng Sawtooth_T_4.png

%% Realistic Sawtooth Wave, Kp: 1,2,4, Ki: 2*kp, Period = 8
close all

file1 = load( 'FS_RST_output_Kp_1_Ki_2_T_8.csv' );
file2 = load( 'FS_RST_output_Kp_2_Ki_4_T_8.csv' );
file3 = load( 'FS_RST_output_Kp_4_Ki_8_T_8.csv' );

figure(8);

plot( file1(:,3), file1(:,1), '--','Color', orange, 'linewidth',1.5 );
hold on
plot( file2(:,3), file2(:,1), '-','Color', red, 'linewidth',1.5 );
hold on
plot( file3(:,3), file3(:,1), ':','Color', blue, 'linewidth',1.5 );
hold on

title( 'Ki = 0, Varying Kp', 'fontweight','bold' );
xlabel( 'Time (Seconds)', 'fontweight','bold'  );
ylabel( 'Output (Radians)', 'fontweight','bold'  );
grid
legend( 'Kp = 1', 'Kp = 2', 'Kp = 4' );
print -dpng Sawtooth_T_8.png

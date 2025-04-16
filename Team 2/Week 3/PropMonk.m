clear; clc; close all;

meow = readmatrix("1000 RPM.xlsx");
% 1000 RPM
% J1000 = [0.0000, 0.0217, 0.0435, 0.0652, 0.0869, 0.1086, 0.1304, 0.1521, 0.1738, 0.1955, 0.2173, 0.2390, 0.2607, 0.2825, 0.3042, 0.3259, 0.3476, 0.3694, 0.3911, 0.4128, 0.4345, 0.4563, 0.4780, 0.4997, 0.5214, 0.5432, 0.5649, 0.5866, 0.6084, 0.6301];
J1000 = meow(:,2);
Pe1000 = meow(:,3);
Ct1000 = meow(:,4);
Cp1000 = meow(:,5);

bruh = readmatrix("13000 RPM.xlsx");
J13000 = bruh(:,2);
Pe13000 = bruh(:,3);
Ct13000 = bruh(:,4);
Cp13000 = bruh(:,5);

K10 = readmatrix("10000 RPM.xlsx");
J10000 = K10(:,2);
Pe10000 = K10(:,3);
Ct10000 = K10(:,4);
Cp10000 = K10(:,5);



figure(1)
hold on
grid on
subplot(1,3,1)
plot(J1000, Pe1000,'-r','LineWidth',1.5)
hold on
plot(J13000, Pe13000,'-b','LineWidth',1.5)
hold on
plot(J10000, Pe10000,'-k','LineWidth',1.5)
legend('1000 RPM','13000 RPM','10000 RPM', 'Location','se')
title('Pe')
xlabel('J')
grid on
box on
set(gca, 'LineWidth', 1)
subplot(1,3,2)
plot(J1000, Ct1000,'-r','LineWidth',1.5)
hold on
plot(J13000, Ct13000,'-b','LineWidth',1.5)
hold on
plot(J10000, Ct10000,'-k','LineWidth',1.5)
legend('1000 RPM','13000 RPM','10000 RPM', 'Location','se')
title('Ct')
xlabel('J')
grid on
box on
set(gca, 'LineWidth', 1)
subplot(1,3,3)
plot(J1000, Cp1000,'-r','LineWidth',1.5)
hold on
plot(J13000, Cp13000,'-b','LineWidth',1.5)
hold on
plot(J10000, Cp10000,'-k','LineWidth',1.5)
legend('1000 RPM','13000 RPM','10000 RPM', 'Location','se')
title('Cp')
xlabel('J')
grid on
box on
set(gca, 'LineWidth', 1)

bruh = readmatrix("13000 RPM.xlsx");
J13000 = bruh(:,2);
Pe13000 = bruh(:,3);
Ct13000 = bruh(:,4);
Cp13000 = bruh(:,5);

% figure(2)
% hold on
% grid on
% subplot(1,3,1)
% plot(J13000, Pe13000,'-r','LineWidth',1.5)
% title('Pe')
% xlabel('J')
% grid on
% subplot(1,3,2)
% plot(J13000, Ct13000,'-b','LineWidth',1.5)
% title('Ct')
% xlabel('J')
% grid on
% subplot(1,3,3)
% plot(J13000, Cp13000,'-g','LineWidth',1.5)
% title('Cp')
% xlabel('J')
% grid on
% 
% 

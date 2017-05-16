% atmo_test
close all
clear
altitude = linspace (0,200e3,1e6);
press= 0*altitude;
tic
for i=1:length(altitude)
    [press(i), dens(i), sspeed(i), temp(i)] = atmo_ext_ISA(altitude(i));
end
toc
figure
subplot(2,2,1)
semilogx(press./1e3, altitude./1e3)
grid on
hold on
xlabel('Press [kPa]')
ylabel('Altitude [km]');

subplot(2,2,2)
semilogx(dens, altitude./1e3)
grid on
hold on
xlabel('Dens [kg/m^3]')
ylabel('Altitude [km]');

subplot(2,2,3)
plot(sspeed, altitude./1e3)
xlabel('Mach = 1 [m/s]')
ylabel('Altitude [km]');

subplot(2,2,4)
plot(temp, altitude./1e3)
xlabel('Temp [K]')
ylabel('Altitude [km]');
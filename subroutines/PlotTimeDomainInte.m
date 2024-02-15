function [] = PlotTimeDomainInte(acc,dt,VT,VF,DT,DF)
%PLOTINTE Summary of this function goes here
%   Detailed explanation goes here
subplot(311);
hold off;
plot((1:length(acc))*dt,acc);
xlabel('Time (s)')
ylabel('Acceleration (g)')
grid on

subplot(312);
hold off; 
plot((1:length(acc))*dt,VT,'b');
xlabel('Time (s)')
ylabel('Velocity (cm/s)')
grid on

subplot(313);
hold off; 
plot((1:length(acc))*dt,DT,'b');
xlabel('Time (s)')
ylabel('Displacement (cm)')
grid on
end


function [] = PlotInte(acc,dt,VT,VF,DT,DF)
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
plot((1:length(acc))*dt,VF,'r');
hold on;
plot((1:length(acc))*dt,VT,'b');
xlabel('Time (s)')
ylabel('Velocity (cm/s)')
legend('Frequency Domain Integration','Time Domain Integration','Location','Best')
grid on

subplot(313);
hold off; 
plot((1:length(acc))*dt,DF,'r');
hold on;
plot((1:length(acc))*dt,DT,'b');
xlabel('Time (s)')
ylabel('Displacement (cm)')
legend('Frequency Domain Integration','Time Domain Integration','Location','Best')
grid on
end


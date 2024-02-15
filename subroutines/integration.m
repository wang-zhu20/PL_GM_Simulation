function [VelF,DispF,VelT,DispT] = integration(acc,dt);

VelF = fdinteg(acc,dt)*981;
DispF = fdinteg(VelF,dt);

acc = acc'; % acc is now a row vector
acc2 = [0, acc(1:(length(acc)-1))];
accAvg = (acc+acc2)/2;
VelT = cumsum(accAvg).*dt .*981;  % in unit of cm/s

velAvg= VelT + (acc./3+acc2./6).*dt.*981 ;
DispT = cumsum(velAvg).* dt;

return
 
 
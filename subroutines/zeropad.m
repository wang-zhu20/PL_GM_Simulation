function [acc1,tpad] = zeropad(acc,dt,nroll,HP)
%ZEROPAD Summary of this function goes here
%   Detailed explanation goes here
nroll=2;
tpad = 1.5*nroll/HP;
pad = zeros(round(tpad/dt),1);

acc1 = [acc;pad];

%%%%%%%%% Duruo HUANG
%%%%%%%%% I think we should add pad before event !!!!
%%%acc1 = [pad; acc];

end


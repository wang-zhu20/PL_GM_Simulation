function [acc1] = Removepad(acc,tpad,dt)
%REMOVEPAD Summary of this function goes here
%   Detailed explanation goes here


acc1 = acc(1:length(acc)-round(tpad/dt));

%%%% GW: I do not understand what is the heck of above
%PadLength=round(tpad/dt);
%acc1 = acc(PadLength:length(acc));

end


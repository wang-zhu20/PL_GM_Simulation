function [acc] = LoopBaselinecorrectionUPDATED(acc,displ,dt,option)
 

%baseline correction, written by Gang Wang
np=length(displ);
time = [1:np]'.*dt;

[a,b]=size(acc);
if a<b
    acc=acc'; %%% make sure it is a vector
end
[a,b]=size(displ);
if a<b
    displ=displ'; %%% make sure it is a vector
end

%%%% interpolate the displlacement time history using larger dt

Newdt=1;  %%%% use larger dt to avoid singularity
xx=[Newdt:Newdt:time(np)]';

%%%%% uniform distribution in logspace
%%%tt=logspace(log10(Newdt), log10(time(np)), floor(np*dt/5));
%%%xx=10.^tt;

Y=interp1(time, displ, xx);

if option==1
    X=[xx.^2, xx.^3, xx.^4, xx.^5, xx.^6];
elseif option==2
    X=[xx.^2, xx.^3, xx.^4, xx.^5];
elseif option==999 %%% extend to stragight line to 50% long
    %%disp('Extend a straight line at end')
    %Factor=input('Input extension factor:::::::');
    %if isempty(Factor)
    Factor=2;
    %end
    timeExtended=[1:1:round(Factor*np)]'.*dt;
    xxExtended=[Newdt:Newdt:timeExtended(length(timeExtended))]';
    ExtendLineSlope=displ(np)/(np*dt);
    displExtended=ExtendLineSlope*timeExtended;  %%% straight line
    displExtended([1:np])=displ; %%% replace with original data
    Y=interp1(timeExtended, displExtended, xxExtended);
    Order=6;
    if isempty(Order)
        Order=6;
    end
    if Order==6
        X=[xxExtended.^2, xxExtended.^3, xxExtended.^4, xxExtended.^5, xxExtended.^6];
    else
        X=[xxExtended.^2, xxExtended.^3, xxExtended.^4, xxExtended.^5];
    end
end
    
%X=[xx.^2, xx.^3];

%%%%% get coefficient by least squre regression
B=inv(X'*X)*X'*Y; 

D_fitline=X*B;


if option==1
    DD=B(1)*time.^2+B(2)*time.^3+B(3)*time.^4+B(4)*time.^5+B(5)*time.^6;
elseif option==2
    DD=B(1)*time.^2+B(2)*time.^3+B(3)*time.^4+B(4)*time.^5;
elseif option==999
    if Order==6
        DD=B(1)*time.^2+B(2)*time.^3+B(3)*time.^4+B(4)*time.^5+B(5)*time.^6;
    else
        DD=B(1)*time.^2+B(2)*time.^3+B(3)*time.^4+B(4)*time.^5;
    end
end
 

if option==1 
    Acc_substract=2*B(1)+6*B(2)*time+12*B(3)*time.^2+20*B(4)*time.^3+30*B(5)*time.^4;
elseif option==999
    if Order==6
        Acc_substract=2*B(1)+6*B(2)*time+12*B(3)*time.^2+20*B(4)*time.^3+30*B(5)*time.^4;
    else
        Acc_substract=2*B(1)+6*B(2)*time+12*B(3)*time.^2+20*B(4)*time.^3;    
    end
else
    Acc_substract=2*B(1)+6*B(2)*time+12*B(3)*time.^2+20*B(4)*time.^3;
end

%B
%size(time)
%figure; plot(acc,'r'); hold on
acc=acc-Acc_substract/981; %%% convert back to g
%plot(acc,'b'); hold on




    

 
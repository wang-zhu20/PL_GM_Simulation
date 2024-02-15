function [acc] = baselinecorrect(acc,disp,dt,q)
%baseline correction
np=length(disp);
xd = [1:np].*dt;


if q ==1
poly =inline('c(1)*xd.^6+c(2)*xd.^5+c(3)*xd.^4+c(4)*xd.^3+c(5)*xd.^2','c','xd');
c0 = [1,1,1,1,1];
c = nlinfit(xd,disp,poly,c0);

for i=1:np;
    acc(i) = acc(i)-(30*c(1)*(i*dt)^4+20*c(2)*(i*dt)^3 ...
        +12*c(3)*(i*dt)^2+6*c(4)*i*dt+2*c(5))/981;
end

else

[c2,s] = polyfit(xd,disp,2);
[v,delta2] = polyval(c2,1:np,s);
[c3,s] = polyfit(xd,disp,3);
[v,delta3] = polyval(c3,1:np,s);
[c4,s] = polyfit(xd,disp,4);
[v,delta4] = polyval(c4,1:np,s);
[c5,s] = polyfit(xd,disp,5);
[v,delta5] = polyval(c5,1:np,s);
[c6,s] = polyfit(xd,disp,6);
[v,delta6] = polyval(c6,1:np,s);
        
m= min([delta2,delta3,delta4,delta5,delta6]);

if m==delta2
    for i=1:np
    acc(i) = acc(i)-(2*c2(1))/981;
    end
    
elseif m==delta3
    for i=1:np
    acc(i) = acc(i)-(6*c3(1)*i*dt+2*c3(3))/981;
    end
    
elseif m==delta4
    for i=1:np
    acc(i) = acc(i)-(12*c4(1)*(i*dt)^2+6*c4(2)*i*dt+2*c4(3))/981;
    end
    
elseif m==delta5
    for i=1:np
    acc(i) = acc(i)-(20*c5(1)*(i*dt)^3 ...
        +12*c5(2)*(i*dt)^2+6*c5(3)*i*dt+2*c5(4))/981;
    end
    
else
for i=1:np
    acc(i) = acc(i)-(30*c6(1)*(i*dt)^4+20*c6(2)*(i*dt)^3 ...
        +12*c6(3)*(i*dt)^2+6*c6(4)*i*dt+2*c6(5))/981;
end
end
end
end
            

    

 
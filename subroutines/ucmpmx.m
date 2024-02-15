function [z]=ucmpmx(kug,ug,time,pr,w,d)
% inputs:
%    kug -- number of time increment
%    ug -- input ground acceleration
%    time -- time sequence
%    pr -- period at which spectra are calculated
%    w -- frequency
%    d -- damping ratio,
% outputs:
%    z - z(1);   maximum relative displacement
%        z(2);   maximum relative velocity
%        z(3);   maximum absolute acceleration
%    x is time history, not outputed yet

wd=sqrt(1.-d*d)*w;
w2=w*w;
w3=w2*w;
for i=1:3
    x(1,i)=0.0;
    z(i)=0.0;
end
f2=1./w2;
f3=d*w;
f4=1./wd;
f5=f3*f4;
f6=2.*f3;
for k=1:kug
    dt=time(k+1)-time(k);
    ns=round(10.*dt/pr-0.01);
    dt=dt/real(ns);  % reduce time step for STIFF system
    f1=2.*d/w3/dt;
    e=exp(-f3*dt);
    g1=e*sin(wd*dt);
    g2=e*cos(wd*dt);
    h1=wd*g2-f3*g1;
    h2=wd*g1+f3*g2;
    dug=(ug(k+1)-ug(k))/real(ns);
    g=ug(k);
    z1=f2*dug;
    z3=f1*dug;
    z4=z1/dt;
    for is=1:ns      % march over reduced substeps
        z2=f2*g;
        b=x(1,1)+z2-z3;
        a=f4*x(1,2)+f5*b+f4*z4;
        x(2,1)=a*g1+b*g2+z3-z2-z1;
        x(2,2)=a*h1-b*h2-z4;
        x(2,3)=-f6*x(2,2)-w2*x(2,1);
        for l=1:3
            c(l)=abs(x(2,l));
            if(c(l)>=z(l))
                z(l)=c(l);
                t(l)=time(k)+is*dt;
            else
            end
            x(1,l)=x(2,l);
        end
        g=g+dug;
    end
end



% compute response of single-DOF system under fixed time step
function [z]=cmpmax(kug,ug,pr,w,d,dt)
% inputs:
%    kug -- number of time increment
%    ug -- input ground acceleration
%    pr -- period at which spectra are calculated
%    w- frequency
%    d -- damping ratio,
%    dt -- time step,
% outputs:
%    z - z(1);   maximum relative displacement
%        z(2);   maximum relative velocity
%        z(3);   maximum absolute acceleration
%    x is time history, not outputed yet

wd=sqrt(1.-d*d)*w;
w2=w*w;
w3=w2*w;
for i=1:3
    x(1,i)=0.0;
    z(i)=0.0;
end
f1=2.*d/(w3*dt);
f2=1./w2;
f3=d*w;
f4=1./wd;
f5=f3*f4;
f6=2.*f3;
e=exp(-f3*dt);
g1=e*sin(wd*dt);
g2=e*cos(wd*dt);
h1=wd*g2-f3*g1;
h2=wd*g1+f3*g2;
for k=1:kug
    dug=ug(k+1)-ug(k);
    z1=f2*dug;
    z2=f2*ug(k);
    z3=f1*dug;
    z4=z1/dt;
    b=x(1,1)+z2-z3;
    a=f4*x(1,2)+f5*b+f4*z4;
    x(2,1)=a*g1+b*g2+z3-z2-z1;   % relative disp.
    x(2,2)=a*h1-b*h2-z4;         % relative vel.
    x(2,3)=-f6*x(2,2)-w2*x(2,1); % absolute acc.
    
    % find the maximum of each
    for l=1:3
        c(l)=abs(x(2,l));
        if(c(l)>=z(l))
            z(l)=c(l);
            t(l)=dt*real(k);
        else
        end
        x(1,l)=x(2,l);
    end
end
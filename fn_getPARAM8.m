function [outprm]=fn_getPARAM8(Mw,Rrup,Rhyp,Vs30,filereg)
% generating parameters based on the regression equations
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu


if (nargin < 5)
 disp('need more input');
 stop;
end

D=Rrup;
% initialize outprm
outprm = fn_initDBs(2);

% generate parameter with corrrelation of mulivariate normal distribution

% Mw: Moment Magnitude
% Dhyp: Hypocentral Distance(km)
% Drup: Closest Distance(km)
% Vs30: Average of shear wave velocity(m/s)
% filereg: name of file that contains regression equations


%     median prediction of each parameter
    ind=1;
    hcoef=10;
    hcoef1=1;
    
%     Coeff = [-1.58	0.73	0.12	-0.21
% 0.28	0.59	-0.03	-0.35
% 0.15	-0.02	0.03	0.26
% -0.93	0.22	0.05	0.20
% -0.73	0.06	-0.04	0.05
% -2.38	0.77	0.15	-0.17
% -0.71	0.59	-0.01	-0.27
% 0.43	-0.14	0.08	0.25
% -1.33	0.17	0.17	0.16
% -0.79	0.02	-0.06	0.11
% -3.63	0.46	-0.89	0.03
% -4.86	-0.01	-0.98	0.01];


  Coeff = [-2.342426488	0.68799677	0.227414135	-0.097600939
-1.657789454	0.596876915	0.286321993	-0.148944489
2.205870259	0.035744168	-0.098301376	0.018698363
1.277279461	0.149512253	0.014687628	0.050559948
-0.544966685	0.043923012	-0.01824002	0.021124879
-3.055804437	0.747193353	0.267777643	-0.09716879
-2.269888287	0.586818353	0.363700077	-0.143378833
2.986360677	-0.004213658	-0.109479299	-0.106823394
2.262305805	0.097405858	-0.004620375	-0.122844183
-0.533706594	0.034498168	-0.025133858	0.025108618
-6.904231234	1.251257109	-1.372244069	-0.323508914
-9.847664715	0.700461509	-1.470451386	-0.093534688];




base_x = [1, Mw, log(sqrt(Rrup^2+36)),log(Vs30)];
Pred = Coeff*base_x';
% totalEnergyM = Pred(1);
% majorEaM = Pred(2);
% majorExM= Pred(3);
% majorSxM    = Pred(4);
% minorExM    = Pred(5);
% minorSxM= Pred(6);
% majorEyM= Pred(7);
% majorSyM    = Pred(8);
% minorEyM    = Pred(9);
% minorSyM= Pred(10);
% majorRxyM= Pred(11);
% minorRxyM= Pred(12);
 majorExM  = Pred(1);
 majorEyM  = Pred(2);
    majorSxM  = Pred(3);
    majorSyM  = Pred(4);
    majorRxyM = Pred(5);
    
    minorExM  =Pred(6);
    minorEyM  = Pred(7);
    minorSxM  = Pred(8);
    minorSyM  = Pred(9);
    minorRxyM = Pred(10);
    
    totalEnergyM = Pred(11);
    majorEaM     = Pred(12);
    

    load('covrho_H_V.mat')

        while true
	%               inter-event residuals
    si = mvnrnd(Ei,covi);
    %             intra-event residuals
    sm=mvnrnd(Em,covm);
    %             total residuals
    epst=si+sm; % 24*1 vector
    epst = epst(13:end);

%             random factor of wavelet packets of minor group
                mn13=1.09;
    sn13=0.1926;
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
            
%     prediction of each parameter with residuals
            ind=1;
            majorEx = exp(majorExM + epst(ind)); ind=ind+1;
            majorEy = exp(majorEyM + epst(ind)); ind=ind+1;
            majorSx = exp(majorSxM + epst(ind)); ind=ind+1;
            majorSy = exp(majorSyM + epst(ind)); ind=ind+1;
            majorRxy = 2*normcdf(majorRxyM + epst(ind),0,1)-1; ind=ind+1;

            minorEx = exp(minorExM + epst(ind)); ind=ind+1;
            minorEy = exp(minorEyM + epst(ind)); ind=ind+1;
            minorSx = exp(minorSxM + epst(ind)); ind=ind+1;
            minorSy = exp(minorSyM + epst(ind)); ind=ind+1;
            minorRxy = 2*normcdf(minorRxyM + epst(ind),0,1)-1; ind=ind+1;
            
            totalEnergy = exp(totalEnergyM + epst(ind)); ind=ind+1;
            majorEa = exp(majorEaM + epst(ind)); ind=ind+1;
            
            
            majorVx = majorSx^2;
            majorVy = majorSy^2;
            majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;

            majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
            majorVlx = log(majorVx/majorEx^2 + 1);
            majorSlx = sqrt(majorVlx);
            majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
            majorVly = log(majorVy/majorEy^2 + 1);
            majorSly = sqrt(majorVly);

            majorCovlxly = log(majorExy/majorEx/majorEy);
            majorRlxly = majorCovlxly/majorSlx/majorSly;

            minorVx = minorSx^2;
            minorVy = minorSy^2;
            minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
            minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
            minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
            minorVlx = minorSlx^2;
            minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
            minorSly = sqrt(log(minorVy/minorEy^2 + 1));
            minorVly = minorSly^2;

            minorCovlxly = log(minorExy/minorEx/minorEy);
            minorRlxly = minorCovlxly/minorSlx/minorSly;
            
            majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
            [T,err] = cholcov(majorLLCov);
            if err == 0;
                break;
            end;
        end;

%         output parameters
        outprm(1) = struct('M',Mw,'hdist',D,'vs30',Vs30, ...
                'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
                'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
                'totalEnergy',totalEnergy, 'majorEa',majorEa,'minorRnd',minorRnd ...
            );

%         for median prediction
    majorExM = exp(majorExM);
    majorSxM = exp(majorSxM);
    majorEyM = exp(majorEyM);
    majorSyM = exp(majorSyM);
    majorRxyM = 2*normcdf(majorRxyM,0,1)-1;

    minorExM = exp(minorExM);
    minorSxM = exp(minorSxM);
    minorEyM = exp(minorEyM);
    minorSyM = exp(minorSyM);
    minorRxyM = 2*normcdf(minorRxyM,0,1)-1;

    majorEaM = exp(majorEaM);
    totalEnergyM = exp(totalEnergyM);

%     about major
    majorVxM = majorSxM^2;
    majorVyM = majorSyM^2;
    majorExyM = majorRxyM*majorSxM*majorSyM + majorExM*majorEyM;
    
    majorElxM = log(majorExM^2 / sqrt(majorVxM+majorExM^2));
    majorSlxM = sqrt(log(majorSxM/majorExM^2 + 1));
    majorVlxM = majorSlxM^2;
    majorElyM = log(majorEyM^2 / sqrt(majorVyM+majorEyM^2));
    majorSlyM = sqrt(log(majorSyM/majorEyM^2 + 1));
    majorVlyM = majorSlyM^2;
    
    majorCovlxlyM = log(majorExyM/majorExM/majorEyM);
    majorRlxlyM = majorCovlxlyM/majorSlxM/majorSlyM;
    
%     about minor
    minorVxM = minorSxM^2;
    minorVyM = minorSyM^2;
    minorExyM = minorRxyM*minorSxM*minorSyM + minorExM*minorEyM;

    minorElxM = log(minorExM^2 / sqrt(minorVxM+minorExM^2));
    minorSlxM = sqrt(log(minorSxM/minorExM^2 + 1));
    minorVlxM = minorSlxM^2;
    minorElyM = log(minorEyM^2 / sqrt(minorVyM+minorEyM^2));
    minorSlyM = sqrt(log(minorSyM/minorEyM^2 + 1));
    minorVlyM = minorSlyM^2;
    
    minorCovlxlyM = log(minorExyM/minorExM/minorEyM);
    minorRlxlyM = minorCovlxlyM/minorSlxM/minorSlyM;
    
        outprm(2) = struct('M',Mw,'hdist',D,'vs30',Vs30, ...
                'minorElx',minorElxM, 'minorSlx',minorSlxM, 'minorVlx',minorVlxM, 'minorEly',minorElyM, 'minorSly',minorSlyM, 'minorVly',minorVlyM,'minorRlxly',minorRlxlyM, ...
                'majorElx',majorElxM, 'majorSlx',majorSlxM, 'majorVlx',majorVlxM, 'majorEly',majorElyM, 'majorSly',majorSlyM, 'majorVly',majorVlyM,'majorRlxly',majorRlxlyM, ...
                'totalEnergy',totalEnergyM, 'majorEa',majorEaM, 'minorRnd',minorRnd ...
            );

    

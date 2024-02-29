% Program for generating artificial ground motions by stochastic ground
% motion model using wavelet packets
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% load('E:\PulseResearch_NGA2\PredictPara.mat')
% PredictPara = PredictPara(10,:);
% M = PredictPara(:,4);R = PredictPara(:,6);Vs30 = PredictPara(:,7);
% SorD = PredictPara(:,8);


clear;
close all
clc;

addpath('./subroutines');
% Rrup:     Rupture distance (km)
% inVs30:   Average shear wave velocity with 30m in surface
% Mw:       Minimum Moment Magnitude
% SorD:     the length of the portion of the rupture between the hypocenter and the site
% nsmpl:    # of samples to generate
% iflg:     =1  1 realization                %%% with variability G.W.
%           =2  wave from median parameters  %%% no variability   G.W.




iflg = 1;
offset=0;

% % conditions for each case

ind=1;


Mw(ind)=7.0;	Rrup(ind)=10;	Rhyp(ind)= 10;	Vs30(ind)= 760; 
SorD(ind) = 53;    nsmpl(ind)= 10;  	ind=ind+1;

% 
% Mw(ind)=7.5;	Rrup(ind)=10;	Rhyp(ind)= 10;	Vs30(ind)= 760; 
% SorD(ind) = 108;    nsmpl(ind)= 100;  	ind=ind+1;

ncase=ind-1;

% initializing random number generater
rs.TotalSign = RandStream('mt19937ar');
rs.MinorRand = RandStream('mt19937ar');
rs.MajorLoc  = RandStream('mt19937ar');
rs.MajorAmp  = RandStream('mt19937ar');
inrand=ceil(rand(1)*100);
for i=1:1:inrand
    rand(rs.TotalSign);
    rand(rs.MinorRand);
    rand(rs.MajorLoc);
    rand(rs.MajorAmp);
end

% loop for each case
for j=1:1:ncase
    % making directory for output files
    filename=sprintf('M%03.1f_Rr%08.4f_Vs%06.1f',Mw(j),Rrup(j),Vs30(j));
    system(sprintf('mkdir %s',filename));

    % loop for each sample
    for i=1:1:nsmpl(j)
        disp(['computing... smpl' num2str(i,'% 6d') '/' num2str(nsmpl(j),'% 6d') '| case' num2str(j,'% 6d') '/' num2str(ncase,'% 6d')]);

        [prmcoef,plcoef]=fn_PredictPara_PL(Mw(j),Rrup(j),Vs30(j),SorD(j),iflg);
        
        [th dt] = fn_get1Sim_PL(prmcoef,rs);
        
        [Pulse_record, dt,Vp,Tp,Et,Eacc] = sim_single(plcoef);
        
        acc_raw = extract_major(th);

        acc_raw = acc_raw'; 
        High = 0.1; Low = 30;nrl = 2;
        [acc2,dt,HPused,LPused,nrolused,baselineOption] = Process(acc_raw,dt,High,Low,nrl);
         acc2 = acc2'; 
         
         
        acc_temp = extract_major(acc2);
        acc = [zeros(1,756),acc_temp]; % 7.56 sec
        vel = cumsum(acc) .* dt .* 981;
        
        
        a = 0.05;
        Pulse_temp = taper(Pulse_record',a);
        Pulse_record = Pulse_temp';
        % same length
%         [acc_total,vel_total] = add_pulse_res(vel,Pulse_record);
        sim1 = Pulse_record;dt = 0.01;
        
                
        if length(sim1)<length(vel)
            sim1=[sim1,zeros(1,length(vel)-length(sim1))];
        else
            sim1 = sim1(1:length(vel));
        end
        vel_total = vel + sim1;
        acc_total = [0,diff(vel_total)./981./dt];
        
        
        
figure(999)
t = tiledlayout(3,1) ;t.Padding = 'compact';t.TileSpacing = 'compact';

nexttile(1);
plot((1:length(vel))*dt,sim1,'k','Linewidth',0.5);hold on

set(gca,'xtick',[]);

nexttile(2);
plot((1:length(vel))*dt,vel,'k','Linewidth',0.5);hold on
set(gca,'xtick',[]);

nexttile(3);
plot((1:length(vel))*dt,vel_total,'k','Linewidth',0.5);hold on
ly = ylim;
xlabel( 'Time (s)');
ylabel(t,'Vel (cm/s)');
for iii = 1:3
    nexttile(iii);set(gca,'fontsize',12);xlim([0 length(vel_total)*dt]);
    %     xlim([0 50])
    if iii ~= 2
        ylim(ly);
    end
end
set(gcf,'Units','centimeters','Position',[12 8 12 9]); % 图片大小

% exportgraphics(gcf,[filename,'/Figure/figure',num2str(i),'.jpg'],'Resolution',300)

Simdata.acc = acc;
Simdata.vel = vel;
Simdata.pulse = sim1;
Simdata.vel_total = vel_total;
Simdata.acc_total = acc_total;


Simdata.dt = 0.01;
Simdata.pulsepara = [Vp,Tp,Et,Eacc];
DataName = [filename,'/Simdata',num2str(i),'.mat'];
save(DataName,'Simdata');

    
    end
    

end

rmpath('./subroutines');

function [th] = extract_major(acc)
Ecum=cumsum(acc.^2);
total=Ecum(length(Ecum));


th = acc((Ecum>=0.0001*total)&(Ecum<=0.9999*total));

end


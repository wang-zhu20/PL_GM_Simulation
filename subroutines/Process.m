function [acc2,dt,HP,LP,nroll,option] = Process(acc,dt,H,L,n)
%DR Huang PROCESS Summary of this function goes here %%% August 15, 2016
%
%  input:
%     acc: unprocessed acceleration time history
%     H: high-pass filter; L: low-pass filter; n: nroll; all from NGA
%     flatfiles

%     [VelF,DispF,VelT,DispT] = integration(acc,dt);
%     figure(1)
%     PlotInte(acc,dt,VelT,VelF,DispT,DispF);
%     title('Time Histories of Unprocessed Data')

%     figure(2)
%     hold off;
%     [f,h] = FourierSpc(acc,dt); %plot Fourier Spectrum to review
%     PlotH0=loglog(f,h);
%     grid on;
%     xlabel('Frequency (Hz)')
%     ylabel('FAS of Acceleration')

%%%disp(sprintf('PEER High pass filter freq=%f',H));
% HP = input('HP= (enter to use default): ');
HP=[];
if isempty(HP)
    HP=H;
    %%%disp(sprintf('Use default HP=%f',HP))
end

%%%disp(sprintf('PEER Low pass filter freq=%f',L));

%LP = input('LP= (enter to use default): ');
LP=[];
if isempty(LP)
    LP=L;
    %%%disp(sprintf('Use default HP=%f',LP))
end
nroll = n;

RedoFlag=1;

acc_save=acc;

while RedoFlag==1
    RedoFlag = 0;
    %%%%%%%%%% set initial range to ZERO!!!
    acc=acc_save;
    % SetZeroSecondRange = input('Set zero acc. range (sec):');
    SetZeroSecondRange=[];
    if isempty(SetZeroSecondRange)
    else
        acc(1:SetZeroSecondRange/dt)=[];
        %%%disp('######## Set range of acc to zeros!!!!! ')
    end
    
    acc=acc-mean(acc);
    
    a = 0.05;
    acc2 = taper(acc,a);
    [acc2,tpad] = zeropad(acc2,dt,nroll,HP);
    acc2 = acausal(LP,HP,nroll,acc2,dt); %filter
    
    %         figure(8888);
    %         subplot(3,1,1)
    %         hold off; plot([dt:dt:length(acc2)*dt],acc2);
    %         xlabel('Time')
    %         ylabel('Acc (g)')
    %         title('Tapered, Padded Acc afer Acausal Filter')
    
    acc2 = Removepad(acc2,tpad,dt);
    [VelF,dispF,VelT,DispT] = integration(acc2,dt);
    
    %         figure(2)
    %         [f,h] = FourierSpc(acc2,dt); %
    %         hold on; PlotH1=loglog(f,h,'g');
    %         hold on; PlotH3=plot([LP,LP],ylim,'b');
    %         hold on; PlotH4=plot([HP,HP],ylim,'b');
    
    %         figure(3)
    %         hold off;
    %         PlotTimeDomainInte(acc2,dt,VelT,VelF,DispT,DispF);
    %         figure(3); subplot(3,1,2); hold on
    %         plot([1,length(acc2)].*dt, [avgPGV,avgPGV],'r'); hold on
    %         plot([1,length(acc2)].*dt, [-avgPGV,-avgPGV],'r'); hold on
    %         ylim(1.1*[-max(max(abs(VelT)),avgPGV),max(max(abs(VelT)), avgPGV)])
    %         figure(3); subplot(3,1,3); hold on
    %         plot([1,length(acc2)].*dt, [avgPGD,avgPGD],'r'); hold on
    %         plot([1,length(acc2)].*dt, [-avgPGD,-avgPGD],'r'); hold on
    %         ylim(1.1*[-max(max(abs(DispT)),avgPGD),max(max(abs(DispT)), avgPGD)])
    %  pause;
    
    option=1;
    % keyboard = input('Need Baseline Correction? 0/[1]: ');
    keyboard=1;
    
    if isempty(keyboard)|keyboard==1
        %option=input('Input baseline correction option (1/2/999): ');
        option=[];
        if isempty(option)
            option=1;
            disp(sprintf('BaselineOption=%f',1))
        else
            %%%disp(sprintf('BaselineOption=%f',option))
        end
        %%%%% changed by Gang Wang
        option=1;
        acc2 = baselinecorrectionUPDATED(acc2,DispT,dt,option); % baseline
    end
    
    
    if isempty(SetZeroSecondRange)
    else
        acc2=[zeros(SetZeroSecondRange/dt,1); acc2];
    end
    
    
    
    [VelF,DispF,VelT,DispT] = integration(acc2,dt); %integration
    
    
    
end
%%% end while do



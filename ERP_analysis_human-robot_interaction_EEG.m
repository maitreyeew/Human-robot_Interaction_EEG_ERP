%% ERP analysis for IoT paper 
% Author - Maitreyee Wairagkar 2021

eeg = [];
emo = [];
cnt = 1;
    
% Load EEG data, segment it into individual trials, collect all time-locked trials
for p=1:10
    load(strcat('P',num2str(p)));
    P=eval(strcat('P',num2str(p)));

    Fs = P.Fs;
    
    % Segmented EEG into trials -200ms to +800ms from the stimulus
    segment_start = round((Fs/10)*2);   % 200 ms
    segment_end = round((Fs/10)*8);     % 400 ms
    
    for r = 1:length(P.EEG) %collect all trials from all runs
        for tr = 1:length(P.EEG(r).stimuli.index)
            ind = P.EEG(r).stimuli.index(tr);
            emo(cnt) = P.EEG(r).stimuli.emotion(tr);

            eeg(cnt,:,:)= P.EEG(r).clean(:,ind-segment_start:ind+segment_end);
            cnt = cnt+1;
        end %tr
    end %r
end %p

time=(-segment_start:segment_end)*(1000/Fs);

%% Pre-process EEG trials

% 1. laplacian filter for channels Cz and Fz
%'Fp1','Fpz','Fp2', 'F3','Fz', 'F4', 'C3', 'Cz', 'C4','T7','T8', 'P3','Pz','P4','POz','Oz'
%  1     2     3      4    5     6     7     8     9   10   11    12   13   14   15   16
% Cz = Cz - 1/4(C3 + C4 + Fz + Pz) i.e. channels 8 = 8 -1/4(7+9+5+13);
% Fz = Fz - 1/4(F3 + F4 + Fpz + Cz) i.e. channels 5 = 5 - 1/4(4+6+2+8);
%%{
for tr = 1:size(eeg,1)
    eeg(tr,8,:) = eeg(tr,8,:) - ((eeg(tr,7,:)+eeg(tr,9,:)+eeg(tr,5,:)+eeg(tr,13,:))/4);%laplacian filter on Cz
    eeg(tr,5,:) = eeg(tr,5,:) - ((eeg(tr,4,:)+eeg(tr,6,:)+eeg(tr,2,:)+eeg(tr,8,:))/4);%laplacian filter on Fz
end
%}

EEG = eeg;
[~,zero_s] = find(time==0);      % Index of stimulus onset
[~,ms170]=min(abs(time - 170));  % Index of 170ms for N170 or VPP


% 2. lowpass and highpass filter each trial
[b_l,a_l]=butter(4, 5/(Fs*0.5),'low'); 
[b_h,a_h]=butter(4, 1/(Fs*0.5),'high'); 
for tr = 1:size(EEG,1)
    for ch = 1:size(EEG,2)
        EEG(tr,ch,:) = filtfilt(b_l,a_l,squeeze(double(EEG(tr,ch,:))));%low pass filter
        EEG(tr,ch,:) = filtfilt(b_h,a_h,squeeze(double(EEG(tr,ch,:))));%high pass filter
    end
end

% 3. remove baseline from each trial ie remove the mean of prestimulus 200ms from all samples
for tr = 1:size(EEG,1)
    for ch = 1:size(EEG,2)
        baseline = mean(squeeze(EEG(tr,ch,1:zero_s)));%mean of prestimulus onset 
        EEG(tr,ch,:) = EEG(tr,ch,:)-baseline;
        
        %EEG(tr,ch,:) = movmean(EEG(tr,ch,:),round(Fs/10)); % moving average filter on EEG trials
    end
end

%% For average ERP of all expressions together 
% 5. Compute ERP
ERP = squeeze(mean(EEG,1));

figure,
% 6. Plot ERP 
chs ={'Fp1','Fpz','Fp2', 'F3','Fz', 'F4', 'C3', 'Cz', 'C4','T7','T8', 'P3','Pz','P4','POz','Oz'};

for i = 1:16
   subplot(4,4,i), plot(time, ERP(i,:)','Linewidth',1.5);hold on; 
   axis tight;
   plot(zeros(1:2),ylim,'k','Linewidth',1);
   plot(xlim,zeros(1:2),'k','Linewidth',1);
   plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1); %170 ms is 49th value in time
   title(chs(i));
end

%% ERP of individual expressions 

% 5. separate processed trials of each expression (1-Angry, 2-Happy, 3-Sad, 4-Surprised)
cnt_ang = 1; cnt_hap = 1; cnt_sad = 1; cnt_sur = 1;
for tr = 1:size(EEG,1)
    % Angry
    if emo(tr) == 1 
        EEG_ang(cnt_ang,:,:) = EEG(tr,:,:);
        cnt_ang = cnt_ang+1;
    end
    % Happy 
    if emo(tr) == 2 
        EEG_hap(cnt_hap,:,:) = EEG(tr,:,:);
        cnt_hap = cnt_hap+1;
    end
    % Sad 
    if emo(tr) == 3 
        EEG_sad(cnt_sad,:,:) = EEG(tr,:,:);
        cnt_sad = cnt_sad+1;
    end
    % Surprised
    if emo(tr) == 4 
        EEG_sur(cnt_sur,:,:) = EEG(tr,:,:);
        cnt_sur = cnt_sur+1;
    end
        
end

% 6. Compute ERP for each expression
ERP_ang = squeeze(mean(EEG_ang,1));
ERP_hap = squeeze(mean(EEG_hap,1));
ERP_sad = squeeze(mean(EEG_sad,1));
ERP_sur = squeeze(mean(EEG_sur,1));

% 7. Plot ERPs
chs ={'Fp1','Fpz','Fp2', 'F3','Fz', 'F4', 'C3', 'Cz', 'C4','T7','T8', 'P3','Pz','P4','POz','Oz'};

for i = 1:16
   subplot(4,4,i), 
   p1 = plot(time, ERP_ang(i,:)','Linewidth',1.5,'color',[0.8 0 0]);hold on; 
   p2 = plot(time, ERP_hap(i,:)','Linewidth',1.5,'color',[0 0.5 0]);hold on; 
   p3 = plot(time, ERP_sad(i,:)','Linewidth',1.5,'color',[0.5 0 0.5]);hold on; 
   p4 = plot(time, ERP_sur(i,:)','Linewidth',1.5,'color',[1 0.7 0]);hold on; 
   axis tight;
   plot(zeros(1:2),ylim,'k','Linewidth',1);
   plot(xlim,zeros(1:2),'k','Linewidth',1);
   plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1); %170 ms is 49th value in time
   title(chs(i));
end

legend([p1 p2 p3 p4],'Angry', 'Happy', 'Sad', 'Surprised','Orientation','horizontal', 'FontSize', 18)

%% additionally filter ERP of expressions
ERP_ang = double(ERP_ang);
ERP_hap = double(ERP_hap);
ERP_sad = double(ERP_sad);
ERP_sur = double(ERP_sur);

[b,a]=butter(4, 4/(Fs*0.5),'low'); 
ch = 8; %cz
ERP_ang(ch,:) = filtfilt(b,a,double(ERP_ang(ch,:)));
ERP_hap(ch,:) = filtfilt(b,a,double(ERP_hap(ch,:)));
ERP_sad(ch,:) = filtfilt(b,a,double(ERP_sad(ch,:)));
ERP_sur(ch,:) = filtfilt(b,a,double(ERP_sur(ch,:)));

%% Figure 7 for Journal %%%%%%%%%
[~,zero_s] = find(time==0);    %Index of stimulus onset
[~,ms170]=min(abs(time - 170));  %Index of 170ms for N170 or VPP
minus100ms = 1;

scale = 1.5; %because laplacian filter reduces amplitude so restore it to original 

chs ={'Pz'};
channel = 13;

subplot(1,2,1), 
p1 = plot(time, ERP(channel,:).*scale','Linewidth',1.5);hold on;
axis tight;
plot(zeros(1:2),ylim,'k','Linewidth',1);
plot(xlim,zeros(1:2),'k','Linewidth',1);
plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1);
title('Cz','FontSize', 16); xlabel({'Time (ms)';'.'}); ylabel('Amplitude (\muV)');
set(gca,'FontSize', 16)


subplot(1,2,2), 
p2 = plot(time, ERP_ang(channel,:).*scale','Linewidth',1.5, 'color',[0.8 0 0]);hold on;
p3 = plot(time, ERP_hap(channel,:).*scale','Linewidth',1.5, 'color',[0 0.5 0]);hold on;
p4 = plot(time, ERP_sad(channel,:).*scale','Linewidth',1.5, 'color',[0.5 0 0.5]);hold on;
p5 = plot(time, ERP_sur(channel,:).*scale','Linewidth',1.5, 'color',[1 0.7 0]);hold on;
axis tight;
plot(zeros(1:2),ylim,'k','Linewidth',1);
plot(xlim,zeros(1:2),'k','Linewidth',1);
plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1);
title('Cz'); 
set(gca,'FontSize', 16)


legend([p1 p2 p3 p4 p5],'All emotions', 'Angry', 'Happy', 'Sad', 'Surprised','Orientation','horizontal', 'FontSize', 18)


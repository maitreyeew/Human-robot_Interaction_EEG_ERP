%% Final ERP analysis for IoT paper - 21/05/21
% Copyright 2021, Maitreyee Wairagkar

% Load EEG segmented data -200ms to +800ms from the stimulus
eeg = EEG_all(:,:,1:1230);
time = ((0:size(eeg,3)-1)/2048).*1000; time = time - 200;

% 1. laplacian filter for channels Cz and Fz
%'Fp1','Fpz','Fp2', 'F3','Fz', 'F4', 'C3', 'Cz', 'C4','T7','T8', 'P3','Pz','P4','POz','Oz'
%  1     2     3      4    5     6     7     8     9   10   11    12   13   14   15   16
% Cz = Cz - 1/4(C3 + C4 + Fz + Pz) i.e. channels 8 = 8 -1/4(7+9+5+13);
% Fz = Fz - 1/4(F3 + F4 + Fpz + Cz) i.e. channels 5 = 5 - 1/4(4+6+2+8);
for tr = 1:size(eeg,1)
    eeg(tr,8,:) = eeg(tr,8,:) - ((eeg(tr,7,:)+eeg(tr,9,:)+eeg(tr,5,:)+eeg(tr,13,:))/4);%laplacian filter on Cz
    eeg(tr,5,:) = eeg(tr,5,:) - ((eeg(tr,4,:)+eeg(tr,6,:)+eeg(tr,2,:)+eeg(tr,8,:))/4);%laplacian filter on Fz
end


e = eeg;
% 2. Downsample to 128 Hz
for tr = 1:size(e,1)
    for ch = 1:size(e,2)
        EEG(tr,ch,:) = downsample(e(tr,ch,:), 16);
    end
end
time = downsample(time,16);
Fs = 128; 
zero_s = 26;
ms170 = 49;

% 3. lowpass and highpass filter each trial
[b_l,a_l]=butter(4, 5/(Fs*0.5),'low'); 
[b_h,a_h]=butter(4, 1/(Fs*0.5),'high'); 
for tr = 1:size(EEG,1)
    for ch = 1:size(EEG,2)
        EEG(tr,ch,:) = filtfilt(b_l,a_l,squeeze(double(EEG(tr,ch,:))));%low pass filter
        EEG(tr,ch,:) = filtfilt(b_h,a_h,squeeze(double(EEG(tr,ch,:))));%high pass filter
    end
end

% 4. remove baseline from each trial ie remove the mean of prestimulus 200ms from all samples
for tr = 1:size(EEG,1)
    for ch = 1:size(EEG,2)
        baseline = mean(squeeze(EEG(tr,ch,1:zero_s)));%mean of first 200 ms 
        EEG(tr,ch,:) = EEG(tr,ch,:)-baseline;
    end
end

%% For average ERP of all expressions together 
% 5. Compute ERP
ERP = squeeze(mean(EEG,1));

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

%% ERP of individual expressions separately

% 5. separate processed trials of each expression (1-Angry, 2-Happy, 3-Sad, 4-Surprised)
cnt_ang = 1; cnt_hap = 1; cnt_sad = 1; cnt_sur = 1;
for tr = 1:size(EEG,1)
    % Angry
    if emo_all(tr) == 1 
        EEG_ang(cnt_ang,:,:) = EEG(tr,:,:);
        cnt_ang = cnt_ang+1;
    end
    % Happy 
    if emo_all(tr) == 2 
        EEG_hap(cnt_hap,:,:) = EEG(tr,:,:);
        cnt_hap = cnt_hap+1;
    end
    % Sad 
    if emo_all(tr) == 3 
        EEG_sad(cnt_sad,:,:) = EEG(tr,:,:);
        cnt_sad = cnt_sad+1;
    end
    % Surprised
    if emo_all(tr) == 4 
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
   plot(time, ERP_ang(i,:)','Linewidth',1.5,'color',[0.8 0 0]);hold on; 
   plot(time, ERP_hap(i,:)','Linewidth',1.5,'color',[0 0.5 0]);hold on; 
   plot(time, ERP_sad(i,:)','Linewidth',1.5,'color',[0.5 0 0.5]);hold on; 
   plot(time, ERP_sur(i,:)','Linewidth',1.5,'color',[1 0.7 0]);hold on; 
   axis tight;
   plot(zeros(1:2),ylim,'k','Linewidth',1);
   plot(xlim,zeros(1:2),'k','Linewidth',1);
   plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1); %170 ms is 49th value in time
   title(chs(i));
end

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

%% Figure for Journal Paper %%%%%%%%%
zero_s = 26;
ms170 = 49;
minus100ms = 13;
scale = 1.5; %because laplacian filter reduces amplitude so restore it to original 

chs ={'Cz'};
channel = 8;

subplot(1,2,1), 
p1 = plot(time(minus100ms:end), ERP(channel,minus100ms:end).*scale','Linewidth',1.5);hold on;
axis tight;
plot(zeros(1:2),ylim,'k','Linewidth',1);
plot(xlim,zeros(1:2),'k','Linewidth',1);
plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1);
title('Cz','FontSize', 16); xlabel({'Time (ms)';'.'}); ylabel('Amplitude (\muV)');
set(gca,'FontSize', 16)


subplot(1,2,2), 
p2 = plot(time(minus100ms:end), ERP_ang(channel,minus100ms:end).*scale','Linewidth',1.5, 'color',[0.8 0 0]);hold on;
p3 = plot(time(minus100ms:end), ERP_hap(channel,minus100ms:end).*scale','Linewidth',1.5, 'color',[0 0.5 0]);hold on;
p4 = plot(time(minus100ms:end), ERP_sad(channel,minus100ms:end).*scale','Linewidth',1.5, 'color',[0.5 0 0.5]);hold on;
p5 = plot(time(minus100ms:end), ERP_sur(channel,minus100ms:end).*scale','Linewidth',1.5, 'color',[1 0.7 0]);hold on;
axis tight;
plot(zeros(1:2),ylim,'k','Linewidth',1);
plot(xlim,zeros(1:2),'k','Linewidth',1);
plot([time(ms170),time(ms170)],ylim,'k--','Linewidth',1);
title('Cz'); 
set(gca,'FontSize', 16)


legend([p1 p2 p3 p4 p5],'All emotions', 'Angry', 'Happy', 'Sad', 'Surprised','Orientation','horizontal', 'FontSize', 18)


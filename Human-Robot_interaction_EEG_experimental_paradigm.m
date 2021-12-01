%% Miko Experimental Paradigm
% Maitreyee Wairagkar (c) 2019
% Last Update: 27/03/2019
% The experiment contains three components:
% 1. EEG Acquisition
% 2. Webcam video to record stimuli on Miko
% 3. Stimulus cue generaator for the Miko operator

%Turn off the sound on the phone to avoid Miko making sound
%Keep Miko on Charging to avoid the flashing lights on Miko
%Keep Miko on an elevated platform such that it won't move
%% Check the framerate of camera
webcamlist
web_cam_no=2;
cam = webcam(web_cam_no);
preview(cam);
%% Set some variables
ParticipantID = 'P10';
run_number = 'run_6'; %%% REMEMBER TO CHANGE RUN NUMBER

filename = strcat(ParticipantID,'_',run_number);
vid_filename = strcat(filename,'.avi');

webcamlist
web_cam_no=2;

num_trials = 10; %number of trials for each condition in a single run
num_expressions = 4; %there are 4 conditions
cue =[];
for i = 1:num_trials  
    cue = [cue,randperm(num_expressions)];
end
cue_prep_time=4; %4sec
cue_send_time=4; %4sec

%emotion={'Angry','Disgust','stern','surprised','happy', 'sad','afraid','tired'};
emotion={'Angry','Happy','Sad','Surprised'};

close all
%% Set up everything

%Set up the EEG device TMSi

%Example of plotting and saving to Poly5 file at the same time.
%

% Subset of channels to show
channel_subset = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; % only show 16 bipolar channels 

% Step 1: Setup library and choose connection type ( 'usb', 'bluetooth', 'network' or 'wifi' ).
library = TMSi.Library('usb');

% Step 2: Find device to connect to. Keep on trying every second.
while numel(library.devices) == 0
    library.refreshDevices();
    pause(1);
end

try
    % Step 3: Get first device an retrieve information about device.
    device = library.getFirstDevice();

    % Step 4: Create a sampler with which we are going to retrieve samples.
    sampler = device.createSampler();

    try
        % Step 5: Set settings for sampler.
        % sampler.setSampleRate(2048);
        % sampler.setReferenceCalculation(true);

        % Step 6: Create a RealTimePlot.
        realTimePlot = TMSi.RealTimePlot('RealTimePlot Example', sampler.sample_rate, device.channels(channel_subset));
        realTimePlot.setWindowSize(10);
        realTimePlot.show();

        % Step 7: Create a Poly5 file to save samples to.
        poly5 = TMSi.Poly5(strcat(filename,'.Poly5'), filename, sampler.sample_rate, device.channels(channel_subset)); 

        % Step 8: Connect to device through the sampler.
        sampler.connect();
        
    catch matlabException
        warning(sprintf('Error while trying to connect sampler to device.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
    end
    
catch matlabException
    warning(sprintf('Error while trying to create device/sampler.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
end

% Set up webcam
cam = webcam(web_cam_no);
%preview(cam);
vidWriter = VideoWriter(vid_filename);
open(vidWriter);



%% Start the experiment
index = 1;
cnt = 1;
cue_on_flag = 1;
disp(emotion{cue(cnt)})
old_time_cue_prep = now;

% Sample the EEG from device
    try
        % Step 9: Start sampling.
        sampler.start();
        
        % Step 10: Sample as long as plot is visible. 
        %   Can be closed with 'q' or cross.
        %   All ranges can be set by 'r', a dialog allows to put in a range.
        %   Press 'a' for autoscale.
        samples_with_time =[];
        while realTimePlot.is_visible
            
            %Show the cue to send the emotion
        if (cnt <=length(cue))    
            if (etime(datevec(now),datevec(old_time_cue_prep))>cue_prep_time) && cue_on_flag==1
                disp([emotion{cue(cnt)},' send']);
                cue_time(cnt,:)= [cue(cnt), now];

                cnt = cnt+1;
                cue_on_flag = 0;
                old_time_cue_prep = now;
            elseif (etime(datevec(now),datevec(old_time_cue_prep))>cue_send_time) && cue_on_flag==0
                disp(emotion{cue(cnt)});
                cue_on_flag = 1;
                old_time_cue_prep = now;
            end
        else
            disp('End of Run')
        end 
            % Acquire frame for processing
            img = snapshot(cam);
            
            samples = sampler.sample();

            % add timestamp to the samples as additional row %MNW
            samples_with_time = [samples_with_time, [samples(channel_subset,:); ones(1,size(samples(channel_subset,:),2))*now]];
            % add timestamp to video
            video_time(index)=now; %save timestamp
            index = index+1;
            % Step 11: Append samples to plot.
            realTimePlot.append(samples(channel_subset,:));

            % Step 12: Append and draw samples to poly5 file.
            poly5.append(samples(channel_subset,:));
            realTimePlot.draw();
            
            % Write frame to video
            writeVideo(vidWriter,img);
    
        end

        % Step 13: Close file.
        poly5.close();
        
    catch matlabException
        warning(sprintf('Error while trying to sample.\n\tIdentifier: %s\n\tMessage: %s\n', matlabException.identifier, matlabException.message));
    end
    
%correct the time stamps on EEG
cnt2 = 1;
while cnt2<(size(samples_with_time,2))
    
start = cnt2; 
while (samples_with_time(end,cnt2)==samples_with_time(end,cnt2+1)) && ((cnt2+1)<=size(samples_with_time,2))
    cnt2=cnt2+1;
    if (cnt2+1)>=size(samples_with_time,2)
        break;
    end
end

diff = cnt2+1-start;
incriment = (samples_with_time(end,cnt2+1)-samples_with_time(end,1))/diff;
%%{
for i=start+1:cnt2
samples_with_time(end,i)=samples_with_time(end,i-1)+incriment;
end
%}
cnt2 = cnt2+1;
end 
    
%% Disconnect everything
% Step 14: Stop sampler.
sampler.stop();

% Step 15: Disconnect with device.
sampler.disconnect();

% Step 16: Cleanup library.
library.destroy();

clear cam vidWriter

%% save 
variable_name = strcat(filename,'.mat');
save(variable_name,'samples_with_time','video_time','cue_time');

%% read video to check
obj = VideoReader(vid_filename);
k=1;
while hasFrame(obj)
  frame(k).fr = readFrame(obj);
  k=k+1;
  %image(this_frame, 'Parent', thisax);
end
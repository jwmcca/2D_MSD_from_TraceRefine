%% Calculate 2D MSD from TraceRefine Files.
% Section 1. 
% Script written by Joshua McCausland in the Jacobs-Wagner Laboratory,
% 2022.
%
% Input: A directory of TraceRefine files, named however you wish. 
% Note that you must define four major parameters before you can execute
% this script.
%   - Px_size: This is the pixel size in your experiment, calculated in
%       micrometers.
%
%   - ExpT: This is the "exposure time." Define the time per frame in
%       seconds.
%
%   - traj_length: This is the minimum length of the trajectories to
%       consider. I usually only go with 10 for my experiments, as I 
%       feel I have good resolution with that number to determine if 
%       a trajectory is real.
%
%   - TraceRefine_File_Folder: This is the source folder where you keep the
%       TraceRefine files. This will be a prompt where Matlab will ask you
%       to select your folder with TraceRefine files.
%
% The output is a single matrix termed "Raw_MSDs." This has the calculated
% MSDs of all the trajectories concatenated on one another, so it will be a
% very large matrix. The code then will take these Raw MSDs and
% find the average. The average MSD is used to compute the diffusion
% coefficient, D, as well as the confinement term, alpha. This structure is
% termed "MSD." The fourth column is the number of trajectories used to
% compute each step.
%
% This script does not compute the D for individual MSD traces from every
% trajectory. Instead, the SEM for the average D and alpha is estimated by 
% extrapolating from the residuals.

clear; clc;

%%%%%%%% Experiment Parameters %%%%%%%%%%%%%%%%
px_size = 0.16; %micrometers.
ExpT = 0.03; %seconds.

% Minimum trajectory length to consider. 
traj_length = 10; 

% Source subfolder with all of the trace refine files.
TraceRefine_File_Folder = uigetdir('Please select the TraceRefine folder.');

%MSD Equation for fitting. 
msd_eqn2D = @(P,x) 4*P(1)*x.^(P(2)) + P(3); %2D MSD Equation.

%%%%%%% Search through the TR directory and make MSDs for each particle. %%%%%%%%
Raw_MSDs = [];
MSD_params = [];
opts = optimset('Display','off');

currentPath = pwd;
TR_files = dir([TraceRefine_File_Folder '/*.mat']);
for idxa = 1:length(TR_files)
    load([TR_files(idxa).folder '/' TR_files(idxa).name]);
    TracksROI = tracksRefine.TracksROI;
    for idxb = 1:length(TracksROI)
        MSD_temp = [];
        traj = TracksROI(idxb).Coordinates; %Defines the trajectory to plot.
        if length(traj) >= traj_length
            MSD_temp = MSD_Calc(traj,ExpT,px_size); %Create the MSD matrix.

            %Fit the individual raw MSD
            if sum(isnan(MSD_temp(:))) == 0
                
                MSD_params(end+1,:) = lsqcurvefit(msd_eqn2D,[0 1 0],MSD_temp(:,1),MSD_temp(:,2),[0 0 0],[1e7 1e7 1e7],opts);
           
                % Concatenate the raw MSDs
                Raw_MSDs = [Raw_MSDs; MSD_temp];
            end
        end
    end
end

%%%%%%% Construct an average MSD from all MSDs %%%%%%%%
MSD = [];
for idxd = 1:max(unique(Raw_MSDs(:,1)))
    MSD_temp = [];
    MSD_temp = Raw_MSDs(find(Raw_MSDs(:,1) == idxd),2);
    MSD(idxd,1) = idxd; % Frame
    MSD(idxd,2) = nanmean(MSD_temp); % MSD
    MSD(idxd,3) = nanstd(MSD_temp)/sqrt(length(MSD_temp)); % Standard Err
    MSD(idxd,4) = length(MSD_temp); % N
end
    
MSD(:,1) = ExpT * MSD(:,1); % Time (s)

output = [MSD_params];
VariableNames = {'D (um2/s)','alpha','D0'};
outputparameters = array2table(output,'VariableNames',VariableNames);
writetable(outputparameters,['IndMSD_Fit_Params.xlsx'],'sheet','FitResults');
save('IndMSD_Fit_Params.mat','MSD_params')

%% Section 2. Plot MSD
% Requires the Matlab package "rgb," included in this directory. 
% rgb citation: Ben Mitch, https://www.mathworks.com/matlabcentral/fileexchange/1805-rgb-m
%
%
% For this section, you have inputs that specify the output Matlab plot.
% However, all of the MSD data (including the fit results) are exported as
% an excel file, so this can be replotted however you wish later. The most
% important parameter to note here is "index." You should pick an index
% that captures the most number of trajectories in your MSD data set. 
% For example, if your smallest trajectory has only 10 steps, then you 
% should plot an MSD of up to 10 time lags. If you plot 11 or more, then 
% you are estimating an MSD that does not include these smaller 
% trajectories, so your estimate is not as accurate. That's also why the 
% error bars get bigger the farther you go out.


close all;

% Set up parameters
index = 11; %Farthest point to plot. Do not exceed smallest trajectory length.
SaveName = 'Avg_MSD'
figure_size = [5,4] % [width x height], inches
y_limits = [0, 0.06]
SaveFormat = 'PDF' % Input either PDF or PNG
save_plot = 0 % 0 if you don't want to save the figure, 1 if you do.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the MSD
msd_eqn2D = @(P,x) 4*P(1)*x.^(P(2)) + P(3); %2D MSD Equation.
[param,~,residual,~,~,~,jacobian] = lsqcurvefit(msd_eqn2D,[0 1 0],MSD(1:index,1),MSD(1:index,2));

% Use the residuals and jacobian to estimate the error for each parameter.
ci = nlparci(param,residual,jacobian);
t = tinv(1-0.05/2,index-3); %Degrees of freedom must be calculated. ('INDEX' MSD observations - 3 parameters, see index variable above.)
se_msd = (ci(:,2)-ci(:,1)) ./ (2*t); %SEM for fitted parameters D, alpha, and D0.

xfit = linspace(0,MSD(index,1),500);
yfit = msd_eqn2D(param,xfit); 

% Plot the MSD. 
f1 = figure('Units','inches','Position',[0 0 figure_size],'PaperUnits','inches','PaperPosition',[0 0 figure_size],'PaperSize',[figure_size],'CreateFcn','movegui center');
e1 = errorbar(MSD(1:index,1),MSD(1:index,2),MSD(1:index,3),'ok');
hold on;
p1 = plot(xfit,yfit,'-r','LineWidth',1.5,'Color',rgb('FireBrick'));
leg = legend([p1],{['D = ' num2str(round(param(1),3)) ' \pm ' num2str(round(se_msd(1),3)) '\mum^2/s' 10,...
    '\alpha = ' num2str(round(param(2),2)) ' \pm ' num2str(round(se_msd(2),3)) 10,...
    'N = ' num2str(MSD(index,4))]},'Box','off','location','southeast');
ylim(y_limits)
xlabel('Time (s)');
ylabel('MSD (\mum^2)');
set([e1],'MarkerSize',5,'LineWidth',1,'MarkerEdgeColor',rgb('Black'),'MarkerFaceColor',rgb('LightGray'));
set(gca,'LineWidth',1.75,'FontSize',14,'TickDir','Out','Box','Off','XColor',[0 0 0],'YColor',[0 0 0]);
set(leg,'FontSize',10)

if strcmp(SaveFormat,'PDF') & save_plot
    print(f1,[SaveName '.pdf'],'-dpdf');
elseif strcmp(SaveFormat,'PNG') & save_plot
    print(f1,[SaveName '.png'],'-dpng','-r600');
end

% Save the MSD. Fourth column is the number of trajectories that went into
% each time lag.

output = [MSD(1:index,:)];
VariableNames = {'Time Lag (s)','MSD (um2)','StdEr','N'};
outputparameters = array2table(output,'VariableNames',VariableNames);
writetable(outputparameters,[SaveName '.xlsx'],'sheet','MSD');

% Save the fits for plotting.

output = [xfit',yfit'];
VariableNames = {'xfit','yfit'};
outputparameters = array2table(output,'VariableNames',VariableNames);
writetable(outputparameters,[SaveName '.xlsx'],'sheet','Fit');

% Save all fit parameters as well as "N," the total trajectories used in
% MSD. Note the three parameters:
%   - D: Diffusion Coefficient
%   - alpha: measure of confinement (alpha < 1) or super diffusion (alpha >
%       1)
%   - D0: Uncertainty. Used for adjusting starting y value and
%       approximating best fit.

output = [param(1), se_msd(1), param(2), se_msd(2), param(3), se_msd(3), MSD(index,4)];
VariableNames = {'D (um2/s)','Derror','alpha','aerror','D0','D0error','N'};
outputparameters = array2table(output,'VariableNames',VariableNames);
writetable(outputparameters,[SaveName '.xlsx'],'sheet','FitResults');

%% Section 3. Plot distribution of D and alpha. 
% This is just a plotting section for quick visualization of the diffusion
% coefficient and alpha values from the MSD fit. They are also separately
% saved as an excel file for your own plotting purposes.
% Individual D and alpha determined from fitting each trajectory's MSD.

close all;

%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%
SaveName = 'MSD_Params' %Save name for the final plot.
SaveFormat = 'PDF' % Input either PDF or PNG
save_plot = 0 % 0 if you don't want to save the figure, 1 if you do.

figure_size = [10,3] % [width x height], inches
num_bins_diffusion = 20; % number of bins for diffusion coefficient histogram.
log_diffusion = 1; % 0 or 1. Input "1" if the histogram should be in log scale.

num_bins_alpha = 20; % number of bins for the alpha coefficient histogram.


%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%
if log_diffusion
    [~,edges] = histcounts(log10(MSD_params(:,1)),num_bins_diffusion);
    edges = 10.^edges;
else
    [~,edges] = histcounts(MSD_params(:,1),num_bins_diffusion);
end

f = figure('Units','inches','Position',[0 0 figure_size],'PaperUnits','inches','PaperPosition',[0 0 figure_size],'PaperSize',[figure_size],'CreateFcn','movegui center');
ax1 = subplot(1,2,1)
histogram(MSD_params(:,1),edges,'Normalization','Probability','FaceColor',rgb('LightGray'),'EdgeColor','k','LineWidth',1.5);
title('Diffusion Coefficient')
xlabel('D (\mum^2/s)')
ylabel('Fraction')

if log_diffusion
    set(ax1,'XScale','Log')
end

ax2 = subplot(1,2,2)
histogram(MSD_params(:,2),num_bins_alpha,'Normalization','Probability','FaceColor',rgb('LightGray'),'EdgeColor','k','LineWidth',1.5);
title('Alpha')
xlabel('\alpha')

set(f.Children,'Box','Off','TickDir','Out','LineWidth',1.5,'FontSize',14,'YColor','k');
if strcmp(SaveFormat,'PDF') & save_plot
    print(f,[SaveName '.pdf'],'-dpdf');
elseif strcmp(SaveFormat,'PNG') & save_plot
    print(f,[SaveName '.png'],'-dpng','-r600');
end

%% Functions for this script.
function MSD = MSD_Calc(traj,ExpT,pixelsize);
    frames = traj(:,1)-traj(1,1)+1; %Makes the frame count starting from one. 
    r=[]; %creates a blank matrix named "r."
    MSD=[]; %creates a blank matrix named "MSD."
    N = (max(traj(:,1))-min(traj(:,1)));
    for j = 1:N;
        r = [];
        for i = 1:min(frames(end)-j,length(frames));
            Endf=find(frames==(frames(i)+j)); %End frame for calculating endpoint.
            if Endf > 0;
            Endp=traj(Endf,2:3); %defines the endpoint for timelag.
            r(i,1) = sum((((traj(i,2:3)-Endp)*(pixelsize)).^2)); %squared displacement in micrometers
            end;
            %r = [r; squared]; %Writes the matrix
        end;
        MSD(j, 1) = j; %First column, time lag.
        MSD(j, 2) = mean(r); %Second column, MSD value.
        MSD(j, 3) = std(r)/sqrt(length(r)); %Third column, standard error.
    end;
    MSD(:,1) = MSD(:,1);
end
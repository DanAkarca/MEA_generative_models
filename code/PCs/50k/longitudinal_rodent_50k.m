%% Simulate developmental trajectories of 50k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data
clear; clc;
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/rodent_50k_001Hz_trajectory_data');
%% set data from previous run
integration_comparison = rodent_50k_001Hz_trajectory_data.efficiency_comparison;
integrationi = rodent_50k_001Hz_trajectory_data.efficiency_trajectory.simulated;
plotdatai = rodent_50k_001Hz_trajectory_data.efficiency_trajectory.observed;
segregation_comparison = rodent_50k_001Hz_trajectory_data.modularity_comparison;
segregationi = rodent_50k_001Hz_trajectory_data.modularity_trajectory.simulated;
plotdatas = rodent_50k_001Hz_trajectory_data.modularity_trajectory.observed;
%% scale simulated segregation against observed segregation
% visualise the simulations of these scaled simulations
h = figure;
h.Position = [100 100 400 400];
for i = 1:length(segregationi);
    % set the simulated scale here
    plot(linspace(0,14,length(segregationi{i})),segregationi{i},...
        'color',[169 169 169]./256,...
        'linewidth',1);
    % plot networks on top
    hold on;
end
% visualise the observed segregation data on top
hold on;
iosr.statistics.boxPlot(plotdatas,...
    'theme','colorall',...
    'boxcolor',[202 143 66]./256);
xticks([0 7 14]);
xlim([-1 15]);
xticklabels({'0%','50%','100%'});
ylabel('Modularity Q'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
box off;
%% scale simulated integration against observed integration
% visualise the simulations of these scaled simulations
h = figure;
h.Position = [100 100 400 400];
for i = 1:length(integrationi);
    % set the simulated scale here
    plot(linspace(0,14,length(integrationi{i})),integrationi{i},...
        'color',[169 169 169]./256,...
        'linewidth',1);
    % plot networks on top
    hold on;
end
% visualise the observed integration data on top
hold on;
iosr.statistics.boxPlot(plotdatai,...
    'theme','colorall',...
    'boxcolor',[102 167 197]./256);
xticks([0 7 14]);
xlim([-1 15]);
xticklabels({'0%','50%','100%'});
ylabel('Global Efficiency'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
box off;
%% compute variance explained
% compute correlations
a = squeeze(integration_comparison(:,:,1));
b = squeeze(integration_comparison(:,:,2));
c = squeeze(segregation_comparison(:,:,1));
d = squeeze(segregation_comparison(:,:,2));
% integration
h = figure; h.Position = [100 100 400 400];
[r p] = corr(a(:),b(:));
facecolor = 1-copper(4);
for div = 1:4;
scatter(a(:,div),b(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
sgtitle(sprintf('R^2=%.3g, r=%.3g, p=%.3d',r^2*100,r,p));
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0 0.4]); ylim([0 0.4]); xticks([0:0.1:0.4]);
hold on;
end
% segregation
h = figure; h.Position = [100 100 400 400];
[r p] = corr(c(:),d(:));
facecolor = copper(4);
for div = 1:4;
scatter(c(:,div),d(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
sgtitle(sprintf('R^2=%.3g, r=%.3g, p=%.3d',r^2*100,r,p));
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0.1 0.9]); ylim([0.1 0.9]); 
hold on;
end
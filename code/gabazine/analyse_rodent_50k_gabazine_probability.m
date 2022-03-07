%% Analyse experiemental probability distribiutions
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load experimental data
clear; clc;
% load 0.01Hz gabazine networks
a=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03200_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag5_jitter5_prob.mat');
b=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03200_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
c=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03200_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag20_jitter20_prob.mat');
d=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03212_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag5_jitter5_prob.mat');
e=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03212_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
f=load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03212_rec_1800s_min_rate_0.01hz_alpha_0.001_1800s_min_rate_001hz_lag20_jitter20_prob.mat');
gabazine = {a.gm_data b.gm_data c.gm_data};
control = {d.gm_data e.gm_data f.gm_data};
data = {gabazine control};
% addpaths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master');
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
addpath('/imaging/astle/users/da04/PhD/toolboxes/bluewhitered');
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/functions/');
% load generative model data
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/rodent_gabazine_50k_001Hz_generative_models.mat');
% load matching probability distributions
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/rodent_gabazine_50k_001Hz_matching_probability.mat');
% addpaths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master');
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
addpath('/imaging/astle/users/da04/PhD/toolboxes/bluewhitered');
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/functions/');
% equation
eqn = string({'KS_k','KS_c','KS_b','KS_d'});
%% plot probability distributions at half way through the simulation
% take data
matching_P = rodent_gabazine_50k_001Hz_matching_probability;
% visualise
u = figure; u.Position = [100 100 600 400];
% set which group to plot
plotg = 1;
% set the smoothing factor
sm = 10;
% set cultures from each group to plot
gabazine_cultures = [1 2 3];
control_cultures = [4 5 6 7 8 9];
% set proportion
proportions = [.1:.01:1];
% col
gabacol = 1-copper(2*length(proportions)); gabacol = gabacol(4:length(proportions)+3,:);
contcol = pink(2*length(proportions)); contcol = contcol(4:length(proportions)+3,:);
cols = {gabacol contcol};
% initialise
gabazine_keep = []; control_keep = [];
pprop = []; dprop = [];
% loop over proportions
for j = 1:length(proportions);
    proportion = proportions(j);
    % set scale limits
    lim = 1000;
    % take the mid point distributions
    distribution = {}; distribution_scaled = {}; distribution_scaled_sized = zeros(9,lim);
    for culture = 1:9;
        A = matching_P{culture};
        t = size(A,1);
        n = size(A,2);
        B = squeeze(A(round(t*proportion),:,:));
        C = B.*triu(ones(n),1);
        %distribution{culture} = C(C>0); % removed non zeros
        distribution{culture} = C;
        distribution_scaled{culture} = imresize(distribution{culture},[lim 1]);
        distribution_scaled_sized(culture,:) = rescale(distribution_scaled{culture},0,1);
    end
    % average across groups - pick group here
    gabazine_mean = mean(distribution_scaled_sized(gabazine_cultures,:));
    control_mean = mean(distribution_scaled_sized(control_cultures,:));
    % compute statistics
    [~,pprop(j)] = ttest(gabazine_mean',control_mean');
    dprop(j) = computeCohen_d(gabazine_mean,control_mean);
    % keep these
    gabazine_keep(j,:) = gabazine_mean;
    control_keep(j,:) = control_mean;
    % plot ksdensity 
    if plotg == 1;
        [f x] = ksdensity(gabazine_mean);
        h = plot(x,smooth(f,sm),'linewidth',3); h.Color = gabacol(j,:);
    else if plotg == 2;
            [fy y] = ksdensity(control_mean);
            h = plot(y,smooth(fy,sm),'linewidth',3); h.Color = contcol(j,:);
        end
    end
    xlabel('Probability score (P_i_j)'); ylabel('Frequency'); 
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
    xlim([0 1.1]); ylim([0 5]); yticks([]); box off;
    hold on;
end
c=colorbar; c.Ticks = [0:0.1:1]; c.TickLabels = {'0%','','','','','50%','','','','','100%'}; 
caxis([0 1]); 
c.Label.String = 'Simulation time'; c.Label.FontSize = 25; c.Label.FontName = 'Arial';
colormap(cols{plotg});
%% compute a visualisation between controls and gabazine
% average over all
mean_control_mean = mean(control_keep);
mean_gabazine_mean = mean(gabazine_keep);
% take std
std_control_mean = std(control_keep);
std_gabazine_mean = std(gabazine_keep);
% new colours
gabacol = 1-copper(20); gabacol = gabacol(4:13,:);
contcol = pink(20); contcol = contcol(4:13,:);
cols = {gabacol contcol};
% statistical test
[h p] = ttest(mean_control_mean,mean_gabazine_mean);
d = computeCohen_d(mean_control_mean,mean_gabazine_mean); 

% set figure density
u = figure; u.Position = [100 100 600 400];
% plot ksdensity 
[f x] = ksdensity(mean_gabazine_mean);
h = plot(x,smooth(f,sm),'linewidth',6); h.Color = gabacol(4,:);
hold on;
[fy y] = ksdensity(mean_control_mean);
h = plot(y,smooth(fy,sm),'linewidth',6); h.Color = contcol(4,:);
xlabel('Probability score (P_i_j)'); ylabel('Frequency'); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlim([0 inf]); 
yticks([]); box off;
% plot legend
legend(sprintf('p=%g, \nd=%g',p,d));

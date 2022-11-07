%% supplementary figure 5
% written by danyal akarca
%% load data
clear; clc;
% set directory
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/';
% load generative models
load(strcat(gendir,'50k_rodent_div14_randomization_exp/randomized_D/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/combined/results_generative_models')); random_d = results_generative_models;
load(strcat(gendir,'50k_rodent_div14_randomization_exp/randmio/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/combined/results_generative_models')); random = results_generative_models;
load(strcat(gendir,'50k_rodent_div14_randomization_exp/full_random_network/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/combined/results_generative_models')); random_all = results_generative_models;
load(strcat(gendir,'50k_rodent_div14_te/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_TE_allin/combined/results_generative_models')); te = results_generative_models;
load(strcat(gendir,'50k_rodent_div14_te_tril/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_TE_tril_allin/combined/results_generative_models')); te_l = results_generative_models;
load(strcat(gendir,'50k_rodent_div14_te_triu/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_TE_triu_allin/combined/results_generative_models')); te_u = results_generative_models;
% load transfer entropy
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/TE_results.mat');
spk_data = {};
for well = 1:6;
    spk_data{well} = load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/50k_rodent_div14_te/spk_data_%g.mat',well));
end
% collect the data together
data = {random_d random_all random te te_l te_u};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
% addpath of roundsd
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
% addpath of MIT giant
addpath('/imaging/astle/users/da04/PhD/toolboxes/MIT_code_2011/');
%% take energy values
% set number of models
nmodels = 13;
% set the sample size for the 50k and 100k cultures
ngroup = [6 6 6 6 6 6];
% define the top n parameters
nset = [1 10 50 100];
% compute the minimum
for group = 1:length(ngroup);
    % get the group
    energy = data{group}.energy;
    parameters = data{group}.parameters;
    % initialise
    top_energy = cell(length(nset),1);
    top_energy_mean = zeros(length(nset),13,ngroup(group));
    top_parameters = cell(length(nset),1);
    top_parameters_mean = zeros(length(nset),13,ngroup(group),2);
    % ouput
    null_generative{group} = struct;
    % loop over sets
    for no = 1:length(nset)
        % take the actual amount of top performing parameters
        nnet_100k = nset(no);
        % loop over models
        for model = 1:nmodels;
            % take energies for this model
            pipeline_d = squeeze(energy(:,model,:))';
            % rank them for each subject
            [v i] = sort(pipeline_d);
            % take top n energies and their indices
            n_e = v(1:nnet_100k,:);
            n_i = i(1:nnet_100k,:);
            % take the corresponding parameters
            u = zeros(ngroup(group),nnet_100k,2);
            for s = 1:ngroup(group);
                % keep parameters
                u(s,:,:) = squeeze(parameters(s,model,n_i(:,s),:));
            end
            % if top parameter only
            if nnet_100k == 1
                % squeeze the matrices
                u = squeeze(u);
                % assign
                top_energy{no}(model,:) = n_e';
                top_parameters{no}(model,:,:) = u;
                % and assign it to the mean
                top_energy_mean(no,model,:) = n_e';
                top_parameters_mean(no,model,:,:) = u;
            else
                top_energy{no}(model,:,:) = n_e';
                top_parameters{no}(model,:,:,:) = u;
                % keep a mean value too
                top_energy_mean(no,model,:) = squeeze(mean(n_e',2));
                top_parameters_mean(no,model,:,:) = squeeze(mean(u,2));
            end
        end
    end
    null_generative{group}.top_energy = top_energy;
    null_generative{group}.top_energy_mean = top_energy_mean;
    null_generative{group}.top_parameters = top_parameters;
    null_generative{group}.top_parameters_mean = top_parameters_mean;
end
%% visualise randomness models
% take the top n energy values
ncomb = 1;
% nulls to test
nnull = [1];
% set null labels
dataplot = [];
for null = 1:length(nnull);
    dataplot(:,:,null) = squeeze(null_generative{nnull(null)}.top_energy_mean(ncomb,:,:));
end
i = [2:3 9:13 4:8 1];
data_plot = dataplot(i,:,:);
% permute the order
data_plot = permute(data_plot,[2 3 1]);
% iosr boxplot - all models
% visualise
h = figure; h.Position = [100 100 400 450];
u = iosr.statistics.boxPlot(data_plot,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticks([]); xlabel('{\itK_{ij}} only');
ylim([0 1]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% visualise randomness models
% take the top n energy values
ncomb = 1;
% nulls to test
nnull = [2 3];
% set null labels
nulllabels = {'Random','{\itrandmio(5)}'};
dataplot = [];
for null = 1:length(nnull);
    dataplot(:,:,null) = squeeze(null_generative{nnull(null)}.top_energy_mean(ncomb,:,:));
end
i = [2:3 9:13 4:8 1];
data_plot = dataplot(i,:,:);
% permute the order
data_plot = permute(data_plot,[2 3 1]);
% iosr boxplot - all models
% visualise
h = figure; h.Position = [100 100 800 450];
u = iosr.statistics.boxPlot(data_plot,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticklabels(nulllabels); xlabel('Random null models');
ylim([0 1]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% visualise transfer entropy models
% take the top n energy values
ncomb = 1;
% nulls to test
nnull = [4 5 6];
% set null labels
nulllabels = {'TE-symmetrical','TE-out connections','TE-in connections'};
dataplot = [];
for null = 1:length(nnull);
    dataplot(:,:,null) = squeeze(null_generative{nnull(null)}.top_energy_mean(ncomb,:,:));
end
i = [2:3 9:13 4:8 1];
data_plot = dataplot(i,:,:);
% permute the order
data_plot = permute(data_plot,[2 3 1]);
% iosr boxplot - all models
% visualise
h = figure; h.Position = [100 100 1200 450];
u = iosr.statistics.boxPlot(data_plot,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticklabels(nulllabels); xlabel('Transfer entropy models');
ylim([0 1]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% visualise transfer entropy
% take empirical data
A = TE_results.well_1;
% set coordinates
coordinates = spk_data{1}.spk_data.coordinates;
% remove nan
A(isnan(A)) = 0;
% visualise
h = figure; h.Position = [100 100 1200 450];
% matrix
subplot(1,2,1);
imagesc(A);
b = gca; b.TickDir = 'out'; xlabel('Node'); ylabel('Node'); 
box off; c = colorbar; c.Label.String = 'Transfer Entropy';
b.FontName = 'Arial'; b.FontSize = 25;
subplot(1,2,2);
% remove non-nodes
ids = find(degrees_und(A));
Ag = A(ids,ids);
coords = coordinates(ids,:);
d = squareform(pdist(coords));
% get backbone of each upper triangle
k = 2;
upper = triu(mst,1); upper = upper + upper'; umst = backbone_wu(upper,k);
lower = tril(mst,1); lower = lower + lower'; lmst = backbone_wu(lower,k);
mst = zeros(size(Ag,1)); mst(find(triu(umst,1))) = 1; mst(find(tril(lmst,1))) = 1;
% plot
g = digraph(mst);
u = plot(g,...
    'XData',coords(:,1),...
    'YData',coords(:,2),...
    'markersize',.1*degrees_und(Ag),...
    'nodelabel',[],...
    'edgecolor',[.45 .45 .45],...
    'edgealpha',.4,...
    'linewidth',1,...
    'nodecolor',[237 136 152]./256,...
    'arrowsize',12);
axis off;
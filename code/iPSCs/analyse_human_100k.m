%% generative modelling findings of 100k human ipsc hd gnm data
% written by danyal akarca
%% load data for rat_primary_cortex matching
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_ipsc_DIV28_30min');
% load sttc data
gluta = load('GN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
motor = load('MN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
dopa = load('DN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
% collect the data together
human_data = {gluta.gm_data motor.gm_data dopa.gm_data};
% load generative models
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/human_100k_001Hz_generative_models');
% addpaths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master');
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
addpath('/imaging/astle/users/da04/PhD/toolboxes/bluewhitered');
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/functions/');
% define model types
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
eqn = string({'KS_k','KS_c','KS_b','KS_d'});
%% form labels
% 3 types
type_labels = string({'GN','MN','DN'});
% 3 pipelines
pipeline_labels = string({'10'});
% 4 time points
div_labels = string({'28'});
%% set hyperparameters of the dataset
% number of types
ntypes = 3;
% number of cultures
ncultures = [8 7 6];
% number of samples
nsamp = ncultures(1)+ncultures(2)+ncultures(3);
% nuber of models
nmodels = 13;
%% get all the cultures
index = size(nsamp,2);
all_o_networks = {};
all_d_networks = {};
step = 1;
% loop through and visualise
for type = 1:ntypes;
    for culture = 1:ncultures(type);
            % form the index
            index(step,1) = type;
            index(step,2) = culture;
            % get the data as a cell array
            all_o_networks{step} = human_data{type}.table_bu{culture};
            all_d_networks{step} = squareform(pdist(human_data{type}.table_xy{culture}));
            % update the step
            step = step + 1;
    end
end 
%% load 0.01Hz generative model data
%{
% initialise
energy_sample = zeros(nsamp,13,20000);
ks_sample = zeros(nsamp,13,20000,4);
networks_sample = cell(nsamp,1);
parameters_sample = zeros(nsamp,13,20000,2);
errors = zeros(nsamp,1);
% loop over parameters and load
step = 1;
for type = 1:ntypes;
    for culture = 1:ncultures(type);
        for model = 1:nmodels;
            % form the string
            str = sprintf(...
                '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_ipsc_DIV28_30min/%s_GNM/1800s_min_rate_001hz_lag10_jitter10_prob/human_ipsc_1_%g_1_generative_model_%g.mat',...
                type_labels(type),culture,model);
            % load
            load(str);
            % keep the data
            energy_sample(step,model,:) = output.energy;
            ks_sample(step,model,:,:) = output.ks;
            networks_sample{step}{model,:,:} = output.networks;
            parameters_sample(step,model,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('%s network_1_%g_1_model %s loaded',type_labels(type),culture,models(model)));
        end
        % update step
        step = step + 1;
    end
end
human_100k_001Hz_generative_models = struct;
human_100k_001Hz_generative_models.energy = energy_sample;
human_100k_001Hz_generative_models.ks = ks_sample;
human_100k_001Hz_generative_models.networks = networks_sample;
human_100k_001Hz_generative_models.parameters = parameters_sample;
human_100k_001Hz_generative_models.models = models;
human_100k_001Hz_generative_models.index = index;
human_100k_001Hz_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 6 cultures - index column 3:4 days in vitro (7 days, 10 days, 12 days, 14 days)'});
human_100k_001Hz_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
human_100k_001Hz_generative_models.author = string('written by dr danyal akarca, 23/01/22');
%}
%% place all data
energy_sample = human_100k_001Hz_generative_models.energy;
ks_sample = human_100k_001Hz_generative_models.ks;
parameters_sample = human_100k_001Hz_generative_models.parameters;
nexist = nsamp;
exist = index;
%% look at energy landscape for a selected rule and network
% select model and network
model = 3;
net = 1;
% take the measure
e = squeeze(energy_sample(net,model,:));
pipeline_p = squeeze(parameters_sample(net,model,:,:));
% visualise
h = figure;
if model == 1
    scatter(pipeline_p(:,1),e,100,e,'.');
    xlabel('eta'); ylabel('energy'); 
    ylim([0 1]);
    caxis([0 1]); c = colorbar; c.Label.String = 'energy';
else
    % plot the energy landscape
    scatter(pipeline_p(:,1),pipeline_p(:,2),100,e,'.'); 
    xlabel('eta'); ylabel('gamma'); 
    caxis([0 1]); c = colorbar; c.Label.String = 'energy';
end
b = gca; b.TickDir = 'out';
%% look at energy landscape for all networks by group
% select model
model = 3;
% select the criteria for selection
criteria = index(:,1)==2;
criteria_ind = find(criteria);
% take the measure
e = squeeze(energy_sample(:,model,:));
p = squeeze(parameters_sample(:,model,:,:));
% breakdown by group to visualise, if wished
e_select = e(:,:);
p_select = p(:,:,:);
nsamp_select = size(e_select,1);
% visualise
if model == 1
    h = figure;
    h.Position = [100 100 600 300];
    eta = squeeze(p(:,:,1));
    m = parula(10000);
    for net = 1:nsamp_select;
        % get colours
        pipeline_d = round(e_select(net,:) .* 10000);
        col = m(pipeline_d,:);
        % plot
        scatter(eta(net,:),e_select(net,:),20,col,'.'); ylabel('energy'); xlabel('eta'); ylim([0 1]);
        b = gca; b.TickDir = 'out';
        b.FontName = 'Arial';
        yticks([]); xticks([]); hold on;
    end
else
% plot the energy landscape
h = figure;
h.Position = [100 100 500 400];
for net = 1:nsamp_select
    scatter(squeeze(p_select(net,:,1)),squeeze(p_select(net,:,2)),20,e_select(net,:),'.'); hold on;
    %xlabel('eta'); ylabel('gamma');
end
caxis([0 1]);
c = colorbar; c.Label.String = 'energy';
b = gca; b.TickDir = 'out';
xticks([]); yticks([]);
end
% set any limits
xlim([-10 10]); ylim([-10 10]);
%% look at ks landscape for all networks
% select model
model = 3;
% take the measure
ks = squeeze(ks_sample(:,model,:,:));
pipeline_p = squeeze(parameters_sample(:,model,:,:));
% plot the energy landscape
h = figure;
h.Position = [100 100 1400 220];
for net = 1:nsamp
    for j = 1:4;
        subplot(1,4,j);
        scatter(squeeze(pipeline_p(net,:,1)),squeeze(pipeline_p(net,:,2)),20,ks(net,:,j),'.'); hold on;
        xlabel('eta'); ylabel('gamma'); title(eqn(j));
        caxis([0 1]); c = colorbar;
    end
end
%% compute summary statistics over the sample
% define the top n parameters
nset = [1 10 50];
% initialise a matrix to keep the data, which will be subject by n by parameters
top_e = cell(length(nset),1);
top_e_mean = zeros(length(nset),13,nexist);
top_p = cell(length(nset),1);
top_p_mean = zeros(length(nset),13,nexist,2);
% compute the minimum
% run just top 2 sets for memory
for no = 1:3;
    % take the actual amount of top performing parameters
    n = nset(no);
    % loop over models
    for model = 1:13;
        % take energies for this model
        pipeline_d = squeeze(energy_sample(:,model,:))';
        % rank them for each subject
        [v i] = sort(pipeline_d);
        % take top n energies and their indices
        n_e = v(1:n,:);
        n_i = i(1:n,:);
        % take the corresponding parameters
        u = zeros(nexist,n,2);
        for s = 1:nexist;
            % keep parameters
            u(s,:,:) = squeeze(parameters_sample(s,model,n_i(:,s),:));
        end
        % if top parameter only
        if n == 1
            % squeeze the matrices
            u = squeeze(u);
            % assign
            top_e{no}(model,:) = n_e';
            top_p{no}(model,:,:) = u;
            % and assign it to the mean
            top_e_mean(no,model,:) = n_e';
            top_p_mean(no,model,:,:) = u;
            
        else
            top_e{no}(model,:,:) = n_e';
            top_p{no}(model,:,:,:) = u;
            % keep a mean value too
            top_e_mean(no,model,:) = squeeze(mean(n_e',2));
            top_p_mean(no,model,:,:) = squeeze(mean(u,2));
        end
    end
    disp(sprintf('set %g of %g complete',no,length(nset)));
end
%% count what drives the energy 
% specify the model
model = 3;
% take the data
ks_data = squeeze(ks_sample(:,model,:,:));
% initialise
driver = [];
% loop over networks
for net = 1:nexist;
    % find the max ks statistics for this network
    [v i] = max(squeeze(ks_data(net,:,:))');
    % group data
    driver(net,:) = [sum(i==1),sum(i==2),sum(i==3),sum(i==4)];,
end
% form a percentage
driver = driver ./ 20000 * 100;
% visualise
figure;
iosr.statistics.boxPlot(driver,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'showScatter',logical(1),...
    'scatterColor',[.5 .5 .5],...
    'scatterAlpha',0.5,...
    'symbolColor',[.5 .5 .5],...
    'boxColor','c',...
    'boxAlpha',0.15); 
ylim([0 100]);
b = gca; b.TickDir = 'out';
xticklabels({'degree','clustering','betweenness','edge length'});
ylabel('max(KS)'); yticklabels({'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'});
%% visualise summary energy for each model
% select which set of nset to view
set = 1;
% set the new model order
i = [2:3 9:13 4:8 1];
% seperate by type
e_select = nan(nmodels,nsamp,ntypes);
for type = 1:ntypes;
    x = squeeze(top_e_mean(set,i,exist(:,1)==type));
    e_select(:,1:size(x,2),type) = x;
end
% permute
data = permute(e_select,[2 3 1]);
% iosr boxplot
h = figure;
h.Position = [100 100 1200 600];
iosr.statistics.boxPlot(data,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylim([0 1]); ylabel('Energy'); xlabel('Human iPSC line');
xticklabels({'Glutamatergic neurons','Motor neurons','Dopaminergic neurons'});
yticks([0:0.1:1]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% group by generative rule
% select which set of nset to view
set = 1;
% set the new model order
i = [2:3 9:13 4:8 1];
% seperate by type
e_select = nan(nmodels,nsamp,ntypes);
for type = 1:ntypes;
    x = squeeze(top_e_mean(set,i,exist(:,1)==type));
    e_select(:,1:size(x,2),type) = x;
end
% permute
data = permute(e_select,[2 3 1]);
% group by meaning
e_mean_rules = [];
e_mean_rules(:,:,1) = squeeze(mean(data(:,:,[1:2]),3));
e_mean_rules(:,:,2) = squeeze(mean(data(:,:,[3:7]),3));
e_mean_rules(:,:,3) = squeeze(mean(data(:,:,[8:12]),3));
e_mean_rules(:,:,4) = squeeze(data(:,:,13));
% group by binning
e_rules = nan(nsamp*5,ntypes+1,4);
a = squeeze(data(:,:,[1:2])); a = permute(a,[2 1 3]); e_rules(1:nsamp*2,1:3,1) = a(:,:)';
b = squeeze(data(:,:,[3:7])); b = permute(b,[2 1 3]); e_rules(1:nsamp*5,1:3,2) = b(:,:)';
c = squeeze(data(:,:,[8:12])); c = permute(c,[2 1 3]); e_rules(1:nsamp*5,1:3,3) = c(:,:)';
d = squeeze(data(:,:,[13])); d = permute(d,[2 1 3]); e_rules(1:nsamp,1:3,4) = d(:,:)';
% add in cerebral organoid data here
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/hCO_energy_rules.mat');
% place in here
e_rules(1:30,4,:) = hCO_energy_rules;
% set xlimits
xlims = [0 0.85];
% visualise
h = figure;
h.Position = [100 100 3000 700];
iosr.statistics.boxPlot(e_rules,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
ylim(xlims); ylabel('Energy'); %xlabel('Human iPSC line');
xticklabels({'GNs','MNs','DNs','hCOs'});
yticks([0:0.1:xlims(end)]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.TickLength = [.01 .01];
b.FontName = 'Arial';
b.FontSize = 25;
%% compute statistics between rules for each cell-type
% set the data
data = e_rules;
% logical for saving
saveTable = 1;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/human_100k_001Hz_energy_rules_all.csv';
% initialise
p = []; d = []; stats = {}; anovatab = {}; compare = {};
% label comparisons
rules = string({'Homophily','Degree','Clustering','Spatial'});
% loop over divs
for type = 1:ntypes;
    % take the data for this div
    datad = squeeze(data(:,type,:));
    % run an anova
    [p(type),anovatab{type},stats{type}] = anova1(datad,{'homophily','degree','clustering','spatial'},'off');
    % run a tukey-kramer
    compare{type} = multcompare(stats{type},'display','off');
    % compute pairwise cohen d
    step = 1;
    for i = 1:4;
        for j = 1:4;
            d(step,type) = computeCohen_d(datad(:,i),datad(:,j),'paired');
            comparison_labels(step,:) = [rules(i) rules(j)];
            step = step + 1;
        end
    end
end
% filter by singular comparisons
i = [2 3 4 7 8 12];
% take
comp_rules = comparison_labels(i,:);
d_rules = -d(i,:); % minus
p_rules = [compare{1}(:,6) compare{2}(:,6) compare{3}(:,6)];
% set significant figures
d_rules = roundsd(d_rules,3);
p_rules = roundsd(p_rules,3);
% form table
t = table(...
    comp_rules(:,1),...
    comp_rules(:,2),...
    p_rules(:,1),...
    d_rules(:,1),...
    p_rules(:,2),...
    d_rules(:,2),...
    p_rules(:,3),...
    d_rules(:,3),...
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('GN p=%.3g',p(1)),'GN Cohen d',...
    sprintf('MN p=%.3g',p(2)),'MN Cohen d',...
    sprintf('DN p=%.3g',p(3)),'DN Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
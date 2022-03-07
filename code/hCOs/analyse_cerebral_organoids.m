%% Analyse generative modeling outcomes from human cerebral organoid networks
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load cerebral organoid data
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_cerebral_organoids_30min/')
% load sttc data
a = load('M03912_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
b = load('M03912_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
organoid_data = {a.gm_data b.gm_data};
% load generative model data
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_cerebral_organoids_30min/M03921_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/org_generative_model.mat');
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
% 3 pipelines
pipeline_labels = string({'10ms'})';
% 12 hd cultures
culture_labels = string(1:6);
%% set hyperparameters of the dataset
% number of pipelines
npipelines = 1;
% number of cultures
ncultures = 6;
% total number of networks
nsamp = npipelines * ncultures;
% number of models
nmodels = 13;
%% form a loop index of the cultures and keep the data
% initialise the indices: columns for pipeline, then culture, then div
index = zeros(nsamp,3); 
all_o_networks = {};
all_d_networks = {};
step = 1;
% loop through and visualise
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        % form the index
        index(step,1) = pipeline;
        index(step,2) = culture;
        % get the data as a cell array
        all_o_networks{step} = organoid_data{pipeline}.table_bu{culture};
        all_d_networks{step} = squareform(pdist(organoid_data{pipeline}.table_xy{culture}));
        % update the step
        step = step + 1;
    end
end 
%% load generative model data
%{
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_cerebral_organoids_30min/M03921_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob');
% initialise
energy_sample = zeros(nsamp,13,20000);
ks_sample = zeros(nsamp,13,20000,4);
networks_sample = cell(nsamp,1);
parameters_sample = zeros(nsamp,13,20000,2);
errors = zeros(nsamp,1);
% loop over pipelines, cultures and divs
step = 1;
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            try % load this network's generative model output
                load(sprintf('rodent_50k_qc_%g_%g_%g_generative_model.mat',pipeline,culture,div));
            catch
                % keep if it doesn't load
                errors(step) = 1;
                % display
                disp(sprintf('network_%g_%g_%g non-existant',pipeline,culture,div));
            end
            % assign
            energy_sample(step,:,:) = output.energy;
            ks_sample(step,:,:,:) = output.ks;
            networks_sample{step} = output.networks;
            parameters_sample(step,:,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('network_%g_%g_%g loaded',pipeline,culture,div));
            % upate step
            step = step + 1;
        end
    end
end
%}
%% replace all data
energy_sample = org_generative_model.energy;
ks_sample = org_generative_model.ks;
parameters_sample = org_generative_model.parameters;
%networks_sample = rodent_50k_generative_models.networks;
%% look at energy landscape for a selected rule and network
% select model and network
model = 3;
net = 6;
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
criteria = index(:,1)==2 & index(:,3)==4;
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
        u = scatter(eta(net,:),e_select(net,:),20,col,'.'); ylabel('energy'); xlabel('eta'); ylim([0 1]);
        b = gca; b.TickDir = 'out';
        b.FontName = 'Arial';
        yticks([]); xticks([]); hold on;
    end
else
% plot the energy landscape
h = figure;
h.Position = [100 100 500 400];
for net = 1:nsamp_select
    u = scatter(squeeze(p_select(net,:,1)),squeeze(p_select(net,:,2)),20,e_select(net,:),'.'); hold on;
    %xlabel('eta'); ylabel('gamma');
end
caxis([0 1]); %c = colorbar; c.Label.String = 'energy';
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
top_e_mean = zeros(length(nset),13,nsamp);
top_p = cell(length(nset),1);
top_p_mean = zeros(length(nset),13,nsamp,2);
% compute the minimum
for no = 1:length(nset);
    % take the actual amount of top performing parameters
    n = nset(no);
    for model = 1:13;
        % take energies for this model
        pipeline_d = squeeze(energy_sample(:,model,:))';
        % rank them for each subject
        [v i] = sort(pipeline_d);
        % take top n energies and their indices
        n_e = v(1:n,:);
        n_i = i(1:n,:);
        % take the corresponding parameters
        u = zeros(nsamp,n,2);
        for s = 1:nsamp;
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
for net = 1:nsamp;
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
    'boxColor','r',...
    'boxAlpha',0.15); 
ylim([0 100]);
b = gca; b.TickDir = 'out';
xticklabels({'degree','clustering','betweenness','edge length'});
ylabel('max(KS)'); yticklabels({'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'});
%% visualise energy for each model
% select which set of nset to view
set = 1;
% set the new model order
i = [2:3 9:13 4:8 1];
% take the top set
e_select = squeeze(top_e_mean(set,i,:));
% place into a new dimension
data = nan(1,nmodels,nsamp);
data(1,:,:) = e_select;
% permute
data = permute(data,[3 1 2]);
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
ylim([0 1]); ylabel('Energy'); xlabel('Human Cerebral Organoids');
xticklabels([]);
yticks([0:0.1:1]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% simply rank
ranked = median(squeeze(data));
[u x] = sort(ranked);
disp(models(i(x)));
%% compute statistics between each model for each time point
% initialise
p = []; d = []; stats = {}; anovatab = {}; compare = {};
% take the data
array = squeeze(data(:,1,:));
% run an anova
[p,anovatab,stats] = anova1(array,models,'off');
% run a tukey-kramer
compare = multcompare(stats,'display','off');
% compute pairwise cohen d
p = []; d = [];
step = 1;
for i = 1:13;
    for j = 1:13;
        d(i,j) = abs(computeCohen_d(array(:,i),array(:,j)));
    end
end
p = compare(:,6);
%% visualise summary statistics by generative rule and take the take the energy values
% take the set
set = 1;
% set the new model order
i = [2:13 1];
% take the top set
e_select = squeeze(top_e_mean(set,i,:));
% permute the model order
e_select = permute(e_select,[2 3 1]);
% averaged
e_mean_rules = nan(nsamp,1,4);
e_mean_rules(:,1,1) = squeeze(mean(e_select(:,:,[1:2]),3));
e_mean_rules(:,1,2) = squeeze(mean(e_select(:,:,[3:7]),3));
e_mean_rules(:,1,3) = squeeze(mean(e_select(:,:,[8:12]),3));
e_mean_rules(:,1,4) = squeeze(e_select(:,:,13));
% all together
e_rules = nan(nsamp*5,1,4);
a = squeeze(e_select(:,:,[1:2])); e_rules(1:2*nsamp,1,1) = a(:);
b = squeeze(e_select(:,:,[3:7])); e_rules(1:5*nsamp,1,2) = b(:);
c = squeeze(e_select(:,:,[8:12])); e_rules(1:5*nsamp,1,3) = c(:);
d = squeeze(e_select(:,:,[13])); e_rules(1:nsamp,1,4) = d(:);
% visualise
h = figure; h.Position = [100 100 700 600];
u = iosr.statistics.boxPlot(e_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
ylim([0 0.85]); ylabel('Energy'); %xlabel('Human Cerebral Organoids');
xticklabels('hCOs');
yticks([0:0.1:0.8]);
b = gca;
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.TickLength = [.01 .01];
b.FontName = 'Arial';
b.FontSize = 25;
%% compute statistics between rules
% set the data
data = e_rules;
% logical for saving
saveTable = 1;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/cerebral_organoids_001Hz_energy_rules_all.csv';
% initialise
p = []; d = []; stats = {}; anovatab = {}; compare = {};
% label comparisons
rules = string({'Homophily','Degree','Clustering','Spatial'});
% take the data
datad = data;
% run an anova
[p,anovatab,stats] = anova1(datad,{'homophily','degree','clustering','spatial'},'off');
% run a tukey-kramer
compare = multcompare(stats,'display','off');
% compute pairwise cohen d
step = 1;
for i = 1:4;
    for j = 1:4
        d(step) = computeCohen_d(datad(:,i),datad(:,j),'paired');
        comparison_labels(step,:) = [rules(i) rules(j)];
        step = step + 1;
    end
end
% filter by singular comparisons
i = [2 3 4 7 8 12];
% take
comp_rules = comparison_labels(i,:);
d_rules = -d(i); % minus
p_rules = compare(:,6);
% set significant figures
d_rules = roundsd(d_rules,3);
p_rules = roundsd(p_rules,3);
% form table
t = table(...
    comp_rules(:,1),...
    comp_rules(:,2),...
    p_rules,...
    d_rules',...
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('p=%.3g',p),'Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
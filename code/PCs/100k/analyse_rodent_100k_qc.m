%% Analyse generative modeling outcomes from 100k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortical neurons
clear; clc;
% change directory
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_100k_tracking/');
% load sttc data
a = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
rodent_data = {a.gm_data b.gm_data c.gm_data};
% load generative model data
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/rodent_100k_001Hz_generative_models');
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
nmodels = 13;
%% form labels
% 3 pipelines
pipeline_labels = string({'5','10','20'})';
% 4 time points
div_labels = string({'14','21','28'});
% 12 hd cultures
culture_labels = string(1:12);
%% set hyperparameters of the dataset
% number of pipelines
npipelines = 1;
% number of cultures
ncultures = 12;
% number of divs
ndivs = 3;
% total number of networks
nsamp = npipelines*ncultures*ndivs;
%% form a loop index of the cultures
% initialise the indices: columns for pipeline, then culture, then div
index = zeros(nsamp,3); 
step = 1;
% loop through and visualise
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;;
            index(step,1) = pipeline;
            index(step,2) = culture;
            index(step,3) = div;
            step = step + 1;
        end
    end
end 
% form an array of all networks
for k = 1:size(index,1);
    all_o_networks{k} = rodent_data{index(k,1)}.table_bu_shared{index(k,2)}{index(k,3)};
    all_c_networks{k} = squareform(pdist(rodent_data{index(k,1)}.table_xy_shared{index(k,2)}{index(k,3)}));
end
%% load 0.1Hz generative model data
%{
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/rodent_100k_qc');
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
                load(sprintf('rodent_qc_%g_%g_%g_generative_model.mat',pipeline,culture,div));
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
% remove networks that are erroneous
if sum(errors) > 0;
    energy_sample(errors,:,:) = [];
    ks_sample(errors,:,:) = [];
    networks_sample(errors) = [];
    parameters_sample(errors,:,:,:) = [];
    index(errors,:) = [];
end
% save data as a struct
rodent_100k_01Hz_generative_models = struct;
rodent_100k_01Hz_generative_models.energy = energy_sample;
rodent_100k_01Hz_generative_models.ks = ks_sample;
rodent_100k_01Hz_generative_models.networks = networks_sample;
rodent_100k_01Hz_generative_models.parameters = parameters_sample;
rodent_100k_01Hz_generative_models.models = models;
rodent_100k_01Hz_generative_models.index = index;
rodent_100k_01Hz_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 12 cultures - index column 3:3 days in vitro (14 days, 21 days, 28 days)'});
rodent_100k_01Hz_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
rodent_100k_01Hz_generative_models.author = string('written by dr danyal akarca, 01/10/21');
save('rodent_100k_01Hz_generative_models.mat','rodent_100k_generative_models','-v7.3');
%}
%% load 0.01Hz generative model data
%{
% nsamp
nsamp = ncultures*ndivs;
% initialise
energy_sample = zeros(nsamp,13,20000);
ks_sample = zeros(nsamp,13,20000,4);
networks_sample = cell(nsamp,1);
parameters_sample = zeros(nsamp,13,20000,2);
errors = zeros(nsamp,1);
% loop over pipelines, cultures and divs
step = 1;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            for model = 1:nmodels;
                % load this network's generative model output
                load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_100k_tracking/GNM_tracking_rec%g/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/rodent_1_%g_%g_generative_model_%g.mat',...
                    div,culture,div,model));
                % assign
                energy_sample(step,model,:) = output.energy;
                ks_sample(step,model,:,:) = output.ks;
                networks_sample{step}(model,:,:) = output.networks;
                parameters_sample(step,model,:,:) = output.parameters;
                % clear the variable
                clear output
                % display
                disp(sprintf('network_1_%g_%g %s loaded',culture,div,models(model)));
            end
        step = step + 1;
        end
    end
% save data as a struct
rodent_100k_001Hz_generative_models = struct;
rodent_100k_001Hz_generative_models.energy = energy_sample;
rodent_100k_001Hz_generative_models.ks = ks_sample;
rodent_100k_001Hz_generative_models.networks = networks_sample;
rodent_100k_001Hz_generative_models.parameters = parameters_sample;
rodent_100k_001Hz_generative_models.models = models;
rodent_100k_001Hz_generative_models.index = index;
rodent_100k_001Hz_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 12 cultures - index column 3:3 days in vitro (14 days, 21 days, 28 days)'});
rodent_100k_001Hz_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
rodent_100k_001Hz_generative_models.author = string('written by dr danyal akarca, 21/01/22');
%}
%% replace all data
energy_sample = rodent_100k_001Hz_generative_models.energy;
ks_sample = rodent_100k_001Hz_generative_models.ks;
parameters_sample = rodent_100k_001Hz_generative_models.parameters;
networks_sample = rodent_100k_001Hz_generative_models.networks;
%% look at energy landscape for a selected rule and network
% select model and network
model = 3;
net = 30;
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
criteria = index(:,1)==1 & index(:,3)==1; % note we now only have a single 10ms pipeline
criteria_ind = find(criteria);
% take the measure
e = squeeze(energy_sample(criteria_ind,model,:));
p = squeeze(parameters_sample(criteria_ind,model,:,:));
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
h.Position = [100 100 550 500];
for net = 1:nsamp_select
    u = scatter(squeeze(p_select(net,:,1)),squeeze(p_select(net,:,2)),20,e_select(net,:),'.'); hold on;
    xlabel('\eta'); ylabel('\gamma');
end
xticks([-10 0 10]); yticks([-10 0 10]);
caxis([0 1]); %c = colorbar; c.Label.String = 'Energy';
u = sgtitle(sprintf('Matching generative model',models(model))); u.FontSize = 25; u.FontName = 'Arial';
b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
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
top_network = {};
% loop
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
            % keep top network for top parameter combination
        if n == 1;
        A = [];
        for s = 1:nsamp;
            % get the best performing networks for this model and subject
            aind = squeeze(networks_sample{no}(model,:,n_i(1)));
            % get number of connections
            m = size(aind,1);
            % reshape into a matrix for this subject and model
            b = zeros(size(all_o_networks{no},1));
            b(aind) = 1;
            b = b + b';
            A(s,:,:) = b;
        end
        A = squeeze(A);
        top_n = A;
            % squeeze the matrices
            u = squeeze(u);
            % a
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
% specify the model/s
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
    'boxColor','c',...
    'boxAlpha',0.15); 
ylim([0 100]);
b = gca; b.TickDir = 'out';
xticklabels({'degree','clustering','betweenness','edge length'});
ylabel('max(KS)'); yticklabels({'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'});
%% get statistics of top simulated networks
% set model to consider
model = 3;
% set network criteria
criteria = index(:,1)==1;
% get the top n network
n = 1;      
% take energy
pipeline_d = squeeze(energy_sample(criteria,model,:))';
% rank them for each network
[v i] = sort(pipeline_d);
% take top n energies and their indices
n_i = i(1,:);
% find indices
aind = find(criteria);
% take the corresponding parameters
A = [];
for s = 1:length(aind);
    % take the network
    uu = aind(s);
    % get the best performing networks for this model and subject
    u = squeeze(networks_sample{uu}(model,:,n_i(1)));
    % reshape into a matrix for this subject and model
    b = zeros(size(all_o_networks{uu},1));
    b(u) = 1;
    b = b + b';
    % take top network
    top_networks{s}=b;
end
% take the distribuitions of each simulation
simulated_statistics = {};
empirical_statistics = {};
for net = 1:length(aind);
    % take observed
    observed = all_o_networks{aind(net)};
    d = all_c_networks{aind(net)};
    % take simulated
    simulated = top_networks{net};
    % take n
    nnode = size(simulated,1);
    % initialise
    empirical_statistics{net}=zeros(4,nnode);
    simulated_statistics{net}=zeros(4,nnode);
    % compute observed
    empirical_statistics{net}(1,:) = degrees_und(observed);
    empirical_statistics{net}(2,:) = clustering_coef_bu(observed)';
    empirical_statistics{net}(3,:) = betweenness_bin(observed);
    empirical_statistics{net}(4,:) = sum(triu(observed.*d,1));
    % computesimulated
    simulated_statistics{net}(1,:) = degrees_und(simulated);
    simulated_statistics{net}(2,:) = clustering_coef_bu(simulated)';
    simulated_statistics{net}(3,:) = betweenness_bin(simulated);
    simulated_statistics{net}(4,:) = sum(triu(simulated.*d,1));
end
%% visualise summary statistics across a specific group
% select which set of nset to view
set = 1;
% select which pipelines to view over time
pview = 1;
% plot over pipelines
e_select = [];
for pipeline = 1:npipelines;
    criteria = index(:,1)==pipeline;
    % take specific energy values
    e_select(:,:,pipeline) = squeeze(top_e_mean(set,:,criteria));
end
% set the new model order
i = [2:3 9:13 4:8 1];
e_select = e_select(i,:,:);
% permute the order
e_select = permute(e_select,[2 3 1]);
% iosr boxplot
h = figure;
h.Position = [100 100 1200 600];
iosr.statistics.boxPlot(e_select,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylim([0 1]); ylabel('Energy'); xlabel('STTC \Deltat');
xticklabels(pipeline_labels);
yticks([0:0.1:1]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% plot over divs * note, set any other setting here
e_select = [];
for div = 1:ndivs;
criteria = index(:,1)==pview&index(:,3)==div;
% take specific energy values
e_select(:,:,div) = squeeze(top_e_mean(set,:,criteria));
end
% set the new order
i = [2:3 8:12 4:8 1];
e_select = e_select(i,:,:);
% permute the model order
e_select = permute(e_select,[2 3 1]);
% iosr boxplot
h = figure;
h.Position = [100 100 1200 600];
iosr.statistics.boxPlot(e_select,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylim([0 1]); ylabel('Energy'); xlabel('Days {\itin vitro} (DIV)');
xticklabels(div_labels);
yticks([0:0.1:1]);
ylim([0 1]); % added to scale
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute energy over time for specific models or all
% model
models = [1:13];
e_pick = e_select(:,:,models);
% squeeze over all rules
e_time = permute(e_pick,[2 3 1]);
e_time = e_time(:,:)';
% anova
[p anovatab stats] = anova1(e_time);
comp = multcompare(stats);
%% group by generative rule and look at energy values over time
% set pipeline to view
pview = 1;
% plot over divs
e_select = [];
for div = 1:ndivs;
    % set the crieria
    criteria = index(:,1)==pview & index(:,3)==div;
    % take specific energy values
    e_select(:,:,div) = squeeze(top_e_mean(set,:,criteria));
end
% set the new order
i = [2:13 1];
e_select = e_select(i,:,:);
% permute the model order
e_select = permute(e_select,[2 3 1]);
% mean
e_div_mean_rules = [];
e_div_mean_rules(:,:,1) = squeeze(mean(e_select(:,:,[1:2]),3));
e_div_mean_rules(:,:,2) = squeeze(mean(e_select(:,:,[3:7]),3));
e_div_mean_rules(:,:,3) = squeeze(mean(e_select(:,:,[8:12]),3));
e_div_mean_rules(:,:,4) = squeeze(e_select(:,:,13));
% all
e_div_rules = nan(12*5,ndivs,4);
a = squeeze(e_select(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); e_div_rules(1:ncultures*2,:,1) = a';
b = squeeze(e_select(:,:,[8:12])); b = permute(b,[2 1 3]); b = b(:,:); e_div_rules(1:ncultures*5,:,2) = b';
c = squeeze(e_select(:,:,[3:7])); c = permute(c,[2 1 3]); c = c(:,:); e_div_rules(1:ncultures*5,:,3) = c';
d = squeeze(e_select(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); e_div_rules(1:ncultures,:,4) = d';
% visualise
h = figure;
h.Position = [100 100 1200 600];
u = iosr.statistics.boxPlot(e_div_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]}';
ylim([0 1]); ylabel('Energy'); xlabel('Days {\itin vitro} (DIV)');
yticks([0:0.1:1]);
xticklabels({'14','21','28'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute statistics between rules for each time point
% set the data
data = e_div_rules;
% logical for saving
saveTable = 1;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/rodent_100k_001Hz_energy_rules_all.csv';
% initialise
p = []; d = []; stats = {}; anovatab = {}; compare = {};
% label comparisons
rules = string({'Homophily','Degree','Clustering','Spatial'});
% loop over divs
for div = 1:ndivs;
    % take the data for this div
    datad = squeeze(data(:,div,:));
    % run an anova
    [p(div),anovatab{div},stats{div}] = anova1(datad,{'homophily','degree','clustering','spatial'},'off');
    % run a tukey-kramer
    compare{div} = multcompare(stats{div},'display','off');
    % compute pairwise cohen d
    step = 1;
    for i = 1:4;
        for j = 1:4;
            d(step,div) = computeCohen_d(datad(:,i),datad(:,j),'paired');
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
    sprintf('DIV14 p=%.3g',p(1)),'DIV14 Cohen d',...
    sprintf('DIV21 p=%.3g',p(2)),'DIV21 Cohen d',...
    sprintf('DIV28 p=%.3g',p(3)),'DIV28 Cohen d'})
if saveTable==1;
    writetable(t,tablename);
end
%% take local measures and compute the correlation matrix but loop over time points and models
% set pipeline 
pipeline = 1;
% set how many measures are computed
nmeasures = 6;
% set how many models are computed
nmodels = 13;
% set up to how many divs are considered 
ndivs = 3;
% initialise correlation matrices
corr_local_observed = nan(ndivs,12,nmeasures,nmeasures);
corr_local_simulated = nan(ndivs,12,nmodels,nmeasures,nmeasures);
corr_local_together = nan(ndivs,12,nmodels,nmeasures*2,nmeasures*2);
% initalise the topological organization dissimilairty
todissim = nan(ndivs,12,nmodels);
% labels
var_labels = string({...
    'obs degree','obs clustering','obs betweenn','obs length','obs eff','obs match',...
    'sim degree','sim clustering','sim between','sim length','sim eff','sim match'});
% loop over divs
for div = 1:ndivs
    % index the networks
    index_loop = find(index(:,1)==pipeline&index(:,3)==div);
    % compute the sample size
    nsamp = size(index_loop,1);
    % take the observed networks
    o_networks = {};
    c_networks = {};
    for culture = 1:ncultures;
        o_networks{culture} = rodent_data{pipeline}.table_bu_shared{culture}{div};
        c_networks{culture} = squareform(pdist(rodent_data{pipeline}.table_xy_shared{culture}{div}));
    end
    % take the simulation data for these networks
    E = energy_sample(index_loop,:,:);
    Ps = parameters_sample(index_loop,:,:,:);
    N = networks_sample(index_loop);    
    % compute observed local statistics for each network
    for i = 1:nsamp;
        % take the observed network
        w = o_networks{i};
        % take the cost
        d = c_networks{i};
        % compute the number of nodes
        nnode = size(w,1);
        % initalise array
        observed = zeros(nnode,nmeasures);
        % compute local observed statistics
        observed(:,1) = degrees_und(w)';
        observed(:,2) = clustering_coef_bu(w);
        observed(:,3) = betweenness_bin(w)';
        observed(:,4) = sum(w.*d)';
        observed(:,5) = efficiency_bin(w,1);
        observed(:,6) = mean(matching_ind(w)+matching_ind(w)')';
        % keep
        corr_local_observed(div,i,:,:) = corr(observed);
        for model = 1:nmodels
            % loop over models
            disp(sprintf('evaluating day %s, network %g, generative model %g...',div_labels(div),i,model)); 
            % take the optimal network simulation parameters, given the set model, for this network
            e_net = squeeze(E(i,model,:));
            [~,i_low] = min(e_net);
            % take the end simulated network
            si = squeeze(N{i}(model,:,i_low));
            s = zeros(nnode);
            s(si) = 1;
            s = s + s';
            % initalise the array
            simulated = zeros(nnode,nmeasures);
            % compute local simulated statistics
            simulated(:,1) = degrees_und(s)';
            simulated(:,2) = clustering_coef_bu(s);
            simulated(:,3) = betweenness_bin(s)';
            simulated(:,4) = sum(s.*d)';
            simulated(:,5) = efficiency_bin(s,1);
            simulated(:,6) = mean(matching_ind(s)+matching_ind(s)')';
            % compute the correlation of the simulated
            corr_local_simulated(div,i,model,:,:) = corr(simulated);
            % form a matrix together and correlate these
            corr_local_together(div,i,model,:,:) = corr([observed simulated]);
            % compute the topological organization dissimilarity
            % todissim(div,i,model) = topological_organization_dissimilarity(corr(observed),corr(simulated));
            todissim(div,i,model) = norm(corr(observed)-corr(simulated));
        end
    end
end
%% visualise the tfd
%  alter the order
todissim_permute = permute(todissim,[2 1 3]);
% make zeros into nans
todissim_permute(todissim_permute==0)=NaN;
% model order
i = [2:3 9:13 4:8 1];
% plot
h = figure; h.Position = [100 100 1200 600];
h = iosr.statistics.boxPlot(todissim_permute(:,:,i),...
    'showViolin',logical(0),...
    'theme','colorall',...
    'showScatter',logical(0),...
    'scatterMarker','x',...
    'scatterColor',[.5 .5 .5],...
    'symbolColor',[.5 .5 .5],...
    'symbolMarker','x',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylabel('{\itTFdissimilarity}'); xlabel('Days {\itin vitro} (DIV)');
xticklabels({'14','21','28'});
%ylim([0 0.25]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% group todissimilaity into rules
to_div_mean_rules = [];
to_div_mean_rules(:,:,1) = squeeze(mean(todissim_permute(:,:,[2:3]),3));
to_div_mean_rules(:,:,2) = squeeze(mean(todissim_permute(:,:,[4:8]),3));
to_div_mean_rules(:,:,3) = squeeze(mean(todissim_permute(:,:,[9:13]),3));
to_div_mean_rules(:,:,4) = squeeze(todissim_permute(:,:,1));
% all
to_div_rules = nan(ncultures*5,ndivs,4);
a = squeeze(todissim_permute(:,:,[2:3])); a = permute(a,[2 1 3]); a = a(:,:); to_div_rules(1:ncultures*2,:,1) = a';
b = squeeze(todissim_permute(:,:,[9:13])); b = permute(b,[2 1 3]); b = b(:,:); to_div_rules(1:ncultures*5,:,2) = b';
c = squeeze(todissim_permute(:,:,[4:8])); c = permute(c,[2 1 3]); c = c(:,:); to_div_rules(1:ncultures*5,:,3) = c';
d = squeeze(todissim_permute(:,:,[1])); d = permute(d,[2 1 3]); d = d(:,:); to_div_rules(1:ncultures,:,4) = d';
%% compute statistics between rules for each time point
% set the data
data = to_div_rules;
% logical for saving
saveTable = 0;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/rodent_100k_001Hz_tfd_rules_all.csv';
% initialise
p = []; d = []; stats = {}; anovatab = {}; compare = {};
% label comparisons
rules = string({'Homophily','Degree','Clustering','Spatial'});
% loop over divs
for div = 1:ndivs;
    % take the data for this div
    datad = squeeze(data(:,div,:));
    % run an anova
    [p(div),anovatab{div},stats{div}] = anova1(datad,{'homophily','degree','clustering','spatial'},'off');
    % run a tukey-kramer
    compare{div} = multcompare(stats{div},'display','off');
    % compute pairwise cohen d
    step = 1;
    for i = 1:4;
        for j = 1:4;
            d(step,div) = computeCohen_d(datad(:,i),datad(:,j),'paired');
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
    sprintf('DIV14 p=%.3g',p(1)),'DIV14 Cohen d',...
    sprintf('DIV21 p=%.3g',p(2)),'DIV21 Cohen d',...
    sprintf('DIV28 p=%.3g',p(3)),'DIV28 Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
%% visualise a specific comparison in terms of a culture and averaged over cultures
div = 1;
culture = 5;
model = 1;
% take the specific data
specific_obs = squeeze(corr_local_observed(div,culture,:,:));
specific_sim = squeeze(corr_local_simulated(div,culture,model,:,:));
% visualise the specific culture
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); imagesc(specific_obs); title('observed'); caxis([-1 1]); colorbar; xticks([]); yticks([]);
subplot(1,2,2); imagesc(specific_sim); title('simulated'); caxis([-1 1]); colorbar; xticks([]); yticks([]);
sgtitle(sprintf('culture %g, div %g, model %g, tod %g',culture,div,model,squeeze(squeeze(todissim(div,culture,model)))));
% take the averaged data
averaged_obs = squeeze(mean(corr_local_observed(div,:,:,:),2));
averaged_sim = squeeze(mean(corr_local_simulated(div,:,model,:,:),2));
% visualise the averaged data
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); imagesc(averaged_obs); title('observed'); caxis([-1 1]); colorbar; xticks([]); yticks([]);
subplot(1,2,2); imagesc(averaged_sim); title('simulated'); caxis([-1 1]); colorbar; xticks([]); yticks([]);
sgtitle(sprintf('averaged div %g, model %g, tod %g',div,model,squeeze(mean(todissim(div,:,model),2))));
xticks([]); yticks([]);
%% plot a specific plot
g = figure; g.Position = [100 100 900 500];
plotind = [1 6 11; 2 7 12; 3 8 13; 4 9 14; 5 10 15];
for div = 1:ndivs
    % observed
    subplot(3,5,plotind(1,div)); imagesc(squeeze(mean(corr_local_observed(div,:,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % matching
    subplot(3,5,plotind(2,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,3,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % clu-avg
    subplot(3,5,plotind(3,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,4,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % deg-avg
    subplot(3,5,plotind(4,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,9,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % sptl
    subplot(3,5,plotind(5,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,1,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
end
%% correlate energy and the to matrix
% pick the set
set = 1; 
% take the energy
ep = squeeze(top_e_mean(set,:,:))';
% take only those we have calculated the todissim of
epto = zeros(3,12,13);
% form into the same size to form correlations
for div = 1:3;
    epto(div,:,:) = ep(index(:,1)==pipeline&index(:,3)==div,:);
end
% form indices for rules colours
ind = [ones(3*12*1,1);2*ones(3*12*2,1);3*ones(3*12*5,1);4*ones(3*12*5,1)];
col = {'k','r','g','b'};
% across all models
[rall pall] = corr(epto(:),todissim(:));
% take data
x = epto(:); y = todissim(:);
% keep individual correlations based on rules
rrule  = []; prule = [];
% visualise
figure;
for rule = 1:4;
    [rrule(rule) prule(rule)] = corr(x(ind==rule),y(ind==rule));
    h = scatter(x(ind==rule),y(ind==rule),70,'o',...
        'markerfacecolor',col{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
end
sgtitle('all models'); xlabel('energy'); ylabel('TFdissimilarity');
xlim([0 0.6]);
%ylim([0 0.3]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial';
% visualise scatter hist
u = figure; u.Position = [100 100 900 600];
scatterhist(x,y,'group',ind,...
    'marker','.','markersize',15,...
    'color',{'k','r','g','b'},'linestyle','-',...
    'kernel','on','direction','out');
b = gca; b.FontName = 'Arial'; b.FontSize = 14; b.TickDir = 'out';
legend({'Spatial','Homophily','Clustering','Degree'},'box','off');
xlabel('Energy'); ylabel('TFdissimilarity'); box off;
xlim([0 0.65]); b = gca; b.TickDir = 'out'; b.FontSize = 25;
% plot the average at each time point
eptd = []; totd = [];
for div = 1:3;
    % energy
    eptd(div,:,1) = squeeze(epto(div,:,1));
    eptd(div,:,2) = mean(squeeze(epto(div,:,[2 3])),2);
    eptd(div,:,3) = mean(squeeze(epto(div,:,[4:8])),2);
    eptd(div,:,4) = mean(squeeze(epto(div,:,[9:13])),2);
    % todissim
    totd(div,:,1) = squeeze(todissim(div,:,1));
    totd(div,:,2) = mean(squeeze(todissim(div,:,[2 3])),2);
    totd(div,:,3) = mean(squeeze(todissim(div,:,[4:8])),2);
    totd(div,:,4) = mean(squeeze(todissim(div,:,[9:13])),2);
end
eptd = squeeze(mean(eptd,2)); totd = squeeze(mean(totd,2));
% visualise
figure;
shape = {'o','square','pentagram'};
for rule = 1:4;
    for div = 1:3;
    scatter(eptd(div,rule),totd(div,rule),180,...
        'marker',shape{div},...
        'markerfacecolor',col{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
    end
end
xlim([0 0.6]); ylim([0 0.1]);
xlabel('Energy'); ylabel('TFdissimilarity');
b = gca; b.FontName = 'Arial'; b.FontSize = 14; b.TickDir = 'out';
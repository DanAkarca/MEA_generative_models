%% Analyse generative modeling outcomes from 50k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortical neurons
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/');
% load sttc data
a = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
rodent_data = {a.gm_data b.gm_data c.gm_data};
% load generative model data
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/rodent_50k_001Hz_generative_models.mat');
% addpaths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master');
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
addpath('/imaging/astle/users/da04/PhD/toolboxes/bluewhitered');
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
addpath('/imaging/astle/users/da04/PhD/toolboxes/voronoi');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MIT_code_2011')
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/functions/');
% define model types
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
eqn = string({'KS_k','KS_c','KS_b','KS_d'});
nmodels = length(models);
%% form labels
% 3 pipelines
pipeline_labels = string({'5','10','20'})';
% 4 time points
div_labels = string({'7','10','12','14'});
% 12 hd cultures
culture_labels = string(1:12);
%% set hyperparameters of the dataset
% number of pipelines
npipelines = 3;
% number of cultures
ncultures = 6;
% number of divs
ndivs = 4;
% total number of networks
nsamp = npipelines * ncultures * ndivs;
%% form a loop index of the cultures and keep the data
% initialise the indices: columns for pipeline, then culture, then div
index = zeros(nsamp,3); 
all_o_networks = {};
all_d_networks = {};
step = 1;
% loop through and visualise
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            % form the index
            index(step,1) = pipeline;
            index(step,2) = culture;
            index(step,3) = div;
            % get the data as a cell array
            all_o_networks{step} = rodent_data{pipeline}.table_bu_shared{culture}{div};
            all_d_networks{step} = squareform(pdist(rodent_data{pipeline}.table_xy_shared{culture}{div}));
            % update the step
            step = step + 1;
        end
    end
end 
%% load 0.1Hz generative model data
%{
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/rodent_50k_qc');
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
% remove networks that are erroneous
if sum(errors) > 0;
    energy_sample(errors,:,:) = [];
    ks_sample(errors,:,:) = [];
    networks_sample(errors) = [];
    parameters_sample(errors,:,:,:) = [];
    index(errors,:) = [];
end
%% save data as a struct
rodent_50k_generative_models = struct;
rodent_50k_generative_models.energy = energy_sample;
rodent_50k_generative_models.ks = ks_sample;
rodent_50k_generative_models.networks = networks_sample;
rodent_50k_generative_models.parameters = parameters_sample;
rodent_50k_generative_models.models = models;
rodent_50k_generative_models.index = index;
rodent_50k_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 6 cultures - index column 3:4 days in vitro (7 days, 10 days, 12 days, 14 days)'});
rodent_50k_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
rodent_50k_generative_models.author = string('written by dr danyal akarca, 01/10/21');
%}
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
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            for model = 1:nmodels;
                % form the string
                str = sprintf(...
                    '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/GNM_tracking_rec%g/min_rate_0.01hz_alpha_0.001_lag%s_jitter%s_prob/rodent_1_%g_%g_generative_model_%g.mat',...
                    div,pipeline_labels(pipeline),pipeline_labels(pipeline),culture,div,model);
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
                disp(sprintf('network_%g_%g_%g model %s loaded',pipeline,culture,div,models(model)));
            end
            % update step 
            step = step + 1;
        end
    end
end
% keep in struct
rodent_50k_001Hz_generative_models = struct;
rodent_50k_001Hz_generative_models.energy = energy_sample;
rodent_50k_001Hz_generative_models.ks = ks_sample;
rodent_50k_001Hz_generative_models.networks = networks_sample;
rodent_50k_001Hz_generative_models.parameters = parameters_sample;
rodent_50k_001Hz_generative_models.models = models;
rodent_50k_001Hz_generative_models.index = index;
rodent_50k_001Hz_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 6 cultures - index column 3:4 days in vitro (7 days, 10 days, 12 days, 14 days)'});
rodent_50k_001Hz_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
rodent_50k_001Hz_generative_models.author = string('written by dr danyal akarca, 09/01/22');
%}
%% replace all data
energy_sample = rodent_50k_001Hz_generative_models.energy;
ks_sample = rodent_50k_001Hz_generative_models.ks;
parameters_sample = rodent_50k_001Hz_generative_models.parameters;
networks_sample = rodent_50k_001Hz_generative_models.networks;
%% look at energy landscape for a selected rule and network
% select model and network
model = 3;
net = 36;
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
        u = scatter(eta(net,:),e_select(net,:),20,col,'.'); ylabel('energy'); xlabel('eta'); ylim([0 1]);
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
nset = [1 10 50 100];
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
% specify further criteria
criteria = index(:,3)==4;
ind = find(criteria);
% take the data
ks_data = squeeze(ks_sample(criteria,model,:,:));
% new size 
nnsamp = sum(criteria);
% initialise
driver = [];
% loop over cultures
for net = 1:nnsamp;
    % find the max ks statistics for this network
    [v i] = max(squeeze(ks_data(net,:,:))');
    % group data
    driver(net,:) = [sum(i==1),sum(i==2),sum(i==3),sum(i==4)];,
end
% form a percentage
driver = driver ./ 20000 * 100;
% visualise
h = figure; h.Position = [100 100 1000 400];
iosr.statistics.boxPlot(driver,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'showScatter',logical(1),...
    'scatterColor',[.5 .5 .5],...
    'scatterAlpha',0.5,...
    'symbolColor',[.5 .5 .5],...
    'boxColor','k',...
    'boxAlpha',0.15); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xticklabels({'degree','clustering','betweenness','edge length'});
ylabel('max(KS)'); yticklabels({'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'});
%% visualise summary statistics for specifically defined groups 
% select which set of nset to view
set = 1;
% select which pipeline to view over time
pview = 2;
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
h.Position = [100 100 600 300];
iosr.statistics.boxPlot(e_select,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylim([0 1]); ylabel('Energy'); xlabel('\Deltat');
xticklabels(pipeline_labels);
yticks([0:0.1:1]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 14;
% plot over divs
e_select = [];
for div = 1:ndivs;
criteria = index(:,1)==pview&index(:,3)==div;
% take specific energy values
e_select(:,:,div) = squeeze(top_e_mean(set,:,criteria));
end
% set the new order
i = [2:3 9:13 4:8 1];
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
ylabel('Energy'); xlabel('Days {\itin vitro}');
xticklabels(div_labels);
yticks([0:0.1:1]);
ylim([0 1]); % added to scale
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% model fits from bootstrapping
% set the div
div = 4;
% set the pipeline
pipeline = 2;
% set number of stats computed
nstat = 6;
% get the criteria 
criteria = index(:,1)==pipeline&index(:,3)==div;
% number selected
nsel = sum(criteria);
% get indices
inda = find(criteria);
% set the model and title
model = 9;
title_place = 'Degree Average';
% set the number to permute
nperm = 99;
% get the energies of these cultures
e_net = squeeze(energy_sample(criteria,model,:));
% get their indices
[~,i_low] = sort(e_net','ascend');
% get only the top nperms
ind = i_low(1:nperm,:);
% Determine if there is a slgiht mismatch in number of edges and node
num_edges = []; 
num_nodes = [];
for culture = 1:nsel;
    % get culture;
    A = all_o_networks{inda(culture)};
    % get the giant component
    Ag = giant_component(A);
    % compute number of edges
    num_edges(culture,1) = nnz(Ag)/2;
    % get a simulation
    indi = ind(1,culture);
    % find this network
    net = squeeze(networks_sample{inda(culture)}{model}(:,indi));
    % reconstruct the simulation network
    nnode = size(all_o_networks{inda(culture)},1);
    Asynth = zeros(nnode,nnode);
    Asynth(net) = 1;
    Asynth = Asynth+Asynth';
    % compute number of edges
    num_edges(culture,2) = nnz(Asynth)/2;
    % compute number of edges in the actual observed
    num_edges(culture,3) = nnz(A)/2;
    % do the same but for nnodes
    num_nodes(culture,1) = size(Ag,1); 
    num_nodes(culture,2) = size(Asynth,1);
    num_nodes(culture,3) = size(A,1);
end

% the generative models have been made of lenght of the giant component,
% but compared to the original network.
% this explains why the energy is so high, so early on.
% thus, we have the same number of edges
% but different number of nodes

% initialise permutations
observed_statistics = {}; perm_statistics = {}; test = {};
% take the simulated networks for each culture
for culture = 1:nsel;
    % get the observed networks
    A = all_o_networks{inda(culture)};
    D = all_d_networks{inda(culture)};
    % get the giant component
    %[Ag di] = giant_component(A);
    % compute the giant costs
    %Dg = D(di,di);
    % get the statistics
    observed_statistics{culture}{1} = degrees_und(A);
    observed_statistics{culture}{2} = clustering_coef_bu(A)';
    observed_statistics{culture}{3} = betweenness_bin(A);
    %observed_statistics{culture}{4} = D(triu(A,1)>0)';
    observed_statistics{culture}{4} = sum(D.*triu(A,1))./10^3; %for scale
    % outside of energy
    observed_statistics{culture}{5} = efficiency_bin(A,1)';
    [ci,q] = community_louvain(A);
    observed_statistics{culture}{6} = participation_coef(A,ci)';
    % test
    test{culture}{1} = matching_ind(A);
    test{culture}{2} = degrees_und(A);
    % loop over top nperm networks
    for i = 1:nperm;
        % index the network
        indi = ind(i,culture);
        % find this network
        net = squeeze(networks_sample{inda(culture)}{model}(:,indi));
        % reconstruct the simulation network
        nnode = size(A,1);
        Asynth = zeros(nnode,nnode);
        Asynth(net) = 1;
        Asynth = Asynth+Asynth';
        % get the number of edges
        num_edges_sim(culture) = nnz(Asynth)/2;
        % get the statistics
        perm_statistics{culture}{1}(i,:) = degrees_und(Asynth)';
        perm_statistics{culture}{2}(i,:) = clustering_coef_bu(Asynth);
        perm_statistics{culture}{3}(i,:) = betweenness_bin(Asynth)';
        %perm_statistics{culture}{4}(i,:) = D(triu(Asynth,1)>0);
        perm_statistics{culture}{4}(i,:) = sum(D.*triu(Asynth,1))./10^3;  %for scale;
        % outside of energy
        perm_statistics{culture}{5}(i,:) = efficiency_bin(Asynth,1);
        [ci,q] = community_louvain(Asynth);
        perm_statistics{culture}{6}(i,:) = participation_coef(Asynth,ci);
        % display
        disp(sprintf('culture %g permutataion %g of %g complete',culture,i,nperm));
    end
end
%% degree average and matching seem to have no relationship in observed networks
figure;
for culture = 1:nsel
    % get the unparametrised homophily and deg-avg
    n = size(test{culture}{1});
    ind = find(triu(ones(n),1));
    match = test{culture}{1}(ind);
    nc = nchoosek(1:n,2);
    degavv = (test{culture}{2}(nc(:,1)) + test{culture}{2}(nc(:,2)))/2;
    degavm = zeros(n); degavm(ind) = degavv;
    degavg = degavm(ind);
    % correlate them
    [r(culture) p(culture)] = corr(match,degavg);
    % plot
    scatter(match,degavg,'.'); hold on;
end
ylabel('Degree average'); xlabel('Matching');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
%% compute ks statistics for culture and measures
% initialise
all_perms = 1:nperm; obs_ks = []; boot_ks = []; p_rank = []; close = [];
% initialise
observed_cdf_all = {}; simulated_cdf_all = {}; mvalue = [];
% loop over
for culture = 1:nsel;
    for stat = 1:nstat;
        % take staitstics
        null = perm_statistics{culture}{stat};
        obs = observed_statistics{culture}{stat}';
        % keep the maximum values
        mvalue(culture,stat) = max(max(max(null)),max(max(obs)));
        % compute the observed ks with the mean null
        obs_ks(culture,stat) = fcn_ks(mean(null,1)',obs);
        for i = 1:nperm;
            boot = all_perms;
            boot(i) = [];
            nnull = [null(boot,:)' obs];
            perm = null(i,:)';
            % compute the observed ks with the mean null
            boot_ks(culture,stat,i) = fcn_ks(mean(nnull,2),perm);
        end
        % compute how the rank of how close it is to the median of the distribution, and double it for a two-tailed test
        a = obs_ks(culture,stat);
        b = squeeze(boot_ks(culture,stat,:));
        % compute how close obs is to the mean
        obsclose = abs(mean(b)-a);
        % compute how close each is to the mean and rank
        for j = 1:length(b);
            close(j) = abs(mean(b)-b(j));
        end
        % now we do a one tailed test on this distribution
        p_rank(culture,stat) = sum(obsclose<close)/(nperm+1);
        % display
        disp(sprintf('Culture %g, Statistic %g model fits evaluated.',culture,stat));
    end
end
% round
p_rank = roundsd(p_rank,3);
% compute a median p
median_p_rank = median(p_rank);
% number under alpha
alpha = 0.05;
countunder = sum(p_rank<alpha);
% form the table
t = table([1:nsel]',...
    p_rank(:,1),...
    p_rank(:,2),...
    p_rank(:,3),...
    p_rank(:,4),...
    p_rank(:,5),...
    p_rank(:,6),...
    'VariableNames',{'Culture',...
    sprintf('Degree (%g/6 p<%g)',countunder(1),alpha),...
    sprintf('Clustering (%g/6 p<%g)',countunder(2),alpha),...
    sprintf('Betweenness (%g/6 p<%g)',countunder(3),alpha),...
    sprintf('Edge length (%g/6 p<%g)',countunder(4),alpha),...
    sprintf('Efficiency (%g/6 p<%g)',countunder(5),alpha),...
    sprintf('Participation (%g/6 p<%g)',countunder(6),alpha)});
disp(t);
%% visualise bootstraps for a single culture
% set the labels
xlabels = {'Degree, k','Clustering, c','Betweenness, b','Edge length, d','Efficiency, e','Participation, p'};
ylabels = {'F(k)','F(c)','F(b)','F(d)','F(e)','F(p)'};
% colormap 
% colormaps = [173,216,230; 104 179 179; 255 179 179; 241 208 118; 192 192 192; 255 228 181]./256; % old
colormaps = [83, 62, 133; 72 143 177; 79 211 196; 193 248 207; 192 192 192; 255 228 181]./256; %new
% set the culture
culture = 2;
% initialise figure
h = figure; h.Position = [100 100 1500 600];
% set placing inside and outside the equation
place = [1 2 3 4 6 7];
for stat = 1:nstat;
    % take place on subplot
    put = place(stat);
    % customise placing
    subplot(2,4,put);
    % plot the distribution of cdfs
    for perm = 1:nperm;
        u = cdfplot(perm_statistics{culture}{stat}(perm,:));
        u.Color = colormaps(stat,:); u.LineWidth=1; hold on;
    end
    % plot the observed culture on top
    v = cdfplot(observed_statistics{culture}{stat});
    v.Color=[.25 .25 .25]; v.LineWidth=3; hold on;
    % labels and settings
    xlabel(xlabels{stat}); ylabel(ylabels{stat});
    grid off; b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
    % place p values 
    title(sprintf('Median p_{rank}=%.3g',median_p_rank(stat)),...
        'FontName','Helvetica','FontSize',25);
    box off;
end
% settings
sgtitle(title_place,'FontSize',25);
%% visualise only energy equation bookstraps
% initialise figure
h = figure; h.Position = [100 100 2200 350];
% set placing inside and outside the equation
for stat = 1:4;
    % customise placing
    subplot(1,4,stat);
    % plot the distribution of cdfs
    for perm = 1:nperm;
        u = cdfplot(perm_statistics{culture}{stat}(perm,:));
        u.Color = colormaps(stat,:); u.LineWidth=1; hold on;
    end
    % plot the observed culture on top
    v = cdfplot(observed_statistics{culture}{stat});
    v.Color=[.25 .25 .25]; v.LineWidth=3; hold on;
    % labels and settings
    xlabel(xlabels{stat}); ylabel(ylabels{stat});
    grid off; b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
    box off;
    % no title
    title([]);
end
%% visualise summary statistics for just the best rules
% plot over divs
e_select = [];
for div = 1:ndivs;
criteria = index(:,1)==pview&index(:,3)==div;
% take specific energy values
e_select(:,:,div) = squeeze(top_e_mean(set,:,criteria));
end
% set the new order
i = [2:13 1];
e_select = e_select(i,:,:);
% permute the model order
e_select = permute(e_select,[2 3 1]);
% iosr boxplot
h = figure;
h.Position = [100 100 600 300];
iosr.statistics.boxPlot(e_select(:,:,[2 8 3 13]),...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
ylabel('Energy'); xlabel('Days {\itin vitro} (DIV)');
sgtitle('Best rules');
xticklabels(div_labels);
yticks([0:0.1:1]);
ylim([0 1]); % added to scale
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 14;
%% visualise summary statistics by generative rule and take the take the energy values over time
% set the pipeline to view
pview = 2;
% plot over divs
e_select = [];
for div = 1:ndivs;
    % set the crieria
    criteria = index(:,1)==pview&index(:,3)==div;
    % take specific energy values
    e_select(:,:,div) = squeeze(top_e_mean(set,:,criteria));
end
% set the new order
i = [2:13 1];
e_select = e_select(i,:,:);
% permute the model order
e_select = permute(e_select,[2 3 1]);
% visualise mean summary energy across divs across rule types: note, they've already been reordered above
e_div_mean_rules = [];
e_div_mean_rules(:,:,1) = squeeze(mean(e_select(:,:,[1:2]),3));
e_div_mean_rules(:,:,2) = squeeze(mean(e_select(:,:,[8:12]),3));
e_div_mean_rules(:,:,3) = squeeze(mean(e_select(:,:,[3:7]),3));
e_div_mean_rules(:,:,4) = squeeze(e_select(:,:,13));
% same but all together
to_div_rules = nan(ncultures*5,ndivs,4);
a = squeeze(e_select(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); to_div_rules(1:ncultures*2,:,1) = a';
b = squeeze(e_select(:,:,[8:12])); b = permute(b,[2 1 3]); b = b(:,:); to_div_rules(1:ncultures*5,:,2) = b';
c = squeeze(e_select(:,:,[3:7])); c = permute(c,[2 1 3]); c = c(:,:); to_div_rules(1:ncultures*5,:,3) = c';
d = squeeze(e_select(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); to_div_rules(1:ncultures,:,4) = d';
% visualise
h = figure; h.Position = [100 100 1200 600];
u = iosr.statistics.boxPlot(to_div_rules,...
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
xticklabels({'7','10','12','14'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute statistics between rules for each time point
% set the data
data = to_div_rules;
% logical for saving
saveTable = 0;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/rodent_50k_001Hz_energy_rules_all.csv';
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
p_rules = [compare{1}(:,6) compare{2}(:,6) compare{3}(:,6) compare{4}(:,6)];
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
    p_rules(:,4),...
    d_rules(:,4),...
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('DIV7 p=%.3g',p(1)),'DIV7 Cohen d',...
    sprintf('DIV10 p=%.3g',p(2)),'DIV10 Cohen d',...
    sprintf('DIV12 p=%.3g',p(3)),'DIV12 Cohen d',...
    sprintf('DIV14 p=%.3g',p(4)),'DIV14 Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
%% compute the tfd for all networks within a pipeline
% set pipeline 
pipeline = 2;
% set how many measures are computed
nmeasures = 6;
% set how many models are computed
nmodels = 13;
% set up to how many divs are considered 
ndivs = 4;
% initialise correlation matrices, empty because sample sizes may alter
corr_local_observed = [];
corr_local_simulated = [];
corr_local_together = [];
% initalise the topological organization dissimilairty
todissim = [];
% labels
var_labels = string({...
    'obs degree','obs clustering','obs betweenn','obs length','obs eff','obs match',...
    'sim degree','sim clustering','sim between','sim length','sim eff','sim match'});
% loop over divs
for div = 1:ndivs
    % index the networks
    index_loop = find(index(:,1)==pipeline&index(:,3)==div);
    % compute the sample size for this selection
    nsamp_i = size(index_loop,1);
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
    for i = 1:nsamp_i;
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
            si = squeeze(N{i}{model}(:,i_low));
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
% set model order for visualisation
i = [2:3 9:13 4:8 1];
%  alter the order
todissim_permute = permute(todissim,[2 1 3]);
% make zeros into nans
todissim_permute(todissim_permute==0)=NaN;
% plot
h = figure; h.Position = [100 100 1200 600];
h = iosr.statistics.boxPlot(todissim_permute(:,:,i),...
    'showViolin',logical(0),...
    'theme','colorall',...
    'showScatter',logical(0),...
    'scatterMarker','x',...
    'scatterColor',[.5 .5 .5],...
    'symbolColor','k',...
    'symbolMarker','x',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylabel('{\itTFdissimilarity}'); xlabel('Days {\itin vitro} (DIV)');
xticklabels({'7','10','12','14'});
%ylim([0 0.2]);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% group meaned todissimilaity into rules
to_mean_div_rules = [];
to_mean_div_rules(:,:,1) = squeeze(mean(todissim_permute(:,:,[2:3]),3));
to_mean_div_rules(:,:,2) = squeeze(mean(todissim_permute(:,:,[4:8]),3));
to_mean_div_rules(:,:,3) = squeeze(mean(todissim_permute(:,:,[9:13]),3));
to_mean_div_rules(:,:,4) = squeeze(todissim_permute(:,:,1));
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
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/0.01Hz/rodent_50k_001Hz_tfd_rules_all.csv';
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
p_rules = [compare{1}(:,6) compare{2}(:,6) compare{3}(:,6) compare{4}(:,6)];
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
    p_rules(:,4),...
    d_rules(:,4),...
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('DIV7 p=%.3g',p(1)),'DIV7 Cohen d',...
    sprintf('DIV10 p=%.3g',p(2)),'DIV10 Cohen d',...
    sprintf('DIV12 p=%.3g',p(3)),'DIV12 Cohen d',...
    sprintf('DIV14 p=%.3g',p(4)),'DIV14 Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
%% visualise a specific comparison in terms of a culture and averaged over cultures
% settings
div = 4;
culture = 5;
model = 3;
% set colours
nstep = 500;
lin = linspace(0,1,nstep)'; 
c = ones(nstep,3); c(:,1)=lin; c(:,2)=lin;
d = ones(nstep,3); d(:,2)=lin; d(:,3)=lin;
col = [d;flip(c)];
% take the specific data
specific_obs = squeeze(corr_local_observed(div,culture,:,:));
specific_sim = squeeze(corr_local_simulated(div,culture,model,:,:));
% visualise the specific culture
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); imagesc(specific_obs); title('observed'); caxis([-1 1]); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r';  xticks([]); yticks([]);
subplot(1,2,2); imagesc(specific_sim); title('simulated'); caxis([-1 1]); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r';  xticks([]); yticks([]);
sgtitle(sprintf('culture %g, div %g, model %g, tod %g',culture,div,model,squeeze(squeeze(todissim(div,culture,model)))));
% take the averaged data
averaged_obs = squeeze(mean(corr_local_observed(div,:,:,:),2));
averaged_sim = squeeze(mean(corr_local_simulated(div,:,model,:,:),2));
% visualise the averaged data
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); imagesc(averaged_obs); title('observed'); colormap(col); caxis([-1 1]); colormap(col); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r'; xticks([]); yticks([]);
subplot(1,2,2); imagesc(averaged_sim); title('simulated'); colormap(col); caxis([-1 1]); colormap(col); u = colorbar;  u.TickDirection = 'out'; u.Label.String = 'r'; xticks([]); yticks([]);
sgtitle(sprintf('averaged div %g, model %g, tod %g',div,model,squeeze(mean(todissim(div,:,model),2))));
xticks([]); yticks([]);
%% plot a specific plot
g = figure; g.Position = [100 100 900 700];
plotind = [1 6 11 16; 2 7 12 17; 3 8 13 18; 4 9 14 19; 5 10 15 20];
for div = 1:ndivs
    % observed
    subplot(4,5,plotind(1,div)); imagesc(squeeze(mean(corr_local_observed(div,:,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % matching
    subplot(4,5,plotind(2,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,3,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % clu-avg
    subplot(4,5,plotind(3,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,9,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % deg-avg
    subplot(4,5,plotind(4,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,4,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % sptl
    subplot(4,5,plotind(5,div)); imagesc(squeeze(mean(corr_local_simulated(div,:,1,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
end
%% correlate energy and the to matrix
% pick the set
set = 1; 
% take the energy
ep = squeeze(top_e_mean(set,:,:))';
% take only those we have calculated the todissim of
epto = zeros(ndivs,ncultures,13);
% form into the same size to form correlations
for div = 1:ndivs;
    epto(div,:,:) = ep(index(:,1)==pipeline&index(:,3)==div,:);
end
% form indices for rules colours
ind = [ones(ndivs*ncultures*1,1);2*ones(ndivs*ncultures*2,1);3*ones(ndivs*ncultures*5,1);4*ones(ndivs*ncultures*5,1)];
colsim = {'y','r','g','b'};
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
        'markerfacecolor',colsim{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
end
sgtitle('all models'); xlabel('energy'); ylabel('TFdissimilarity');
xlim([0 0.6]);
%ylim([0 0.3]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial';
% set colours (ordered to make work)
sptl_c = [255 255 102];
homo_c = [255 202 205];
clu_c = [206 253 205];
deg_c = [206 204 255];
scatcol = [sptl_c; homo_c; clu_c; deg_c]./256;
scatcol_o = flip(scatcol);
% visualise scatter hist for all
u = figure; u.Position = [100 100 850 700];
k = scatterhist(x,y,'group',ind,...
    'marker','o','markersize',12,...
    'color',scatcol,'linestyle','-',...
    'kernel','on','direction','out','linewidth',5);
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
legend({'Spatial','Homophily','Clustering','Degree'},'box','off');
xlabel('Energy'); ylabel('{\itTFdissimilarity}');
xlim([0 0.65]);
box off;
%{
for col = 1:4;
    k(1).Children(col).MarkerFaceColor = scatcol_o(col,:);
end
%}
% plot the average at each time point
eptd = []; totd = [];
for div = 1:4;
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
        'markerfacecolor',colsim{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
    end
end
xlim([0 0.6]); ylim([0 0.1]);
xlabel('Energy'); ylabel('TFdissimilarity');
b = gca; b.FontName = 'Arial'; b.FontSize = 14; b.TickDir = 'out';
%% take optimal parameters for each network
% set model to evaluate
model = 3;
% set top n parameters
n = 2;
% sort the parameter combinations by energy
[e i] = sort(energy_sample,2);
% order the parameters by this
pselect = [];
for net = 1:nsamp;
    e = squeeze(energy_sample(net,model,:));
    [v,i] = sort(e);
    p = squeeze(parameters_sample(net,model,i(1:n),:));
    if n == 1;
    k = squeeze(p(:,1,:));
    else
    k = squeeze(mean(p,1));
    end
    pselect(net,:) = k;
end
% correlations to check
[r p] = corr(pselect);
%% correlate energy and the to matrix
% pick the set
set = 1; 
% take the energy
ep = squeeze(top_e_mean(set,:,:))';
% set the div crtieria
divs = 2:4;
nndivs = length(divs);
% update the old todissim
todissim_divs = todissim(divs,:,:);
% take only those we have calculated the todissim of
epto = zeros(length(divs),6,13);
% form into the same size to form correlations
for div = divs;
    i = div-1;
    epto(i,:,:) = ep(index(:,1)==pipeline&index(:,3)==div,:);
end
% form indices for rules colours
ind = [ones(nndivs*6*1,1);2*ones(nndivs*6*2,1);3*ones(nndivs*6*5,1);4*ones(nndivs*6*5,1)];
col = {'k','r','g','b'};
% across all models
[rall pall] = corr(epto(:),todissim_divs(:));
% take data
x = epto(:); y = todissim_divs(:);
% keep individual correlations based on rules
rrule  = []; prule = [];
% visualise
figure;
for rule = 1:4;
    [rrule(rule) prule(rule)] = corr(x(ind==rule&~isnan(y)),y(ind==rule&~isnan(y)));
    h = scatter(x(ind==rule),y(ind==rule),70,'o',...
        'markerfacecolor',col{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
end
sgtitle('all models'); xlabel('energy'); ylabel('TOdissimilarity');
xlim([0 0.6]);
ylim([0 0.3]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial';
% visualise scatter hist
u = figure;
scatterhist(x,y,'group',ind,...
    'marker','*','markersize',8,...
    'color',{'k','r','g','b'},'linestyle','-',...
    'kernel','on','direction','in');
b = gca; b.FontName = 'Arial'; b.FontSize = 14; b.TickDir = 'out';
legend({'spatial','homophily','clustering','degree'});
xlabel('Energy'); ylabel('TFdissimilarity');
xlim([0 0.6]);
% plot the average at each time point
eptd = []; totd = [];
for div = 1:nndivs;
    % energy
    eptd(div,:,1) = squeeze(epto(div,:,1));
    eptd(div,:,2) = mean(squeeze(epto(div,:,[2 3])),2,'omitnan');
    eptd(div,:,3) = mean(squeeze(epto(div,:,[4:8])),2,'omitnan');
    eptd(div,:,4) = mean(squeeze(epto(div,:,[9:13])),2,'omitnan');
    % todissim
    totd(div,:,1) = squeeze(todissim_divs(div,:,1));
    totd(div,:,2) = mean(squeeze(todissim_divs(div,:,[2 3])),2,'omitnan');
    totd(div,:,3) = mean(squeeze(todissim_divs(div,:,[4:8])),2,'omitnan');
    totd(div,:,4) = mean(squeeze(todissim_divs(div,:,[9:13])),2,'omitnan');
end
eptd = squeeze(mean(eptd,2,'omitnan')); totd = squeeze(mean(totd,2,'omitnan'));
% visualise
figure;
shape = {'o','square','pentagram','d'};
for rule = 1:4;
    for div = 1:nndivs;
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
%% run generative simulations capturing generative, integration and segregation nodal measures
% initialise the nodal measures
matching_D     = {}; % edge-wise D
matching_Dsum  = {}; % node-wise D
matching_K     = {}; % edge-wise K
matching_Ksum  = {}; % node-wise K
matching_Fk    = {}; % edge-wise parameterised K
matching_Fksum = {}; % node-wise parameterised K
matching_P     = {}; % edge-wise wiring probability
matching_Psum  = {}; % node-wise wiring probability
matching_A     = {}; % networks
segregation    = {}; % segregation: modularity q statistic
integration    = {}; % integration: global efficiency
% loop through subjecst
for net = 1:nsamp;
    tic  
    % set target networks
    Atgt_set    = all_o_networks;
    % seed network is empty
    A           = zeros(size(Atgt_set{net}));
    % compute the connections in the seed
    mseed       = nnz(A)/2;
    % euclidean distance
    D           = all_d_networks{net};
    % set target
    Atgt        = Atgt_set{net};
    % set the parameter for this subject
    params      = pselect(net,:);
    % number of bi-directional connections
    m           = nnz(Atgt)/2;
    % model var
    modelvar    = [{'powerlaw'},{'powerlaw'}];
    % minimum edge
    epsilon     = 1e-5;
    % run initial matching
    n           = length(D);
    nparams     = size(params,1);
    b           = zeros(m,nparams);
    K           = matching_ind(A);
    K           = K + K';
    % keep the K, Fk, P, A and statistics at each iteration, under the current parameter
    Kall        = zeros(m-mseed,n,n);
    Fkall       = zeros(m-mseed,n,n);
    Pall        = zeros(m-mseed,n,n);
    Aall        = zeros(m-mseed,n,n);
    Seg         = zeros(m-mseed,1);
    Int         = zeros(m-mseed,1);
    % save the first K
    Kall(1,:,:) = K;
    % save the first A
    Aall(1,:,:) = A;
    % compute the initial segregation: classic modularity
    try % in case the seed is empty
    [~,q]  = modularity_und(A,1);
    Seg(1) = q;
    end
    % compute the initial integration: global efficiency
    try
    Int(1) = efficiency_bin(A);
    end
    % loop thorugh the subsequent connections
    for iparam = 1:nparams
        eta = params(iparam,1);
        gam = params(iparam,2);
        K = K + epsilon;
        n = length(D);
        mseed = nnz(A)/2;
        mv1 = modelvar{1};
        mv2 = modelvar{2};
        switch mv1
            case 'powerlaw'
                Fd = D.^eta;
            case 'exponential'
                Fd = exp(eta*D);
        end
        switch mv2
            case 'powerlaw'
                Fk = K.^gam;
            case 'exponential'
                Fk = exp(gam*K);
        end
        Ff = Fd.*Fk.*~A;
        [u,v] = find(triu(ones(n),1));
        indx = (v - 1)*n + u;
        P = Ff(indx);
        % save the first parameterised K in Fk
        Fkall(1,:,:)  = Fk;
        % save the first probabilities
        Ff(isinf(Ff)) = 0;
        Pall(1,:,:)   = Ff;
        step = 2; % added in
        for ii = (mseed + 1):m
            C = [0; cumsum(P)];
            r = sum(rand*C(end) >= C);
            uu = u(r);
            vv = v(r);
            A(uu,vv) = 1;
            A(vv,uu) = 1;
            updateuu = find(A*A(:,uu));
            updateuu(updateuu == uu) = [];
            updateuu(updateuu == vv) = [];
            updatevv = find(A*A(:,vv));
            updatevv(updatevv == uu) = [];
            updatevv(updatevv == vv) = [];
            c1 = [A(:,uu)', A(uu,:)];
            for i = 1:length(updateuu)
                j = updateuu(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(uu) = 0;  use(uu+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(uu,j) = epsilon;
                    K(j,uu) = epsilon;
                else
                    K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,uu) = K(uu,j);
                end
            end
            c1 = [A(:,vv)', A(vv,:)];
            for i = 1:length(updatevv)
                j = updatevv(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(vv) = 0;  use(vv+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(vv,j) = epsilon;
                    K(j,vv) = epsilon;
                else
                    K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,vv) = K(vv,j);
                end
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            Kall(step,:,:)  = K;                     % added in (K)
            Fkall(step,:,:) = Fk;                    % added in (Fk)
            Aall(step,:,:)  = A;                     % added in (A)
            % compute the segregation: classic modularity
            [~,q]   = modularity_und(A,1);
            Seg(ii) = q;
            % compute the integration: global efficiency
            Int(ii) = efficiency_bin(A);
            % compute probabilites
            Ff = Fd.*Fk.*~A;
            P = Ff(indx);
            % remove infinite values (self connections)
            Ff(isinf(Ff))   = 0;
            Pall(step,:,:)  = Ff;                    % added in (P)
            % change the step
            step = step+1;
            % display progress
            disp(sprintf('network %g connection %g of %g formed',net,ii,m));
        end
        % keep the data
        b(:,iparam)         = find(triu(A,1));
        matching_K{net}     = Kall;                  
        matching_Fk{net}    = Fkall;                
        matching_Ksum{net}  = squeeze(sum(Kall,2)); 
        matching_Fksum{net} = squeeze(sum(Fkall,2));
        matching_P{net}     = Pall;                  
        matching_Psum{net}  = squeeze(sum(Pall,2)); 
        matching_A{net}     = Aall;                  
        segregation{net}    = Seg;
        integration{net}    = Int;
    end
    time = toc;
    disp(sprintf('generative model complete: network %g took %g seconds',net,time));
end
%% compute summary statistics over time
% number of measures
nmeasures = 7;
% initialise
statistics = zeros(npipelines,ncultures,ndivs,nmeasures);
% statistics labels
statistics_labels = string({...
    'Density',...
    'Degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Edge length',...
    'Modularity'});
% loop over networks
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            % take the binary network
            network = rodent_data{pipeline}.table_bu_shared{culture}{div};
            % take the coordinates
            coords = rodent_data{pipeline}.table_xy_shared{culture}{div};
            % compute the euclidean distance
            euclidean = squareform(pdist(coords));
            % compute global statistics
            % density
            statistics(pipeline,culture,div,1) = density_und(network);
            % total number of connections
            statistics(pipeline,culture,div,2) = sum(degrees_und(network));
            % global efficiency
            statistics(pipeline,culture,div,3) = efficiency_bin(network,0);
            % betweenness
            statistics(pipeline,culture,div,4) = mean(betweenness_bin(network));
            % clustering
            statistics(pipeline,culture,div,5) = mean(clustering_coef_bu(network));
            % average weighted length
            statistics(pipeline,culture,div,6) = mean(mean(euclidean.*network));
            % modularity
            [~,statistics(pipeline,culture,div,7)] = modularity_und(network);
        end
    end
end
%% plot segregation and integation as a function of connections
% visualise segregation
h = figure; h.Position = [100 100 600 300];
% loop over networks across all networks
for net = 1:nsamp;
    plot(segregation{net},'linewidth',1,'color',[238 50 51]./256);
    hold on;
end
ylabel('measure'); xlabel('simulated time');
b = gca; b.TickDir = 'out';
b.FontName = 'Arial';
sgtitle('all networks');
% visualise integation across all networks
for net = 1:nsamp;
    plot(integration{net},'linewidth',1,'color',[102 167 197]./256);
    hold on;
end
ylabel('measure'); xlabel('simulated time');
b = gca; b.TickDir = 'out';
b.FontName = 'Arial';
sgtitle('all networks');
%% scale simulated segregation against observed segregation
% pick the pipeline to compute
pipeline = 2;
% take the observed modularity
stat = 7;
% visualise the number of connections
statv = squeeze(statistics(pipeline,:,:,stat));
% alter the obseved statistics to be at specific points
plotdatas = nan(size(statv,1),16);
plotdatas(:,7)=statv(:,1);
plotdatas(:,10)=statv(:,2);
plotdatas(:,12)=statv(:,3);
plotdatas(:,14)=statv(:,4);
% get the 14 day networks of a particular pipeline
neti = index(:,1)==pipeline&index(:,3)==4;
% get the simulated segregation data of these networks
segregationi = segregation(find(neti));
% form into a matrix
u = [];
for k = 1:length(segregationi);
    u(k) = size(segregationi{k},1);
end
[v,~] = max(u);
j = nan(sum(neti),v);
for k = 1:length(segregationi);
    j(k,1:length(segregationi{k})) = segregationi{k};
end
% take the average but omitting nans
y = mean(j,'omitnan');
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
xticks([-1:15]); xlim([-1 15]);
%xticklabels({'','0','','','','','','','7','','','10','','12','','14',''});
xticklabels({'','0%','','','','','','','50%','','','','','','','100%',''});
ylabel('Modularity Q'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
%% get the segregation 
% take the observed segregation
ind = find(sum(isnan(plotdata)));
observed_segregation = plotdata;
observed_segregation(:,ind) = [];
% take the simualated segregation at the 7/14, 10/14, 12/14, 14/14 time points
simulated_segregation = [];
for k = 1:length(segregationi);
    % take the segregation data
    u = segregationi{k};
    % get total length
    m = length(u);
    % get the four timepoints
    a = round(1/2*m); b = round(10/14*m); c = round(12/14*m); d = m;
    % index those four time points
    simulated_segregation(k,:) = [u(a),u(b),u(c),u(d)];
end
% put together
segregation_comparison = [];
segregation_comparison(:,:,1) = observed_segregation;
segregation_comparison(:,:,2) = simulated_segregation;
% visualise
figure;
h = iosr.statistics.boxPlot(segregation_comparison,...
    'theme','colorall');
b = gca; b.TickDir = 'out';
xticklabels(div_labels); 
ylabel('Modularity');
%% scale simulated integration against observed integration
% pick the pipeline to compute
pipeline = 2;
% take the observed modularity
stat = 3;
% visualise the number of connections
statv = squeeze(statistics(pipeline,:,:,stat));
% alter the obseved statistics to be at specific points
plotdatai = nan(size(statv,1),16);
plotdatai(:,7)=statv(:,1);
plotdatai(:,10)=statv(:,2);
plotdatai(:,12)=statv(:,3);
plotdatai(:,14)=statv(:,4);
% get the 14 day networks of a particular pipeline
neti = index(:,1)==pipeline&index(:,3)==4;
% get the simulated segregation data of these networks
integrationi = integration(find(neti));
% form into a matrix
u = [];
for k = 1:length(integrationi);
    u(k) = size(integrationi{k},1);
end
[v,~] = max(u);
j = nan(sum(neti),v);
for k = 1:length(integrationi);
    j(k,1:length(integrationi{k})) = integrationi{k};
end
% take the average but omitting nans
y = mean(j,'omitnan');
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
xticks([-1:15]); xlim([-1 15]);
% xticklabels({'','0','','','','','','','7','','','10','','12','','14',''});
xticklabels({'','0%','','','','','','','50%','','','','','','','100%',''});
% xticks([]);
ylabel('Global Efficiency'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
%% get the observed and simulated integration data
% take the observed integration
ind = find(sum(isnan(plotdata)));
observed_integration = plotdata;
observed_integration(:,ind) = [];
% take the simualated integration at the 7/14, 10/14, 12/14, 14/14 time points
simulated_integration = [];
for k = 1:length(integrationi);
    % take the integration data
    u = integrationi{k};
    % get total length
    m = length(u);
    % get the four timepoints
    a = round(1/2*m); b = round(10/14*m); c = round(12/14*m); d = m;
    % index those four time points
    simulated_integration(k,:) = [u(a),u(b),u(c),u(d)];
end
% put together
integration_comparison = [];
integration_comparison(:,:,1) = observed_integration;
integration_comparison(:,:,2) = simulated_integration;
% visualise
figure;
h = iosr.statistics.boxPlot(integration_comparison,...
    'theme','colorall');
b = gca; b.TickDir = 'out';
xticklabels(div_labels); 
ylabel('Global Efficiency');
%% compute variance explained
%{
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/integration_comparison');
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/statistics/segregation_comparison');
%}
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
%% keep data
rodent_50k_001Hz_trajectory_data = struct;
rodent_50k_001Hz_trajectory_data.efficiency_comparison = integration_comparison;
rodent_50k_001Hz_trajectory_data.efficiency_trajectory.simulated = integrationi;
rodent_50k_001Hz_trajectory_data.efficiency_trajectory.observed = plotdatai;
rodent_50k_001Hz_trajectory_data.modularity_comparison = segregation_comparison;
rodent_50k_001Hz_trajectory_data.modularity_trajectory.simulated = segregationi;
rodent_50k_001Hz_trajectory_data.modularity_trajectory.observed = plotdatas;
rodent_50k_001Hz_trajectory_data.info = string({'6 cultures x 4 divs x 2 conditions (observed v simulated). The simulated statistics were taken from the best fitting 14 day simulated networks at ~50%, 71%, 86% and 100% to correspond to 7, 10, 12 and 14 days.'});
rodent_50k_001Hz_trajectory_data.author = string({'computed by danyal akarca, university of cambridge, 25/01/22'});
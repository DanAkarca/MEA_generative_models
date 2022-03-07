%% Exploration and generative model analysis of gabazine experimentation
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load experimental data
clear; clc;
% load gabazine networks
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
% equation
eqn = string({'KS_k','KS_c','KS_b','KS_d'});
% define model types
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
% group label
group_label = {'Gabazine','Control'};
% pipeline
pipeline_labels = {'5ms','10ms','20ms'};
%% set hyperparameters of the dataset
% number of cultures
ncultures = [3 6];
% culture ids
cultureid{1} = [2 3 6]; cultureid{2} = 1:6;
% total number of networks
nsamp = sum(ncultures);
% number of models
nmodels = 13;
%% load 0.01Hz generative model data
%{
% initialise
energy_sample = zeros(nsamp,13,20000);
ks_sample = zeros(nsamp,13,20000,4);
networks_sample = cell(nsamp,1);
parameters_sample = zeros(nsamp,13,20000,2);
errors = zeros(nsamp,1);
% loop over pipelines, cultures and divs
step = 1;
for group = 1:2;
    for culture = 1:ncultures(group);;
        for model = 1:nmodels;
            if group==1;
                % load this network's generative model output
                load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03200_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/human_org_1_%g_1_generative_model_%g.mat',...
                    cultureid{group}(culture),model));
            else if group == 2;
                    % load this network's generative model output
                    load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/human_org_1_%g_1_generative_model_%g.mat',...
                    culture,model));
                end
            end
            % keep index
            index(step,1) = group; index(step,2) = cultureid{group}(culture);
            % assign
            energy_sample(step,model,:) = output.energy;
            ks_sample(step,model,:,:) = output.ks;
            networks_sample{step}(model,:,:) = output.networks;
            parameters_sample(step,model,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('network_%g_%g %s loaded',group,culture,models(model)));
        end
        step = step + 1;
    end
end
% save data as a struct
rodent_gabazine_50k_001Hz_generative_models = struct;
rodent_gabazine_50k_001Hz_generative_models.energy = energy_sample;
rodent_gabazine_50k_001Hz_generative_models.ks = ks_sample;
rodent_gabazine_50k_001Hz_generative_models.networks = networks_sample;
rodent_gabazine_50k_001Hz_generative_models.parameters = parameters_sample;
rodent_gabazine_50k_001Hz_generative_models.models = models;
rodent_gabazine_50k_001Hz_generative_models.index = index;
rodent_gabazine_50k_001Hz_generative_models.info = string({'index column 1: 3 pipelines (STTC 5ms, 10ms, 20ms) - index column 2: 12 cultures - index column 3:3 days in vitro (14 days, 21 days, 28 days)'});
rodent_gabazine_50k_001Hz_generative_models.procedure = string({'voronoi tesselation n=20,000 parameters (5 steps of 4000 parameters, alpha = 2), eta and gamma with limits [-10 10]'});
rodent_gabazine_50k_001Hz_generative_models.author = string('written by dr danyal akarca, 21/01/22');
%}
%% place all data
index = rodent_gabazine_50k_001Hz_generative_models.index;
models = rodent_gabazine_50k_001Hz_generative_models.models;
energy_sample = rodent_gabazine_50k_001Hz_generative_models.energy;
ks_sample = rodent_gabazine_50k_001Hz_generative_models.ks;
parameters_sample = rodent_gabazine_50k_001Hz_generative_models.parameters;
nsamp = 9;
% get all networks
all_o_networks = {}; all_c_networks = {};
for k = 1:nsamp;
    all_o_networks{k} = data{index(k,1)}{2}.table_bu{index(k,2)};
    all_d_networks{k} = squareform(pdist(data{index(k,1)}{2}.table_xy{index(k,2)}));
end
%% firing rate matrix visualisation
% set the group
group = 2;
% set the group labels
glabels = {'M03200','M03212'};
% colormaps
colorset = [177 206 220; 140 108 108]./256;
% set xlims
xlims = [0 300];
ylims = [0 600];
% set the well to plot
well = 2;
% load the recording
load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/%s/210811/well%g/spk_data',...
    glabels{group},well));
% get spike times
spike_times = spk_data.spike_times;
% set bin
bin = 1;
% initialise
array = [];
% for each electrode bin the spikes
for electrode = 1:length(spike_times);
    % get spike times for the electrode
    x = spk_data.spike_times{electrode};
    % keep the number of spikes
    m(electrode) = length(x);
    % define the bins
    binlims = [0:bin:(ceil(max(x/10))*10)];
    [count,idx] = histc(x,binlims);
    % keep count
    array(electrode,1:length(count)) = count;
end
% visualise
h = figure; h.Position = [100 100 600 320];
imagesc(array);
xlabel('Time (s)');
ylabel('Putative neurons');
c = colorbar; c.Label.String = 'Firing Rate (Hz)'; 
% set(gca,'ColorScale','log');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
% plot the summed spike vector
summed_spike = sum(array);
h = figure; h.Position = [100 100 600 320];
plot(summed_spike,'color',colorset(group,:),'linewidth',2)
xlim(xlims); ylim(ylims);
xlabel('Time (s)');
ylabel('Co-active units');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
box off;
%% comptue group differences in firing rates
% initialsie step
step = 1;
% initialse
cat_fr = {}; ind = [];
for group = 1:2;
    for culture = 1:ncultures(group);
        well = cultureid{group}(culture);
        load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/%s/210811/well%g/spk_data',...
            glabels{group},well));
        % set numbers
        nmins = 30;
        % number 
        mean_firing_rate = [];
        nelectrode = length(spk_data.spike_times);
        for electrode = 1:nelectrode;
            mean_firing_rate(electrode) = length(spk_data.spike_times{electrode})/(nmins*60);
        end
        % mean the mean firing rate
        cat_fr{step} = mean_firing_rate; 
        mmfr(step) = mean(mean_firing_rate);
        smfr(step) = std(mean_firing_rate);
        % update
        step = step + 1;
    end
end
% statistics
p = ranksum(mmfr(1:3),mmfr(4:9));
disp(median(mmfr(1:3))); disp(iqr(mmfr(1:3)));
disp(median(mmfr(4:9))); disp(iqr(mmfr(4:9)));
% set up for visualise
uu = nan(6,2);
uu(:,1) = mmfr(4:9);
uu(1:3,2) = mmfr(1:3);
% visualise
figure; 
iosr.statistics.boxPlot(uu,...
    'theme','colorall','boxColor','w','showScatter',logical(1),'scatterColor','k','scatterMarker','x');
b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
ylabel('Mean firing rate (Hz)'); xticklabels({'Control','Gabazine'});
%% computing burstiness isi distributions
% calculate interspike interval distributions
% set the group labels
glabels = {'M03200','M03212'};
% colormaps
colorset = [140 180 180; 140 108 108]./256;
% wells
wells = {[2 3 6] [1:6]};
% set the timeframe
divide = 1000;
unit = 'ms';
% initialise isi
cat_isi_gabazine = [];
cat_isi_control = [];
isi_all = cell(1,2);
isi_p_gaba = [];
isi_p_control = [];
% loop and load
for group = 1:length(glabels);
    % wells
    for well = length(wells{group});
        % get well
        j = wells{group}(well);
        % load the recording
        load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_DIV14_Gabazine_30min/%s/210811/well%g/spk_data',...
            glabels{group},j));
        % spiketimes
        h = spk_data.spike_times;
        % nelectrodes
        nelectrodes = length(h);
        step = 1;
        % loop over electrodes
        keep_isi = [];
        for i = 1:nelectrodes;
            % take the interspike intervals
            isi = diff(h{i});
            % keep this
            keep_isi = [keep_isi; isi];
            % concatenate across
            if group == 1;
                cat_isi_gabazine = [cat_isi_gabazine;isi];
                isi_all{group}{well} = isi;
            end
            if group == 2;
                cat_isi_control = [cat_isi_control;isi];
            end
        end
    end
end
% visualise
h = figure; h.Position = [100 100 1000 600];
histogram(log(cat_isi_control/divide),'edgecolor','w','facecolor','r','facealpha',0.2);
hold on;
histogram(log(cat_isi_gabazine/divide),'edgecolor','w','facecolor','b','facealpha',0.2);
ylabel('Frequency'); xlabel(sprintf('log(ISI) (%s)',unit));
u = get(gca,'XTick');
xticklabels(10.^u);
b = gca; b.FontSize = 25; b.FontName = 'Arial'; b.TickDir = 'out';
box off;
% legend
legend({'Control','Gabazine'},'box','off');
%% visualise the networks
% select group
group = 2;
% select a pipeline
pipeline = 2;
% select a culture (2,3,6 for gabazine)
culture = 6;
% set the colours
gabacol = 1-copper(20); gabacol = gabacol(6,:);
contcol = pink(20); contcol = contcol(7,:);
cols = {gabacol contcol};
% compute the cultures
networks = data{group}{pipeline}.table_bu{culture};
coords = data{group}{pipeline}.table_xy{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 400 325]);
% sgtitle(sprintf('%s culture %g %s: binary networks',group_label{group},culture,pipeline_labels{pipeline}));
% loop through divs
g = graph(networks);
h = plot(g,...
    'XData',coords(:,1),...
    'YData',coords(:,2),...
    'EdgeColor',[120 120 120]./256,...
    'LineWidth',1.5,...
    'NodeColor',cols{group},...
    'NodeLabel',[],...
    'MarkerSize',0.75*degrees_und(networks)+0.01);
xticks([]); yticks([]);
set(gca,'visible','off');
%% compute group differences in statistics
% set pipeline 
pipeline = 2;
% number of measures
nmeasures = 10;
% initialise
global_statistics = zeros(nsamp,nmeasures);
% set number of permutations for small worldness calculation
nperm = 1000;
% statistics labels
global_statistics_labels = string({...
    'Network density (%)',...
    'Strength',...
    'Total degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Small-worldness',...
    'Edge length',...
    'Modularity',...
    'Matching'});
% initialise
step = 1;
% loop over 100k networks
for group = 1:2;
    for culture = 1:ncultures(group);
    % take the binary network
    network = data{group}{pipeline}.table_bu{cultureid{group}(culture)};
    % take the weighted network
    wnetwork = data{group}{pipeline}.table_wu{cultureid{group}(culture)};
    % take the euclidean distances
    euclidean = squareform(pdist(data{group}{pipeline}.table_xy{cultureid{group}(culture)}));
    % compute global statistics
    % density
    statistics(step,1) = density_und(network);
    % total strength
    statistics(step,2) = sum(strengths_und(wnetwork(~isnan(wnetwork))));
    % total number of connections
    statistics(step,3) = sum(degrees_und(network));
    % global efficiency
    statistics(step,4) = efficiency_bin(network,0);
    % betweenness
    statistics(step,5) = mean(betweenness_bin(network));
    % clustering
    statistics(step,6) = mean(clustering_coef_bu(network));
    % small worldness
        % compute nnode
        nnode = size(network,1);
        % compute number of edges
        m = nnz(network)/2;
        % compute observed
        clu = mean(clustering_coef_bu(network));
        cpl = charpath(network);
        % initialise
        clu_perm = [];
        cpl_perm = [];
        % compute nulls
        for j = 1:nperm;
            % form a random
            Arand = makerandCIJ_und(nnode,m);
            % form a lattic
            clu_perm(j) = mean(clustering_coef_bu(Arand));
            cpl_perm(j) = charpath(Arand);
        end
        % calclate means
        mean_clu_perm = mean(clu_perm);
        mean_cpl_perm = mean(cpl_perm);
        % calculate smw
        statistics(step,7) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    % average weighted length
    statistics(step,8) = mean(mean(euclidean.*network));
    % modularity
    [~,statistics(step,9)] = modularity_und(network);
    % homophily
    statistics(step,10) = mean(matching_ind(network),'all');
    % display
    disp(sprintf('network %g statistics computed',step));
    % update stp
    step = step + 1;
    end
end
%% group comparisons
% set stat
stat = 9;
% compute groups
gaba = statistics(1:3,stat);
con  = statistics(4:9,stat);
% staitistcs
p = ranksum(gaba,con);
% group
k = nan(6,2);
k(:,1) = con;
k(1:3,2) = gaba;
% visualise
h = figure; h.Position = [100 100 500 400];
iosr.statistics.boxPlot(k);
b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
xticklabels({'Control','Gabazine'});
ylabel(global_statistics_labels(stat));
% compute a range
% 1 - density, 3 - degree, 5 - betweenness, 7 - smw, 9 - edge length, 10 - matching
range = [1:10];
h = figure; h.Position = [100 100 1600 550];
for stat = 1:length(range);
    subplot(2,ceil(length(range)/2),stat);
    % compute groups
    gaba = statistics(1:3,range(stat));
    con  = statistics(4:9,range(stat));
    % staitistcs
    p = ranksum(gaba,con);
    % group
    k = nan(6,2);
    k(:,1) = con;
    k(1:3,2) = gaba;
    % boxplot
    iosr.statistics.boxPlot(k,...
        'theme','colorall',...
        'boxColor','w',...
        'showScatter',logical(1),...
        'scatterMarker','x',...
        'scatterColor','k');
    xticklabels({'C','G'});
    ylabel(global_statistics_labels(range(stat)));
    title(sprintf('p=%.3g',p),'FontSize',15);
    b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
end
%% look at energy landscape for a selected rule and network
% select model and network
model = 3;
net = 4;
% take the measure
e = squeeze(energy_sample(net,model,:));
pipeline_p = squeeze(parameters_sample(net,model,:,:));
% visualise
if model == 1
    h = figure;
    scatter(pipeline_p(:,1),e,100,e,'.');
    xlabel('eta'); ylabel('energy'); 
    ylim([0 1]);
    caxis([0 1]); c = colorbar; c.Label.String = 'energy';
else
    % plot the energy landscape in 3D
    h = figure;
    col = parula(20000);
    colu = col(round(e*20000),:);
    h = scatter3(pipeline_p(:,1),pipeline_p(:,2),e,'.'); 
    h.CData = colu;
    xlabel('eta'); ylabel('gamma'); 
    caxis([0 1]); c = colorbar; c.Label.String = 'energy';
    grid off;
    % plot the energy landscape in 2D
    h = figure;
    scatter(pipeline_p(:,1),pipeline_p(:,2),100,e,'.'); 
    xlabel('eta'); ylabel('gamma'); 
    caxis([0 1]); c = colorbar; c.Label.String = 'energy';
    grid off;
end
b = gca; b.TickDir = 'out';
%% look at energy landscape for all networks by group
% select model
model = 3;
% select the criteria for selection
criteria = index(:,1)==2;
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
% plot the energy landscape in 2D
h = figure;
h.Position = [100 100 500 400];
for net = 1:nsamp_select
    scatter(squeeze(p_select(net,:,1)),squeeze(p_select(net,:,2)),20,e_select(net,:),'.'); hold on;
    xlabel('\eta'); ylabel('\gamma');
end
caxis([0 1]); %c = colorbar; c.Label.String = 'energy';
b = gca; b.TickDir = 'out'; b.FontSize = 25;
xticks([]); yticks([]);
% set any limits
xlim([-10 10]); ylim([-10 10]);
% plot the energy landscape in 3D
h = figure;
h.Position = [100 100 500 400];
for net = 1:nsamp_select
    col = parula(20000);
    colu = col(round(e_select(net,:)*20000),:);
    u = scatter3(squeeze(p_select(net,:,1)),squeeze(p_select(net,:,2)),e_select(net,:),'.');
    u.CData = colu; xticks([]); yticks([]); zticks([]);
    hold on;
    %xlabel('eta'); ylabel('gamma');
end
caxis([0 1]); %c = colorbar; c.Label.String = 'energy';
b = gca; b.TickDir = 'out';
xticks([]); yticks([]);
% set any limits
xlim([-10 10]); ylim([-10 10]);
end
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
% specify the models
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
h = figure; h.Position = [100 100 1000 500];
iosr.statistics.boxPlot(driver,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'showScatter',logical(1),...
    'scatterColor',[.5 .5 .5],...
    'scatterAlpha',0.5,...
    'symbolColor',[.5 .5 .5],...
    'boxColor','k',...
    'boxAlpha',0.2); 
ylim([0 100]);
b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
xticklabels({'degree','clustering','betweenness','edge length'});
ylabel('max(KS)'); yticklabels({'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'});
%% visualise summary statistics across a specific group
% select which set of nset to view
set = 1;
% visualise both sets
% set the new order
i = [2:3 9:13 4:8 1];
% set the criteria
array = nan([6 2 13]);
for group = 1:2;
    order = [2 1];
    gview = order(group);
    criteria = index(:,1)==gview;
    e_select = top_e_mean(set,i,criteria);
    % form the array
    array(1:sum(criteria),group,:) = squeeze(e_select)';
end
% iosr boxplot
h = figure;
h.Position = [100 100 1000 500];
u = iosr.statistics.boxPlot(array,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
ylim([0 1]); ylabel('Energy'); 
xticklabels({'Control','Gabazine'});
xlabel('Generative models');
yticks([0:0.1:1]);
ylim([0 0.7]); % added to scale
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% group by generative rule
% visualise summary energy across divs across rule types: note, they've already been reordered above
% note that clustering and degree have already been swapped
e_rules_mean = [];
for group = 1:2;
    e_rules_mean(:,group,1) = squeeze(mean(array(:,group,[1:2]),3));
    e_rules_mean(:,group,2) = squeeze(mean(array(:,group,[3:7]),3));
    e_rules_mean(:,group,3) = squeeze(mean(array(:,group,[8:12]),3));
    e_rules_mean(:,group,4) = squeeze(array(:,group,13));
end
% all 
e_rules = nan(30,2,4);
a = squeeze(array(:,:,[1:2])); a = reshape(a,[6*2 2]); e_rules(1:6*2,:,1) = a;
b = squeeze(array(:,:,[3:7])); b = reshape(b,[6*5 2]); e_rules(1:6*5,:,2) = b;
c = squeeze(array(:,:,[8:12])); c = reshape(c,[6*5 2]); e_rules(1:6*5,:,3) = c;
d = squeeze(array(:,:,[13])); d = reshape(d,[6 2]); e_rules(1:6,:,4) = d;
% visualise
h = figure;
h.Position = [100 100 800 650];
u = iosr.statistics.boxPlot(e_rules,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k',...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
ylim([0 0.7]); ylabel('Energy');
xticklabels({'Control','Gabazine'});
xlabel('Generative models');
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% run the statistics on group differences in model fit for homophily
p = ranksum(squeeze(e_rules(:,1,1)),squeeze(e_rules(:,2,1)));
%% take optimal parameters for each network
% set model to evaluate
model = 3;
% set top n parameters
n = 1;
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
%% plot the parameters and the statistics
% seperate by group
k = nan(sum(index(:,1)==2),2,2);
k(4:9,:,1) = pselect(index(:,1)==2,:);
k(1:3,:,2) = pselect(index(:,1)==1,:);
% statistics on wiring parmeters
p = []; d = [];
eta = squeeze(k(:,1,:)); gamma = squeeze(k(:,2,:));
p(1) = ranksum(eta(:,1),eta(:,2));
d(1) = computeCohen_d(eta(:,1),eta(:,2));
p(2) = ranksum(gamma(:,1),gamma(:,2));
d(2) = computeCohen_d(gamma(:,1),gamma(:,2));
p = roundsd(p,3); d = roundsd(d,3);
% form the colours
gabacol = 1-copper(20); gabacol = gabacol(10,:);
contcol = pink(20); contcol = contcol(5,:);
cols{1} = gabacol; cols{2} = contcol;
% visualise
h = figure;
h.Position = [100 100 800 500];
for group = 1:2;
    subplot(1,2,group);
    u = iosr.statistics.boxPlot(k(:,:,group),...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor',cols{group},...
        'boxAlpha',1);
    ylabel('Magnitude'); xticks([1,2]);
    %ylim([-1.5 1.5]);
    xticklabels({'\eta','\gamma'});
    xlabel('Wiring parameter'); 
    b = gca; 
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 25;
end
% same but based on cultures
h = figure;
h.Position = [100 100 1000 400];
ylimitsg = [-1.2 1.2;-0.1 1.1];
ylabelsg = {'\eta','\gamma'};
for parameter = 1:2;
    subplot(1,2,parameter);
    u = iosr.statistics.boxPlot(permute(k(:,parameter,:),[1 3 2]),...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor','w',...
        'boxAlpha',1);
    ylabel(ylabelsg{parameter}); 
    xticks([1,2]);
    %ylim(ylimitsg(parameter,:));
    title(sprintf('p=%g, d=%g',p(parameter),d(parameter)));
    xticklabels({'Control','Gabazine'});
    b = gca;
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 25;
end
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
for net = 1:nsamp; %nsamp
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
%% change name of matching probability distributions and keep
rodent_gabazine_50k_001Hz_matching_probability = matching_P;
%% compute the variability in probability distributions
% seperate out the groups
gabazine_p = matching_P(1:3);
control_p = matching_P(4:9);
% plot the coefficient of variation of probability over time
cvg = {}; cvc = {}; cvg_scaled = []; cvc_scaled = [];
% compute gabazine
for culture = 1:3;
    % take network
    Ag = gabazine_p{culture};
    % initailise
    cvg{culture} = zeros(size(Ag,1),1);
    nnode = size(Ag,2);
    ind = triu(ones(nnode),1);
    pg{culture} = nan(size(Ag,1),size(find(ind),1));
    % loop over time
    for t = 1:size(Ag,1);
        % take the top triu only
        vector = squeeze(Ag(t,:,:)).*ind;
        vector = vector(vector>0);
        % keep this vector
        pc{culture}(1:length(vector),t)=vector;
        % take coefficient of variation
        cvg{culture}(t) = std(vector)/mean(vector);
        % display
        disp(sprintf('computed gabazine network %g time point %g',culture,t));
    end
    % keep a scaled version
    cvg_scaled(:,culture) = imresize(cvg{culture},[1000 1]);
    % diplay
    disp(sprintf('computed gabazine network %g',culture));
end
% compute controls
for culture = 1:6;
    % take network
    Ac = control_p{culture};
    % initailise
    cvc{culture} = zeros(size(Ac,1),1);
    nnode = size(Ac,2);
    ind = triu(ones(nnode),1);
    pc{culture} = nan(size(Ac,1),size(find(ind),1));
    % loop over time
    for t = 1:size(Ac,1);
        % take the top triu only
        vector = squeeze(Ac(t,:,:)).*ind;
        vector = vector(vector>0);
        % keep this vector
        pc{culture}(1:length(vector),t)=vector;
        % take coefficient of variation
        cvc{culture}(t) = std(vector)/mean(vector);
        % display
        disp(sprintf('computed control network %g time point %g',culture,t));
    end
    % keep a scaled version
    cvc_scaled(:,culture) = imresize(cvc{culture},[1000 1]);
    % diplay
    disp(sprintf('computed control network %g',culture));
end
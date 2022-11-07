%% figure 6
% written by danyal akarca
%% load data
clear; clc;
% load sttc
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/50k_rodent_div14_gabazine_exp';
load(strcat(gendir,'/M03200_gm_data_min_rate_0.01hz_0.001_alpha_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); gbz = gm_data;
load(strcat(gendir,'/M03212_gm_data_min_rate_0.01hz_0.001_alpha_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); ctrl = gm_data;
gendir_2 = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/50k_rodent_div14_Gbz_2';
load(strcat(gendir_2,'/M05565_gm_data_min_rate_0.01hz_0.001_alpha_Gbz_1800s_min_rate_001hz_lag10_jitter10_prob_gbz.mat')); gbz_2 = gm_data;
load(strcat(gendir_2,'/M05565_gm_data_min_rate_0.01hz_0.001_alpha_Washout_1800s_min_rate_001hz_lag10_jitter10_prob_washout.mat')); wash = gm_data;
rodent_data_all = {ctrl gbz gbz_2 wash}; clear gm_data;
% load generative model data
load(strcat(gendir,'/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); 
ctrl = results_generative_models;
load(strcat(gendir,'/M03200_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); 
gbz = results_generative_models;
load(strcat(gendir_2,'/Gbz_M05565_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); 
gbz_2 = results_generative_models;
load(strcat(gendir_2,'/Washout_M05565_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); 
wash = results_generative_models;
% collect the data together
gbz_data = {ctrl gbz gbz_2 wash};
% addpath of brewerplot
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer/');
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
% addpath of Cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen/');
% addpath of round
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd/');
%% set hyperparameters
% number of groups
ngroups = 4;
% number of cultures
ncultures = [6 3 6 6];
% culture ids
cultureid = {1:6 [2 3 6] 1:6 1:6};
% total number of networks
nsamp = sum(ncultures);
% number of models
nmodels = 13;
% initailise
noconnect = {};
% loop over data
for group = 1:ngroups;
    % loop over cultures
    networks = {};
    for culture = 1:ncultures(group);
        % initialise
        networks{culture} = [];
        % keep the data
        networks{culture} = [];
        % loop over models
        for model = 1:nmodels;
            if group == 1;
                load(strcat(gendir,sprintf(...
                    '/M03212_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_1_generative_model_%g.mat',...
                    cultureid{1}(culture),model)));
                networks{culture}(:,:,model) = output.networks;
                noconnect{group}(culture) = size(output.networks,1);
            end
            if group == 2;
                load(strcat(gendir,sprintf(...
                    '/M03200_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_1_generative_model_%g.mat',...
                    cultureid{2}(culture),model)));
                networks{culture}(:,:,model) = output.networks;
                noconnect{group}(culture) = size(output.networks,1);
            end
            if group == 3;
                load(strcat(gendir_2,sprintf(...
                    '/Gbz_M05565_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_1_generative_model_%g.mat',...
                    cultureid{3}(culture),model)));
                networks{culture}(:,:,model) = output.networks;
                noconnect{group}(culture) = size(output.networks,1);
            end
            if group == 4;
                load(strcat(gendir_2,sprintf(...
                    '/Washout_M05565_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_1_generative_model_%g.mat',...
                    cultureid{4}(culture),model)));
                networks{culture}(:,:,model) = output.networks;
                noconnect{group}(culture) = size(output.networks,1);
            end
        end
    end
    % assign
    gbz_data{group}.networks = networks;
    % display
    disp(sprintf('Group %g, %g simulated networks loaded',group,culture));
end
%% remove nodes with no connections from empirical data
% initialise
rodent_data = rodent_data_all;
% loop over networks and remove nodes with no connections
for group = 1:ngroups;
    for culture = 1:ncultures;
        % get network
        Atgt = rodent_data_all{group}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = rodent_data_all{group}.table_wu{culture};
        xy = rodent_data_all{group}.table_xy{culture};
        % keep only nodes with edges
        rodent_data{group}.table_bu{culture} = Atgt(ids,ids);
        rodent_data{group}.table_wu{culture} = wu(ids,ids);
        rodent_data{group}.table_xy{culture} = xy(ids,:);
    end
end
%% pool together the gabazine cultures (run once at the beginning)
ix = [7 8 9];
jx = [2 3 6];
% pool the empirical data
for i = 1:3;
    rodent_data_all{3}.table_bu{ix(i)} = rodent_data_all{2}.table_bu{jx(i)};
    rodent_data_all{3}.table_wu{ix(i)} = rodent_data_all{2}.table_wu{jx(i)};
    rodent_data_all{3}.table_xy{ix(i)} = rodent_data_all{2}.table_xy{jx(i)};
end
rodent_data_all = {rodent_data_all{1} rodent_data_all{3} rodent_data_all{4}};
% pool the generative data
for i = 1:3;
    gbz_data{3}.energy(ix(i),:,:) = gbz_data{2}.energy(i,:,:);
    gbz_data{3}.ks(ix(i),:,:,:) = gbz_data{2}.ks(i,:,:,:);
    gbz_data{3}.parameters(ix(i),:,:,:) = gbz_data{2}.parameters(i,:,:,:);
    gbz_data{3}.networks{ix(i)} = gbz_data{2}.networks{i};
end
gbz_data = {gbz_data{1} gbz_data{3} gbz_data{4}};
% update groups
ngroups = 3;
ncultures = [6 9 6];
cultureid = {1:6 1:9 1:6};
%% sanity check on the control versus gabazine config files
dc = {}; dg = {};
for i = 1:6;
    dc{i} = squareform(pdist(rodent_data_all{1}.table_xy{i}));
end
for i = 1:9;
    dg{i} = squareform(pdist(rodent_data_all{2}.table_xy{i}));
end
for i = 1:6;
    dw{i} = squareform(pdist(rodent_data_all{3}.table_xy{i}));
end
h = figure; h.Position = [100 100 1300 400];
subplot(1,3,1); 
for i = 1:6;
    histogram(dc{i},'normalization','countdensity','facealpha',.5,'edgecolor','w');
    hold on; box off;
    xlim([0 4500]);
    ylim([0 19]);
    title('Control (DIV14)');
    ylabel('Count Density');
    xlabel('Distance');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
end
subplot(1,3,2);
for i = 1:9;
    histogram(dg{i},'normalization','countdensity','facealpha',.5,'edgecolor','w');
    hold on; box off;
    xlim([0 4500]);
    ylim([0 19]);
    title('Gabazine (DIV14)');
    ylabel('Count Density');
    xlabel('Distance');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
end
for i = 1:6;
    subplot(1,3,3);
    histogram(dw{i},'normalization','countdensity','facealpha',.5,'edgecolor','w');
    hold on; box off;
    xlim([0 4500]);
    ylim([0 19]);
    title('Washout (DIV15)');
    ylabel('Count Density');
    xlabel('Distance');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
end
%% cellular properties
%% visualise the networks
% set which data to plot
dataplot = rodent_data_all;
% select group
group = 3;
% select a culture (2,3,6 for gabazine)
culture = 2;
% set the colours 
contcol = pink(20); contcol = contcol(7,:);
gabacol = 1-copper(20); gabacol = gabacol(6,:);
washcol = parula(20); washcol = washcol(1,:);
cols = {contcol gabacol washcol};
% compute the cultures
networks = dataplot{group}.table_bu{culture};
coords = dataplot{group}.table_xy{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 400 325]);
% sgtitle(sprintf('%s culture %g %s: binary networks',group_label{group},culture,pipeline_slabels{pipeline}));
% loop through divs
g = graph(networks);
h = plot(g,...
    'XData',dataplot{group}.table_xy{culture}(:,2),...
    'YData',dataplot{group}.table_xy{culture}(:,1),...
    'markersize',degrees_und(networks)*.75+0.01,...
    'nodelabel',[],...
    'edgecolor',[.6 .6 .6],...
    'nodecolor',cols{group});
xticks([]); yticks([]);
set(gca,'visible','off');
%% compute global statistics
% set which data to plot
dataplot = rodent_data_all;
% number of measures
nmeasures = 11;
% initialise
global_statistics = zeros(nsamp,nmeasures);
% set number of permutations for small worldness calculation
nperm = 1000;
% statistics labels
global_statistics_labels = string({...
    'Network density (%)',...
    'Number of nodes',...
    'Strength',...
    'Total degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Small-worldness',...
    'Edge length (\mum)',...
    'Modularity',...
    'Matching'});
% initialise
step = 1;
% loop over 100k networks
for group = 1:ngroups;
    for culture = 1:ncultures(group);
    % take the binary network
    network = dataplot{group}.table_bu{cultureid{group}(culture)};
    % take the weighted network
    wnetwork = dataplot{group}.table_wu{cultureid{group}(culture)};
    % take the euclidean distances
    euclidean = squareform(pdist(dataplot{group}.table_xy{cultureid{group}(culture)}));
    % compute global statistics
    % density
    statistics(step,1) = density_und(network)*100;
    % number of nodes
    statistics(step,2) = nnz(network)/2;
    % total strength
    statistics(step,3) = sum(strengths_und(wnetwork(~isnan(wnetwork))));
    % total number of connections
    statistics(step,4) = sum(degrees_und(network));
    % global efficiency
    statistics(step,5) = efficiency_bin(network,0);
    % betweenness
    statistics(step,6) = mean(betweenness_bin(network));
    % clustering
    statistics(step,7) = mean(clustering_coef_bu(network));
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
        statistics(step,8) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    % average weighted length
    statistics(step,9) = mean(mean(euclidean.*network));
    % modularity
    [~,statistics(step,10)] = modularity_und(network);
    % homophily
    statistics(step,11) = mean(matching_ind(network),'all');
    % display
    disp(sprintf('network %g statistics computed',step));
    % update stp
    step = step + 1;
    end
end
%% group comparisons in network statistics
% compute a range
% 1 - density, 3 - degree, 5 - betweenness, 7 - smw, 9 - edge length, 10 - matching
range = [1:11];
h = figure; h.Position = [100 100 2000 450];
for stat = 1:length(range);
    subplot(2,ceil(length(range)/2),stat);
    % compute groups
    con  = statistics(1:6,range(stat));
    gaba = statistics(7:15,range(stat));
    wash = statistics(16:21,range(stat));
    % staitistcs
    p = ranksum(gaba,con);
    % group
    k = nan(9,3);
    k(1:6,1) = con;
    k(:,2) = gaba;
    k(1:6,3) = wash;
    % boxplot
    a = iosr.statistics.boxPlot(k,...
        'theme','colorall',...
        'boxColor','w',...
        'showScatter',logical(1),...
        'scatterMarker','x',...
        'scatterColor','k',...
        'symbolMarker','x',...
        'symbolColor','k');
    % set colour
    a.handles.box(1).FaceColor = contcol;
    a.handles.box(1).FaceAlpha = .6;
    a.handles.box(2).FaceColor = gabacol;
    a.handles.box(2).FaceAlpha = .6;
    xticklabels({'C','G','W'});
    ylabel(global_statistics_labels(range(stat)));
    title(sprintf('p=%.3g',p),'FontSize',15);
    b = gca; b.TickDir = 'out'; b.FontSize = 15; b.FontName = 'Arial';
end
%% compute top performing simulations
% set number of models
nmodels = 13;
% set the sample size for the 50k and 100k cultures
ngroup = [6 9 6];
% define the top n parameters
nset = [1 10 50 100];
% compute the minimum
for group = 1:length(ngroup);
    % get the group
    energy = gbz_data{group}.energy;
    parameters = gbz_data{group}.parameters;
    % initialise
    top_energy = cell(length(nset),1);
    top_energy_mean = zeros(length(nset),13,ngroup(group));
    top_parameters = cell(length(nset),1);
    top_parameters_mean = zeros(length(nset),13,ngroup(group),2);
    % ouput
    gbz_generative{group} = struct;
    % loop over sets
    for no = 1:length(nset)
        % take the actual amount of top performing parameters
        n = nset(no);
        % loop over models
        for model = 1:nmodels;
            % take energies for this model
            pipeline_d = squeeze(energy(:,model,:))';
            % rank them for each subject
            [v i] = sort(pipeline_d);
            % take top n energies and their indices
            n_e = v(1:n,:);
            n_i = i(1:n,:);
            % take the corresponding parameters
            u = zeros(ngroup(group),n,2);
            for s = 1:ngroup(group);
                % keep parameters
                u(s,:,:) = squeeze(parameters(s,model,n_i(:,s),:));
            end
            % if top parameter only
            if n == 1
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
    gbz_generative{group}.top_energy = top_energy;
    gbz_generative{group}.top_energy_mean = top_energy_mean;
    gbz_generative{group}.top_parameters = top_parameters;
    gbz_generative{group}.top_parameters_mean = top_parameters_mean;
end
%% visualise
% take the top n energy values
n = 1;
% combine
data = nan(nmodels,ngroup(2),3);
data(:,1:ngroup(1),1) = squeeze(gbz_generative{1}.top_energy_mean(n,:,:));
data(:,1:ngroup(2),2) = squeeze(gbz_generative{2}.top_energy_mean(n,:,:));
data(:,1:ngroup(3),3) = squeeze(gbz_generative{3}.top_energy_mean(n,:,:));
% set order
i = [2:3 9:13 4:8 1];
data = data(i,:,:);
% permute the order
data = permute(data,[2 3 1]);
% iosr boxplot - all models
% visualise
h = figure; h.Position = [100 100 1000 450];
u = iosr.statistics.boxPlot(data,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticklabels({'Control','Gabazine','Washout'}); xlabel('Condition');
ylim([0 0.6]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% iosr boxplot - grouped models
energy_gbz_rules = nan(ngroup(2)*5,length(ngroup),4);
a = squeeze(data(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); energy_gbz_rules(1:ngroup(2)*2,:,1) = a';
b = squeeze(data(:,:,[3:7])); b = permute(b,[2 1 3]); b = b(:,:); energy_gbz_rules(1:ngroup(2)*5,:,2) = b';
c = squeeze(data(:,:,[8:12])); c = permute(c,[2 1 3]); c = c(:,:); energy_gbz_rules(1:ngroup(2)*5,:,3) = c';
d = squeeze(data(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); energy_gbz_rules(1:ngroup(2),:,4) = d';
% visualise
h = figure; h.Position = [100 100 1000 450];
u = iosr.statistics.boxPlot(energy_gbz_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticklabels({'Control','Gabazine','Washout'}); xlabel('Condition');
ylim([0 0.6]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% statistical comparisons between homophily
rule = [1];
a = squeeze(energy_gbz_rules(:,1,rule));
b = squeeze(energy_gbz_rules(:,2,rule));
c = squeeze(energy_gbz_rules(:,3,rule));
p1 = ranksum(a(:),b(:));
p2 = ranksum(b(:),c(:));
p3 = ranksum(a(:),c(:));
disp(p1); disp(p2); disp(p3);
% plot homophily energy difference
gabacol = [0.6053 0.7533 0.8429];
contcol = [0.7375 0.4588 0.4588];
washcol = [0.5514 0.5393 0.7533];
u = figure; u.Position = [100 100 800 650];
z = iosr.statistics.boxPlot(energy_gbz_rules(:,:,rule),...
    'theme','colorall',...
    'showScatter',logical(1),...
    'scatterMarker','x',...
    'scatterColor','k',...
    'boxColor','r');
xlabel('Condition');
ylabel('Homophily energy');
xticklabels({'C','G','W'});
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
z.handles.box(1).FaceColor = contcol;
z.handles.box(2).FaceColor = gabacol;
z.handles.box(3).FaceColor = washcol;

%% plot the energy landscape
% form the modeltype
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','diff-prod'});
% select model
model = 3;
% plot the energy landscape
h = figure;
h.Position = [100 100 1600 450];
for group = 1:length(ngroup);
    % form the subplot
    subplot(1,3,group);
    % take the measure
    e = squeeze(gbz_data{group}.energy(:,model,:));
    p = squeeze(gbz_data{group}.parameters(:,model,:,:));
    % visualise
    if model == 1
        eta = squeeze(p(:,:,1));
        m = parula(10000);
        for net = 1:ngroup(group);
            % get colours
            pipeline_d = round(e(net,:) .* 10000);
            col = m(pipeline_d,:);
            % plot
            u = scatter(eta(net,:),e(net,:),20,col,'.'); ylabel('energy'); xlabel('/eta'); ylim([0 0.6]);
            b = gca; b.TickDir = 'out';
            b.FontName = 'Arial';
            yticks([]); xticks([]); hold on;
        end
    else
    for net = 1:ngroup(group)
        u = scatter(squeeze(p(net,:,1)),squeeze(p(net,:,2)),20,e(net,:),'.'); hold on;
        xlabel('\eta'); ylabel('\gamma');
    end
    xticks([-10 0 10]); yticks([-10 0 10]);
    caxis([0 1]); %c = colorbar; c.Label.String = 'Energy';
    % u = sgtitle(sprintf('Matching generative model',models(model))); u.FontSize = 25; u.FontName = 'Arial';
    b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
    end
    % set any limits
    xlim([-10 10]); ylim([-10 10]);
end
% plot of a single culture
group = 3;
net = 2;
e = squeeze(gbz_data{group}.energy(:,model,:));
p = squeeze(gbz_data{group}.parameters(:,model,:,:));
h = figure;
h.Position = [100 100 550 450];
if model == 1
    eta = squeeze(p(:,:,1));
    m = parula(10000);
    pipeline_d = round(e(net,:) .* 10000);
    col = m(pipeline_d,:);
    u = scatter(eta(net,:),e(net,:),20,col,'.'); ylabel('energy'); xlabel('/eta'); ylim([0 0.6]);
    b = gca; b.TickDir = 'out';
    b.FontName = 'Arial';
    yticks([]); xticks([]); hold on;
else
    u = scatter(squeeze(p(net,:,1)),squeeze(p(net,:,2)),20,e(net,:),'.'); hold on;
    xlabel('\eta'); ylabel('\gamma');
    xticks([-10 0 10]); yticks([-10 0 10]);
    caxis([0 1]); 
    b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
end
xlim([-10 10]); ylim([-10 10]);
%% compare parameters
% take top n parameters
n = 3;
% select the model
model = 9;
% seperate
ctrl_parameters = squeeze(gbz_generative{1}.top_parameters_mean(n,model,:,:));
gbz_parameters = squeeze(gbz_generative{2}.top_parameters_mean(n,model,:,:));
wash_parameters = squeeze(gbz_generative{3}.top_parameters_mean(n,model,:,:));
% collect
k = nan(ncultures(2),2,3); 
k(1:ncultures(1),:,1) = ctrl_parameters; 
k(1:ncultures(2),:,2) = gbz_parameters;
k(1:ncultures(3),:,3) = wash_parameters;
eta = nan(ncultures(2),3); 
eta(1:ncultures(1),1) = ctrl_parameters(:,1);
eta(1:ncultures(2),2) = gbz_parameters(:,1);
eta(1:ncultures(3),3) = wash_parameters(:,1);
gam = nan(ncultures(2),3); 
gam(1:ncultures(1),1) = ctrl_parameters(:,2);
gam(1:ncultures(2),2) = gbz_parameters(:,2);
gam(1:ncultures(3),3) = wash_parameters(:,2);
% statistics on wiring parmeters combinations
p = []; d = []; step = 1;
for i = 1:3;
    for j = 1:3;
        p(i,j,1) = ranksum(eta(:,i),eta(:,j));
        d(i,j,1) = computeCohen_d(eta(:,i),eta(:,j));
        p(i,j,2) = ranksum(gam(:,i),gam(:,j));
        d(i,j,2) = computeCohen_d(gam(:,i),gam(:,j));
        p = roundsd(p,3); d = roundsd(d,3);
    end
end
% form the colours
gabacol = 1-copper(20); gabacol = gabacol(10,:);
contcol = pink(20); contcol = contcol(5,:);
washcol = parula(20); washcol = washcol(4,:);
cols{1} = gabacol; cols{2} = contcol; cols{3} = washcol;
% remove cultures not to plot (group, culture) ****** <<<<<< here <<<<<< ******
remindex = [3 2];
for i = 1:size(remindex,1);
    k(remindex(2),:,remindex(1)) = [NaN NaN];
end
% visualise
h = figure;
h.Position = [100 100 1000 650];
for group = 1:ngroups;
    subplot(1,3,group);
    u = iosr.statistics.boxPlot(k(:,:,group),...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor',cols{group},...
        'boxAlpha',.6,...
        'showScatter',logical(1),...
        'scatterMarker','x',...
        'scatterColor','k');
    ylabel('Magnitude'); xticks([1,2]);
    ylim([-1.5 1.5]);
    xticklabels({'\eta','\gamma'});
    xlabel('Wiring parameter'); 
    b = gca; 
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 25;
    yline(0);
    if group == 3;
        u.handles.box(1).FaceColor = washcol;
        u.handles.box(1).FaceAlpha = .3;
        u.handles.box(2).FaceColor = washcol;
        u.handles.box(2).FaceAlpha = .3;
    end
end
% same but based on cultures
h = figure;
h.Position = [100 100 1000 650];
ylimitsg = [-1.2 1.2;-0.1 1.1];
ylabelsg = {'\eta','\gamma'};
for parameter = 1:2;
    subplot(1,2,parameter);
    uu = squeeze(k(:,parameter,:));
    a = iosr.statistics.boxPlot(uu,...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor','w',...
        'boxAlpha',1,...
        'showScatter',logical(1),...
        'scatterMarker','x',...
        'scatterColor','k');
    ylabel(ylabelsg{parameter}); 
    xticks([1:3]);
    % set colour
    a.handles.box(1).FaceColor = contcol;
    a.handles.box(1).FaceAlpha = .6;
    a.handles.box(2).FaceColor = gabacol;
    a.handles.box(2).FaceAlpha = .6;
    a.handles.box(3).FaceColor = washcol;
    a.handles.box(3).FaceAlpha = .3;
    %ylim(ylimitsg(parameter,:));
    %title(sprintf('p=%g, d=%g',p(parameter),d(parameter)));
    xticklabels({'C','G','W'});
    xlabel('Condition');
    b = gca;
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 25;
    if parameter == 1
        %ylim([-1.5 0.25]);
        ylim([-1.5 1]);
    end
    if parameter == 2;
        %ylim([0 1]);
        ylim([0 1.5]);
    end
end
%% %% run generative simulations capturing generative, integration and segregation nodal measures
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
step           = 1;
% set data to use
dataplot = rodent_data_all;
% loop between the two groups
for group = 1:ngroups;
    for culture = 1:ncultures(group); %nsamp
        tic  
        % set target networks
        Atgt_set    = dataplot{group}.table_bu;
        % seed network is empty
        A           = zeros(size(Atgt_set{cultureid{group}(culture)}));
        % compute the connections in the seed
        mseed       = nnz(A)/2;
        % euclidean distance
        D           = squareform(pdist(dataplot{group}.table_xy{cultureid{group}(culture)}));
        % set target
        Atgt        = Atgt_set{cultureid{group}(culture)};
        % set the parameter for this subject
        params      = squeeze(gbz_generative{group}.top_parameters_mean(1,model,culture,:))';
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
                disp(sprintf('Network %g connection %g of %g formed',culture,ii,m));
            end
            % keep the data
            b(:,iparam) = find(triu(A,1));
            matching_K{group,culture}     = Kall;                  
            matching_Fk{group,culture}    = Fkall;                
            matching_Ksum{group,culture}  = squeeze(sum(Kall,2)); 
            matching_Fksum{group,culture} = squeeze(sum(Fkall,2));
            matching_P{group,culture}     = Pall;                  
            matching_Psum{group,culture}  = squeeze(sum(Pall,2)); 
            matching_A{group,culture}     = Aall;                  
            segregation{group,culture}    = Seg;
            integration{group,culture}    = Int;
        end
        time = toc;
        disp(sprintf('Generative model complete: network %g took %g seconds',culture,time));
    end
end
% change name of matching probability distributions and keep
rodent_gabazine_50k_001Hz_matching_probability = matching_P;
%% plot probability distributions at half way through the simulation
% take data
matching_P = rodent_gabazine_50k_001Hz_matching_probability;
% place the cell into a list
matching_P = ...
    [matching_P(1,1:ncultures(1))';...
    matching_P(2,1:ncultures(2))';...
    matching_P(3,1:ncultures(3))'];
% visualise
u = figure; u.Position = [100 100 600 400];
% set the smoothing factor (e.g. 10)
sm = 10;
% set which to plot
plotg = 2;
% set cultures from each group to plot
control_cultures = [1:6];
gabazine_cultures = [7:15];
wash_cultures = [16:21];
% set proportion
proportions = [.01:.01:1];
% set colours
contcol = pink(2*length(proportions)); contcol = flip(contcol(round(linspace(65,65+length(proportions),length(proportions))),:));
gabacol = 1-copper(2*length(proportions)); gabacol = gabacol(4:length(proportions)+100,:);
washcol = brewermap(2*length(proportions),'purples'); washcol = washcol(round(linspace(30,30+length(proportions),length(proportions))),:);
cols = {contcol gabacol washcol};
% initialise
control_keep = []; gabazine_keep = []; wash_keep = [];
pprop = []; dprop = [];
% loop over proportions
for j = 1:length(proportions);
    proportion = proportions(j);
    % set scale limits
    lim = 1000;
    % take the mid point distributions
    distribution = {}; distribution_scaled = {}; distribution_scaled_sized = [];
    for culture = 1:sum(ncultures);
        A = matching_P{culture};
        t = size(A,1);
        n = size(A,2);
        B = squeeze(A(round(t*proportion),:,:));
        C = B.*triu(ones(n),1);
        %distribution{culture} = C(C>0); % removed non zeros
        distribution{culture} = C;
        distribution_scaled{culture} = imresize(distribution{culture},[lim 1]);
        distribution_scaled_array(culture,:) = distribution_scaled{culture}; % form an array
        distribution_scaled_sized(culture,:) = rescale(distribution_scaled{culture},0,1); % note this may artificially elevate probabilities as zeros elevate depending on node size
    end
    % state which to plot 
    dataplot = distribution_scaled_sized;
    % remove specific cultures here *******<<<<<<<****** manual
    % average across groups
    control_mean = mean(dataplot(control_cultures,:));
    gabazine_mean = mean(dataplot(gabazine_cultures,:));
    wash_mean = mean(dataplot([16 18:21],:));
    % keep these
    control_keep(j,:) = control_mean;
    gabazine_keep(j,:) = gabazine_mean;
    wash_keep(j,:) = wash_mean;
    % plot ksdensity
    % control
    if plotg == 1;
        [fy y] = ksdensity(control_mean);
        u = plot(y,smooth(fy,sm),'linewidth',3); u.Color = contcol(j,:);
    end
    if plotg == 2;
        [f x] = ksdensity(gabazine_mean);
        u = plot(x,smooth(f,sm),'linewidth',3); u.Color = gabacol(j,:);
    end
    if plotg == 3;
        [f z] = ksdensity(wash_mean);
        u = plot(z,smooth(f,sm),'linewidth',3); u.Color = washcol(j,:);
    end
    hold on;
    xlabel('Probability score ({\itP_i_j})'); ylabel('Frequency');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
    xlim([0 1]); yticks([]); box off; ylim([0 12]);
    c = colorbar; c.Ticks = [0:0.1:1]; c.TickLabels = {'0%','','','','','50%','','','','','100%'};  caxis([0 1]); 
    c.Label.String = 'Simulation time'; c.Label.FontSize = 25; c.Label.FontName = 'Arial'; colormap(cols{plotg});
    set(gcf,'Position',[100 100 1000 450]);
end
%% compute a visualisation between controls, gabazine and washout
% average over all
mean_control_mean = mean(control_keep);
mean_gabazine_mean = mean(gabazine_keep);
mean_wash_mean = mean(wash_keep);
% take std
std_control_mean = std(control_keep);
std_gabazine_mean = std(gabazine_keep);
std_wash_mean = std(wash_keep);
% new colours
gabacol = [0.6053 0.7533 0.8429];
contcol = [0.7375 0.4588 0.4588];
washcol = [0.5514 0.5393 0.7533];
% statistical test
stats = [mean_control_mean' mean_gabazine_mean' mean_wash_mean'];
p = []; d = [];
for i = 1:3;
    for j = 1:3;
        [~,p(i,j)] = ttest(stats(:,i),stats(:,j));
        d(i,j) = computeCohen_d(stats(:,i),stats(:,j));
    end
end
% set figure density
u = figure; u.Position = [100 100 800 650];
% plot ksdensity 
[f x] = ksdensity(mean_gabazine_mean);
h = plot(x,smooth(f,sm),'linewidth',6); h.Color = gabacol;
hold on;
[fy y] = ksdensity(mean_control_mean);
h = plot(y,smooth(fy,sm),'linewidth',6); h.Color = contcol;
hold on;
[fy z] = ksdensity(mean_wash_mean);
h = plot(z,smooth(fy,sm),'linewidth',6); h.Color = washcol;
xlabel('Probability score ({\itP_i_j})'); ylabel('Frequency'); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlim([0 inf]); 
yticks([]); box off;
% plot legend
% legend(sprintf('p=%.3g, \nd=%.3g',p,d),'box','off');
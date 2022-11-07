%% figure 7
% written by danyal akarca
%% load data
clear; clc;
% set directory
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022/';
% load sttc
load(strcat(gendir,'human_ipsc_gn/GN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); gn = gm_data;
load(strcat(gendir,'human_ipsc_mn/MN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat'));  mn = gm_data;
load(strcat(gendir,'human_ipsc_dn/DN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat'));  dn = gm_data;
ipsc_data_all = {gn mn dn};
% human COs
load(strcat(gendir,'human_organoids/M03912_GNM_data_min_rate_0.01hz_0.001_alpha_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); hco_data_all = gm_data;
% load generative model data
load(strcat(gendir,'/human_ipsc_gn/1800s_min_rate_001hz_lag10_jitter10_prob_re-analysis/combined/results_generative_models.mat')); gn = results_generative_models;
load(strcat(gendir,'/human_ipsc_mn/1800s_min_rate_001hz_lag10_jitter10_prob_re-analysis/combined/results_generative_models.mat')); mn = results_generative_models;
load(strcat(gendir,'/human_ipsc_dn/1800s_min_rate_001hz_lag10_jitter10_prob_re-analysis/combined/results_generative_models.mat')); dn = results_generative_models;
load(strcat(gendir,'/human_organoids/211231_t075/M03912_GNM_rec_1800s_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/combined/results_generative_models.mat')); hco = results_generative_models;
data = {gn mn dn hco};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
% addpath of Cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
% addpath of round
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
%% remove nodes with no connections from empirical data
%%% ipsc %%%
% set number of cultures
ncultures = [8 7 6];
% set number of dgroups
ngroups = 3;
% set number of models
nmodels = 13;
% set number of parameters
nparams = 20000;
% label 3 groups
group_labels = string({'GN','MN','DN'});
% initialise
ipsc_data = ipsc_data_all;
% loop over networks and remove nodes with no connections
for group = 1:ngroups;
    for culture = 1:ncultures(group);
        % get network
        Atgt = ipsc_data_all{group}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = ipsc_data_all{group}.table_wu{culture};
        xy = ipsc_data_all{group}.table_xy{culture};
        % keep only nodes with edges
        ipsc_data{group}.table_bu{culture} = Atgt(ids,ids);
        ipsc_data{group}.table_wu{culture} = wu(ids,ids);
        ipsc_data{group}.table_xy{culture} = xy(ids,:);
    end
end
%%% hco %%%
% set number of cultures
ncultures = 6;
% set number of models
nmodels = 13;
% set number of parameters
nparams = 20000;
% initialise
hco_data = hco_data_all;
% loop over networks and remove nodes with no connections
for culture = 1:ncultures;
    % get network
    Atgt = hco_data_all.table_bu{culture};
    % find nodes with edges
    ids = find(degrees_und(Atgt));
    % get other data
    wu = hco_data_all.table_wu{culture};
    xy = hco_data_all.table_xy{culture};
    % keep only nodes with edges
    hco_data.table_bu{culture} = Atgt(ids,ids);
    hco_data.table_wu{culture} = wu(ids,ids);
    hco_data.table_xy{culture} = xy(ids,:);
end
%% visualise the networks
% ipsc
% set data to plot
dataplot = ipsc_data;
% set group
group = 3;
% select a culture e.g. GN=6, MN=3, DN=5
culture = 5;
% compute the cultures
network = dataplot{group}.table_bu{culture};
coords = dataplot{group}.table_xy{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 400 400]);
% form the graph
g = graph(network);
% other nice colors
colbar = 1-summer(4);
% plot
h = plot(g,...
    'XData',coords(:,1),...
    'YData',coords(:,2),...
    'EdgeColor',[130 130 130]./256,...
    'LineWidth',1,...
    'NodeColor',colbar(group,:),...
    'NodeLabel',[],...
    'MarkerSize',.6*degrees_und(network)+0.01);
xticks([]); yticks([]);
set(gca,'visible','off');
%%% hco %%%
% set data to plot
dataplot = hco_data;
% select a culture e.g. hCO 2/6
culture = 6;
% compute the cultures
network = dataplot.table_bu{culture};
coords = dataplot.table_xy{culture};
% calculate the distance
d = squareform(pdist(coords));
% visualise over divs
figure;
set(gcf,'Position',[100 100 400 400]);
% form the graph
g = graph(network);
% plot
h = plot(g,...
    'XData',coords(:,1),...
    'YData',coords(:,2),...
    'EdgeColor',[130 130 130]./256,...
    'LineWidth',1,...
    'NodeColor',[10 153 153]./256,...
    'NodeLabel',[],...
    'MarkerSize',.6*degrees_und(network)+0.01);
xticks([]); yticks([]);
set(gca,'visible','off');
%% plot topological comparisons between human ipscs
% set which data to plot
dataplot = ipsc_data_all;
% set number of cultures
ncultures = [8 7 6];
% set number of dgroups
ngroups = 3;
% number of measures
nmeasures = 11;
% initialise
statistics = nan(ncultures(1),ngroups,nmeasures);
% statistics labels
statistics_labels = string({...
    'Number of nodes',...
    'Mean STTC',...
    'Network density (%)',...
    'Total degree',...
    'Efficiency',...
    'Betweenness',...
    'Global clustering',...
    'Edge length',...
    'Modularity',...
    'Matching',...
    'Small-worldness'});
% loop over cultures
for group = 1:ngroups;
    for culture = 1:ncultures(group);
        % take the weighted network
        wnetwork = dataplot{group}.table_wu{culture};
        % take the binary network
        network = dataplot{group}.table_bu{culture};
        % take the coordinates
        coords = dataplot{group}.table_xy{culture};
        % compute the euclidean distance
        euclidean = squareform(pdist(coords));
        % compute global statistics
        % nnodes
        statistics(culture,group,1) = size(network,1);
        % total STTC
        sttc = triu(wnetwork,1);
        statistics(culture,group,2) = mean(sttc(:),'all');
        % density
        statistics(culture,group,3) = density_und(network)*100;
        % total number of connections
        statistics(culture,group,4) = sum(degrees_und(network));
        % global efficiency
        statistics(culture,group,5) = efficiency_bin(network,0);
        % betweenness
        statistics(culture,group,6) = mean(betweenness_bin(network));
        % clustering
        statistics(culture,group,7) = mean(clustering_coef_bu(network));
        % average weighted length
        statistics(culture,group,8) = mean(mean(euclidean.*network));
        % modularity
        [~,statistics(culture,group,9)] = modularity_und(network);
        % homophily
        statistics(culture,group,10) = mean(matching_ind(network),'all');
        % small-worldness
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
        for j = 1:1000;
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
        statistics(culture,group,11) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
        % display
        disp(sprintf('computed statistics for culture %g div %s',...
            culture,group_labels(group)));
    end
end
%%% plot statistics over time %%%
% stat, e.g. mean sttc 2, density 3, small-worldness 11
stat = 3;
% visualise 
k = squeeze(statistics(:,:,stat));
% plot
h = figure; h.Position = [100 100 400 400];
iosr.statistics.boxPlot(k,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'boxColor',[.5 .5 .5],...
    'boxAlpha',.5,...
    'symbolColor','k',...
    'showScatter',logical(1),...
    'scatterColor','k');
xticklabels(group_labels);
xlabel('Human iPSCs');
ylabel(statistics_labels(stat));
b = gca; b.TickDir = 'out';
b.FontName = 'Arial'; 
b.FontSize = 25;
box off;
% tabulate
med_tab = median(k);
iqr_tab = iqr(k);
stat_tab = table(med_tab,iqr_tab);
% run a statistical test
p = anova1(k,[1 2 3],'off');
% columns are divs, rows are pipelines
disp(stat_tab);
%%% plot proportion of edge lengths over time %%%
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binsum = []; binprop = [];
% compute
for group = 1:ngroups;
    for culture = 1:ncultures(group);
    % network
    a = dataplot{group}.table_bu{culture}; 
    % euclidean
    b = squareform(pdist(dataplot{group}.table_xy{culture}));
    % nnode
    nnode = size(a,1);
    % all possible connections
    pos = (nnode*nnode-1)/2;
    % sum across bin limits
    for k = 2:length(binvec);
        % filter by length
        u = a.*(b<binvec(k)&b>binvec(k-1));
        % compute proportion
        prop = sum(u,'all')/pos;
        % keep
        binprop(culture,group,k-1) = prop;
    end
    end
end
% compute mean and std
x = 10*squeeze(mean(binprop))*100;
y = 10*squeeze(std(binprop))*100;
% form a jitter
j = 50;
% set colorbars
colbar = 1-summer(4);
% visualise
h = figure; h.Position = [100 100 400 400];
for group = 1:ngroups;
    uu = errorbar(binvec(2:end)+j*group,x(group,:),y(group,:)); 
    uu.Color = colbar(group,:); uu.Color(4) = .2;
    uu.Marker = '.'; 
    uu.MarkerSize = 60;
    uu.LineWidth = 5;
    hold on;
end
xlim([0 binlim+2*j]);
ylim([0 inf]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Distance (\mum)'); ylabel('Prop. connections (%)');
legend(group_labels,'box','off'); box off;
%% plot topology of human cerebral organoids
% set data to plot
dataplot = hco_data_all;
% set number of cultures
ncultures = 6;
% number of measures
nmeasures = 10;
% initialise
statistics = nan(ncultures,nmeasures);
% statistics labels
statistics_labels = string({...
    'Number of nodes',...
    'Mean STTC',...
    'Network density (%)',...
    'Total degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Edge length',...
    'Modularity',...
    'Matching',...
    'Small-worldness'});
% loop over cultures
for culture = 1:ncultures;
     % take weighted network
     wnetwork = dataplot.table_wu{culture};
     % take the binary network
     network = dataplot.table_bu{culture};
     % take the coordinates
     coords = dataplot.table_xy{culture};
     % compute the euclidean distance
     euclidean = squareform(pdist(coords));
     % number of nodes
     statistics(culture,1) = size(network,1);
     % compute global statistics
     sttc = triu(wnetwork,1);
     statistics(culture,2) = mean(sttc(:),'all');
     % density
     statistics(culture,3) = density_und(network)*100;
     % total number of connections
     statistics(culture,4) = sum(degrees_und(network));
     % global efficiency
     statistics(culture,5) = efficiency_bin(network,0);
     % betweenness
     statistics(culture,6) = mean(betweenness_bin(network));
     % clustering
     statistics(culture,7) = mean(clustering_coef_bu(network));
     % average weighted length
     statistics(culture,8) = mean(mean(euclidean.*network));
     % modularity
     [~,statistics(culture,9)] = modularity_und(network);
     % homophily
     statistics(culture,10) = mean(matching_ind(network),'all');
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
     for j = 1:1000;
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
     statistics(culture,11) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
     % display
     disp(sprintf('computed statistics for culture %g hCO',...
         culture));
end
%%% plot statistics %%%
% stat
stat = 2;
% visualise the number of connections
k = squeeze(statistics(:,stat));
% plot
h = figure; h.Position = [100 100 400 400];
iosr.statistics.boxPlot(k,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'boxColor',[.5 .5 .5],...
    'boxAlpha',.5,...
    'symbolColor','k',...
    'showScatter',logical(1),...
    'scatterColor','k');
xticks([]);
ylabel(statistics_labels(stat));
b = gca; b.TickDir = 'out';
b.FontName = 'Arial'; 
b.FontSize = 25;
b = gca; b.TickDir = 'out';
b.FontName = 'Arial';
y = b.YLim;
ylim([y(1)-0.05*y(1) y(2)]);
xlabel('hCOs','FontSize',25);
%%% plot connections over distance binned %%%
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binsum = []; binprop = [];
% compute
for culture = 1:ncultures;
    % network
    a = dataplot.table_bu{culture}; 
    % euclidean
    b = squareform(pdist(dataplot.table_xy{culture}));
    % nnode
    nnode = size(a,1);
    % all possible connections
    pos = (nnode*nnode-1)/2;
    % sum across bin limits
    for k = 2:length(binvec);
        % filter by length
        u = a.*(b<binvec(k)&b>binvec(k-1));
        % compute proportion
        prop = sum(u,'all')/pos;
        % keep
        binprop(culture,k-1) = prop;
    end
end
% compute mean and std
x = squeeze(mean(binprop))*100;
y = squeeze(std(binprop))*100;
% form a jitter
j = 50;
% set colorbars
colbar = [10 153 153]./256;
% visualise
h = figure; h.Position = [100 100 400 400];
uu = errorbar(binvec(2:end)+j,x,y); 
uu.Color = colbar; uu.Color(4) = .2;
uu.Marker = '.'; 
uu.MarkerSize = 50;
uu.LineWidth = 5;
hold on;
xlim([0 binlim+2*j]);
ylim([0 inf]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
xlabel('Distance (\mum)'); ylabel('Prop. connections (%)'); xticks([0:2000:4000]);
%% cellular properties
%% generative models
% number of groups
ngroups = 4;
% ncultures
ncultures = [8 7 6 6];
% number of models
nmodels = 13;
% number of params
nparams = 20000;
% form energy matrix
energy = nan(ncultures(1),nmodels,nparams,ngroups);
parameters = nan(ncultures(1),nmodels,nparams,2,ngroups);
% place in the divs
for group = 1:ngroups;
    energy(1:ncultures(group),:,:,group) = data{group}.energy;
    parameters(1:ncultures(group),:,:,:,group) = data{group}.parameters;
end
%% compute top performing simulations
% define the top n parameters
nset = [1 10 50 100];
% initialise
top_energy = cell(length(nset),ngroups);
top_energy_mean = zeros(length(nset),13,ncultures(1),ngroups);
top_parameters = cell(length(nset),ngroups);
top_parameters_mean = zeros(length(nset),13,ncultures(1),2,ngroups);
% compute the minimum
for no = 1:length(nset);
    % take the actual amount of top performing parameters
    n = nset(no);
    % loop over divs
    for group = 1:ngroups;
        % loop over models
        for model = 1:nmodels;
            % take energies for this model
            pipeline_d = squeeze(energy(:,model,:,group))';
            % rank them for each subject
            [v i] = sort(pipeline_d);
            % take top n energies and their indices
            n_e = v(1:n,:);
            n_i = i(1:n,:);
            % take the corresponding parameters
            u = zeros(ncultures(group),n,2);
            for s = 1:ncultures(group);
                % keep parameters
                u(s,:,:) = squeeze(parameters(s,model,n_i(:,s),:,group));
            end
            % if top parameter only
            if n == 1
                % squeeze the matrices
                u = squeeze(u);
                % assign
                top_energy{no,group}(model,:) = n_e';
                top_parameters{no,group}(model,:,:) = u;
                % and assign it to the mean
                top_energy_mean(no,model,:,group) = n_e';
                top_parameters_mean(no,model,1:ncultures(group),:,group) = u;
            else
                top_energy{no,group}(model,:,:) = n_e';
                top_parameters{no,group}(model,:,:,:) = u;
                % keep a mean value too
                top_energy_mean(no,model,:,group) = squeeze(mean(n_e',2));
                top_parameters_mean(no,model,1:ncultures(group),:,group) = squeeze(mean(u,2));
            end
        end
    end
    disp(sprintf('set %g of %g complete',no,length(nset)));
end
%% visualise energy over groups
% set the top n parameters
ncomb = 1;
% visualsie all models
data_plot = squeeze(top_energy_mean(ncomb,:,:,:));
% set order
i = [2:3 9:13 4:8 1];
data_plot = data_plot(i,:,:);
% permute the order
data_plot = permute(data_plot,[2 3 1]);
% iosr boxplot - all models
h = figure; h.Position = [100 100 1800 600];
u = iosr.statistics.boxPlot(data_plot,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','r','b','b','b','b','b','g','g','g','g','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]}';
ylim([0 0.75]); ylabel('Energy'); xlabel('');
xticklabels({'GNs','MNs','DNs','hCOs'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% iosr boxplot - grouped models
energy_group_rules = nan(ncultures(1)*5,ngroups,4);
a = squeeze(data_plot(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); energy_group_rules(1:ncultures*2,:,1) = a';
b = squeeze(data_plot(:,:,[3:7])); b = permute(b,[2 1 3]); b = b(:,:); energy_group_rules(1:ncultures*5,:,2) = b';
c = squeeze(data_plot(:,:,[8:12])); c = permute(c,[2 1 3]); c = c(:,:); energy_group_rules(1:ncultures*5,:,3) = c';
d = squeeze(data_plot(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); energy_group_rules(1:ncultures,:,4) = d';
% visualise
h = figure; h.Position = [100 100 1800 600];
u = iosr.statistics.boxPlot(energy_group_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]}';
ylim([0 0.75]); ylabel('Energy'); xlabel('');
xticklabels({'GNs','MNs','DNs','hCOs'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute energy statistics between rules for each time point
% set the data
data = energy_group_rules;
% ndivs
ndivs = 4;
% logical for saving
saveTable = 1;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/human_100k_energy_rules_all.csv';
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
for div = 1:ndivs;
    p_rules(:,div) = compare{div}(:,6);
end
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
    sprintf('GNs p=%.3g',p(1)),'GNs Cohen d',...
    sprintf('MNs p=%.3g',p(2)),'MNs Cohen d',...
    sprintf('DNs p=%.3g',p(3)),'DNs Cohen d',...
    sprintf('hCOs p=%.3g',p(4)),'hCOs Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
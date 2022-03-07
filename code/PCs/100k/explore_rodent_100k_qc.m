%% Exploration of 100k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortical neurons
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_100k_tracking');
a = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
data = {a.gm_data b.gm_data c.gm_data};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
%% form labels
% 3 pipelines
pipeline_labels = string({'5ms','10ms','20ms'})';
% 4 time points
div_labels = string({'14','21','28'});
% 12 hd cultures
culture_labels = string(1:12);
%% set hyperparameters of the dataset
% number of pipelines
npipelines = 3;
% number of cultures
ncultures = 12;
% number of divs
ndivs = 3;
%% visualise the networks
% select a pipeline
pipeline = 2;
% select a culture
culture = 2;
% set divs to explore
divs = [1 2 3];
% compute the cultures
networks = data{pipeline}.table_bu_shared{culture};
coords = data{pipeline}.table_xy_shared{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 1200 300]);
sgtitle(sprintf('culture %s: binary networks',culture_labels(culture)));
% loop through divs
for i = 1:length(divs);
    div = divs(i);
    subplot(1,length(divs),i);
    g = graph(networks{div});
    h = plot(g,...
        'XData',coords{div}(:,1),...
        'YData',coords{div}(:,2),...
        'EdgeColor',[130 130 130]./256,...
        'LineWidth',1.5,...
        'NodeColor',[239 152 170]./256,...
        'NodeLabel',[],...
        'MarkerSize',0.5*degrees_und(networks{div})+0.01);
    title(div_labels(i));
    xticks([]); yticks([]);
    grid on;
end
% visualise euclidean distances over divs
figure;
set(gcf,'Position',[100 100 1200 300]);
sgtitle(sprintf('culture %s: euclidean distance',culture_labels(culture)));
% loop through divs
for i = 1:length(divs);
    subplot(1,length(divs),i);
    % compute euclidean distances
    imagesc(squareform(pdist(coords{div}))); caxis([0 5e3]);
    title(div_labels(i));
    xticks([]); yticks([]);
end
%% visualise across all cultures 
% set hyperparameters 
pipeline = 2;
divs = [1 2 3];
% loop over cultures
h = figure;
h.Position = [100 100 600 300];
m = [];
for culture = 1:12;
    network = data{pipeline}.table_bu_shared{culture}(divs);
    coord = data{pipeline}.table_xy_shared{culture}(divs);
    % loop over divs. Note, this is not really needed as euclidean is the same over divs
    for div = 1:length(divs);
        % get coordinates
        a = coord{div};
        % make eucliean
        a = squareform(pdist(a));
        % get size
        n = size(a,1);
        % take single edges
        u = a(find(triu(ones(n),1)));
        % make into a ks density plot
        [f xi] = ksdensity(u);
        % plot and shift to right
        h = plot(xi + 500,f,'linewidth',2);
        m(culture,div) = max(xi) + 500;
        % keep plotting
        hold on;
    end
    hold on;
end
% keep max euclidean for later
ma = max(m,[],'all');
% set limits
xlim([0 5000]);
yticks([]); u = gca; u.TickDir = 'out'; u.FontName = 'Arial';
xlabel('Euclidean \mum'); ylabel('Frequency');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; box off; b.TickLength = [0 0];
% Plot the eucldiean distances of each culture
h = figure;
h.Position = [100 100 1300 400];
n = [];
for culture = 1:12;
    % subplot
    subplot(2,6,culture);
    % take the coordinates
    coords = data{pipeline}.table_xy_shared{culture}{1};
    % take the euclidean
    d = squareform(pdist(coords));
    % keep size
    n(culture) = size(d,1);
    % plot the euclidean matrix
    imagesc(d);
    % plot all the same caxis to the max
    caxis([0 5000]);
    xticks([]); yticks([]);
end
%% visualise a example network
% take the top 20 nodes of 10 rodent network
culture = 10;
% take network and trim it
A = data001Hz{2}.table_bu_shared{culture}{1};
C = data001Hz{2}.table_xy_shared{culture}{1};
D = squareform(pdist(C));
% indexs
ind = 1:size(A,1);
% remove nodes 
rem = [33 34 1 2 14 17 11 20 21 22 24 28];
% take new ones
ind(rem)=[];
% update
A = A(ind,ind);
C = C(ind,:);
D = D(ind,ind);
% calculate statistics
k = degrees_und(A);
c = clustering_coef_bu(A);
b = betweenness_bin(A);
d = sum(triu(A.*D,1));
stats = [k;c';b;d];
% cols = [.5 .8 .9;.01 .5 .3;1 .5 .5; .9 .7 .1]; % old one
cols = [83, 62, 133; 72 143 177; 79 211 196; 193 248 207]./256;
scalar = [4 25 .15 .005];
% plot the euclidean matrix
h = figure; h.Position = [100 100 300 300];
imagesc(D); xticks([]); yticks([]); axis off;
% plot graphs
g = graph(A);
h = figure; h.Position = [100 100 300 300];
h = plot(g,...
    'XData',C(:,1),...
    'YData',C(:,2),...
    'markersize',10,...
    'nodecolor',[.3 .3 .3],...
    'nodelabel',[],...
    'edgecolor',[.7 .7 .7],...
    'linewidth',2);
xticks([]); yticks([]); axis off;
% now with the statistics
h = figure; h.Position = [100 100 1200 200];
for stat = 1:4;
subplot(1,4,stat);
h = plot(g,...
    'XData',C(:,1),...
    'YData',C(:,2),...
    'markersize',scalar(stat)*stats(stat,:)+0.01,...
    'nodecolor',cols(stat,:),...
    'nodelabel',[],...
    'edgecolor',[.7 .7 .7],...
    'linewidth',2);
xticks([]); yticks([]); axis off;
end
% set nbins
nbins = 8;
% plot distributions
h = figure; h.Position = [100 100 1200 200];
subplot(1,4,1); y=histogram(k,nbins,'facecolor',cols(1,:)); 
u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14; ylabel('Frequency'); xlabel('Degree, k'); box off; ylim([0 16]);
subplot(1,4,2); histogram(c,nbins,'facecolor',cols(2,:)); 
u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14; xlabel('Clustering, c'); ylim([0 25]); box off;
subplot(1,4,3); histogram(b,nbins,'facecolor',cols(3,:)); 
u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14; xlabel('Betweenness, b'); ylim([0 25]); box off;
subplot(1,4,4); histogram(d,nbins,'facecolor',cols(4,:)); 
u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14; xlabel('Edge length, d'); ylim([0 16]); box off;
% form a random network
Asynth = randmio_und(A,3);
% compute examples based on this
ind = find(Asynth);
Azero = zeros(size(A));
Aone  = Azero; Aone(ind(1))=1; Aone = Aone+Aone';
Atwo  = Azero; Atwo(ind([1 2]))=1; Atwo = Atwo+Atwo';
Athree = Azero; Athree(ind([1 2 4]))=1; Athree = Athree+Athree';
% set and plot
g = graph(Asynth);
h = figure; h.Position = [100 100 300 300];
h = plot(g,...
    'XData',C(:,1),...
    'YData',C(:,2),...
    'markersize',10,...
    'nodecolor',[.3 .3 .3],...
    'nodelabel',[],...
    'edgecolor',[.4 .1 .8],...
    'linewidth',2);
xticks([]); yticks([]); axis off;
%% compute summary statistics over time
% number of measures
nmeasures = 9;
% initialise
statistics = zeros(3,12,3,nmeasures);
% statistics labels
statistics_labels = string({...
    'Number of nodes',...
    'Density',...
    'Degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Edge length',...
    'Modularity',...
    'Matching'});
% loop over networks
for pipeline = 1:3;
    for culture = 1:12;
        for div = 1:3;
            % take the binary network
            network = data{pipeline}.table_bu_shared{culture}{div};
            % take the coordinates
            coords = data{pipeline}.table_xy_shared{culture}{div};
            % compute the euclidean distance
            euclidean = squareform(pdist(coords));
            % compute global statistics
            % number of nodes
            statistics(pipeline,culture,div,1) = size(network,1);
            % density
            statistics(pipeline,culture,div,2) = density_und(network);
            % total number of connections
            statistics(pipeline,culture,div,3) = sum(degrees_und(network));
            % global efficiency
            statistics(pipeline,culture,div,4) = efficiency_bin(network,0);
            % betweenness
            statistics(pipeline,culture,div,5) = mean(betweenness_bin(network));
            % clustering
            statistics(pipeline,culture,div,6) = mean(clustering_coef_bu(network));
            % average weighted length
            statistics(pipeline,culture,div,7) = mean(mean(euclidean.*network));
            % modularity
            [~,statistics(pipeline,culture,div,8)] = modularity_und(network);
            % homophily
            statistics(pipeline,culture,div,9) = mean(matching_ind(network),'all');
        end
    end
end
%% plot number of connections and densities by pipeline and div
% stat
stat = 1;
% visualise the number of connections
k = squeeze(statistics(:,:,:,stat));
% plot
figure;
iosr.statistics.boxPlot(permute(k,[2 3 1]),...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k');
xticklabels({'14 days','21 days','28 days'});
xlabel('Days in vitro');
ylabel(statistics_labels(stat));
u = gca; u.TickDir = 'out';
u.FontName = 'Arial';
% tabulate
med_tab = squeeze(median(k,2));
iqr_tab = squeeze(iqr(k,2));
stat_tab = table(med_tab,iqr_tab);
% columns are divs, rows are pipelines
%% visualise statistics over time
% set statistics
stat = 6;
% set limits for each
limits = [0,0.35;0,35;0.15,0.65;20,180;0.15,0.7;0,450];
% compute for the statistic
single_statistic = squeeze(statistics(:,:,:,stat));
% set colours
colour = 1-jet(400);
col = [colour(100,:); colour(200,:); colour(300,:)];
% visualise 
figure;
set(gcf,'Position',[100 100 900 300]);
% xaxis
xaxis = [14 21 28];
for pipeline = 1:3;
    % all cultures
    subplot(1,2,1);
    plot(xaxis,squeeze(single_statistic(pipeline,:,:))',...
        'LineWidth',1,...
        'Color',col(pipeline,:));
    hold on; grid on;
    xticks([14:28]); ylim(limits(stat,:));
    ylabel(statistics_labels(stat)); xlabel('days'); title('all cultures');
    % averaged over cultures
    subplot(1,2,2);
    plot(xaxis,squeeze(mean(squeeze(single_statistic(pipeline,:,:)),1))',...
   'LineWidth',2,...
    'Color',col(pipeline,:)); 
    hold on; grid on;
    xticks([14:28]); ylim(limits(stat,:));
    ylabel(statistics_labels(stat)); xlabel('days'); title('averaged over cultures');
end
% on top of each other
figure;
set(gcf,'Position',[100 100 400 300]);
% make minor colours
minor_col = col;
minor_col(:,4) = [0.15;0.15;0.15];
for pipeline = 1:3;
    % all cultures
    h = plot(xaxis,squeeze(single_statistic(pipeline,:,:))',...
        'LineWidth',1,...
        'Color',minor_col(pipeline,:));
    hold on; grid on;
    % averaged over cultures
    plot(xaxis,squeeze(mean(squeeze(single_statistic(pipeline,:,:)),1))',...
   'LineWidth',3.5,...
    'Color',col(pipeline,:)); 
    hold on; grid on;
    xticks([14:28]); ylim(limits(stat,:));
    ylabel(statistics_labels(stat),'FontSize',15); xlabel('Days','FontSize',15);
end
legend(pipeline_labels);
set(gca, 'FontName', 'Arial')
u = gca; u.TickDir = 'out';
%% correlate statistics between pipelines
% visualise correlation matrices across sttc parameters
h = figure; h.Position = [100 100 900 400];
for stat = 1:nmeasures;
    subplot(2,4,stat);
    x = statistics(:,:,:,stat);
    y = x(:,:);
    r = corr(y');
    imagesc(r); caxis([-1 1]);
    xticklabels(pipeline_labels);
    yticks([]);
    title(statistics_labels(stat));
    u = gca; u.TickDir = 'out';
    u.FontName = 'Arial';
end
% visualise a specific statistic
stat = 2;
% visualise
x = statistics(:,:,:,stat);
y = x(:,:);
h = figure; h.Position = [100 100 1000 250];
sgtitle(statistics_labels(stat));
subplot(1,3,1); scatter(y(1,:),y(2,:)); xlabel('5ms'); ylabel('10ms'); grid on; r = corr(y(1,:)',y(2,:)'); title(sprintf('r=%.3g',r));
subplot(1,3,2); scatter(y(2,:),y(3,:)); xlabel('10ms'); ylabel('20ms'); grid on; r = corr(y(2,:)',y(3,:)'); title(sprintf('r=%.3g',r));
subplot(1,3,3); scatter(y(1,:),y(3,:)); xlabel('5ms'); ylabel('20ms'); grid on; r = corr(y(1,:)',y(3,:)'); title(sprintf('r=%.3g',r));
%% histogram of statistics for each div
% set statistics
stat = 1;
% compute for the statistic
single_statistic = squeeze(statistics(:,:,:,stat));
% set colours
colour = 1-jet(400);
col = [colour(100,:); colour(200,:); colour(300,:)];
% visualise 
figure;
set(gcf,'Position',[100 100 600 300]);
for pipeline = 1:3;
    % all cultures
    h = histogram(squeeze(single_statistic(pipeline,:,:))',...
        'LineWidth',2,...
        'EdgeColor',col(pipeline,:));
    hold on; grid on;
    ylabel(sprintf('F(%s)',statistics_labels(stat)));
    xlabel(statistics_labels(stat));
end
%% compute the weight versus distance
% select a pipeline
pipeline = 1;
% select a culture
culture = 1;
% compute the cultures
networks_bu = data{pipeline}.table_bu_shared{culture};
networks_wu = data{pipeline}.table_wu_shared{culture};
coords = data{pipeline}.table_xy_shared{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 1400 300]);
sgtitle(sprintf('pipeline %s, culture %s: weight and distance relationships',...
    pipeline_labels(pipeline),culture_labels(culture)));
% loop through divs
for i = 1:4;
    % take weights of "extant" edges according to this pipeline
    u = networks_wu{i}.*networks_bu{i};
    ext_wu = u(u>0);
    % compute euclidean distance
    eu = squareform(pdist(coords{i}));
    ext_eu = eu(u>0);
    % plot the scatter
    subplot(1,4,i);
    scatter(ext_wu,ext_eu,30,'o','markerfacecolor',[50 50 50]./256,'markerfacealpha',0.2);
    % quote the correlation
    [r p] = corr(ext_wu,ext_eu);
    title(sprintf('%s: r=%.2g, p=%.2g',div_labels(i),r,p));
end
figure;
set(gcf,'Position',[100 100 1400 300]);
sgtitle(sprintf('pipeline %s: weight and distance relationships',...
    pipeline_labels(pipeline)));
% across all cultures
ext_wu_keep = [];
ext_eu_keep = [];
for i = 1:4
    for culture = 1:12
        % compute the cultures
        networks_bu = data{pipeline}.table_bu_shared{culture};
        networks_wu = data{pipeline}.table_wu_shared{culture};
        coords = data{pipeline}.table_xy_shared{culture};
        % take weights of "extant" edges according to this pipeline
        u = networks_wu{i}.*networks_bu{i};
        ext_wu = u(u>0);
        % compute euclidean distance
        eu = squareform(pdist(coords{i}));
        ext_eu = eu(u>0);
        % keep the data
        ext_wu_keep = [ext_wu_keep; ext_wu];
        ext_eu_keep = [ext_eu_keep; ext_eu];
    end
    ext_wu_keep_div{i} = ext_wu_keep;
    ext_eu_keep_div{i} = ext_eu_keep;
    % plot the scatter
    subplot(1,4,i);
    scatter(ext_wu_keep_div{i},ext_eu_keep_div{i},30,...
        '.','markerfacecolor',[50 50 50]./256,'markerfacealpha',0.2);
    % quote the correlation
    [r p] = corr(ext_wu_keep_div{i},ext_eu_keep_div{i});
    title(sprintf('%s: r=%.2g, p=%.2g',div_labels(i),r,p));
end
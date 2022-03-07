%% Exploration of human cerebral organoid networks
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load cerebral organoid data
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_cerebral_organoids_30min/')
% load sttc data
a = load('M03912_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
b = load('M03912_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
data = {a.gm_data b.gm_data};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
%% hyperparameters of the dataset
% number of pipelines
npipelines = 2;
% pipeline labels
pipeline_labels = string({'10ms','20ms'});
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
        all_o_networks{step} = data{pipeline}.table_bu{culture};
        all_d_networks{step} = squareform(pdist(data{pipeline}.table_xy{culture}));
        % update the step
        step = step + 1;
    end
end 
%% visualise the network
% select a pipeline
pipeline = 1;
% select a culture - culture 2 or 6 for visualisation
culture = 6;
% compute the cultures
network = data{pipeline}.table_bu{culture};
coords = data{pipeline}.table_xy{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 400 400]);
% form the graph
g = graph(network);
h = plot(g,...
    'XData',coords(:,1),...
    'YData',coords(:,2),...
    'EdgeColor',[130 130 130]./256,...
    'LineWidth',1,...
    'NodeColor',[10 153 153]./256,...
    'NodeLabel',[],...
    'MarkerSize',0.5*degrees_und(network)+0.01);
xticks([]); yticks([]);
set(gca,'visible','off');
%% visualise across all cultures 
% set hyperparameters 
pipeline = 1;
ncultures;
% initialise
h = figure;
h.Position = [100 100  600 300];
m = [];
% loop over cultures
for culture = 1:ncultures;
    % take the network
    network = data{pipeline}.table_bu{culture};
    % take the coordinates of this network
    coord = data{pipeline}.table_xy{culture};
    a = squareform(pdist(coord));
    % get size
    n = size(a,1);
    % take single edges
    b = a(find(triu(ones(n),1)));
    % make into a ks density plot
    [f xi] = ksdensity(b);
    % plot and shift to right
    h = plot(xi + 500,f,'linewidth',2);
    m(culture) = max(xi) + 500;
    % keep plotting
    hold on;
end
% keep max euclidean for later
ma = max(m,[],'all');
% set limits
xlim([0 5000]);
yticks([]); b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Euclidean \mum'); ylabel('Frequency'); box off; b.TickLength = [0 0];
% Plot the eucldiean distances of each culture
h = figure;
h.Position = [100 100 1300 160];
n = [];
for culture = 1:ncultures;
    % subplot
    subplot(1,6,culture);
    % take the coordinates
    coords = data{pipeline}.table_xy{culture};
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
%% compute summary statistics over time
% number of measures
nmeasures = 10;
% initialise
statistics = nan(npipelines,ncultures,nmeasures);
% statistics labels
statistics_labels = string({...
    'Mean STTC',...
    'Number of nodes',...
    'Network density (%)',...
    'Total degree',...
    'Efficiency',...
    'Betweenness',...
    'Clustering',...
    'Edge length',...
    'Modularity',...
    'Matching',...
    'Small-worldness'});
for pipeline = 1:npipelines;
    % display pipeline
    disp(sprintf('computing %s...',pipeline_labels(pipeline)));
    for culture = 1:ncultures;
        % display culture and div
        disp(sprintf('computing culture %g...',culture));
         % take weighted network
         wnetwork = data{pipeline}.table_wu{culture};
         % take the binary network
         network = data{pipeline}.table_bu{culture};
         % take the coordinates
         coords = data{pipeline}.table_xy{culture};
         % compute the euclidean distance
         euclidean = squareform(pdist(coords));
         % compute global statistics
         statistics(pipeline,culture,1) = mean(triu(wnetwork,1),'all');
         % number of nodes
         statistics(pipeline,culture,2) = size(network,1);
         % density
         statistics(pipeline,culture,3) = density_und(network);
         % total number of connections
         statistics(pipeline,culture,4) = sum(degrees_und(network));
         % global efficiency
         statistics(pipeline,culture,5) = efficiency_bin(network,0);
         % betweenness
         statistics(pipeline,culture,6) = mean(betweenness_bin(network));
         % clustering
         statistics(pipeline,culture,7) = mean(clustering_coef_bu(network));
         % average weighted length
         statistics(pipeline,culture,8) = mean(mean(euclidean.*network));
         % modularity
         [~,statistics(pipeline,culture,9)] = modularity_und(network);
         % homophily
         statistics(pipeline,culture,10) = mean(matching_ind(network),'all');
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
         statistics(pipeline,culture,11) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    end
end
%% plot number of connections and densities by pipeline and div
% pipeline
pipeline = 1;
% stat
stat = 3;
% visualise the number of connections
k = squeeze(statistics(pipeline,:,stat));
% plot
h = figure; h.Position = [100 100 400 400];
iosr.statistics.boxPlot(100*k',...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'boxColor',[.5 .5 .5],...
    'boxAlpha',.5,...
    'symbolColor','k');
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
%% plot connections over distance binned
% set pipeline 
pipeline = 1;
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binsum = []; binprop = [];
% compute
for culture = 1:ncultures;
    % network
    a = data{pipeline}.table_bu{culture}; 
    % euclidean
    b = squareform(pdist(data{pipeline}.table_xy{culture}));
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
xlabel('Distance (\mum)'); ylabel('Prop. connections'); xticks([0:1000:4000]);
%% firing rate matrix
% set the well to plot
well = 2;
% load the recording
load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_cerebral_organoids_30min/M03912/211231/well%g/spk_data',...
    well));
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
h = figure; h.Position = [100 100 600 280];
imagesc(array);
xlabel('Time (s)');
ylabel('Putative neurons');
c = colorbar; c.Label.String = 'Firing Rate (Hz)'; 
set(gca,'ColorScale','log');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
% plot the summed spike vector
summed_spike = sum(array);
h = figure; h.Position = [100 100 600 280];
plot(summed_spike,'color',[10 153 153]./256,'linewidth',2)
xlim([0 300]); ylim([0 2400]);
xlabel('Time (s)');
ylabel('Co-active units');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out'; box off;
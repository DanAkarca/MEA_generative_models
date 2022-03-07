%% descriptives of human cortical hd gnm data
% written by danyal akarca
%% load data
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_ipsc_DIV28_30min');
% load sttc data
gluta = load('GN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
motor = load('MN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
dopa = load('DN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat');
% collect the data together
data = {gluta.gm_data motor.gm_data dopa.gm_data};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
%% form labels
% 3 cell types
type_labels = string({'GN','MN','DN'});
% 3 pipelines
pipeline_labels = string({'10ms'});
% 4 time points
div_labels = string({'28'});
%% set hyperparameters of the dataset
% number of types
ntypes = 3;
% number of cultures
ncultures = [8 7 6];
% find missing networks
missing = [];
step = 1;
% form the array
for type = 1:ntypes;
    for culture = 1:ncultures(type);
        % take if empty (i.e. not existant) or an empty network itself (dn pipeline 2 culture 6 div 14)
        if isempty(data{type}.table_bu{culture}) | sum(data{type}.table_bu{culture},'all')==0;
            missing(step,:) = [type pipeline culture div];
            step = step + 1;
        end
    end
end
%% visualise the networks
% cell type
type = 3;
% select a culture - GN=6, MN=3, DN=5
culture = 5;
% compute the cultures
network = data{type}.table_bu{culture};
coords = data{type}.table_xy{culture};
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
    'NodeColor',colbar(type,:),...
    'NodeLabel',[],...
    'MarkerSize',.6*degrees_und(network)+0.01);
xticks([]); yticks([]);
set(gca,'visible','off');
%% visualise across all cultures 
% set hyperparameters 
type = 3;
% compute n
n = ncultures(type);
% initialise
h = figure;
h.Position = [100 100 600 300];
m = [];
% loop over cultures
for culture = 1:n;
    % take the network
    network = data{type}.table_bu{culture};
    % take the coordinates of this network
    coord = data{type}.table_xy{culture};
    % make eucliean
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
yticks([]); b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Euclidean \mum'); ylabel('Frequency'); box off; b.TickLength = [0 0];
% Plot the eucldiean distances of each culture
h = figure;
h.Position = [100 100 1300 400];
n = [];
for culture = 1:ncultures(type);
    % subplot
    subplot(2,6,culture);
    % take the coordinates
    coords = data{type}.table_xy{culture};
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
nmeasures = 11;
% initialise
statistics = nan(ntypes,max(ncultures),nmeasures);
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
% loop over networks
for type = 1:ntypes;
    % display cell type
    disp(sprintf('computing %s networks...',type_labels(type)));
    for culture = 1:ncultures(type);
        % display culture and div
        disp(sprintf('computing culture %g ...',culture));
        % take weighted network
        wnetwork = data{type}.table_wu{culture};
        % take the binary network
        network = data{type}.table_bu{culture};
        % take the coordinates
        coords = data{type}.table_xy{culture};
        % compute the euclidean distance
        euclidean = squareform(pdist(coords));
        % compute global statistics
        statistics(type,culture,1) = mean(triu(wnetwork,1),'all');
        % number of nodes
        statistics(type,culture,2) = size(network,1);
        % density
        statistics(type,culture,3) = density_und(network);
        % total number of connection
        statistics(type,culture,4) = sum(degrees_und(network));
        % global efficiency
        statistics(type,culture,5) = efficiency_bin(network,0);
        % betweenness
        statistics(type,culture,6) = mean(betweenness_bin(network));
        % clustering
        statistics(type,culture,7) = mean(clustering_coef_bu(network));
        % average weighted length
        statistics(type,culture,8) = mean(mean(euclidean.*network));
        % modularity
        [~,statistics(type,culture,9)] = modularity_und(network);
        % homophily
        statistics(type,culture,10) = mean(matching_ind(network),'all');
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
        statistics(type,culture,11) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    end
end
%% plot a speific time and pipeline over types
% set stat
stat = 11;
% visualise 
k = squeeze(statistics(:,:,stat))';
% compute statistics
[p anovatab stats] = anova1(k);
[comparison means h gnames] = multcompare(stats);
% plot
h = figure; h.Position = [100 100 400 400];
iosr.statistics.boxPlot(k,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'boxColor',[.5 .5 .5],...
    'boxAlpha',.5,...
    'symbolColor','k');
xticklabels({'GN','MN','DN'});
xlabel('Human iPSC line');
ylabel(statistics_labels(stat));
b = gca; b.TickDir = 'out';
b.FontName = 'Arial'; 
b.FontSize = 25;
% tabulate
med_tab = median(k);
iqr_tab = iqr(k);
stat_tab = table(med_tab,iqr_tab);
% columns are divs, rows arepipelines
disp(stat_tab);
%% run statistics
% compte statistics
culture_id = 1:size(k,1);
ndivs = 4;
rtbl = table(culture_id',k(:,1),k(:,2),k(:,3),k(:,4),...
    'VariableNames',{'Culture','DIV7','DIV10','DIV12','DIV14'});
stats = table([1:ndivs]','VariableNames',{'DIVs'});
% fit rm
rm = fitrm(rtbl,'DIV7-DIV14~Culture','WithinDesign',stats);
% run ranova
[ranovatbl a b c] = ranova(rm);
disp(ranovatbl);
%% plot connections over distance binned
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binprop = nan(3,8,length(binvec)-1);
% compute
for type = 1:ntypes;
    for culture = 1:ncultures(type);
    % network
    a = data{type}.table_bu{culture}; 
    % euclidean
    b = squareform(pdist(data{type}.table_xy{culture}));
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
        binprop(type,culture,k-1) = prop;
    end
    end
end
% compute mean and std
x = squeeze(mean(binprop,2,'omitnan'))*100;
y = squeeze(std(binprop,[],2,'omitnan'))*100;
% form a jitter
j = 50;
% set colorbars
colbar = 1-summer(4);
% visualise
h = figure; h.Position = [100 100 400 400];
for type = 1:ntypes;
uu = errorbar(binvec(2:end)+j,x(type,:),y(type,:)); 
uu.Color = colbar(type,:); uu.Color(4) = .2;
uu.Marker = '.'; 
uu.MarkerSize = 50;
uu.LineWidth = 5;
hold on;
end
xlim([0 binlim+2*j]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Distance (\mum)'); ylabel('Prop. connections'); xticks([0:1000:4000]);
legend({'GN','MN','DN'});
box off;
%% firing rate matrix
% set the type and culture
type = 1;
culture = 6; % 6, 3, 5
% set axis limits for number of spikes
yaxisl = [0 2400];
xaxisl = [0 300];
% load the recording
load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/human_ipsc_DIV28_30min/%g/210316/spk_data',...
    data{type}.chip_id(culture)));
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
    % define the bins
    binlims = [0:bin:(ceil(max(x/10))*10)];
    [count,idx] = histc(x,binlims);
    % keep count
    array(electrode,1:length(count)) = count;
end
% set colours
colbar = 1-summer(4);
% visualise
h = figure; h.Position = [100 100 600 280];
imagesc(array);
xlabel('Time (s)');
ylabel('Putative neuron');
c = colorbar; c.Label.String = 'Firing Rate (Hz)';
set(gca,'ColorScale','log');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
% plot the summed spike vector
summed_spike = sum(array);
h = figure; h.Position = [100 100 600 280];
plot(summed_spike,'color',colbar(type,:),'linewidth',2.5);
ylim(yaxisl);
xlim(xaxisl);
xlabel('Time (s)');
%ylabel('Co-active units');
ylabel({'# spikes';'per 1 s bin'});
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
box off;
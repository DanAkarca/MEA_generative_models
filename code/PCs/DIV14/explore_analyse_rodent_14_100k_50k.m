%% Analyse DIV14 50k versus 100k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortcal neurons
% written by danyal akarca
clear; clc;
%{
% change directory to the 100k rodent project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ForAlex_updated/Oct_2021/Euler_100k_rodent/tracking/');
% load 100k sttc data
a = load('gm_data_shared_div14_21_28_min_r_0.1_alpha_0.001_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div14_21_28_min_r_0.1_alpha_0.001_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div14_21_28_min_r_0.1_alpha_0.001_lag20_jitter20_prob.mat');
% collect the 100k together
k100 = {a.gm_data b.gm_data c.gm_data};
% change directory to the 50k rodent project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ForAlex_updated/Oct_2021/Euler_50k_rodent/tracking/');
% load 50k sttc data
a = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag5_jitter5.mat');
b = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag10_jitter10.mat');
c = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag20_jitter20.mat');
% collect the 50k data together
k50 = {a.gm_data b.gm_data c.gm_data};
%}
% change directory to the 100k rodent project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_100k_tracking/');
% load sttc data
a = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div14_21_28_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
k100 = {a.gm_data b.gm_data c.gm_data};
% change directory to the 50k rodent project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/');
% load sttc data
a = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
k50 = {a.gm_data b.gm_data c.gm_data};
% addpaths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master');
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
addpath('/imaging/astle/users/da04/PhD/toolboxes/bluewhitered');
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/toolboxes/small_world/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/functions/');
% define model types
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
eqn = string({'KS_k','KS_c','KS_b','KS_d'});
nmodels = 13;
% define pipeline labels
pipeline_labels = string({'5ms','10ms','20ms'});
%% take 100 14 day data
% take 50k data
npipelines = 3;
ncultures = 12;
ndivs = 3;
% initialise the indices: columns for pipeline, then culture, then div
k100_index = []; 
all_100k_networks = {};
all_100k_wnetworks = {};
all_100k_coordinates = {};
all_100k_costs = {};
step = 1;
% loop through and visualise
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            % form the index
            k100_index(step,1) = pipeline;
            k100_index(step,2) = culture;
            k100_index(step,3) = div;
            % get the data as a cell array
            all_100k_networks{step} = k100{pipeline}.table_bu_shared{culture}{div};
            all_100k_wnetworks{step} = k100{pipeline}.table_wu_shared{culture}{div};
            all_100k_coordinates{step} = k100{pipeline}.table_xy_shared{culture}{div};
            all_100k_costs{step} = squareform(pdist(k100{pipeline}.table_xy_shared{culture}{div}));
            % update the step
            step = step + 1;
        end
    end
end 
% keep 14 day 100k
i = find(k100_index(:,3)==3);
k100_14_index = k100_index(i,:);
k100networks = all_100k_networks(i);
k100wnetworks = all_100k_wnetworks(i);
k100coordinates = all_100k_coordinates(i);
k100costs = all_100k_costs(i);
k100n = length(k100networks);
%% take 50k 14 day data
% take 50k data
npipelines = 3;
ncultures = 6;
ndivs = 4;
% initialise the indices: columns for pipeline, then culture, then div
k50_index = []; 
all_50k_networks = {};
all_50k_wnetworks = {};
all_50k_coordinates = {};
all_50k_costs = {};
step = 1;
% loop through and visualise
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            % form the index
            k50_index(step,1) = pipeline;
            k50_index(step,2) = culture;
            k50_index(step,3) = div;
            % get the data as a cell array
            all_50k_networks{step} = k50{pipeline}.table_bu_shared{culture}{div};
            all_50k_wnetworks{step} = k50{pipeline}.table_wu_shared{culture}{div};
            all_50k_coordinates{step} = k50{pipeline}.table_xy_shared{culture}{div};
            all_50k_costs{step} = squareform(pdist(k50{pipeline}.table_xy_shared{culture}{div}));
            % update the step
            step = step + 1;
        end
    end
end 
% keep 14 day 50k
i = find(k50_index(:,3)==4);
k50_14_index = k50_index(i,:);
k50networks = all_50k_networks(i);
k50wnetworks = all_50k_wnetworks(i);
k50coordinates = all_50k_coordinates(i);
k50costs = all_50k_costs(i);
k50n = length(k50networks);
%% look at statistics side by side
% compute the number of networks
nnets100k = length(k100networks);
nnets50k = length(k50networks);
% number of measures
nmeasures = 10;
% initialise
global_statistics100k = zeros(nnets100k,nmeasures);
global_statistics50k = zeros(nnets50k,nmeasures);
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
    'Edge length (\mum)',...
    'Modularity',...
    'Matching'});
% loop over 100k networks
for net = 1:nnets100k;
    % take the binary network
    network = k100networks{net};
    % take the weighted network
    wnetwork = k100wnetworks{net};
    % take the euclidean distances
    euclidean = k100costs{net};
    % compute global statistics
    % density
    global_statistics100k(net,1) = density_und(network);
    % total strength
    global_statistics100k(net,2) = sum(strengths_und(wnetwork(~isnan(wnetwork))));
    % total number of connections
    global_statistics100k(net,3) = sum(degrees_und(network));
    % global efficiency
    global_statistics100k(net,4) = efficiency_bin(network,0);
    % betweenness
    global_statistics100k(net,5) = mean(betweenness_bin(network));
    % clustering
    global_statistics100k(net,6) = mean(clustering_coef_bu(network));
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
        global_statistics100k(net,7) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    % average weighted length
    global_statistics100k(net,8) = mean(mean(euclidean.*network));
    % modularity
    [~,global_statistics100k(net,9)] = modularity_und(network);
    % homophily
    global_statistics100k(net,10) = mean(matching_ind(network),'all');
    %{
    % weighted global efficiency
    global_statistics100k(net,10) = efficiency_wei(wnetwork(~isnan(wnetwork)),0);
    % weighted betwenness
    global_statistics100k(net,11) = mean(betweenness_wei(wnetwork(~isnan(wnetwork))));
    % weighted clutering
    global_statistics100k(net,12) = mean(clustering_coef_wu(wnetwork(~isnan(wnetwork))));
    %}
    % display
    disp(sprintf('50k network %g statistics computed',net));
end
% loop over 50k networks
for net = 1:nnets50k
    % take the binary network
    network = k50networks{net};
    % take the weighted network
    wnetwork = k50wnetworks{net};
    % take the euclidean distances
    euclidean = k50costs{net};
    % compute global statistics
    % density
    global_statistics50k(net,1) = density_und(network);
    % total strength
    global_statistics50k(net,2) = sum(strengths_und(wnetwork(~isnan(wnetwork))));
    % total number of connections
    global_statistics50k(net,3) = sum(degrees_und(network));
    % global efficiency
    global_statistics50k(net,4) = efficiency_bin(network,0);
    % betweenness
    global_statistics50k(net,5) = mean(betweenness_bin(network));
    % clustering
    global_statistics50k(net,6) = mean(clustering_coef_bu(network));
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
        global_statistics50k(net,7) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
    % average weighted length
    global_statistics50k(net,8) = mean(mean(euclidean.*network));
    % modularity
    [~,global_statistics50k(net,9)] = modularity_und(network);
    % matching
    global_statistics50k(net,10) = mean(matching_ind(network),'all');
    %{
    % weighted global efficiency
    global_statistics50k(net,10) = efficiency_wei(wnetwork(~isnan(wnetwork)),0);
    % weighted betwenness
    global_statistics50k(net,11) = mean(betweenness_wei(wnetwork(~isnan(wnetwork))));
    % weighted clutering
    global_statistics50k(net,12) = mean(clustering_coef_wu(wnetwork(~isnan(wnetwork))));
    %}
    % display
    disp(sprintf('50k network %g statistics computed',net));
end
%% compare the statistics
% set colours
col = 1-jet(10);
% set pipeline
pipeline = 2;
% alter the order
order = [1 4 8 10 5 9 2 3 6 7];
% visualise
h = figure; h.Position = [100 100 1500 600];
% loop over statistics
for stat = 1:nmeasures;
    b = subplot(2,ceil(nmeasures)/2,stat);
    % form array
    a = global_statistics100k(k100_14_index(:,1)==pipeline,order(stat));
    b = global_statistics50k(k50_14_index(:,1)==pipeline,order(stat));
    c = nan(size(a,1),1);
    c(1:size(b,1)) = b;
    array = [c a];
    % statistics
    [p h] = ranksum(array(:,1),array(:,2));
    % view histogram
    u = iosr.statistics.boxPlot(array,...
        'theme','colorall',...
        'showViolin',logical(0),...
        'boxColor','k',...
        'showScatter',logical(1),...
        'scatterColor',[.5 .5 .5],...
        'scatterMarker','x',...
        'symbolMarker','x',...
        'symbolColor','k');
    xticklabels({'50k','100k'}); 
    ylabel(sprintf('%s',global_statistics_labels(order(stat))));
    xlabel('Plating density');
    if p < 0.05
    title(sprintf('p=%.3d',p));
    else
        title(sprintf('p=%.3g',p));
    end
    % update colours
    u.handles.box(1).FaceColor = [142 199 243]./256;
    u.handles.box(2).FaceColor = [180 148 148]./256;
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontName = 'Arial'; b.FontSize = 20;
end
[142 199 243; 180 148 148]./256;

%% plot specific comparisons
% set the specifics to plot
plotv = [1 4 8 10 9 7];
% visualise
h = figure; h.Position = [100 100 1600 800];
% loop over statistics
for stat = 1:length(plotv);
    b = subplot(2,3,stat);
    % form array
    a = global_statistics100k(k100_14_index(:,1)==pipeline,plotv(stat));
    b = global_statistics50k(k50_14_index(:,1)==pipeline,plotv(stat));
    c = nan(size(a,1),1);
    c(1:size(b,1)) = b;
    array = [c a];
    % for density make it a percentage
    if stat == 1;
        array = array*100;
    end
    % statistics
    [p(stat) h] = ranksum(array(:,1),array(:,2));
    % view histogram
    u = iosr.statistics.boxPlot(array,...
        'theme','colorall',...
        'showViolin',logical(0),...
        'boxColor','w',...
        'showScatter',logical(1),...
        'scatterColor',[.5 .5 .5],...
        'scatterMarker','x',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxalpha',1);
    xticklabels({'50k','100k'}); 
    ylabel(sprintf('%s',global_statistics_labels(plotv(stat))));
    xlabel('Plating density');
    % update colours
    u.handles.box(1).FaceColor = [142 199 243]./256;
    u.handles.box(2).FaceColor = [180 148 148]./256;    
    % settings
    b = gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial';
    yz = b.YLim;
    ylim([yz(1)-0.05*yz(1) yz(2)+0.05*yz(2)]);
    % form asterisk
    u = text(1.5,yz(2)-(0.01*yz(2)),'*'); u.FontSize = 25; u.FontName = 'Arial';
    box off;
end
%% show the small-worldness plot
% set stat
stat = 10;
% loop over pipelines
array = []; p = [];
for pipeline = 1:npipelines;
% form array
a = global_statistics100k(k100_14_index(:,1)==pipeline,order(stat));
b = global_statistics50k(k50_14_index(:,1)==pipeline,order(stat));
c = nan(size(a,1),1);
c(1:size(b,1)) = b;
array(:,:,pipeline) = [c a];
[p(pipeline)] = ranksum(c,a);
end
% statistics
h = figure; h.Position = [100 100 300 300];
u = iosr.statistics.boxPlot(array,...
        'theme','colorall',...
        'showViolin',logical(0),...
        'boxcolor',{'r','b','g'},...
        'showScatter',logical(1),...
        'scatterColor',[.5 .5 .5],...
        'scatterMarker','x',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxalpha',0.1);
xticklabels({'50k','100k'}); 
ylabel(sprintf('%s',global_statistics_labels(order(stat))));
% add in the small worldness threshold
ylim([0 inf]);
% add the line
hold on;
o = line('XData',[0 1 2 3],'YData',[1 1 1 1],'linestyle','--','linewidth',3,'color',[.5 .5 .5]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;

% a specific over pipelines
pipeline = 2;
% form array
a = global_statistics100k(k100_14_index(:,1)==pipeline,order(stat));
b = global_statistics50k(k50_14_index(:,1)==pipeline,order(stat));
c = nan(size(a,1),1);
c(1:size(b,1)) = b;
array = [c a];
p = ranksum(c,a);
% statistics
h = figure; h.Position = [100 100 300 300];
u = iosr.statistics.boxPlot(array,...
        'theme','colorall',...
        'showViolin',logical(1),...
        'boxcolor',col(stat,:),...
        'showScatter',logical(1),...
        'scatterColor',[.5 .5 .5],...
        'scatterMarker','x',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxalpha',0.2);
xticklabels({'50k','100k'}); 
ylabel(sprintf('%s',global_statistics_labels(order(stat))));
if p < 0.05;
    title(sprintf('p=%.3d',p));
else
    title(sprintf('p=%.3g',p));
end
% add in the small worldness threshold
ylim([0 7]);
% add the line
hold on;
o = line('XData',[0 1 2 3],'YData',[1 1 1 1],'linestyle','--','linewidth',3,'color',[.5 .5 .5]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
%% compare global correlation matrices within pipelines
% set pipeline
pipeline = 2;
% take statitic from pipeline
a = global_statistics50k(k50_14_index(:,1)==pipeline,:);
b = global_statistics100k(k100_14_index(:,1)==pipeline,:);
% compute correalation matrices
acorr = corr(a);
bcorr = corr(b);
% visualise
h = figure; h.Position = [100 100 750 350];
subplot(1,2,1); imagesc(acorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; 
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
subplot(1,2,2); imagesc(bcorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; 
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
% set comparisons
x = 4; 
y = 7;
% visualise the single example
h = figure; h.Position = [100 100 500 350];
p = scatterhist([a(:,x); b(:,x)],[a(:,y); b(:,y)],...
    'group',[ones(6,1);2*ones(12,1)],...
    'kernel','on',...
    'direction','out',...
    'color',[.1 .6 .9; .5 .3 .3],...
    'linestyle','-',...
    'marker','.',...
    'markersize',20);
legend({'50k','100k'});
xlabel(global_statistics_labels(x)); ylabel(global_statistics_labels(y));
d = gca; d.TickDir = 'out'; d.FontName = 'Arial'; d.FontSize = 14;
%% compare global correlation matrices across all pipelines
% compute correalation matrices
acorr = corr(global_statistics50k);
bcorr = corr(global_statistics100k);
% visualise
h = figure; h.Position = [100 100 750 350];
subplot(1,2,1); imagesc(acorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; 
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
subplot(1,2,2); imagesc(bcorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; 
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
% reorder matrix
order = [1 3 8 2 4 6 10 5 7 9];
% visualise a single matrix for both
tricorr = [tril(acorr,-1) + triu(bcorr,1)]; 
h = figure; h.Position = [100 100 800 750];
u = imagesc(tricorr(order,order));
caxis([-1 1]); xticks(1:nmeasures);
global_statistics_labels{1} = 'Network Density'; % change back for this
xticklabels(global_statistics_labels(order));
xtickangle(45);
ylabel('Global network statistics'); yticks(1:nmeasures);
b=gca; b.TickDir = 'out'; b.FontSize = 25; b.FontName = 'Arial'; b.TickLength = [0 0];
nstep = 500;
lin = linspace(0,1,nstep)'; 
c = ones(nstep,3); c(:,1)=lin; c(:,2)=lin;
d = ones(nstep,3); d(:,2)=lin; d(:,3)=lin;
col = [d;flip(c)]; colormap(col); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r';
brighten(.6);
% set comparisons
x = 4; 
y = 7;
% visualise the single example
h = figure; h.Position = [100 100 900 700];
p = scatterhist([global_statistics50k(:,x); global_statistics100k(:,x)],[global_statistics50k(:,y); global_statistics100k(:,y)],...
    'group',[ones(18,1);2*ones(36,1)],...
    'kernel','on',...
    'direction','out',...
    'color',[142 199 243; 180 148 148]./256,...
    'linestyle','-',...
    'linewidth',7,...
    'marker','.',...
    'markersize',50);
legend({'50k plating density','100k plating density'},'box','off'); box off;
xlabel(global_statistics_labels(x)); ylabel(global_statistics_labels(y));
d = gca; d.TickDir = 'out'; d.FontName = 'Arial'; d.FontSize = 25;
box off;
%% view example networks side by side
% set which 50k network to load
k50net = 12;
% set which 100k network to load
k100net = 19;
% form the networks
A = abs(k50networks{k50net}); A(find(eye(size(A))))=0;
B = abs(k100networks{k100net}); B(find(eye(size(B))))=0;
% get edge widths
Aw = triu(A,1); Aw = Aw(:) + 1e-10;
% visualise
h = figure; h.Position = [100 100 800 300]; 
subplot(1,2,1); plot(graph(A),...
    'XData',k50coordinates{k50net}(:,1),...
    'YData',k50coordinates{k50net}(:,2),...
    'nodelabel',[],...
    'nodecolor',[142 199 243]./256,...
    'markersize',1.5*degrees_und(A)+0.001,...
    'edgecolor',[120 120 120]./256,...
    'edgealpha',.1);
set(gca,'Visible','off');
b = gca; b.TickDir = 'out';
subplot(1,2,2); plot(graph(B),...
    'XData',k100coordinates{k100net}(:,1),...
    'YData',k100coordinates{k100net}(:,2),...
    'nodelabel',[],...
    'nodecolor',[180 148 148]./256,...
    'markersize',1.5*degrees_und(B)+0.001,...
    'edgecolor',[120 120 120]./256,...
    'edgealpha',.1);
b = gca; b.TickDir = 'out';
set(gca,'Visible','off');
%{
%% load the 100k 14 days generative models 0.1Hz
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/rodent_100k_qc');
% 14 day only
div = 1;
% set hyperparameters
npipelines = 3;
ncultures = 12;
% initialise
k100_energy_sample = zeros(k100n,13,20000);
k100_ks_sample = zeros(k100n,13,20000,4);
k100_networks_sample = cell(k100n,1);
k100_parameters_sample = zeros(k100n,13,20000,2);
errors = zeros(k100n,1);
% loop over pipelines, cultures and divs
step = 1;
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        try % load this network's generative model output
            load(sprintf('rodent_qc_%g_%g_%g_generative_model.mat',pipeline,culture,div));
        catch
            % keep if it doesn't load
            errors(step) = 1;
            % display
            disp(sprintf('network_%g_%g_%g non-existant',pipeline,culture,div));
        end
        % assign
        k100_energy_sample(step,:,:) = output.energy;
        k100_ks_sample(step,:,:,:) = output.ks;
        k100_networks_sample{step} = output.networks;
        k100_parameters_sample(step,:,:,:) = output.parameters;
        % clear the variable
        clear output
        % display
        disp(sprintf('network_%g_%g_%g loaded',pipeline,culture,div));
        % upate step
        step = step + 1;
    end
end
%% load the 50k 14 days generative models 0.1Hz
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/rodent_50k_qc');
% 14 day only
div = 4;
% set hyperparameters
npipelines = 3;
ncultures = 6;
% initialise
k50_energy_sample = zeros(k50n,13,20000);
k50_ks_sample = zeros(k50n,13,20000,4);
k50_networks_sample = cell(k50n,1);
k50_parameters_sample = zeros(k50n,13,20000,2);
errors = zeros(k50n,1);
% loop over pipelines, cultures and divs
step = 1;
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        try % load this network's generative model output
            load(sprintf('rodent_50k_qc_%g_%g_%g_generative_model.mat',pipeline,culture,div));
            % assign
            k50_energy_sample(step,:,:) = output.energy;
            k50_ks_sample(step,:,:,:) = output.ks;
            k50_networks_sample{step} = output.networks;
            k50_parameters_sample(step,:,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('network_%g_%g_%g loaded',pipeline,culture,div));
            % upate step
            step = step + 1;
        catch
            % keep if it doesn't load
            errors(step) = 1;
            % display
            disp(sprintf('network_%g_%g_%g non-existant',pipeline,culture,div));
        end
    end
end
%}
%% load the 100k 14 days generative models 0.01Hz
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_100k_tracking/GNM_tracking_rec1/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/');
ncultures = 12;
% initialise
k100_energy_sample = zeros(ncultures,nmodels,20000);
k100_ks_sample = zeros(ncultures,nmodels,20000,4);
k100_networks_sample = cell(ncultures,1);
k100_parameters_sample = zeros(ncultures,nmodels,20000,2);
% loop over pipelines, cultures and divs
step = 1;
    for culture = 1:ncultures;
        for model = 1:nmodels
            % load
            load(sprintf('rodent_1_%g_1_generative_model_%g.mat',culture,model));
            % assign
            k100_energy_sample(step,model,:) = output.energy;
            k100_ks_sample(step,model,:,:) = output.ks;
            k100_networks_sample{step}(model,:,:) = output.networks;
            k100_parameters_sample(step,model,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('rodent_100k_1_%g_1_generative_model %s loaded',culture,models(model)));
        end
        step = step + 1;
    end
%% load the 50k 14 days generative models 0.01Hz
% change directory to generative model data
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/GNM_tracking_rec4/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob/');
ncultures = 6;
% initialise
k50_energy_sample = zeros(ncultures,nmodels,20000);
k50_ks_sample = zeros(ncultures,nmodels,20000,4);
k50_networks_sample = cell(ncultures,1);
k50_parameters_sample = zeros(ncultures,nmodels,20000,2);
% loop over pipelines, cultures and divs
step = 1;
    for culture = 1:ncultures;
        for model = 1:nmodels
            % load
            load(sprintf('rodent_1_%g_4_generative_model_%g.mat',culture,model));
            % assign
            k50_energy_sample(step,model,:) = output.energy;
            k50_ks_sample(step,model,:,:) = output.ks;
            k50_networks_sample{step}(model,:,:) = output.networks;
            k50_parameters_sample(step,model,:,:) = output.parameters;
            % clear the variable
            clear output
            % display
            disp(sprintf('rodent_50k_1_%g_1_generative_model %s loaded',culture,models(model)));
        end
        step = step + 1;
    end
%% compare the 14 day generative models 
% set pipeline
pipeline = 2;
% take the top n energy values
n = 1;
% sort rows of energy
x = sort(k100_energy_sample,3); 
x = squeeze(x(:,:,1:n));
y = sort(k50_energy_sample,3); 
y = squeeze(y(:,:,1:n));
% take data
array = nan([size(x) 2]);
array(1:6,:,1) = y;
array(:,:,2) = x;
% alter the order
array = array(:,[2:13 1],:);
% average
array_mean_rules = [];
array_mean_rules(:,1,:) = mean(array(:,1:2,:),2);
array_mean_rules(:,2,:) = mean(array(:,8:12,:),2);
array_mean_rules(:,3,:) = mean(array(:,3:7,:),2);
array_mean_rules(:,4,:) = array(:,13,:);
% all
array_rules = nan(12*5,4,2);
a = squeeze(array(:,[1:2],:)); a = reshape(a,[12*2 2]); array_rules(1:12*2,1,:) = a;
b = squeeze(array(:,[8:12],:)); b = reshape(b,[12*5 2]); array_rules(1:12*5,2,:) = b;
c = squeeze(array(:,[3:7],:)); c = reshape(c,[12*5 2]); array_rules(1:12*5,3,:) = c;
d = squeeze(array(:,[13],:)); d = reshape(d,[12 2]); array_rules(1:12,4,:) = d;
% permute
arraypermute = permute(array_rules,[1 3 2]);
% model order
% visualise
h = figure; h.Position = [100 100 800 650];
iosr.statistics.boxPlot(arraypermute,...
    'theme','colorall',...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.5,...
    'symbolColor','k',...
    'symbolMarker','x',...
    'boxAlpha',0.2);
xticklabels({'50k','100k'}); xlabel('Plating density');
ylim([0 0.7]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
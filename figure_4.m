%% figure 4
% written by danyal akarca
%% load data
clear; clc;
% change directory to the project folder
dir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/';
cd(dir);
% load sttc data
% 50k
% 5ms
load(strcat(dir,'50k_rodent_dev/DIV14_lag5/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag5_jitter5_prob.mat')); div14_50k_5ms = gm_data;
% 10ms
load(strcat(dir,'50k_rodent_dev/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div14_50k_10ms = gm_data;
% 20ms
load(strcat(dir,'50k_rodent_dev/DIV14_lag20/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag20_jitter20_prob.mat')); div14_50k_20ms = gm_data;
% 100k
% 5ms
load(strcat(dir,'100k_rodent_dev/DIV14_lag5/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag5_jitter5_prob.mat')); div14_100k_5ms = gm_data;
% 10ms
load(strcat(dir,'100k_rodent_dev/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div14_100k_10ms = gm_data;
% 20ms
load(strcat(dir,'100k_rodent_dev/DIV14_lag20/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag20_jitter20_prob.mat')); div14_100k_20ms = gm_data;
rodent_data_all = {{div14_50k_5ms div14_50k_10ms div14_50k_20ms},{div14_100k_5ms div14_100k_10ms div14_100k_20ms}};
% load generative model data
load(strcat(dir,'50k_rodent_dev/tracking_rec4/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div14_50k = results_generative_models;
load(strcat(dir,'100k_rodent_dev/tracking_rec1/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div14_100k = results_generative_models;
% collect the data together
data = {div14_50k div14_100k};
clear div14_50k div14_100k;
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
% addpath of roundsd
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
%% compute statistics side by side
% compute the number of networks
nnet = [6 12];
% number of groups
ngroup = 2;
% set number of pipelines
npipelines = 3;
% number of measures
nmeasures = 10;
% remove nodes with no connections from empirical data
rodent_data = rodent_data_all;
% loop over networks and remove nodes with no connections
for group = 1:ngroup;
    for pipeline = 1:npipelines;
    for culture = 1:nnet(group);
        % get network
        Atgt = rodent_data_all{group}{pipeline}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = rodent_data_all{group}{pipeline}.table_wu{culture};
        xy = rodent_data_all{group}{pipeline}.table_xy{culture};
        % keep only nodes with edges
        rodent_data{group}{pipeline}.table_bu{culture} = Atgt(ids,ids);
        rodent_data{group}{pipeline}.table_wu{culture} = wu(ids,ids);
        rodent_data{group}{pipeline}.table_xy{culture} = xy(ids,:);
    end
    end
end
%% compute global statistics
% set which data to plot
dataplot = rodent_data_all;
% initialise
global_statistics = {};
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
% loop over groups
for group = 1:ngroup;
    % initialise
    global_statistics{group} = zeros(nnet(group),nmeasures,npipelines);
    % loop over networks and pipelines
    for pipeline = 1:npipelines;
    for culture = 1:nnet(group);;
        % take the binary network
        network = dataplot{group}{pipeline}.table_bu{culture};
        % take the weighted network
        wnetwork = dataplot{group}{pipeline}.table_wu{culture};
        % take the euclidean distances
        euclidean = squareform(pdist(dataplot{group}{pipeline}.table_xy{culture}));
        % compute global statistics
        % density
        global_statistics{group}(culture,1,pipeline) = density_und(network)*100;
        % total strength
        global_statistics{group}(culture,2,pipeline) = sum(strengths_und(wnetwork(~isnan(wnetwork))));
        % total number of connections
        global_statistics{group}(culture,3,pipeline) = sum(degrees_und(network));
        % global efficiency
        global_statistics{group}(culture,4,pipeline) = efficiency_bin(network,0);
        % betweenness
        global_statistics{group}(culture,5,pipeline) = mean(betweenness_bin(network));
        % clustering
        global_statistics{group}(culture,6,pipeline) = mean(clustering_coef_bu(network));
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
            global_statistics{group}(culture,7,pipeline) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
        % average weighted length
        global_statistics{group}(culture,8,pipeline) = mean(mean(euclidean.*network));
        % modularity
        [~,global_statistics{group}(culture,9,pipeline)] = modularity_und(network);
        % homophily
        global_statistics{group}(culture,10,pipeline) = mean(matching_ind(network),'all');
        % display
        disp(sprintf('pipeline %g, culture %g statistics computed',pipeline,culture));
    end
    end
end
%% compare the statistics
% set pipeline
pipeline = 2;
% set colours
col = 1-jet(10);
% alter the order
order = [1 4 8 10 5 9 2 3 6 7];
% visualise
h = figure; h.Position = [100 100 1500 600];
% loop over statistics
for stat = 1:nmeasures;
    subplot(2,ceil(nmeasures)/2,stat);
    % form array
    a = global_statistics{1}(:,stat,pipeline);
    b = global_statistics{2}(:,stat,pipeline);
    c = nan(size(b,1),1);
    c(1:size(a,1)) = a;
    array = [c b];
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
%% compare global correlation matrices across pipeines
% take statitics and unroll over pipelines
adata = global_statistics{1}; adata = [squeeze(adata(:,:,1)); squeeze(adata(:,:,2)); squeeze(adata(:,:,3))];
bdata = global_statistics{2}; bdata = [squeeze(bdata(:,:,1)); squeeze(bdata(:,:,2)); squeeze(bdata(:,:,3))];
% compute correalation matrices
acorr = corr(adata);
bcorr = corr(bdata);
% visualise seperate correlation matrices
h = figure; h.Position = [100 100 750 350];
subplot(1,2,1); imagesc(acorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14;
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
box off;
subplot(1,2,2); imagesc(bcorr); u = gca; u.TickDir = 'out'; u.FontName = 'Arial'; u.FontSize = 14;
caxis([-1 1]); xticks(1:nmeasures); xticklabels(global_statistics_labels); xtickangle(45);
box off; 
% visualise a single matrix for both
tricorr = [tril(acorr,-1) + triu(bcorr,1)]; 
order = [1 3 8 2 4 6 10 5 7 9];
h = figure; h.Position = [100 100 800 750];
u = imagesc(tricorr(order,order));
caxis([-1 1]); xticks(1:nmeasures);
global_statistics_labels{1} = 'Network density'; % change back for this
global_statistics_labels{8} = 'Edge length';
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
h = figure; h.Position = [100 100 800 750];
p = scatterhist([adata(:,x); bdata(:,x)],[adata(:,y); bdata(:,y)],...
    'group',[ones(6*3,1);2*ones(12*3,1)],...
    'kernel','on',...
    'direction','out',...
    'color',[142 199 243; 180 148 148]./256,...
    'linestyle','-',...
    'marker','.',...
    'markersize',50);
legend({'50k','100k'},'box','off');
xlabel(global_statistics_labels(x)); ylabel(global_statistics_labels(y));
d = gca; d.TickDir = 'out'; d.FontName = 'Arial'; d.FontSize = 25; box off;
%% view example networks side by side
% set pipeline 
pipeline = 2;
% set data to plot
dataplot = rodent_data_all;
% set which 50k network to load
k50net = 6;
% set which 100k network to load
k100net = 3;
% form the networks
A = abs(dataplot{1}{pipeline}.table_bu{k50net}); A(find(eye(size(A))))=0;
B = abs(dataplot{2}{pipeline}.table_bu{k100net}); B(find(eye(size(B))))=0;
% get edge widths
Aw = triu(A,1); Aw = Aw(:) + 1e-10;
% visualise
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); plot(graph(A),...
    'XData',dataplot{1}{pipeline}.table_xy{k50net}(:,1),...
    'YData',dataplot{1}{pipeline}.table_xy{k50net}(:,2),...
    'nodelabel',[],...
    'nodecolor',[142 199 243]./256,...
    'markersize',1.5*degrees_und(A)+0.001,...
    'edgecolor',[120 120 120]./256,...
    'edgealpha',.1);
set(gca,'Visible','off');
b = gca; b.TickDir = 'out';
subplot(1,2,2); plot(graph(B),...
    'XData',dataplot{2}{pipeline}.table_xy{k100net}(:,1),...
    'YData',dataplot{2}{pipeline}.table_xy{k100net}(:,2),...
    'nodelabel',[],...
    'nodecolor',[180 148 148]./256,...
    'markersize',1.5*degrees_und(B)+0.001,...
    'edgecolor',[120 120 120]./256,...
    'edgealpha',.1);
b = gca; b.TickDir = 'out';
set(gca,'Visible','off');
%% organise the generative model findings
% place the div14 data together
div14 = data;
% set number of models
nmodels = 13;
% set the sample size for the 50k and 100k cultures
ngroup = [6 12];
% define the top n parameters
nset = [1 10 50 100];
% compute the minimum
for group = 1:length(ngroup);
    % get the group
    energy = div14{group}.energy;
    parameters = div14{group}.parameters;
    % initialise
    top_energy = cell(length(nset),1);
    top_energy_mean = zeros(length(nset),13,ngroup(group));
    top_parameters = cell(length(nset),1);
    top_parameters_mean = zeros(length(nset),13,ngroup(group),2);
    % ouput
    div14_generative{group} = struct;
    % loop over sets
    for no = 1:length(nset)
        % take the actual amount of top performing parameters
        nnet_100k = nset(no);
        % loop over models
        for model = 1:nmodels;
            % take energies for this model
            pipeline_d = squeeze(energy(:,model,:))';
            % rank them for each subject
            [v i] = sort(pipeline_d);
            % take top n energies and their indices
            n_e = v(1:nnet_100k,:);
            n_i = i(1:nnet_100k,:);
            % take the corresponding parameters
            u = zeros(ngroup(group),nnet_100k,2);
            for s = 1:ngroup(group);
                % keep parameters
                u(s,:,:) = squeeze(parameters(s,model,n_i(:,s),:));
            end
            % if top parameter only
            if nnet_100k == 1
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
    div14_generative{group}.top_energy = top_energy;
    div14_generative{group}.top_energy_mean = top_energy_mean;
    div14_generative{group}.top_parameters = top_parameters;
    div14_generative{group}.top_parameters_mean = top_parameters_mean;
end
%% visualise
% take the top n energy values
nnet_100k = 1;
% combine
data = nan(nmodels,ngroup(2),2);
data(:,1:ngroup(1),1) = squeeze(div14_generative{1}.top_energy_mean(nnet_100k,:,:));
data(:,:,2) = squeeze(div14_generative{2}.top_energy_mean(nnet_100k,:,:));
% set order
i = [2:3 9:13 4:8 1];
data = data(i,:,:);
% permute the order
data = permute(data,[2 3 1]);
% iosr boxplot - all models
% visualise
h = figure; h.Position = [100 100 800 450];
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
xticklabels({'50k','100k'}); xlabel('Plating density');
ylim([0 0.6]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% iosr boxplot - grouped models
energy_div14_rules = nan(ngroup(2)*5,length(ngroup),4);
a = squeeze(data(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); energy_div14_rules(1:ngroup(2)*2,:,1) = a';
b = squeeze(data(:,:,[3:7])); b = permute(b,[2 1 3]); b = b(:,:); energy_div14_rules(1:ngroup(2)*5,:,2) = b';
c = squeeze(data(:,:,[8:12])); c = permute(c,[2 1 3]); c = c(:,:); energy_div14_rules(1:ngroup(2)*5,:,3) = c';
d = squeeze(data(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); energy_div14_rules(1:ngroup(2),:,4) = d';
% visualise
h = figure; h.Position = [100 100 800 450];
u = iosr.statistics.boxPlot(energy_div14_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5]}';
xticklabels({'50k','100k'}); xlabel('Plating density');
ylim([0 0.6]);
ylabel('Energy');
b = gca; 
b.TickDir = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% statistical comparisons between homophily
a = squeeze(energy_div14_rules(:,1,1));
b = squeeze(energy_div14_rules(:,2,1));
[h p] = ttest(a,b);
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
h.Position = [100 100 800 350];
for group = 1:length(ngroup);
    % form the subplot
    subplot(1,2,group);
    % take the measure
    e = squeeze(div14{group}.energy(:,model,:));
    p = squeeze(div14{group}.parameters(:,model,:,:));
    % visualise
    if model == 1
        h = figure;
        h.Position = [100 100 800 350];
        eta = squeeze(p(:,:,1));
        m = parula(10000);
        for culture = 1:ngroup(group);
            % get colours
            pipeline_d = round(e(culture,:) .* 10000);
            col = m(pipeline_d,:);
            % plot
            u = scatter(eta(culture,:),e(culture,:),20,col,'.'); ylabel('energy'); xlabel('/eta'); ylim([0 0.6]);
            b = gca; b.TickDir = 'out';
            b.FontName = 'Arial';
            yticks([]); xticks([]); hold on;
        end
    else
    for culture = 1:ngroup(group)
        u = scatter(squeeze(p(culture,:,1)),squeeze(p(culture,:,2)),20,e(culture,:),'.'); hold on;
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
%% compare parameters
% take top n parameters
n = 1;
% select the model
model = 3;
% seperate
aparameters = squeeze(div14_generative{1}.top_parameters_mean(n,model,:,:));
bparameters = squeeze(div14_generative{2}.top_parameters_mean(n,model,:,:));
k = nan(nnet(2),2,2); k(1:nnet(1),:,1) = aparameters; k(:,:,2) = bparameters;
eta = nan(nnet(2),2); 
eta(1:nnet(1),1) = aparameters(:,1);
eta(:,2) = bparameters(:,1);
gam = nan(nnet(2),2); 
gam(1:nnet(1),1) = aparameters(:,2);
gam(:,2) = bparameters(:,2);
% statistics on wiring parmeters
p = []; d = [];
p(1) = ranksum(eta(:,1),eta(:,2));
d(1) = computeCohen_d(eta(:,1),eta(:,2));
p(2) = ranksum(gam(:,1),gam(:,2));
d(2) = computeCohen_d(gam(:,1),gam(:,2));
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
    ylabel('Magnitude'); 
    xticks([1,2]);
    ylim([-1.5 1.5]);
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
    u = iosr.statistics.boxPlot(k(:,parameter,:),...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor','w',...
        'boxAlpha',1);
    ylabel(ylabelsg{parameter}); 
    title(sprintf('p=%g, d=%g',p(parameter),d(parameter)));
    xticklabels('Comparison');
    b = gca;
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 25;
end
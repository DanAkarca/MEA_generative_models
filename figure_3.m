%% figure 3
% written by danyal akarca
%% load data
clear; clc;
% set directory
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/50k_rodent_dev';
% load sttc
load(strcat(gendir,'/M03212_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div7 = gm_data;
load(strcat(gendir,'/M03212_rec2_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div10 = gm_data;
load(strcat(gendir,'/M03212_rec3_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div12 = gm_data;
load(strcat(gendir,'/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div14 = gm_data;
rodent_data_all = {div7 div10 div12 div14};
% load generative model data
load(strcat(gendir,'/tracking_rec1/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div7 = results_generative_models;
load(strcat(gendir,'/tracking_rec2/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div10 = results_generative_models;
load(strcat(gendir,'/tracking_rec3/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div12 = results_generative_models;
load(strcat(gendir,'/tracking_rec4/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div14 = results_generative_models;
% collect the data together
data = {div7 div10 div12 div14};
clear div7 div10 div12 div14;
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
% addpath of the ks statisics
addpath('/imaging/astle/users/da04/PhD/toolboxes/voronoi');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
% addpath of rounding
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
% addpath of the processing code
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022/Code/Code2share/');
%% load the simulated networks into the genreative model data
% set hyperparameters of the dataset
% set number of cultures
ncultures = 6;
% set number of divs
ndivs = 4;
% set number of models
nmodels = 13;
% set number of parameters
nparams = 20000;
% label 4 time points
div_labels = string({'7','10','12','14'});
% label 6 hd cultures
culture_labels = string(1:6);
% loop over data
for div = 1:ndivs;
    % loop over cultures
    networks = {};
    for culture = 1:ncultures;
        % initialise
        networks{culture} = [];
        % keep the data
        networks{culture} = [];
        % loop over models
        for model = 1:nmodels;
            load(strcat(gendir,sprintf(...
                '/tracking_rec%g/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_%g_generative_model_%g.mat',...
                div,culture,div,model)));
            networks{culture}(:,:,model) = output.networks;
        end
    end
    % assign
    data{div}.networks = networks;
    % display
    disp(sprintf('div %s simulated networks loaded',div_labels{div}));
end
%% remove nodes with no connections from empirical data to keep for later
% initialise
rodent_data = rodent_data_all;
% loop over networks and remove nodes with no connections
for div = 1:ndivs;
    for culture = 1:ncultures;
        % get network
        Atgt = rodent_data_all{div}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = rodent_data_all{div}.table_wu{culture};
        xy = rodent_data_all{div}.table_xy{culture};
        % keep only nodes with edges
        rodent_data{div}.table_bu{culture} = Atgt(ids,ids);
        rodent_data{div}.table_wu{culture} = wu(ids,ids);
        rodent_data{div}.table_xy{culture} = xy(ids,:);
    end
end
%% visualise cultures over time
% set culture
culture = 1;
% set data to plot
dataplot = rodent_data;
% visualise over divs
h = figure; h.Position = [100 100 1200 350];
for div = 1:ndivs;
    subplot(1,4,div);
    A = dataplot{div}.table_bu{culture};
    g = graph(A);
    plot(g,...
        'XData',dataplot{div}.table_xy{culture}(:,2),...
        'YData',dataplot{div}.table_xy{culture}(:,1),...
        'markersize',degrees_und(A)*.75,...
        'nodelabel',[],...
        'edgecolor',[.6 .6 .6],...
        'nodecolor',[.4 .4 .4]);
    b = gca;
    axis off;
end
%% supplementary figure 2
% set which data to plot
dataplot = rodent_data_all;
% number of measures
nmeasures = 10;
% initialise
statistics = zeros(ncultures,ndivs,nmeasures);
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
    'Matching'});
% loop over networks
for culture = 1:ncultures;
    for div = 1:ndivs;
        % take the weighted network
        wnetwork = dataplot{div}.table_wu{culture};
        % take the binary network
        network = dataplot{div}.table_bu{culture};
        % take the coordinates
        coords = dataplot{div}.table_xy{culture};
        % compute the euclidean distance
        euclidean = squareform(pdist(coords));
        % compute global statistics
        % nnodes
        statistics(culture,div,1) = size(network,1);
        % total STTC
        sttc = triu(wnetwork,1);
        statistics(culture,div,2) = mean(sttc(:),'all');
        % density
        statistics(culture,div,3) = density_und(network);
        % total number of connections
        statistics(culture,div,4) = sum(degrees_und(network));
        % global efficiency
        statistics(culture,div,5) = efficiency_bin(network,0);
        % betweenness
        statistics(culture,div,6) = mean(betweenness_bin(network));
        % clustering
        statistics(culture,div,7) = mean(clustering_coef_bu(network));
        % average weighted length
        statistics(culture,div,8) = mean(mean(euclidean.*network));
        % modularity
        [~,statistics(culture,div,9)] = modularity_und(network);
        % homophily
        statistics(culture,div,10) = mean(matching_ind(network),'all');
        % display
        disp(sprintf('computed statistics for culture %g div %s',...
            culture,div_labels(div)));
    end
end

%%% plot statistics over time %%%
% stat, e.g. density 3 (note k * 100), degree 4
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
    'symbolColor','k');
xticklabels(div_labels);
xlabel('Days {\itin vitro} (DIV)');
ylabel(statistics_labels(stat));
b = gca; b.TickDir = 'out';
b.FontName = 'Arial'; 
b.FontSize = 25;
box off;
% tabulate
med_tab = median(k);
iqr_tab = iqr(k);
stat_tab = table(med_tab,iqr_tab);
% columns are divs, rows are pipelines
disp(stat_tab);

%%% plot proportion of edge lengths over time %%%
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binsum = []; binprop = [];
% compute
for div = 1:ndivs;
    for culture = 1:ncultures;
    % network
    a = dataplot{div}.table_bu{culture}; 
    % euclidean
    b = squareform(pdist(dataplot{div}.table_xy{culture}));
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
        binprop(culture,div,k-1) = prop;
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
for div = 1:ndivs;
    uu = errorbar(binvec(2:end)+j*div,x(div,:),y(div,:)); 
    uu.Color = colbar(div,:); uu.Color(4) = .2;
    uu.Marker = '.'; 
    uu.MarkerSize = 60;
    uu.LineWidth = 5;
    hold on;
end
xlim([0 binlim+2*j]);
ylim([0 inf]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Euclidean (\mum)'); ylabel('Prop. connections (%)');
legend({'DIV7','DIV10','DIV12','DIV14'},'box','off'); box off;
%% from the energy and parameter matrices
% form energy matrix
energy = zeros(ncultures,nmodels,nparams,ndivs);
parameters = zeros(ncultures,nmodels,nparams,2,ndivs);
% place in the divs
for div = 1:ndivs;
    energy(:,:,:,div) = data{div}.energy;
    parameters(:,:,:,:,div) = data{div}.parameters;
end
%% compute top performing simulations
% define the top n parameters
nset = [1 10 50 100];
% initialise
top_energy = cell(length(nset),ndivs);
top_energy_mean = zeros(length(nset),13,ncultures,ndivs);
top_parameters = cell(length(nset),ndivs);
top_parameters_mean = zeros(length(nset),13,ncultures,2,ndivs);
% compute the minimum
for no = 1:length(nset);
    % take the actual amount of top performing parameters
    n = nset(no);
    % loop over divs
    for div = 1:ndivs;
        % loop over models
        for model = 1:nmodels;
            % take energies for this model
            pipeline_d = squeeze(energy(:,model,:,div))';
            % rank them for each subject
            [v i] = sort(pipeline_d);
            % take top n energies and their indices
            n_e = v(1:n,:);
            n_i = i(1:n,:);
            % take the corresponding parameters
            u = zeros(ncultures,n,2);
            for s = 1:ncultures;
                % keep parameters
                u(s,:,:) = squeeze(parameters(s,model,n_i(:,s),:,div));
            end
            % if top parameter only
            if n == 1
                % squeeze the matrices
                u = squeeze(u);
                % assign
                top_energy{no,div}(model,:) = n_e';
                top_parameters{no,div}(model,:,:) = u;
                % and assign it to the mean
                top_energy_mean(no,model,:,div) = n_e';
                top_parameters_mean(no,model,:,:,div) = u;
            else
                top_energy{no,div}(model,:,:) = n_e';
                top_parameters{no,div}(model,:,:,:) = u;
                % keep a mean value too
                top_energy_mean(no,model,:,div) = squeeze(mean(n_e',2));
                top_parameters_mean(no,model,:,:,div) = squeeze(mean(u,2));
            end
        end
    end
    disp(sprintf('set %g of %g complete',no,length(nset)));
end
%% visualise energy over time
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
h = figure; h.Position = [100 100 1200 450];
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
ylim([0 0.6]); ylabel('Energy'); xlabel('Days {\itin vitro} (DIV)');
xticklabels({'7','10','12','14'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% iosr boxplot - grouped models
energy_div_rules = nan(ncultures*5,ndivs,4);
a = squeeze(data_plot(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); energy_div_rules(1:ncultures*2,:,1) = a';
b = squeeze(data_plot(:,:,[3:7])); b = permute(b,[2 1 3]); b = b(:,:); energy_div_rules(1:ncultures*5,:,2) = b';
c = squeeze(data_plot(:,:,[8:12])); c = permute(c,[2 1 3]); c = c(:,:); energy_div_rules(1:ncultures*5,:,3) = c';
d = squeeze(data_plot(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); energy_div_rules(1:ncultures,:,4) = d';
% visualise
h = figure; h.Position = [100 100 1200 450];
u = iosr.statistics.boxPlot(energy_div_rules,...
    'showViolin',logical(0),...
    'showScatter',logical(0),...
    'scatterAlpha',.8,...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor',[.5 .5 .5],...
    'boxColor',{'r','b','g','y'},...
    'boxAlpha',0.2);
u.scatterColor = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]}';
ylim([0 0.6]); ylabel('Energy'); xlabel('Days {\itin vitro} (DIV)');
xticklabels({'7','10','12','14'});
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute energy statistics between rules for each time point
% set the data
data = energy_div_rules;
% ndivs
ndivs = 4;
% logical for saving
saveTable = 0;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/rodent_50k_energy_rules_all.csv';
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
    sprintf('DIV7 p=%.3g',p(1)),'DIV7 Cohen d',...
    sprintf('DIV10 p=%.3g',p(2)),'DIV10 Cohen d',...
    sprintf('DIV12 p=%.3g',p(3)),'DIV12 Cohen d',...
    sprintf('DIV14 p=%.3g',p(4)),'DIV14 Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
%% plot the energy landscape
% set div
div = 4;
% form the modeltype
models = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','diff-prod'});
% select model
model = 3;
% plot the energy landscape
h = figure;
h.Position = [100 100 500 450];
% take the measure
e = squeeze(data{div}.energy(:,model,:));
p = squeeze(data{div}.parameters(:,model,:,:));
% visualise
if model == 1
    eta = squeeze(p(:,:,1));
    m = parula(10000);
    for culture = 1:ncultures;
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
    for culture = 1:ncultures;
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
%% compare parameters
% take top n parameters
n = 1;
% select the model
model = 3;
% visualise 
h = figure;
h.Position = [100 100 1200 300];
for div = 1:ndivs;
    % seperate
    params = squeeze(top_parameters_mean(n,model,:,:,div));
    subplot(1,4,div);
    u = iosr.statistics.boxPlot(params,...
        'showViolin',logical(0),...
        'theme','colorall',...
        'symbolMarker','x',...
        'symbolColor','k',...
        'boxColor',[.5 .5 .5],...
        'boxAlpha',1);
    if div == 1;
        ylabel('Magnitude'); 
    end
    xticks([1,2]);
    ylim([-4 4]);
    yline(0);
    xticklabels({'\eta','\gamma'});
    b = gca; 
    b.XAxis.TickDirection = 'out';
    b.YAxis.TickDirection = 'out';
    b.FontName = 'Arial';
    b.FontSize = 18;
    box off;
end
%% form a loop index of the cultures and keep the data
% number of cultures
ncultures = 6;
% number of divs
ndivs = 4;
% set data to use
dataplot = rodent_data_all;
% total number of networks
nsamp = ncultures * ndivs;
% initialise the indices: columns for pipeline, then culture, then div
index = zeros(nsamp,2); 
all_o_networks = {};
all_d_networks = {};
step = 1;
% loop through and visualise
for culture = 1:ncultures;
    for div = 1:ndivs;
        % form the index
        index(step,1) = culture;
        index(step,2) = div;
        % get the data as a cell array
        all_o_networks{step} = dataplot{div}.table_bu{culture};
        all_d_networks{step} = squareform(pdist(dataplot{div}.table_xy{culture}));
        % update the step
        step = step + 1;
    end
end
%% check observations and simlations have the same number of edges
% initailise
noconn = [];
for i = 1:nsamp;
    % empirical
    noconn(i,1) = nnz(all_o_networks{i})/2;
    % simulated
    noconn(i,2) = size(data{index(i,2)}.networks{index(i,1)},1);
end
%% model fits from bootstrapping
% set the div
div = 4;
% set number of stats computed
nstat = 6;
% get the criteria 
criteria = index(:,2)==div;
% number selected
nsel = sum(criteria);
% get indices
inda = find(criteria);
% set the model and title
model = 3;
title_place = 'Matching';
% set the number to permute
nperm = 99;
% get the energies of these cultures
e_net = squeeze(data{div}.energy(:,model,:));
% get their indices
[~,i_low] = sort(e_net','ascend');
% get only the top nperms
ind = i_low(1:nperm,:);
% initialise permutations
observed_statistics = {}; 
perm_statistics = {}; 
test = {};
% take the simulated networks for each culture
for culture = 1:nsel;
    % get the observed networks
    A = all_o_networks{inda(culture)};
    D = all_d_networks{inda(culture)};
    % get the statistics
    observed_statistics{culture}{1} = degrees_und(A);
    observed_statistics{culture}{2} = clustering_coef_bu(A)';
    observed_statistics{culture}{3} = betweenness_bin(A);
    observed_statistics{culture}{4} = D(triu(A,1)>0)';
    %observed_statistics{culture}{4} = sum(D.*triu(A,1))./10^3;
    % outside of energy
    observed_statistics{culture}{5} = efficiency_bin(A,1)';
    [ci,q] = community_louvain(A);
    observed_statistics{culture}{6} = participation_coef(A,ci)';
    % test
    test{culture}{1} = matching_ind(A);
    test{culture}{2} = degrees_und(A);
    % initialise
    perm_statistics{culture}=[];
    % loop over top nperm networks
    for i = 1:nperm;
        % index the network
        indi = ind(i,culture);
        % find this network
        b = squeeze(data{div}.networks{culture}(:,indi,model));
        % reconstruct the simulation network
        nnode = size(A,1);
        Asynth = zeros(nnode,nnode);
        Asynth(b) = 1;
        Asynth = Asynth+Asynth';
        % get the statistics
        perm_statistics{culture}{1}(i,:) = degrees_und(Asynth);
        perm_statistics{culture}{2}(i,:) = clustering_coef_bu(Asynth);
        perm_statistics{culture}{3}(i,:) = betweenness_bin(Asynth)';
        perm_statistics{culture}{4}(i,:) = D(triu(Asynth,1)>0);
        %perm_statistics{culture}{4}(i,:) = sum(D.*triu(Asynth,1))./10^3; % for scale;
        % outside of energy
        perm_statistics{culture}{5}(i,:) = efficiency_bin(Asynth,1);
        [ci,q] = community_louvain(Asynth);
        perm_statistics{culture}{6}(i,:) = participation_coef(Asynth,ci);
        % display
        disp(sprintf('culture %g permutation %g of %g complete',culture,i,nperm));
    end
end
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
h = figure; h.Position = [100 100 1800 350];
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
segregation    = {}; % segregation: modularity q statisti
integration    = {}; % integration: global efficiency
% set data to use
dataplot = rodent_data_all;
% set the top n parameters to average over  from nset
nn = 1;
% loop through networks
for div = 1:ndivs;
    for culture = 1:ncultures;
        tic  
        % set target networks
        Atgt_set    = dataplot{div}.table_bu;
        % seed network is empty
        A           = zeros(size(Atgt_set{culture}));
        % compute the connections in the seed
        mseed       = nnz(A)/2;
        % euclidean distance
        D           = squareform(pdist(dataplot{div}.table_xy{culture}));
        % set target
        Atgt        = Atgt_set{culture};
        % set the parameter for this subject
        params      = squeeze(top_parameters_mean(nn,3,culture,:,div))'; % **** <<< there is an error here somewhere, because the gamma is too high <<< ****
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
                disp(sprintf('network %g connection %g of %g formed',culture,ii,m));
            end
            % keep the data
            b(:,iparam)         = find(triu(A,1));
            matching_K{culture,div}     = Kall;                  
            matching_Fk{culture,div}    = Fkall;                
            matching_Ksum{culture,div}  = squeeze(sum(Kall,2)); 
            matching_Fksum{culture,div} = squeeze(sum(Fkall,2));
            matching_P{culture,div}     = Pall;                  
            matching_Psum{culture,div}  = squeeze(sum(Pall,2)); 
            matching_A{culture,div}     = Aall;                  
            segregation{culture,div}    = Seg;
            integration{culture,div}    = Int;
        end
        time = toc;
        disp(sprintf('generative model complete: network %g took %g seconds',culture,time));
    end
end
%% compute summary statistics over time
% number of measures
nmeasures = 7;
% initialise
statistics = zeros(ncultures,ndivs,nmeasures);
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
for culture = 1:ncultures;
    for div = 1:ndivs;
        % take the binary network
        network = dataplot{div}.table_bu{culture};
        % take the coordinates
        coords = dataplot{div}.table_xy{culture};
        % compute the euclidean distance
        euclidean = squareform(pdist(coords));
        % compute global statistics
        % density
        statistics(culture,div,1) = density_und(network);
        % total number of connections
        statistics(culture,div,2) = sum(degrees_und(network));
        % global efficiency
        statistics(culture,div,3) = efficiency_bin(network);
        % betweenness
        statistics(culture,div,4) = mean(betweenness_bin(network));
        % clustering
        statistics(culture,div,5) = mean(clustering_coef_bu(network));
        % average weighted length
        statistics(culture,div,6) = mean(mean(euclidean.*network));
        % modularity
        [~,statistics(culture,div,7)] = modularity_und(network);
    end
end
%% scale simulated segregation against observed segregation
% take the observed modularity
stat = 7;
% visualise the number of connections
statv = squeeze(statistics(:,:,stat));
% alter the obseved statistics to be at specific points
plotdatas = nan(size(statv,1),16);
plotdatas(:,7)=statv(:,1);
plotdatas(:,10)=statv(:,2);
plotdatas(:,12)=statv(:,3);
plotdatas(:,14)=statv(:,4);
% get the 14 day networks of a particular pipeline
div = 4;
% get the simulated segregation data of these networks
segregationi = segregation(:,div);
% form into a matrix
u = [];
for k = 1:length(segregationi);
    u(k) = size(segregationi{k},1);
end
[v,~] = max(u);
j = nan(length(segregationi),v);
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
box off; axis square;
%% get the segregation 
% take the observed segregation
ind = find(sum(isnan(plotdatas)));
observed_segregation = plotdatas;
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
% take the observed integration
stat = 3;
% visualise the number of connections
statv = squeeze(statistics(:,:,stat));
% alter the obseved statistics to be at specific points
plotdatai = nan(size(statv,1),16);
plotdatai(:,7)=statv(:,1);
plotdatai(:,10)=statv(:,2);
plotdatai(:,12)=statv(:,3);
plotdatai(:,14)=statv(:,4);
% get the 14 day networks of a particular pipeline
div = 4;
% get the simulated segregation data of these networks
integrationi = integration(:,div);
% form into a matrix
u = [];
for k = 1:length(integrationi);
    u(k) = size(integrationi{k},1);
end
[v,~] = max(u);
j = nan(length(integrationi),v);
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
box off; axis square;
%% get the observed and simulated integration data
% take the observed integration
ind = find(sum(isnan(plotdatai)));
observed_integration = plotdatai;
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
% compute correlations
a = squeeze(integration_comparison(:,:,1));
b = squeeze(integration_comparison(:,:,2));
c = squeeze(segregation_comparison(:,:,1));
d = squeeze(segregation_comparison(:,:,2));
% integration
h = figure; h.Position = [100 100 400 400];
[ra pa] = corr(a(:),b(:));
facecolor = 1-copper(4);
for div = 1:4;
scatter(a(:,div),b(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
sgtitle(sprintf('R^2=%.3g, r=%.3g, p=%.3d',ra^2*100,ra,pa));
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0 0.4]); ylim([0 0.4]); xticks([0:0.1:0.4]);
hold on;
end
% keep the variance explained
integration_varex = [ra^2*100 ra pa];
axis square;
% segregation
h = figure; h.Position = [100 100 400 400];
[rb pb] = corr(c(:),d(:));
facecolor = copper(4);
for div = 1:4;
scatter(c(:,div),d(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
sgtitle(sprintf('R^2=%.3g, r=%.3g, p=%.3d',rb^2*100,rb,pb));
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0.1 0.9]); ylim([0.1 0.9]); 
hold on;
end
% keep the variance explained
segregation_varex = [rb^2*100 rb pb];
axis square;
%% keep data
rodent_50k_trajectory_data = struct;
rodent_50k_trajectory_data.efficiency_comparison = integration_comparison;
rodent_50k_trajectory_data.efficiency_trajectory.simulated = integrationi;
rodent_50k_trajectory_data.efficiency_trajectory.observed = plotdatai;
rodent_50k_trajectory_data.efficiency_variance = integration_varex;
rodent_50k_trajectory_data.modularity_comparison = segregation_comparison;
rodent_50k_trajectory_data.modularity_trajectory.simulated = segregationi;
rodent_50k_trajectory_data.modularity_trajectory.observed = plotdatas;
rodent_50k_trajectory_data.modularity_variance = segregation_varex;
rodent_50k_trajectory_data.info = string({'6 cultures x 4 divs x 2 conditions (observed v simulated). The simulated statistics were taken from the best fitting 14 day simulated networks at ~50%, 71%, 86% and 100% to correspond to 7, 10, 12 and 14 days.'});
rodent_50k_trajectory_data.author = string({'computed by danyal akarca, university of cambridge, 08/22'});
%% set data from previous run
integration_comparison = rodent_50k_trajectory_data.efficiency_comparison;
integrationi = rodent_50k_trajectory_data.efficiency_trajectory.simulated;
plotdatai = rodent_50k_trajectory_data.efficiency_trajectory.observed;
segregation_comparison = rodent_50k_trajectory_data.modularity_comparison;
segregationi = rodent_50k_trajectory_data.modularity_trajectory.simulated;
plotdatas = rodent_50k_trajectory_data.modularity_trajectory.observed;
%% visualise all the data
% visualise the simulations of these scaled simulations
h = figure;
h.Position = [100 100 1800 350];
subplot(1,4,1);
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
    'boxcolor',[202 143 66]./256,...
    'symbolMarker','x',...
    'symbolColor','k');
xticks([0 7 14]);
xlim([-1 15]);
xticklabels({'0%','50%','100%'});
ylabel('Modularity Q'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
box off; axis square;
% scale simulated integration against observed integration
% visualise the simulations of these scaled simulations
subplot(1,4,3);
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
    'boxcolor',[102 167 197]./256,...
    'symbolMarker','x',...
    'symbolColor','k');
xticks([0 7 14]);
xlim([-1 15]);
xticklabels({'0%','50%','100%'});
ylabel('Global Efficiency'); xlabel('Simulated time');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
box off; axis square;
% compute variance explained
% compute correlations
a = squeeze(integration_comparison(:,:,1));
b = squeeze(integration_comparison(:,:,2));
c = squeeze(segregation_comparison(:,:,1));
d = squeeze(segregation_comparison(:,:,2));
% integration
subplot(1,4,4);
[r p] = corr(a(:),b(:));
facecolor = 1-copper(4);
for div = 1:4;
scatter(a(:,div),b(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0 0.4]); ylim([0 0.4]); xticks([0:0.1:0.4]);
hold on;
end
axis square;
% segregation
subplot(1,4,2);
[r p] = corr(c(:),d(:));
facecolor = copper(4);
for div = 1:4;
scatter(c(:,div),d(:,div),200,'o','markerfacecolor',facecolor(div,:),'markeredgecolor','k');
xlabel('Observed'); ylabel('Simulated');
u = gca; u.TickDir = 'out'; u.FontSize = 25; u.FontName = 'Arial'; 
xlim([0.1 0.9]); ylim([0.1 0.9]); 
hold on;
end
axis square;
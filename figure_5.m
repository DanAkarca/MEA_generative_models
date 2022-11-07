%% figure 5
% written by danyal akarca
%% load data
clear; clc;
% set directory
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/50k_rodent_dev';
% load 50k sttc
load(strcat(gendir,'/M03212_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div7 = gm_data;
load(strcat(gendir,'/M03212_rec2_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div10 = gm_data;
load(strcat(gendir,'/M03212_rec3_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div12 = gm_data;
load(strcat(gendir,'/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div14 = gm_data;
rodent_data_all = {div7 div10 div12 div14};
% load 50k generative model data
load(strcat(gendir,'/tracking_rec1/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div7 = results_generative_models;
load(strcat(gendir,'/tracking_rec2/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div10 = results_generative_models;
load(strcat(gendir,'/tracking_rec3/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div12 = results_generative_models;
load(strcat(gendir,'/tracking_rec4/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div14 = results_generative_models;
% collect the data together
data = {div7 div10 div12 div14};
clear div7 div10 div12 div14;
gendir100 = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022_updated/100k_rodent_dev';
% load 100k sttc
load(strcat(gendir100,'/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div14 = gm_data;
load(strcat(gendir100,'/100k_combined_rec2_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div21 = gm_data;
load(strcat(gendir100,'/100k_combined_rec3_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div28 = gm_data;
rodent_data_100_all = {div14 div21 div28};
% load 100k generative model data
load(strcat(gendir100,'/tracking_rec1/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div14 = results_generative_models;
load(strcat(gendir100,'/tracking_rec2/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div21 = results_generative_models;
load(strcat(gendir100,'/tracking_rec3/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/combined/results_generative_models.mat')); div28 = results_generative_models;
data_100 = {div14 div21 div28};
% load the topological fingerprint dissimilarity
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/new/tf_dissimilarity.mat');
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the ks statisics
addpath('/imaging/astle/users/da04/PhD/toolboxes/voronoi');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
% addpath of rounding
addpath('/imaging/astle/users/da04/PhD/toolboxes/roundsd');
% addpath of cohen's d
addpath('/imaging/astle/users/da04/PhD/toolboxes/computeCohen');
% addpath of the processing code
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022/Code/Code2share/');
%% load the simulated networks into the genreative model data
% set hyperparameters of the 50k dataset
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
    disp(sprintf('50k div %s simulated networks loaded',div_labels{div}));
end
% set hyperparameters of the 100k dataset
% set number of cultures
ncultures = 12;
% set number of divs
ndivs = 3;
% set number of models
nmodels = 13;
% set number of parameters
nparams = 20000;
% label 4 time points
div_labels = string({'14','21','28'});
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
            load(strcat(gendir100,sprintf(...
                '/tracking_rec%g/min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob_allin/rodent_1_%g_%g_generative_model_%g.mat',...
                div,culture,div,model)));
            networks{culture}(:,:,model) = output.networks;
        end
    end
    % assign
    data_100{div}.networks = networks;
    % display
    disp(sprintf('100k div %s simulated networks loaded',div_labels{div}));
end
%% remove nodes with no connections from empirical data to keep for later
% 50k
% initialise
rodent_data = rodent_data_all;
% settings
ndivs = 4;
ncultures = 6;
% loop over networks and remove nodes with no connections
for div = 1:ndivs;
    for culture = 1:ncultures;
        % get network
        Atgt = rodent_data_all{div}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = rodent_data{div}.table_wu{culture};
        xy = rodent_data{div}.table_xy{culture};
        % keep only nodes with edges
        rodent_data{div}.table_bu{culture} = Atgt(ids,ids);
        rodent_data{div}.table_wu{culture} = wu(ids,ids);
        rodent_data{div}.table_xy{culture} = xy(ids,:);
    end
end
% 100k
% initialis
rodent_data_100 = rodent_data_100_all;
% settings
ndivs = 3;
ncultures = 12;
% loop over networks and remove nodes with no connections
for div = 1:ndivs;
    for culture = 1:ncultures;
        % get network
        Atgt = rodent_data_100_all{div}.table_bu{culture};
        % find nodes with edges
        ids = find(degrees_und(Atgt));
        % get other data
        wu = rodent_data_100{div}.table_wu{culture};
        xy = rodent_data_100{div}.table_xy{culture};
        % keep only nodes with edges
        rodent_data_100{div}.table_bu{culture} = Atgt(ids,ids);
        rodent_data_100{div}.table_wu{culture} = wu(ids,ids);
        rodent_data_100{div}.table_xy{culture} = xy(ids,:);
    end
end
%% set hyperparameters for tfd analysis
% set which data to analyse 
obs_data = {rodent_data_all rodent_data_100_all};
gm_data = {data data_100};
% set number of groups
ngroups = 2;
% set number of cultures
ncultures = [6 12];
% set how many measures are computed
nmeasures = 6;
% set how many models are computed
nmodels = 13;
% set up to how many divs are considered (e.g. 4 or 3)
ndivs = [4 3];
%% run the tf calculations
% initialise
for group = 1:ngroups
    a = zeros(ndivs(group),ncultures(group),nmeasures,nmeasures);
    b = zeros(ndivs(group),ncultures(group),nmodels,nmeasures,nmeasures);
    c = zeros(ndivs(group),ncultures(group),nmodels,nmeasures*2,nmeasures*2);
    d = zeros(ndivs(group),ncultures(group),nmodels);
    corr_local_observed{group} = a;
    corr_local_simulated{group} = b;
    corr_local_together{group} = c;
    todissim{group} = d;
end
% labels
var_labels = string({...
    'obs degree','obs clustering','obs betweenn','obs length','obs eff','obs match',...
    'sim degree','sim clustering','sim between','sim length','sim eff','sim match'});
% loop over groups
for group = 1:ngroups;
    % loop over divs
    for div = 1:ndivs(group);
    % compute observed local statistics for each network
        for culture = 1:ncultures(group);
            % take the observed network
            w = obs_data{group}{div}.table_bu{culture};
            % take the cost
            d = squareform(pdist(obs_data{group}{div}.table_xy{culture}));
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
            corr_local_observed{group}(div,culture,:,:) = corr(observed);
            for model = 1:nmodels
                % loop over models
                disp(sprintf('group %g evaluating day %g, culture %g, generative model %g...',group,div,culture,model)); 
                % take the optimal network simulation parameters, given the set model, for this network
                E = gm_data{group}{div}.energy;
                N = gm_data{group}{div}.networks;
                e_net = squeeze(E(culture,model,:));
                [~,i_low] = min(e_net);
                % take the end simulated network
                si = squeeze(N{culture}(:,i_low,model));
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
                corr_local_simulated{group}(div,culture,model,:,:) = corr(simulated);
                % form a matrix together and correlate these
                corr_local_together{group}(div,culture,model,:,:) = corr([observed simulated]);
                % compute the topological organization dissimilarity
                todissim{group}(div,culture,model) = norm(corr(observed)-corr(simulated));
            end
        end
    end
end
% keep this data
tf_dissimilarity = struct;
tf_dissimilarity.observed_tf = corr_local_observed;
tf_dissimilarity.simulated_tf = corr_local_simulated;
tf_dissimilarity.together_tf = corr_local_together;
tf_dissimilarity.tfdissim = todissim;
tf_dissimilarity.info = '50k (div7,10,12,14), 100k (div14,21,28), div x culture x (model) x measure x measure, written by danyal akarca 08/2022';
%% visualise the tfd
% set the group
group = 1;
% settings
ndivs = 4;
ncultures = 6;
daylabels = {'7','10','12','14'};
% take data
tfd = tf_dissimilarity.tfdissim{group};
% set model order for visualisation
i = [2:3 9:13 4:8 1];
%  alter the order
todissim_permute = permute(tfd,[2 1 3]);
% make zeros into nans
todissim_permute(todissim_permute==0)=NaN;
% plot
h = figure; h.Position = [100 100 1200 450];
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
xticklabels(daylabels);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% group meaned todissimilaity into rules
% all
to_div_rules = nan(ncultures*5,ndivs,4);
a = squeeze(todissim_permute(:,:,[2:3])); a = permute(a,[2 1 3]); a = a(:,:); to_div_rules(1:ncultures*2,:,1) = a';
b = squeeze(todissim_permute(:,:,[9:13])); b = permute(b,[2 1 3]); b = b(:,:); to_div_rules(1:ncultures*5,:,2) = b';
c = squeeze(todissim_permute(:,:,[4:8])); c = permute(c,[2 1 3]); c = c(:,:); to_div_rules(1:ncultures*5,:,3) = c';
d = squeeze(todissim_permute(:,:,[1])); d = permute(d,[2 1 3]); d = d(:,:); to_div_rules(1:ncultures,:,4) = d';
%% compute topological fingerprint statistics between rules for each time point
% set the data
data = to_div_rules;
% ndivs - alligningwith the data
ndivs = 4;
% set group - 1 is 50k, 2 is 100k
group = 1; 
% logical for saving
saveTable = 0;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/rodent_50k_tfd_rules_all.csv';
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
if group == 1;
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
end
% form table
if group == 2;
t = table(...
    comp_rules(:,1),...
    comp_rules(:,2),...
    p_rules(:,1),...
    d_rules(:,1),...
    p_rules(:,2),...
    d_rules(:,2),...
    p_rules(:,3),...
    d_rules(:,3),...
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('DIV14 p=%.3g',p(1)),'DIV14 Cohen d',...
    sprintf('DIV21 p=%.3g',p(2)),'DIV21 Cohen d',...
    sprintf('DIV28 p=%.3g',p(3)),'DIV28 Cohen d'})
    if saveTable==1;
        writetable(t,tablename);
    end
end
%% visualise a specific comparison in terms of a culture and averaged over cultures
% set group
group = 1;
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
specific_obs = squeeze(tf_dissimilarity.observed_tf{group}(div,culture,:,:));
specific_sim = squeeze(tf_dissimilarity.simulated_tf{group}(div,culture,model,:,:));
% visualise the specific culture
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1);
imagesc(specific_obs); title('Observed'); caxis([-1 1]); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r';  xticks([]); yticks([]);
subplot(1,2,2); 
imagesc(specific_sim); title('Simulated'); caxis([-1 1]); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r';  xticks([]); yticks([]);
sgtitle(sprintf('Culture %g, DIV %g, Model %g, TFD %g',culture,div,model,squeeze(tf_dissimilarity.simulated_tf{group}(div,culture,model))));
% take the averaged data
averaged_obs = squeeze(mean(squeeze(tf_dissimilarity.observed_tf{group}(div,:,:,:)),1));
averaged_sim = squeeze(mean(squeeze(tf_dissimilarity.simulated_tf{group}(div,:,model,:,:)),1));
% visualise the averaged data
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1); imagesc(averaged_obs); title('Observed'); caxis([-1 1]); u = colorbar; u.TickDirection = 'out'; u.Label.String = 'r'; xticks([]); yticks([]);
subplot(1,2,2); imagesc(averaged_sim); title('Simulated'); caxis([-1 1]); u = colorbar;  u.TickDirection = 'out'; u.Label.String = 'r'; xticks([]); yticks([]);
sgtitle(sprintf('Averaged DIV %g, Model %g, TFD %g',div,model,squeeze(mean(tf_dissimilarity.simulated_tf{group}(div,:,model),2))));
xticks([]); yticks([]);
%% plot a specific plot
g = figure; g.Position = [100 100 900 700];
plotind = [1 6 11 16; 2 7 12 17; 3 8 13 18; 4 9 14 19; 5 10 15 20];
for div = 1:ndivs
    % observed
    subplot(4,5,plotind(1,div)); imagesc(squeeze(mean(tf_dissimilarity.observed_tf{group}(div,:,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % matching
    subplot(4,5,plotind(2,div)); imagesc(squeeze(mean(tf_dissimilarity.simulated_tf{group}(div,:,3,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % clu-avg
    subplot(4,5,plotind(3,div)); imagesc(squeeze(mean(tf_dissimilarity.simulated_tf{group}(div,:,9,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % deg-avg
    subplot(4,5,plotind(4,div)); imagesc(squeeze(mean(tf_dissimilarity.simulated_tf{group}(div,:,4,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
    % sptl
    subplot(4,5,plotind(5,div)); imagesc(squeeze(mean(tf_dissimilarity.simulated_tf{group}(div,:,1,:,:),2,'omitnan'))); caxis([-1 1]); xticks([]); yticks([]);
end
%% organise the generative model findings
% set number of models
nmodels = 13;
% set the sample size for the 50k and 100k cultures
ngroup = [6 12];
% define the ndiv
ndiv = [4 3];
% define the top n parameters
nset = [1 10 50 100];
% intiailise
density_generative = {};
% compute the minimum
for group = 1:length(ngroup);
    for div = 1:ndiv(group);
        % get the group and div
        energy = gm_data{group}{div}.energy;
        parameters = gm_data{group}{div}.parameters;
        % initialise
        top_energy = cell(length(nset),1);
        top_energy_mean = zeros(length(nset),13,ngroup(group));
        top_parameters = cell(length(nset),1);
        top_parameters_mean = zeros(length(nset),13,ngroup(group),2);
        % ouput
        density_generative{group}{div} = struct;
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
        density_generative{group}{div}.top_energy = top_energy;
        density_generative{group}{div}.top_energy_mean = top_energy_mean;
        density_generative{group}{div}.top_parameters = top_parameters;
        density_generative{group}{div}.top_parameters_mean = top_parameters_mean;
    end
end
%% visualise energy over time
% set the top n parameters
ncomb = 1;
% set group
group = 2;
% set divs
ndivs = 3;
divlabels = {'14','21','28'};
% number of cultures
ncultures = [6 12];
% form all models
data_plot = [];
for div = 1:ndivs;
    data_plot(:,:,div) = squeeze(density_generative{group}{div}.top_energy_mean(ncomb,:,:));
end
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
xticklabels(divlabels);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
% iosr boxplot - grouped models
energy_div_rules = nan(ncultures(group)*5,ndivs,4);
a = squeeze(data_plot(:,:,[1:2])); a = permute(a,[2 1 3]); a = a(:,:); energy_div_rules(1:ncultures(group)*2,:,1) = a';
b = squeeze(data_plot(:,:,[3:7])); b = permute(b,[2 1 3]); b = b(:,:); energy_div_rules(1:ncultures(group)*5,:,2) = b';
c = squeeze(data_plot(:,:,[8:12])); c = permute(c,[2 1 3]); c = c(:,:); energy_div_rules(1:ncultures(group)*5,:,3) = c';
d = squeeze(data_plot(:,:,[13])); d = permute(d,[2 1 3]); d = d(:,:); energy_div_rules(1:ncultures(group),:,4) = d';
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
xticklabels(divlabels);
b = gca; 
b.XAxis.TickDirection = 'out';
b.YAxis.TickDirection = 'out';
b.FontName = 'Arial';
b.FontSize = 25;
%% compute energy statistics between rules for each time point
% set the data
data = energy_div_rules;
% ndivs
ndivs = 3;
% logical for saving
saveTable = 1;
% set name of table if saved
tablename = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/rodent_100k_energy_rules_all.csv';
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
    'VariableNames',...
    {'Rule A','Rule B',...
    sprintf('DIV14 p=%.3g',p(1)),'DIV14 Cohen d',...
    sprintf('DIV21 p=%.3g',p(2)),'DIV21 Cohen d',...
    sprintf('DIV28 p=%.3g',p(3)),'DIV28 Cohen d'});
if saveTable==1;
    writetable(t,tablename);
end
%% correlate energy and the to matrix across divs
% set ndivs and ncultures again
ndiv = [4 3];
ncultures = [6 12];
% set the group
group = 1;
% take the energy
ep = [];
for div = 1:ndiv(group);
    ep(div,:,:) = density_generative{group}{div}.top_energy_mean(1,:,:);
end
% permute ep
ep = permute(ep,[1 3 2]);
% take the tfd
to = squeeze(tf_dissimilarity.tfdissim{group});
% form indices for rules colours
ind = [ones(ncultures(group)*1,1);...
    2*ones(ncultures(group)*2,1);...
    3*ones(ncultures(group)*5,1);...
    4*ones(ncultures(group)*5,1)];
if group == 1;
ind = [ind;ind;ind;ind];
else if group == 2;
        ind = [ind;ind;ind];
    end
end
colsim = {'y','r','g','b'};
% take data
x = ep(:); y = to(:);
% across all models
[r p] = corr(x(y>0),y(y>0));
% keep individual correlations based on rules
rrule  = []; prule = [];
% visualise
figure;
for rule = 1:4;
    [rrule(rule) prule(rule)] = corr(x(ind==rule),y(ind==rule));
    h = scatter(x(ind==rule),y(ind==rule),70,'o',...
        'markerfacecolor',colsim{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor','k');
    hold on;
end
xlabel('Energy'); ylabel('{\itTFdissimilarity}');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
%legend({'Spatial','Homophily','Clustering','Degree'},'box','off');
xlabel('Energy'); ylabel('{\itTFdissimilarity}');
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
for col = 1:4;
    k(1).Children(col).MarkerFaceColor = scatcol_o(col,:);
end
% plot the average at each time point
eptd = []; totd = [];
for div = 1:ndiv(group);
    % energy
    eptd(div,:,1) = squeeze(ep(div,:,1));
    eptd(div,:,2) = mean(squeeze(ep(div,:,[2 3])),2)';
    eptd(div,:,3) = mean(squeeze(ep(div,:,[4:8])),2)';
    eptd(div,:,4) = mean(squeeze(ep(div,:,[9:13])),2)';
    % todissim
    totd(div,:,1) = squeeze(to(div,:,1));
    totd(div,:,2) = mean(squeeze(to(div,:,[2 3])),2)';
    totd(div,:,3) = mean(squeeze(to(div,:,[4:8])),2)';
    totd(div,:,4) = mean(squeeze(to(div,:,[9:13])),2)';
end
eptd = squeeze(mean(eptd,2)); totd = squeeze(mean(totd,2));
% visualise
figure;
shape = {'o','square','pentagram','diamond'};
for rule = 1:4;
    for div = 1:ndiv(group);
    scatter(eptd(div,rule),totd(div,rule),180,...
        'marker',shape{div},...
        'markerfacecolor',colsim{rule},...
        'markerfacealpha',0.5,...
        'markeredgecolor',[.5 .5 .5]);
    hold on;
    end
end
xlabel('Energy'); ylabel('TFdissimilarity');
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
%% correlate the tf dissimilarity and energy only for specific divs
% set div and group
div = 3;
group = 2;
% take the energy
ep = squeeze(density_generative{group}{div}.top_energy_mean(1,:,:))';
% take the tfd
to = squeeze(tf_dissimilarity.tfdissim{group}(div,:,:));
% form indices for rules colours
ind = [ones(ncultures(group)*1,1);...
    2*ones(ncultures(group)*2,1);...
    3*ones(ncultures(group)*5,1);...
    4*ones(ncultures(group)*5,1)];
% set colours (ordered to make work)
sptl_c = [255 255 102];
homo_c = [255 202 205];
clu_c = [206 253 205];
deg_c = [206 204 255];
scatcol = [sptl_c; homo_c; clu_c; deg_c]./256;
scatcol_o = flip(scatcol);
% take data
x = ep(:); y = to(:);
% across all models
[r p] = corr(x,y);
% keep individual correlations based on rules
rrule  = []; prule = [];
% visualise
h = figure; h.Position = [100 100 550 450];
% order to plot
plotorder = [4 3 2 1];
for r = 1:4;
    rule = plotorder(r);
    [rrule(rule) prule(rule)] = corr(x(ind==rule),y(ind==rule));
    h = scatter(x(ind==rule),y(ind==rule),150,'o',...
        'markerfacecolor',scatcol(rule,:),...
        'markerfacealpha',1,...
        'markeredgecolor','k');
    hold on;
end
% get total correlation
[rall pall] = corr(x,y);
b = gca; b.FontName = 'Arial'; b.FontSize = 25; b.TickDir = 'out';
%legend({'Spatial','Homophily','Clustering','Degree'},'box','off');
xlabel('Energy'); ylabel('{\itTFdissimilarity}');
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
for col = 1:4;
    k(1).Children(col).MarkerFaceColor = scatcol_o(col,:);
end
%% Exploration of 50k rodent primary cortical neuronal cultures
% written by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortical neurons
clear; clc;
% change directory to the project folder
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/');
% load sttc data
a = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag5_jitter5_prob.mat');
b = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag10_jitter10_prob.mat');
c = load('gm_data_shared_div7_10_12_14_min_r_0.01_alpha_0.001_min_rate_001hz_lag20_jitter20_prob.mat');
% collect the data together
data = {a.gm_data b.gm_data c.gm_data};
% addpath of the brain connectivity toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
% addpath of the iosr statistics toolbox
addpath('/imaging/astle/users/da04/PhD/toolboxes/MatlabToolbox-master/');
%% form labels
% 3 pipelines
pipeline_labels = string({'5','10','20'})';
% 4 time points
div_labels = string({'7','10','12','14'});
% 12 hd cultures
culture_labels = string(1:12);
%% set hyperparameters of the dataset
% number of pipelines
npipelines = 3;
% number of cultures
ncultures = 6;
% number of divs
ndivs = 4;
%% visualise the networks
% select a pipeline
pipeline = 2;
% select a culture
culture = 4;
% set divs to explore
divs = [1 2 3 4];
% compute the cultures
networks = data{pipeline}.table_bu_shared{culture};
coords = data{pipeline}.table_xy_shared{culture};
% visualise over divs
figure;
set(gcf,'Position',[100 100 1500 300]);
sgtitle(sprintf('culture %s %s: binary networks',culture_labels(culture),pipeline_labels(pipeline)));
% loop through divs
for i = 1:length(divs);
    div = divs(i);
    subplot(1,length(divs),i);
    g = graph(networks{div});
    h = plot(g,...
        'XData',coords{div}(:,1),...
        'YData',coords{div}(:,2),...
        'EdgeColor',[120 120 120]./256,...
        'LineWidth',1.5,...
        'NodeColor',[239 152 170]./256,...
        'NodeLabel',[],...
        'MarkerSize',0.45*degrees_und(networks{div})+0.01);
    title(div_labels(i));
    xticks([]); yticks([]);
    set(gca,'visible','off')
end
% visualise euclidean distances over divs
figure;
set(gcf,'Position',[100 100 1500 300]);
sgtitle(sprintf('culture %s %s: euclidean distance',culture_labels(culture),pipeline_labels(pipeline)));
% loop through divs
for i = 1:length(divs);
    subplot(1,length(divs),i);
    % compute euclidean distances
    imagesc(squareform(pdist(coords{i}))); caxis([0 5e3]);
    title(div_labels(i));
    xticks([]); yticks([]);
end
colorbar;
%% visualise across all cultures 
% set hyperparameters 
pipeline = 2;
divs = [1 2 3 4];
% loop over cultures
h = figure;
h.Position = [100 100 600 300];
m = [];
for culture = 1:ncultures;
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
        b = a(find(triu(ones(n),1)));
        % make into a ks density plot
        [f xi] = ksdensity(b);
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
yticks([]); b = gca; b.TickDir = 'out'; b.FontName = 'Arial';
xlabel('Euclidean \mum'); ylabel('Frequency');
b=gca; b.FontName = 'Arial'; b.FontSize = 25; box off; b.TickLength = [0 0];
% Plot the eucldiean distances of each culture
h = figure;
h.Position = [100 100 650 400];
n = [];
for culture = 1:ncultures;
    % subplot
    subplot(2,3,culture);
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
%% compute summary statistics over time
% number of measures
nmeasures = 10;
% initialise
statistics = zeros(npipelines,ncultures,ndivs,nmeasures);
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
for pipeline = 1:npipelines;
    for culture = 1:ncultures;
        for div = 1:ndivs;
            % take the weighted network
            wnetwork = data{pipeline}.table_wu_shared{culture}{div};
            % take the binary network
            network = data{pipeline}.table_bu_shared{culture}{div};
            % take the coordinates
            coords = data{pipeline}.table_xy_shared{culture}{div};
            % compute the euclidean distance
            euclidean = squareform(pdist(coords));
            % compute global statistics
            % nnodes
            statistics(pipeline,culture,div,1) = size(network,1);
            % total STTC
            sttc = triu(wnetwork,1);
            statistics(pipeline,culture,div,2) = mean(sttc(:),'all');
            % density
            statistics(pipeline,culture,div,3) = density_und(network);
            % total number of connections
            statistics(pipeline,culture,div,4) = sum(degrees_und(network));
            % global efficiency
            statistics(pipeline,culture,div,5) = efficiency_bin(network,0);
            % betweenness
            statistics(pipeline,culture,div,6) = mean(betweenness_bin(network));
            % clustering
            statistics(pipeline,culture,div,7) = mean(clustering_coef_bu(network));
            % average weighted length
            statistics(pipeline,culture,div,8) = mean(mean(euclidean.*network));
            % modularity
            [~,statistics(pipeline,culture,div,9)] = modularity_und(network);
            % homophily
            statistics(pipeline,culture,div,10) = mean(matching_ind(network),'all');
            % display
            disp(sprintf('computed statistics for %s culture %g div %s',...
                pipeline_labels(pipeline),culture,div_labels(div)));
        end
    end
end
%% plot number of connections and densities by pipeline and div
% pipeline
pipeline = 2;
% stat
stat = 3;
% visualise 
k = squeeze(statistics(pipeline,:,:,stat));
% plot
figure;
iosr.statistics.boxPlot(k,...
    'showViolin',logical(0),...
    'theme','colorall',...
    'symbolMarker','x',...
    'symbolColor','k');
xticklabels(div_labels);
xlabel('Days in vitro');
ylabel(statistics_labels(stat));
b = gca; b.TickDir = 'out';
b.FontName = 'Arial'; 
b.FontSize = 14;
% tabulate
med_tab = median(k);
iqr_tab = iqr(k);
stat_tab = table(med_tab,iqr_tab);
% columns are divs, rows are pipelines
disp(stat_tab);
%% plot a specific pipeline
% set pipeline
pipeline = 2;
% set stat
stat = 2;
% set colour
col = [31 43 133; 157 203 231; 107 32 132; 167 219 200; 102 199 167]./256;
colset = 5;
% take data 
k = squeeze(statistics(pipeline,:,:,stat));
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
% visualise
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
% tabulate
med_tab = median(k);
iqr_tab = iqr(k);
stat_tab = table(med_tab,iqr_tab);
% columns are divs, rows arepipelines
disp(stat_tab);
%% load and run the firing rate statistics ranova
% load mean firing rate
load('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/data/0.01Hz/table_rate_mean');
% compte statistics
culture_id = 1:size(k,1);
ndivs = 4;
rtbl = table(culture_id',table_rate_mean(:,1),table_rate_mean(:,2),table_rate_mean(:,3),table_rate_mean(:,4),...
    'VariableNames',{'Culture','DIV7','DIV10','DIV12','DIV14'});
stats = table([1:ndivs]','VariableNames',{'DIVs'});
% fit rm
rm = fitrm(rtbl,'DIV7-DIV14~Culture','WithinDesign',stats);
% run ranova
[ranovatbl a b c] = ranova(rm);
disp(ranovatbl);
%% plot connections over distance binned
% set pipeline 
pipeline = 2;
% loop over cultures and plot the binned distances
binlen = 400;
binlim = 4000;
binvec = [1 binlen:binlen:binlim];
bins = {}; binsum = []; binprop = [];
% compute
for div = 1:ndivs;
    for culture = 1:ncultures;
    % network
    a = data{pipeline}.table_bu_shared{culture}{div}; 
    % euclidean
    b = squareform(pdist(data{pipeline}.table_xy_shared{culture}{div}));
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
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
xlabel('Distance (\mum)'); ylabel('Prop. connections');
legend({'DIV7','DIV10','DIV12','DIV14'}); box off;
%% correlate statistics between pipelines
% visualise correlation matrices across sttc parameters
h = figure; h.Position = [100 100 900 400];
for stat = 1:nmeasures;
    subplot(2,5,stat);
    x = statistics(:,:,:,stat);
    y = x(:,:);
    r = corr(y');
    imagesc(r); caxis([-1 1]);
    xticklabels(pipeline_labels);
    yticks([]);
    title(statistics_labels(stat));
    b = gca; b.TickDir = 'out';
    b.FontName = 'Arial';
end
c = colorbar;
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
%% firing rate matrix
% set the type and culture
culture = 6; % 6, 3, 5
% set the well
well = 1;
% set the div 
div = 4;
% set the pipeline 
pipeline = 2;
% load the recording
load(sprintf('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ToShare_Jan22/rodent_50k_tracking/M03212/210811/well%g/rec%g/min_rate_001hz_lag%s_jitter%s_prob/sttc_data.mat',...
    well,div,pipeline_labels(pipeline),pipeline_labels(pipeline)));
% get spike times
spike_times = sttc_data.spks_data.raster;
% set bin
bin = 1;
% initialise
array = [];
% for each electrode bin the spikes
for electrode = 1:length(spike_times);
    % display
    disp(sprintf('electrode %g computing...',electrode));
    % get spike times for the electrode
    x = spike_times{electrode};
    % define the bins
    binlims = [0:bin:(ceil(max(x/10))*10)];
    [count,idx] = histc(x,binlims);
    % keep count
    array(electrode,1:length(count)) = count;
end
% visualise
h = figure; h.Position = [100 100 1000 600];
imagesc(array);
xlabel('Time (s)');
ylabel('Putative neurons');
c = colorbar; c.Label.String = 'Firing Rate (Hz)';
b = gca; b.FontName = 'Arial'; b.FontSize = 35; b.TickDir = 'out';
caxis([clims]);
xlim([vlim]);
% plot the summed spike vector
summed_spike = sum(array);
h = figure; h.Position = [100 100 1000 600];
plot(summed_spike,'k','linewidth',2)
xlim([0 length(summed_spike)]);
xlabel('Time (s)');
ylabel('Total spikes');
b = gca; b.FontName = 'Arial'; b.FontSize = 35; b.TickDir = 'out';
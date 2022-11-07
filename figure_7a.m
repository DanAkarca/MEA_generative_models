%% tSNE overview of 100k data
% written by Manuel Schroter 
% adapted by Danyal Akarca
%% set pre-requisites
clear all; close all
% set code path
code_path = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/code/';
% add code paths
addpath(genpath(strcat(code_path,'AC_clustering_figure6/')));
addpath(genpath(strcat(code_path,'CellExplorer-master/')));
addpath('/imaging/astle/users/da04/PhD/toolboxes/stdshade');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
% set data path
data_path = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Jan22/';
%% Ensure CellExplorer is compiled
%{
cd(strcat(code_path,'CellExplorer-master/'));
cd('calc_CellMetrics/mex');
% compile
mex -O CCGHeart.c
mex -O FindInInterval.c
%}
%% Rodent 100k neurons
list_rodent100k = [
{strcat(data_path,'rodent_100k_tracking/M02924/well1/')};
{strcat(data_path,'rodent_100k_tracking/M02924/well2/')};
{strcat(data_path,'rodent_100k_tracking/M02924/well3/')};
{strcat(data_path,'rodent_100k_tracking/M02924/well4/')};
{strcat(data_path,'rodent_100k_tracking/M02924/well5/')};
{strcat(data_path,'rodent_100k_tracking/M02924/well6/')};

{strcat(data_path,'rodent_100k_tracking/M03005/well1/')};
{strcat(data_path,'rodent_100k_tracking/M03005/well2/')};
{strcat(data_path,'rodent_100k_tracking/M03005/well3/')};
{strcat(data_path,'rodent_100k_tracking/M03005/well4/')};
{strcat(data_path,'rodent_100k_tracking/M03005/well5/')};
{strcat(data_path,'rodent_100k_tracking/M03005/well6/')};
];
% parameters for CCG
binSize = 1/1000; % 1 ms bins
duration = 1; % 50ms here corresponds to +/-25ms lag
fs = 1/10000;
for i = 1:12
           
    % concat all spk times
    load([list_rodent100k{i} 'spk_data.mat'])
    fprintf(['rodent culture: ' num2str(i),'\n'])
    
    spks_concat = [];
    for unit = 1:length(spk_data.spike_times)
        if unit>1
        spks_concat = [spks_concat; (spk_data.spike_times{unit})+max(spks_concat)];
        else
             spks_concat = [spks_concat; (spk_data.spike_times{unit})+0];
        end
    end
     
    [ccg,r] = CCG(spks_concat,ones(1,length(spks_concat)),'binSize',binSize,'duration',duration,'Fs',fs);
    rodent_table_r(i,:) = r;
    rodent_table_ccg(i,:) = (ccg);
    rodent_table_ccg_zscored(i,:) = zscore(ccg);
    rodent_table_ccg_unity(i,:) = mat2gray(ccg);
    
    clear spk_data spk_max spk_count table_ccg spk_max table_ccg spk_data
end
%% Human Gluta neurons
list_gluta = [
    {strcat(data_path,'human_ipsc_DIV28_30min/3008/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3019/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3020/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3022/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3026/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3062/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3064/210316/')};
];
% parameters for CCG
binSize = 1/1000; % 1 ms bins
duration = 1; % 50ms here corresponds to +/-25ms lag
fs = 1/20000;
for i = 1:7
           
    % concat all spk times
    load([list_gluta{i} 'spk_data.mat'])
    fprintf(['gluta culture: ' num2str(i),'\n'])
    
    spks_concat = [];
    for unit = 1:length(spk_data.spike_times)
        if unit>1
        spks_concat = [spks_concat; (spk_data.spike_times{unit})+max(spks_concat)];
        else
             spks_concat = [spks_concat; (spk_data.spike_times{unit})+0];
        end
    end
     
    [ccg,r] = CCG(spks_concat,ones(1,length(spks_concat)),'binSize',binSize,'duration',duration,'Fs',fs);
    gluta_table_r(i,:) = r;
    gluta_table_ccg(i,:) = (ccg);
    gluta_table_ccg_zscored(i,:) = zscore(ccg);
    gluta_table_ccg_unity(i,:) = mat2gray(ccg);
    
    clear spk_data spk_max spk_count table_ccg spk_max table_ccg
end
%% Human Moto neurons
list_moto = [
    {strcat(data_path,'human_ipsc_DIV28_30min/3009/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3011/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3012/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3017/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3060/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3061/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3063/210316/')};
];
% parameters for CCG
binSize = 1/1000; % 1 ms bins
duration = 1; % 50ms here corresponds to +/-25ms lag
fs = 1/20000;
for i = 1:7
           
    % concat all spk times
    load([list_moto{i} 'spk_data.mat'])
    fprintf(['moto culture: ' num2str(i),'\n'])
    
    spks_concat = [];
    for unit = 1:length(spk_data.spike_times)
        if unit>1
        spks_concat = [spks_concat; (spk_data.spike_times{unit})+max(spks_concat)];
        else
             spks_concat = [spks_concat; (spk_data.spike_times{unit})+0];
        end
    end
     
    [ccg,r] = CCG(spks_concat,ones(1,length(spks_concat)),'binSize',binSize,'duration',duration,'Fs',fs);
    moto_table_r(i,:) = r;
    moto_table_ccg(i,:) = (ccg);
    moto_table_ccg_zscored(i,:) = zscore(ccg);
    moto_table_ccg_unity(i,:) = mat2gray(ccg);
    
    clear spk_data spk_max spk_count table_ccg spk_max table_ccg
end
%% Human Dopa neurons
list_dopa = [
    {strcat(data_path,'human_ipsc_DIV28_30min/3015/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3018/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3023/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3058/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3059/210316/')};
    {strcat(data_path,'human_ipsc_DIV28_30min/3066/210316/')};
];
% parameters for CCG
binSize = 1/1000; % 1 ms bins
duration = 1; % 50ms here corresponds to +/-25ms lag
fs = 1/20000;
for i = 1:6
           
    % concat all spk times
    load([list_dopa{i} 'spk_data.mat'])
    fprintf(['dopa culture: ' num2str(i),'\n'])
    
    spks_concat = [];
    for unit = 1:length(spk_data.spike_times)
        if unit>1
        spks_concat = [spks_concat; (spk_data.spike_times{unit})+max(spks_concat)];
        else
             spks_concat = [spks_concat; (spk_data.spike_times{unit})+0];
        end
    end
     
    [ccg,r] = CCG(spks_concat,ones(1,length(spks_concat)),'binSize',binSize,'duration',duration,'Fs',fs);
    dopa_table_r(i,:) = r;
    dopa_table_ccg(i,:) = (ccg);
    dopa_table_ccg_zscored(i,:) = zscore(ccg);
    dopa_table_ccg_unity(i,:) = mat2gray(ccg);
    
    clear spk_data spk_max spk_count table_ccg spk_max table_ccg
end
%% Human organoids
list_orgs = [
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well1/')};
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well2/')};
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well3/')};
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well4/')};
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well5/')};
    {strcat(data_path,'human_cerebral_organoids_30min/M03912/211231/well6/')};
    ];
% parameters for CCG
binSize = 1/1000; % 1 ms bins
duration = 1; % 50ms here corresponds to +/-25ms lag
fs = 1/10000;
for i = 1:6
           
    % concat all spk times
    load([list_orgs{i} 'spk_data.mat'])
    fprintf(['orgs culture: ' num2str(i),'\n'])
    
    spks_concat = [];
    for unit = 1:length(spk_data.spike_times)
        if unit>1
        spks_concat = [spks_concat; (spk_data.spike_times{unit})+max(spks_concat)];
        else
             spks_concat = [spks_concat; (spk_data.spike_times{unit})+0];
        end
    end
     
    [ccg,r] = CCG(spks_concat,ones(1,length(spks_concat)),'binSize',binSize,'duration',duration,'Fs',fs);
    org_table_r(i,:) = r;
    org_table_ccg(i,:) = (ccg);
    org_table_ccg_zscored(i,:) = zscore(ccg);
    org_table_ccg_unity(i,:) = mat2gray(ccg);
    
    clear spk_data spk_max spk_count table_ccg spk_max table_ccg
end
%% concatencate ACs
concat = [
rodent_table_ccg_unity(:,501:end); % rodent 12
gluta_table_ccg_unity(:,501:end); % gluta 7
moto_table_ccg_unity(:,501:end); % moto 7
dopa_table_ccg_unity(:,501:end); % dopa 6
org_table_ccg_unity(:,501:end);]; % org 6
%% visualisation settings
%{
cmap = [1 1 0;     % 1:yellow - primary
        0 1 0;     % 2:green - retina
        0 0.5 1;   % 3:blue
        1 0.5 1;   % 4:pink
        0 1 1];    % 5:cyan
%}
% set colours
% rodent
rodentcol = [1 1 0];
ipsccol = 1-summer(3); ipsccol = ipsccol(1:3,:);
orgcol = [10 170 170]/256;
cmap = [rodentcol;ipsccol;orgcol];
    
mat = corrcoef(concat'); % concatenated autocorrelations
rng default
input = mat;
Y = tsne(input);
indices = [ones(1,12) ones(1,7)*2 ones(1,6)*3 ones(1,7)*4 ones(1,6)*5];
D=squareform(pdist(Y));
A=D<200;
A(1:length(A)+1:end)=0;
%% visualise
h = figure('Units', 'pixels', 'Position', [100 100 550 450],'color','white');
[grid,cc,cmap,xrng,yrng,himage] = plot_smoothed_cluster(Y,A,indices,cmap);
xlabel('tSNE 1'); ylabel('tSNE 2');
xticks([]); yticks([]);
box off;
set(gca, 'TickDir', 'out','TickLength', [.01 .01],'linew',1,'FontSize',18,'FontName','Arial');
%{
u = legend('Rodent cortex',...
    'Human gluta',...
    'Human moto',...
    'Human dopa',...
    'Human organoid','Location','northeastoutside','FontName','Arial');
legend('boxoff')
%}
% Save
%{
cd '/Users/schroeter/Desktop/AC_clustering_figure6/'
saveas(gcf,['Figure_6_Clustering_' date '.jpg']);
saveas(gcf,['Figure_6_Clustering_' date '.fig']);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure; h.Position = [100 100 550 450]; 
%subplot(121);gscatter(Y(:,1),Y(:,2),indices); b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
imagesc(mat);  b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xticks([]); yticks([]); c = colorbar; c.Label.String = 'r'; caxis([-1 1]); colormap(viridis);
ylabel('Culture');

% figure
h = figure; h.Position = [100 100 1600 250];
subplot(151);plot(rodent_table_ccg_unity','linewidth',1,'color',[.75 .75 .75]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xlabel('Lag (ms)'); ylabel('Autocorrelation');
xticklabels({'-500','0','500'});
hold on; stdshade(rodent_table_ccg_unity,.75,rodentcol);

subplot(152);plot(gluta_table_ccg_unity','linewidth',1,'color',[.75 .75 .75]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xlabel('Lag (ms)'); ylabel('Autocorrelation');
xticklabels({'-500','0','500'});
hold on; stdshade(gluta_table_ccg_unity,.75,ipsccol(1,:));

subplot(153),plot(moto_table_ccg_unity','linewidth',1,'color',[.75 .75 .75]); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xlabel('Lag (ms)'); ylabel('Autocorrelation');
xticklabels({'-500','0','500'});
hold on; stdshade(moto_table_ccg_unity,.75,ipsccol(2,:));

subplot(154),plot(dopa_table_ccg_unity','linewidth',1,'color',[.75 .75 .75]); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xlabel('Lag (ms)'); ylabel('Autocorrelation');
xticklabels({'-500','0','500'}); 
hold on; stdshade(dopa_table_ccg_unity,.75,ipsccol(3,:));

subplot(155); plot(org_table_ccg_unity','linewidth',1,'color',[.75 .75 .75]); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; box off;
xlabel('Lag (ms)'); ylabel('Autocorrelation');
xticklabels({'-500','0','500'});
hold on; stdshade(org_table_ccg_unity,.75,orgcol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

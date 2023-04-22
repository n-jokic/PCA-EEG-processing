close all;
clear all;
clc;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
SAVE = 0;

Fs = 256;                          % frekvencija odabiranja EEG signala
T = 1;                            % trajanje signala u sekundama
time = 0:1/Fs:T-1/Fs;     

% ucitavanje EEG signala

EEG = importdata('EEG_clean.mat');
duration = T;
EEG = EEG(:, 1:Fs*duration);   


% formiranje sintetickog suma pomocu funkcije ECGwaveGen
% amplituda EKG signala će biti reda EEG (ako su jedinice uV)

bpm = 110;
fs = 256;

amp = 1000;
scale = max(max(EEG,[],2))/amp ;
s = ECGwaveGen(bpm,duration,fs,amp)*scale;
%%

A = zeros(size(EEG,1), 1);
num_noised_channels = 5;             % broj kanala koji se zeli zasumiti u EEG signalu (EEG signal ima 20 kanala)

% u noised_channels promenljivu se beleze redni brojevi kanala koji su slucajno odabrani
% za zasumljivanje (sledeca while petlja)

noised_channels = zeros(1,num_noised_channels);

k = 1;
while ~all(noised_channels)
    new_noised_channel = ceil(size(EEG,1).*rand);
    if ~any(noised_channels == new_noised_channel)
        noised_channels(k) = new_noised_channel;
        k = k+1;
    end
end


clear_ch = 0; 
for i  = 1 : length(EEG(:,1))
    if(ismember(noised_channels, i))
    
    else
        clear_ch = i;
        break
    end

end

%podesiti da imaju dve amplitude EKG šuma
alpha =  1;                            %faktor skaliranja šuma
beta  =  alpha;                           %faktor skaliranja šuma


upper_half = floor(num_noised_channels/2);
lower_half = floor(num_noised_channels/2)+1;

A(noised_channels(1:upper_half),1) = alpha;                        % amplituda jedne polovine suma
A(noised_channels(lower_half:num_noised_channels),1) = beta;       % amplituda druge polovine suma


% zasumljivanje EEG signala sa jednim i drugim sumom po random kanalima

noise = A*s;
x = noise + EEG;          % zasumljeni EEG signali na slucajnim kanalima

SNR = snr(x,noise);  
%%
f = figure(1);
f.Name = 'Comp_EEG';
    plot(time, x(clear_ch,:), 'k');
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [mV]');
        title('EEG comparison');
        grid on;
    
    hold on;
    plot(time, x(noised_channels(1),:), 'r--');
        legend({['EEG w/o noise, ch ' num2str(clear_ch)], ...
            ['EEG with noise, ch ' num2str(noised_channels(1))]})

if(SAVE)
    saveas(f, ['.\izvestaj\slike\' f.Name '_' num2str(SNR) '.eps'], 'epsc');
end
%% algoritmi

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROBUST ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_robustICA, H_robustICA] = robustica(x, [], 1e-3, 1000, 1, 'r', 0, [], 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOBI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H_SOBI,y_SOBI] = acsorbiro(x,size(x,1),100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_JADE = jadeR(x,size(x,1)); 
y_JADE = H_JADE*x; 
y_JADE = real(y_JADE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_CCA,H_CCA,r] = ccabss_test_cc(x);
y_CCA = real(y_CCA);

%%
% Pronalazenje odogovarajucih komponenata koje su maksimalno korelisane sa
% dobijenim nezavisnim komponentama

[s, y_robustICA] = matching_components(s,y_robustICA);
[s, y_SOBI] = matching_components(s,y_SOBI);
[s, y_JADE] = matching_components(s,y_JADE);
[s, y_CCA] = matching_components(s,y_CCA);

%% grafici

f = figure(2);
f.Name = 'Comp_Alg';
subplot(5,1,1)
    plot(time, s); title('Originalni signal'); xlim([0 10]);
subplot(5,1,2)
    plot(time, y_robustICA,'r'); title('robustICA'); xlim([0 10]);
subplot(5,1,3)
    plot(time, y_SOBI,'m'); title('SOBI'); xlim([0 10]);
subplot(5,1,4)
    plot(time, y_JADE,'g'); title('JADE'); xlim([0 10]);
subplot(5,1,5)
    plot(time, y_CCA,'c'); title('CCA' ); xlim([0 10]);

if(SAVE)
    saveas(f, ['.\izvestaj\slike\' f.Name '_' num2str(SNR) '.eps'], 'epsc');
end
%% poredjenje algoritama
%    
rr = zeros(size(s,1),4);        % Spearman-ova korelacija
rc = zeros(size(s,1),4);        % Pearson-ova korelacija
r_rms = zeros(size(s,1),4); 

% U GNU Octave se koristi spearman(s(j,:)',y_(j,:)'); 
rr(1) = corr(s',y_robustICA','type','Spearman'); 
rr(2) = corr(s',y_SOBI','type','Spearman'); 
rr(3) = corr(s',y_JADE','type','Spearman'); 
rr(4) = corr(s',y_CCA','type','Spearman'); 

rc(1) = corr(s',y_robustICA'); 
rc(2) = corr(s',y_SOBI'); 
rc(3) = corr(s',y_JADE'); 
rc(4) = corr(s',y_CCA'); 

absCorPearson1 = abs(rc);
absCorSpearman1 = abs(rr); 

r_rms(1) = rmse(s,y_robustICA);    
r_rms(2) = rmse(s,y_SOBI);
r_rms(3) = rmse(s,y_JADE);
r_rms(4) = rmse(s,y_CCA);



Method = {'robustICA';'SOBI';'JADE';'CCA'};             %
                                     %
rms = [r_rms(1);r_rms(2);r_rms(3);r_rms(4)];                                    %
corr = [rr(1);rr(2);rr(3);rr(4)];                                       %                               %
T = table(rms,corr);                              %
T.Properties.RowNames = Method;  

if(SAVE)
    table2latex(T, ['.\izvestaj\tabele\rms_corr' '_' num2str(SNR)]);  
end
%%
time = {};
total = 40;
t = zeros(total, 1);
for i = 1 : total
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% ROBUST ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [y_robustICA, H_robustICA] = robustica(x, [], 1e-3, 1000, 1, 'r', 0, [], 0);
    t(i) = toc;
end
time{1} = t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOBI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : total
    tic;
    [H_SOBI,y_SOBI] = acsorbiro(x,size(x,1),100);
    t(i) = toc;
end
time{2} = t;


for i = 1 : total
    tic;
    H_JADE = jadeR(x,size(x,1)); 
    y_JADE = H_JADE*x; 
    y_JADE = real(y_JADE);
    t(i) = toc;
end
time{3} = t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : total
    tic;
    [y_CCA,H_CCA,r] = ccabss_test_cc(x);
    y_CCA = real(y_CCA);
    t(i) = toc;
end
time{4} = t;


%% 
estimation = {};
step = 0.00001;
tit = {['ROBUST ICA'], ['SOBI'], ['JADE'], ['CCA']};
for i = 1 : length(time)
    f = figure(10);
        subplot(4, 1, i);
        t = 0.95*min(time{i}):step:max(time{i});
        pdSix = fitdist(time{i},'Kernel');
        ySix = pdf(pdSix,t);
        plot(t, ySix, 'k');
            grid('on');
            xlabel('$t$ [s]');
            ylabel('$\hat{f}_{pdf}(t)$');
            xlim([min(t), max(t)]);
            title([tit{i} ' pdf estimation']);
end

if(SAVE)
    f = figure(10);
    f.Name = 'psd';
    saveas(f, ['.\izvestaj\slike\' f.Name '.eps'], 'epsc');
end

Method = {'robustICA';'SOBI';'JADE';'CCA'};             %
                                     %
mean = [mean(time{1}); mean(time{2});mean(time{3});mean(time{4})];                                    %
std = [std(time{1}); std(time{2});std(time{3});std(time{4})];                                      %                               %
T = table(mean,std);                              %
T.Properties.RowNames = Method;  

if(SAVE)
    table2latex(T, ['.\izvestaj\tabele\time' '_' num2str(SNR)]);  
end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% File:	lidiane_movimentoA.m 
% Author:	Tenysson Will de Lemos (c) 2014, All Rights Reserved
%			tenysson@fmrp.usp.br
%         	NAP - DCD 
%         	Av. Bandeirantes, 3900 
%         	Dept. Biomecanica, Medicina e Reabilitação do Aparelho Locomotor
%		  	Faculdade Medicina de Ribeirao Preto - Universidade de Sao Paulo
%         	Ribeirao Preto, SP, Cep: 14049-900
% Date:   07-06-2018
% Version: 1.0
% Objetivo: 
%         Processar o experimento de resistencia com o MCU
% Creditos:

function lidiane_movimentoA
clc
close all

% Carrega o arquivo *.hpf
[filename,pathname] = uigetfile('*.hpf','Selecione o arquivo desejado');
h = actxserver('EW4COMPlus.HPFAccess');
h.invoke('OpenHPFFile',[pathname filename]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Carregar dados EMG
emg_min = h.invoke('GetFirstAvailibleTime');
emg_max = h.invoke('GetChannelLastAvailibleTime',0);
emg_fs =  h.invoke('GetPerChannelSampleRate',0);

emg_time = linspace(emg_min/emg_fs,emg_max/emg_fs,emg_max);
[b,a] = butter(4,[20 500]/(emg_fs/2));

emg_channels = h.invoke('GetChannelCount');
emg_filtered = zeros(emg_channels,emg_max-emg_min);
emg_names = cell(emg_channels,1);

% Filtragem do sinal EMG
for i = 1:emg_channels
    emg_filtered(i,:) = filtfilt(b,a,double(h.invoke('ReadData',i-1,emg_min,emg_max)));
    emg_names{i} = h.invoke('GetChannelName',i-1);
end
h.invoke('CloseHPFFile');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carrega o arquivo *.mat
[filename,pathname] = uigetfile('*.mat','Selecione o arquivo desejado');
mcu = load([pathname filename]);
Fs_mcu = 10;

% Ponto de sincronizacao
sinc_pt = find(diff(sign(mcu.emg_raw-2.5))==2);

% Corte
mcu.mcu_raw(1:sinc_pt) = [];
%mcu.emgdata(1:sinc_pt) = [];

% Filtragem do Mcu Orig
[b,a] = butter(4,0.5/(Fs_mcu/2));  % 10 Hz de suavizacao
mcu.mcudata_filt = filtfilt(b,a,mcu.mcu_raw);
figure('Name','Mcu','visible','on','units','normalized','outerposition',[0 0 1 1])
plot(mcu.mcudata_filt)
hold on
plot(diff(mcu.mcudata_filt),'r')
plot(diff(sign(diff(mcu.mcudata_filt))),'k')
[cut_pts,~] = ginput(2);
cut_pts = round(cut_pts);
print('-djpeg','-r300',[pathname filename(1:end-4) '_MCU'])
close all

%%%% Analise pelo MCU
% Pontos de Inflexao
vel_pts = find(abs(diff(sign(diff(mcu.mcudata_filt))))==2);

p = round((emg_fs/Fs_mcu)*(vel_pts(vel_pts > cut_pts(1) & vel_pts < cut_pts(end))+2));
% Corte dos intervalo util (de maior numero de pontos)
f = figure('Name',filename(1:end-4),'visible','on','units','normalized','outerposition',[0 0 1 1]);
for i=1:emg_channels
    subplot(emg_channels/2,2,i)
    plot(emg_time,emg_filtered(i,:),'k')
    title(emg_names{i})
    y = get(gca,'YLim');
    line([p(1) p(1)]/emg_fs,y,'Color','g','LineWidth',2)
    line([p(2) p(2)]/emg_fs,y,'Color','c','LineWidth',2)
    line([p(3) p(3)]/emg_fs,y,'Color','r','LineWidth',2)
end

emg_pre = cell(emg_channels,1);
emg_ida = cell(emg_channels,1);
emg_vol = cell(emg_channels,1);

% Cortando o sinal em setores
for i=1:emg_channels
    emg_pre{i} = emg_filtered(i,p(1)-(emg_fs*50e-3):p(1));
    emg_ida{i} = emg_filtered(i,p(1):p(2));
    emg_vol{i} = emg_filtered(i,p(2):p(3));
end

% Calculo da Envoltoria
mcu_env_emg_pre = cell(emg_channels,1);
mcu_env_emg_ida = cell(emg_channels,1);
mcu_env_emg_vol = cell(emg_channels,1);

% Tempo da Envoltoria
% env_t_emg_pre = 50e-3;
% env_t_emg_ida = length(p(1):p(2))/emg_fs;
% env_t_emg_vol = length(p(2):p(3))/emg_fs;

% Rms
mcu_rms_emg_pre = cell(emg_channels,1);
mcu_rms_emg_ida = cell(emg_channels,1);
mcu_rms_emg_vol = cell(emg_channels,1);

[b,a] = butter(4,10/(emg_fs/2));
dt = 1/emg_fs; 
for i=1:emg_channels
   mcu_env_emg_pre{i} = dt*trapz(filtfilt(b,a,abs(emg_pre{i})));
   mcu_env_emg_ida{i} = dt*trapz(filtfilt(b,a,abs(emg_ida{i})));
   mcu_env_emg_vol{i} = dt*trapz(filtfilt(b,a,abs(emg_vol{i})));
   mcu_rms_emg_pre{i} = rms(emg_pre{i});
   mcu_rms_emg_ida{i} = fadiga(emg_ida{i},emg_fs);
   mcu_rms_emg_vol{i} = fadiga(emg_vol{i},emg_fs);
end

% Analise On-Off 
% Escolha 1 do Paul Hogdes 25ms/2 SD/ 50Hz
resp = questdlg({'Escolha dos parametros segundo o artigo do Paul Hodges:','1 - Janela de 25 ms e 3 Desvios Padroes',...
    '2 - Janela de 50 ms e 1 Desvios Padrao','3 - Janela de 25 ms e 2 Desvios Padroes',...
    '','Todas opções usam filtro de 50 Hz e baseline de 50 ms',...
    'Os baseline estao no inicio e final da tarefa (entre 0 e 50 ms para on e entre 9.950 e 10s'},'','1','2','3','1');

switch(resp)
    case '1'
        % Escolha 1 do Paul Hogdes 50ms/1 SD/ 50Hz
        wd = 25;  % Janela (ms)
        sd = 3;   % Desvios Padroes
    case '2'
        % Escolha 2 do Paul Hogdes 50ms/1 SD/ 50Hz
        wd = 50;  % Janela (ms)
        sd = 1;   % Desvios Padroes
    case '3'
        % Escolha 3 do Paul Hogdes 50ms/1 SD/ 50Hz
        wd = 25;  % Janela (ms)
        sd = 2;   % Desvios Padroes
end

baseline = 50; % primeiros e ultimos 50 ms
fc = 50;  % Frequencia de Corte

pre = [zeros(1,p(1)-(emg_fs*50e-3)) ones(1,emg_fs*50e-3) zeros(1,length(emg_filtered(1,:))-p(1))];
ida = [zeros(1,p(1)) ones(1,p(2)-p(1)) zeros(1,length(emg_filtered(1,:))-p(2))];
vol = [zeros(1,p(2)) ones(1,p(3)-p(2)) zeros(1,length(emg_filtered(1,:))-p(3))];

flag_ativacao = cell(emg_channels,3);

for i=1:emg_channels
   disp(['Canal ' num2str(i)])
   onoff = emgonoff2(emg_filtered(i,:)',emg_fs,wd,sd,baseline,fc);
   if(isempty(onoff))
       disp('Sem alteracao em relacao ao baseline')
   else
       figure(f)
       subplot(emg_channels/2,2,i)
       y = get(gca,'YLim');
       line([onoff(1) onoff(1)]/emg_fs,y,'Color','b','LineWidth',2)
       line([onoff(2) onoff(2)]/emg_fs,y,'Color','b','LineWidth',2)

       onoff = [zeros(1,onoff(1)-1) ones(1,onoff(2)-onoff(1)+1) zeros(1,length(emg_filtered(i,:))-onoff(2))];

       if(any(onoff&pre))
           disp('ativado no Pre')
           flag_ativacao{i,1} = 'Ativado Pre';
       end
       if(any(onoff&ida))
           disp('ativado na ida')
           flag_ativacao{i,2} = 'Ativado Ida';
       end

       if(any(onoff&vol))
           disp('ativado na volta')
           flag_ativacao{i,3} = 'Ativado Volta';
       end
   end
end

%legend('sinal','pre','inicio','fim','onoff')
legend('sinal','inicio','meio','fim','onoff')
print('-djpeg','-r300',[pathname filename(1:end-4) '_EMG'])

save([pathname filename(1:end-4) '_processados'],'mcu_env_emg_pre','mcu_env_emg_ida',...
    'mcu_env_emg_vol','mcu_rms_emg_pre','mcu_rms_emg_ida','mcu_rms_emg_vol',...
    'emg_names','flag_ativacao');

end



function saida = fadiga(dados,Fs)
% Janela de Tempo
window_time = 0.375;
overlap_time = window_time/2;

% Janela em pontos
window_length = floor(window_time*Fs);
overlap_length = floor(overlap_time*Fs);

Srms = zeros(1,length(1:overlap_length:length(dados)-window_length));

% Calculo da mediana
for i = 1:overlap_length:length(dados)-window_length
	Srms((i-1)/overlap_length+1) = rms(dados(i:i+window_length));
end


% regstats(fmed,fmed_time,'linear','rsquare')
saida = mean(Srms);

end


% function onoff = emgonoff(rawemg, fs, ws, sd,bl)
% % EMGONOFF - Find on/off times and indicies of raw EMG data.
% % Calculate average(mean) value of resting EMG
% % Define "on" EMG as the sample where average value of EMG in a given window
% % range around the sample is a given # of std. dev. above avg. resting EMG.
% %
% % onoff = emgonoff(rawemg, fs, ws, sd)
% %
% % Use the mouse to select two ranges of "resting" EMG from a graph of the
% % full-wave rectified EMG data.  Click four times: start and end of 1st
% % resting range, and start and end of 2nd resting range.  Mouse clicks need
% % to be consecutive and in order of increasing time (i.e. left-to-right on
% % the graph).  The first range should precede the EMG burst associated with
% % the muscle contraction under consideration.  The sedond range should be 
% % the resting EMG data immediately following the EMG burst.
% %
% % rawemg = input file raw emg data (1-column vector)
% % fs = sampling rate of raw EMG data in Hz
% % ws = window size in milliseconds (50ms @ 2400Hz = 120 samples)
% % sd = number of std. deviations above resting rms emg to trigger an "ON"
% % Default values: 
% %   ws = 50ms
% %   sd = 1
% 
% % Algorythm for EMG onset & offset taken from 
% % Hodges, P.W. and B.H. Bui, _A comparison of computer-based methods for 
% % the determination of onset of muscle contraction using electromyography._ 
% % Electroencephalography & Clinical Neurophysiology, 1996. 101(6): p. 511-9
% %
% % Created by: Kieran A. Coghlan, BSME, MSES
% % SUNY at Buffalo, New York
% % <kc_news@sonic.net>
% % Last modified: 1 May, 2006
% 
% %% Check inputs for defaults
% if nargin < 2, error('Not enough inputs. Type "help emgonoff" for help.'); end
% if nargin < 3, ws = 50; end
% if nargin < 4, sd = 1; end;
% %% Full-Wave-Rectify the raw data
% fwlo = abs(rawemg(:,1));
% %% prepare for loop
% % Get two ranges for resting emg (before & after burst) using ginput
% % R = input('\nUse the mouse to select FOUR points to define the begining and \nend of two data ranges that will be used to calculate average\nresting EMG values before and after the EMG burst (muscle contraction)\nPress [RETURN] to begin: ');
% % clear R;
% % f1 = figure;
% % plot(fwlo);
% % title(musculo)
% % close
% % [x,~] = ginput(2); %click four times: two for start/end of resting emg before burst two for resting emg after burst
% % if(isempty(x))
% %     onoff = [];
% %     return;
% % end
% x = round([1 1+bl*(fs/1000) 10*fs-bl*(fs/1000) 10*fs]);
% 
% %close(f1);clear f1; % Leave commented if you want to keep EMG graph up to 
% %                    % do a visual QA of on/off results.
% %% preallocate arrays
% mvgav = zeros(x(4)-x(1),1);
% onoff(1,1) = 0;
% i=0;
% restav = mean(fwlo(x(1):x(2))); %average value of rest EMG before ON
% reststd = std(fwlo(x(1):x(2))); %std. dev. of rest EMG before ON
% restav2 = mean(fwlo(x(3):x(4))); %average value of rest EMG after OFF
% reststd2 = std(fwlo(x(3):x(4))); %std. dev. of rest EMG after OFF
% %% window size (in samples) = ws*fs e.g. 50ms*2400Hz = 120 samples
% sws2 = fs*(0.001*ws);
% sws = 0.5*(sws2);
% sws = round(sws);
% %% find "ON" index:
% % for xi, change from x(1) to x(2) if you want to ignore any "blips"
% % within the resting range.
% %xi = x(1);
% xi = x(2);
% xi = round(xi);
% for n = 2:length(mvgav);
%     mvgav(n,1) = mean(fwlo((xi-sws):(xi+sws)));
%     if mvgav(n) > restav+sd*reststd;
%         i = i+1;
%         onoff(i,1) = xi;
%         break
%     end
%     xi = xi+1;
% end
% 
% %% find "OFF" index:
% clear n xi i
% mvgav2=zeros(x(4)-x(1),1);
% i=0;
% xi=onoff(1,1)+(1/2)*(x(3)-onoff(1,1)); %start OFF search approx. 1/2 way through ON burst.
% %% OFF loop:
% xi=round(xi);
% for n=2:length(mvgav2);
%     mvgav2(n,1)=mean(fwlo((xi-sws):(xi+sws)));
%     if mvgav2(n)<restav2+sd*reststd2;
%         i=i+1;
%         onoff(i,2)=xi;
%         break
%     end
%     xi=xi+1;
% end
% 
% if(onoff(1)==0)
%    onoff = [];    
% end
%   
% 
% end


function onoff = emgonoff2(rawemg, fs, ws, sd,bl,fc)
% EMGONOFF - Find on/off times and indicies of raw EMG data.
% Calculate average(mean) value of resting EMG
% Define "on" EMG as the sample where average value of EMG in a given window
% range around the sample is a given # of std. dev. above avg. resting EMG.
%
% onoff = emgonoff(rawemg, fs, ws, sd)
%
% Use the mouse to select two ranges of "resting" EMG from a graph of the
% full-wave rectified EMG data.  Click four times: start and end of 1st
% resting range, and start and end of 2nd resting range.  Mouse clicks need
% to be consecutive and in order of increasing time (i.e. left-to-right on
% the graph).  The first range should precede the EMG burst associated with
% the muscle contraction under consideration.  The sedond range should be 
% the resting EMG data immediately following the EMG burst.
%
% rawemg = input file raw emg data (1-column vector)
% fs = sampling rate of raw EMG data in Hz
% ws = window size in milliseconds (50ms @ 2400Hz = 120 samples)
% sd = number of std. deviations above resting rms emg to trigger an "ON"
% Default values: 
%   ws = 50ms
%   sd = 1

% Algorythm for EMG onset & offset taken from 
% Hodges, P.W. and B.H. Bui, A comparison of computer-based methods for 
% the determination of onset of muscle contraction using electromyography._ 
% Electroencephalography & Clinical Neurophysiology, 1996. 101(6): p. 511-9
%
% Created by: Kieran A. Coghlan, BSME, MSES
% SUNY at Buffalo, New York
% <kc_news@sonic.net>
% Last modified: 1 May, 2006

%% Full-Wave-Rectify the and filtering the raw data
[b,a] = butter(4,fc/(fs/2));
fwlo = filtfilt(b,a,abs(rawemg(:,1)));

% The base line points
bslnON_srt  = 1;
bslnOFF_end = 10*fs;
x = round([bslnON_srt bslnON_srt+bl*(fs/1000)-1 bslnOFF_end-bl*(fs/1000) bslnOFF_end]);

%close(f1);clear f1; % Leave commented if you want to keep EMG graph up to 
%                    % do a visual QA of on/off results.
%% preallocate arrays
mvgav = zeros(x(4)-x(1),1);
onoff(1,1) = 0;
i=0;
restav = mean(fwlo(x(1):x(2))); %average value of rest EMG before ON
reststd = std(fwlo(x(1):x(2))); %std. dev. of rest EMG before ON
restav2 = mean(fwlo(x(3):x(4))); %average value of rest EMG after OFF
reststd2 = std(fwlo(x(3):x(4))); %std. dev. of rest EMG after OFF
%% window size (in samples) = ws*fs e.g. 50ms*2400Hz = 120 samples
sws2 = fs*(0.001*ws);
sws = 0.5*(sws2);
sws = round(sws);
%% find "ON" index:
% search after x(2)
xi = round(x(2));
for n = 2:length(mvgav)
    mvgav(n,1) = mean(fwlo((xi-sws):(xi+sws)));
    if (mvgav(n) > restav+sd*reststd)
        i = i+1;
        onoff(i,1) = xi;
        break
    end
    xi = xi+1;
end

%% find "OFF" index:
clear n xi i
mvgav2 = zeros(x(4)-x(1),1);
i = 0;
xi=onoff(1,1)+(1/2)*(x(3)-onoff(1,1)); %start OFF search approx. 1/2 way through ON burst.
%% OFF loop:
xi=round(xi);
for n=2:length(mvgav2)
    mvgav2(n,1)=mean(fwlo((xi-sws):(xi+sws)));
    if (mvgav2(n)<restav2+sd*reststd2)
        i=i+1;
        onoff(i,2)=xi;
        break
    end
    xi=xi+1;
end

if(onoff(1)==0)
   onoff = [];    
end
  

end
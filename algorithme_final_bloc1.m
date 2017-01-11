
% var a importer du menu%
%cas_de_figure=6;
%wav_name='chasseurs';
%  debut algo %

% ATTENTION : ALGORITHME TRES LENT
% CONSEIL: LANCER le premier bloc (ctrl+entree) pour execution des
% localisations
% Pui lancer second bloc (ctrl+entree) apres les '%%' pour ecoute et
% affichage simultanee du resultat

% ATTENTION: avoir le fichier circle.m dans le mm repertoire pour le dessin


%   test spatialisation d'un signal qui fait un 360 autour de l'auditeur.

% un signal sonore de 25 sec
% s'ecoute en stereo.
% je voulais:
% Resultat: fonctionne comme attendu !

% INIT SOURCE -> c'est le fichier son non-spatialise
% fs=44100;
% signal=randn(fs,1);

if parole
    if strcmp(wav_name(end-4:end),'.wav')
        file=wav_name;
        wav_name(end-4:end)=[];
    else
        file=[wav_name, '.wav'];
    end
    
    
    [signal_d,fs]=audioread(file);
    if size(signal_d,2)==2
        signal=signal_d(:,1)+signal_d(:,2);
    else
        signal=signal_d;
    end
    %elimination des premiers echantillons nuls (le son ne commence pas
    %encore); si on le fait pas => pb dans le critere qui vaudra Inf
    for i=1:length(signal)
        if signal(i,1)>0
            debut=i; %enregistrement de l'indice ou le son debute
            break;
        end
    end
    clear signal;
    if size(signal_d,2)==2
        signal=signal_d(debut:end,1)+signal_d(debut:end,2);
    else signal=signal_d(debut:end,1);
    end
    
elseif gauss
    fs=44100;
    signal=randn(fs*8,1);
    [b,a]=butter(3,5000/(fs/2),'low');
    signal=filter(b,a,signal);
end
duree_son = length(signal) / fs;% => 25 sec

% o=audioplayer(signal,fs);
% pour jouer le fichier    >>  play(o)



sortie_space=spatialisation(gauss,signal,wav_name);


outputSignal=sortie_space.outputSignal;
D=sortie_space.D;
%Taille_groupe=sortie_space.Taille_groupe;
B=sortie_space.B;
Ntheta=sortie_space.Ntheta;
P=sortie_space.P;
sigma=sortie_space.sigma;
thetaArg=sortie_space.thetaArg;
azimuth=sortie_space.azimuth;



sortie_loca=localisation(gauss,signal,wav_name,cas_de_figure);













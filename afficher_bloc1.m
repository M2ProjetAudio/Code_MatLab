
%       joue le son sur stereo

p=audioplayer(outputSignal,fs);
play(p);

%         AFFICHAGE


for exp=1:length(D)
    
    
    %source
    plot(D(exp)*cos(azimuth(exp)*pi/180+pi/2),...
        D(exp)*sin(azimuth(exp)*pi/180+pi/2),'sr')
    hold on
    grid on
    text(D(exp)*cos(azimuth(exp)*pi/180+pi/2),...
        D(exp)*sin(azimuth(exp)*pi/180+pi/2),'  Source sonore')
    
    
    %Jnorme=(1./J(exp,:))/(1/maxcrit)*D(exp)*.9; %norme pour la figure
    Jnorme=Jpos(exp,:)/maxcrit_exp(exp)*D(exp)*.9;
    
    fill(cos(thetaArg+pi/2).*Jnorme,sin(thetaArg+pi/2).*Jnorme,[.8 .8 .8])
    %plot(0,0,'bo')
    circle([0,0],.5,'b','LineWidth',2)
    circle([-.5,0],.1,'b','LineWidth',2)
    circle([.5,0],.1,'b','LineWidth',2)
    text(0,-.5,'  Tete')
    line([0,0],[0,1],'LineWidth',2)
    quiver(0,0,D(exp)*cos(theta_mle(exp)+pi/2),D(exp)*sin(theta_mle(exp)+pi/2),'r')
    %quiver(0,0,0,1,'b')
    title(sprintf...
        ('Resultat estimation de la Localisation.\n source en  %.0f degres et distance=%.1f m',azimuth(exp),D(exp)))
    xlabel('x')
    ylabel('y')
    axis([-D(exp)-.5 D(exp)+.5 -D(exp)-.5 D(exp)+.5])
    fig=gcf;
    fig.Color='w';
    hold off
    pause(.5)
end

if parole
    audiowrite([wav_name,'360.wav'],outputSignal,fs)
elseif gauss
    audiowrite('gauss360.wav',outputSignal,fs)
end
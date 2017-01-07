
%       joue le son sur stereo

p=audioplayer(outputSignal,fs);

%%
%nbre d'experiences
nb_exp=length(D);
hold off
%carre de la source:
carre_source=line('Visible','off');
hold on
%texte de la source
texte_source=text(0.1,0.314,'source sonore','Visible','off');
%bulles des valeurs des criteres
criteres_J_fill=fill(0,0,[.8 .8 .8],'Visible','off');
%fleche du mle
fleche=quiver(0,0,1,1,'r','Visible','off');

Jnorme{nb_exp}=zeros(size(Jpos));  xsource(nb_exp)=0; ysource(nb_exp)=0;
u(nb_exp)=0; xtext(nb_exp)=0; ytext(nb_exp)=0;
v(nb_exp)=0; xfill{nb_exp}=zeros(1,360); yfill{nb_exp}=zeros(1,360);
for num_exp=1:nb_exp
   Jnorme{num_exp}=Jpos(num_exp,:)/maxcrit_exp(num_exp)*D(num_exp)*.9;
   xfill{num_exp}=cos(thetaArg+pi/2).*Jnorme{num_exp};
   yfill{num_exp}=sin(thetaArg+pi/2).*Jnorme{num_exp};
   u(num_exp)=D(num_exp)*cos(theta_mle(num_exp)+pi/2);
   v(num_exp)=D(num_exp)*sin(theta_mle(num_exp)+pi/2);
   xsource(num_exp)=D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2);
   ysource(num_exp)=D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2);
   xtext(num_exp)=D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2);
   ytext(num_exp)=D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2);
end

circle([0,0],.5,'b','LineWidth',2)
circle([-.5,0],.1,'b','LineWidth',2)
circle([.5,0],.1,'b','LineWidth',2)
text(0,-.5,'  Tete')
line([0,0],[0,1],'LineWidth',2)
title(sprintf...
    ('Resultat estimation de la Localisation.\n source en  %.0f degres et distance=%.1f m',azimuth(num_exp),D(num_exp)))
xlabel('x')
ylabel('y')
axis([-D(1)-.5 D(1)+.5 -D(1)-.5 D(1)+.5])
fig=gcf;
fig.Color='w';

%%
dt=p.TotalSamples/fs/nb_exp;
play(p);

%         AFFICHAGE


for num_exp=1:nb_exp
   
    carre_source.XData=xsource(num_exp);
    carre_source.YData=ysource(num_exp);
    texte_source.Position=[xtext(num_exp) ytext(num_exp) 0];
    fleche.UData=u(num_exp);
    fleche.VData=v(num_exp);
    criteres_J_fill.XData=xfill{num_exp};
    criteres_J_fill.YData=yfill{num_exp};
    
    if num_exp==1
        carre_source=line('Visible','on');
        %texte de la source
        texte_source=text(0,0,'source sonore','Visible','on');
        %bulles des valeurs des criteres
        criteres_J_fill=fill(0,0,[.8 .8 .8],'Visible','on');
        %fleche du mle
        fleche=quiver(0,0,1,1,'r','Visible','on');
    end
    
    pause(dt)
end





%%
% for num_exp=1:nb_exp
%     
%     
%     %source
%     plot(D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2),...
%         D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2),'sr')
%     hold on
%     grid on
%     text(D(num_exp)*cos(azimuth(num_exp)*pi/180+pi/2),...
%         D(num_exp)*sin(azimuth(num_exp)*pi/180+pi/2),'  Source sonore')
%     
%     
%     %Jnorme=(1./J(num_exp,:))/(1/maxcrit)*D(num_exp)*.9; %norme pour la figure
%     Jnorme=Jpos(num_exp,:)/maxcrit_exp(num_exp)*D(num_exp)*.9;
%     
%     fill(cos(thetaArg+pi/2).*Jnorme,sin(thetaArg+pi/2).*Jnorme,[.8 .8 .8])
%     %plot(0,0,'bo')
%     circle([0,0],.5,'b','LineWidth',2)
%     circle([-.5,0],.1,'b','LineWidth',2)
%     circle([.5,0],.1,'b','LineWidth',2)
%     text(0,-.5,'  Tete')
%     line([0,0],[0,1],'LineWidth',2)
%     quiver(0,0,D(num_exp)*cos(theta_mle(num_exp)+pi/2),D(num_exp)*sin(theta_mle(num_exp)+pi/2),'r')
%     %quiver(0,0,0,1,'b')
%     title(sprintf...
%         ('Resultat estimation de la Localisation.\n source en  %.0f degres et distance=%.1f m',azimuth(num_exp),D(num_exp)))
%     xlabel('x')
%     ylabel('y')
%     axis([-D(num_exp)-.5 D(num_exp)+.5 -D(num_exp)-.5 D(num_exp)+.5])
%     fig=gcf;
%     fig.Color='w';
%     hold off
%     pause(.5)
% end

if parole
    audiowrite([wav_name,'360.wav'],outputSignal,fs)
elseif gauss
    audiowrite('gauss360.wav',outputSignal,fs)
end
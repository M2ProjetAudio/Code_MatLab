function src = synchroCallback(src, eventdata)
myStruct = get(src, 'UserData'); %//Unwrap
num_exp=myStruct.num_exp;

myStruct.carre_source.XData=myStruct.xsource(num_exp);
myStruct.carre_source.YData=myStruct.ysource(num_exp);
myStruct.texte_source.Position=[myStruct.xtext(num_exp) myStruct.ytext(num_exp) 0];
myStruct.fleche.UData=myStruct.u(num_exp);
myStruct.fleche.VData=myStruct.v(num_exp);
myStruct.criteres_J_fill.XData=myStruct.xfill{num_exp};
myStruct.criteres_J_fill.YData=myStruct.yfill{num_exp};
myStruct.barre_rouge.XData=[myStruct.xbarre_rouge(num_exp),myStruct.xbarre_rouge(num_exp)];




if myStruct.prendre_video
    myStruct.frame{num_exp}=getframe;
end

myStruct.num_exp=myStruct.num_exp+1;

set(src, 'UserData', myStruct); %//Rewrap

end
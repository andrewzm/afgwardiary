function y=  SetFontSize(fs)

set( get(gca,'YLabel'),'FontSize',fs); 
set( get(gca,'XLabel'),'FontSize',fs); 
set( get(gca,'Title'),'FontSize',fs); 
set(gca,'fontsize',fs)
plot_ics<-function(data,xmin,xmax,title){p1<-ggplot(data[which(data$term=="st.prs"),],aes(y = fct_rev(var))) + 
  geom_point(aes(x=OR,y = fct_rev(var)), shape=15, size=3) +geom_text(aes(x=OR,label=round(OR,2)),hjust=0.5,vjust=-0.7)+
  geom_linerange(aes(xmin=CI_l, xmax=CI_u))+geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="OR (95% IC) por cada SE", y="")+ theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())+xlim(xmin,xmax)

p2<-  ggplot(data[which(data$term=="st.prs"),],aes(y = fct_rev(var)))+ylab("")+ggtitle(title)+
  geom_text(aes(x = 0, label = var), hjust = 0, fontface = "bold",family="sans")+theme_void()+
  theme(plot.title=element_text(hjust=1,family="sans",color="Black",vjust=4))

p3<-ggplot(data[which(data$term=="st.prs"),],aes(y = fct_rev(var)))+geom_text(aes(x = 0, y = fct_rev(var), label = p),
                                                                             hjust = 0,family="sans")+theme_void()

layout <- c(
  area(t = 0, l = 0, b = 32, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 4, b = 32, r = 10),
  area(t = 0, l = 10.5, b = 32, r = 12) # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30    area(t = 0, l = 9, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
plot_eur_basic<-p2+p1+p3+plot_layout(design=layout)+plot_annotation(theme=theme(
  panel.background = element_rect(fill='transparent',colour=NA),
  plot.background = element_rect(fill='transparent', colour=NA))
)
p11<-ggplot(data[which(data$term=="cincovsrest"),],aes(y = fct_rev(var))) + geom_text(aes(x=OR,label=round(OR,2)),hjust=0.5,vjust=-0.7)+
  geom_point(aes(x=OR,y=fct_rev(var)), shape=15, size=3) +
  geom_linerange(aes(xmin=CI_l, xmax=CI_u))+geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="OR (95% IC) Q5 vs resto", y="")+ theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())+xlim(xmin,xmax)
p22<-  ggplot(data[which(data$term=="cincovsrest"),],aes(y = fct_rev(var)))+ylab("")+
  geom_text(aes(x = 0, label = var), hjust = 0, fontface = "bold",family="sans")+theme_void()
p33<-ggplot(data[which(data$term=="cincovsrest"),],aes(y = fct_rev(var)))+geom_text(aes(x = 0, y = fct_rev(var), label = p),   hjust = 0,family="sans")+theme_void()

layout <- c(
  area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 4, b = 30, r = 10),
  area(t = 0, l = 10.5, b = 30, r = 12) )# middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30    area(t = 0, l = 9, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)

plot_eur_terc<-p22+p11+p33+plot_layout(design=layout)+plot_annotation(theme=theme(
  panel.background = element_rect(fill='transparent',colour=NA),
  plot.background = element_rect(fill='transparent', colour=NA))
)

plot_eur_ors<-plot_eur_basic/plot_eur_terc
return(plot_eur_ors)}

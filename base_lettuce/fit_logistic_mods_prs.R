fit_mod<-function(data,covar,outcomes){
mod_prs<-function(p){f1<-as.formula(sprintf('%s ~ %s + st.prs',as.character(p),paste(covar,collapse="+")))
glm(f1,data,family=binomial(link="logit"))}
mod_quintil<-function(p){f2<-as.formula(sprintf('%s ~ %s + cincovsrest',as.character(p),paste(covar,collapse="+")))
glm(f2,data,family=binomial(link="logit"))}
    
p<-outcomes
models_prs<-map(p,~mod_prs(.x))
models_quintil<-map(p,~mod_quintil(.x))
    
results_models<-bind_rows(
map_dfr(models_prs,broom::tidy,.id="index"),
map_dfr(models_quintil,broom::tidy,.id="index"))%>%
filter(term%in%c("st.prs","cincovsrest"))%>%
mutate(pheno=outcomes[as.numeric(index)],OR=exp(estimate),CI_l=exp(estimate-1.96*std.error),CI_u=exp(estimate+1.96*std.error))%>%
select(-index)%>%relocate(pheno,.before=term)
return(results_models)}

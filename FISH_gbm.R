source('FISH_functions.R')



## ---- formulas
formulas = list(
    all = ~REGION + NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
    all1 = ~NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
    Palm = ~ NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
    Magnetic = ~ NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
    Whitsunday = ~ NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
    Keppel = ~ NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP
)

## ----end

## ---- analysis.list
load('data/fish.RData')
load('data/var.lookup.RData')
resp.lookup = var.lookup %>% filter(Type=='Response') %>% droplevels
pred.lookup = var.lookup %>% filter(Type=='Predictor') %>% droplevels

groupings.all = vector('list', nrow(pred.lookup))
groupings.all1 = vector('list', nrow(pred.lookup))
names(groupings.all) = pred.lookup$Abbreviation
groupings.region = vector('list', nrow(pred.lookup))
names(groupings.region) = pred.lookup$Abbreviation
for (i in 1:nrow(pred.lookup)) {
    pred = as.character(pred.lookup[i,'Abbreviation'])
    groupings.all[[pred]] = ifelse(pred=='REGION',NA,'REGION')
    groupings.all1[[pred]] = NA
    groupings.region[[pred]] = ifelse(pred=='NTR.Pooled',NA,'NTR.Pooled')
}
groupings.all = do.call('c',groupings.all)
groupings.all1 = do.call('c',groupings.all1)
groupings.region = do.call('c',groupings.region)

gfun = function(f, form) {
    #print(f)
    #print(form[[f]])
    form = attr(terms(form[[f]]),'term.labels')
    if (f=='all') return(groupings.all[form])
    if (f=='all1') return(groupings.all1[form])
    else return(groupings.region[form])
}

analyses = vector('list',nrow(resp.lookup))
names(analyses) = resp.lookup$Abbreviation
for (i in 1:nrow(resp.lookup)) {
    resp = as.character(resp.lookup[i,'Abbreviation'])
    fun = as.character(resp.lookup[i,'Transform'])
    fam = as.character(resp.lookup[i,'Family'])
    
    analyses[[resp]] = list('formulas' = lapply(formulas, function(f) update(f, paste0(fun,'(',resp,') ~.'))),
                            'family' = fam)
    if (resp %in% c('PI','PTD','PLD')) analyses[[resp]][['formulas']] = lapply(analyses[[resp]][['formulas']], function(f) f=update(f, .~.+PREY.DENSITY))
    if (resp %in% c('PTB','PLB')) analyses[[resp]][['formulas']] = lapply(analyses[[resp]][['formulas']], function(f) f=update(f, .~.+PREY.BIOMASS))
    analyses[[resp]][['groups']] = sapply(c('all','all1','Palm','Magnetic','Whitsunday','Keppel'), gfun, form=analyses[[resp]][['formulas']], USE.NAMES = TRUE,simplify=FALSE)
    if (resp %in% c('PCO1', 'PCO2')) {
      for (ff in 3:6) {  # used to be 2:5 (as in the regions only)
        analyses[[resp]]$formulas[[ff]] = update(analyses[[resp]]$formulas[[ff]],  paste0(fun, '(', resp,'r', ') ~.'))
      }
    }
    analyses[[resp]][['groups']] = sapply(c('all','all1','Palm','Magnetic','Whitsunday','Keppel'), gfun, form=analyses[[resp]]$formulas, USE.NAMES = TRUE,simplify=FALSE)
}
save(analyses, file='data/analyses.RData')
## ----end



## ---- ABT
for (a in 1:length(analyses)) {
    resp=names(analyses)[a]
    print(paste('Response =',resp))
    ## for (f in 1:length(analyses[[a]]$formulas)) {
    for (f in 2)) {  # temporary just so that we can run the new all1 set
        mod.name = names(analyses[[a]]$formulas)[f]
        print(paste('Model =',mod.name))
        MONOTONE = assignMonotone(fish, analyses[[a]]$formulas[[f]])
        if (mod.name %in% c('all','all1')) {fish.sub=fish
        } else {fish.sub = fish %>% filter(REGION==mod.name)}
        if (any(fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))]==0)) {
            val = fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))]
            val=min(val[val>0], na.rm=TRUE)
            fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))] = fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))] + val
        }
        set.seed(123)
        fish.sub = fish.sub %>% mutate_if(is.character,  as.factor)
        mod = abt(analyses[[a]]$formulas[[f]], data=fish.sub, distribution=analyses[[a]]$family,
                  cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                  var.monotone=as.vector(MONOTONE))
        if (f>2) {  # used to be f>1 (so only applies to the regional models)
          var.lookup1 = var.lookup %>%
            mutate(Field.name=ifelse(Field.name %in% c('PCO1', 'PCO2'), paste0(Field.name, 'r'), Field.name),
                   Abbreviation=ifelse(Abbreviation %in% c('PCO1', 'PCO2'), paste0(Abbreviation, 'r'), Abbreviation))
        } else {
          var.lookup1 = var.lookup
        }
        gr <- na.omit(unique(analyses[[a]]$groups[[f]]))
        if (is.logical(gr)) gr=NULL
        p=plot.abts(mod, var.lookup1, center=FALSE, type='response', return.grid=TRUE,
                    groupby=gr,pt.size=14)#
        thresholds=p$thresholds
        save(thresholds, file=paste0('data/thresholds_',mod.name,'_',resp,'.RData'))
        ps = p[['ps']]
        p = p[['p']]
        if ((length(levels(fish$REGION))>1 || length(levels(fish$NTR.Pooled))>1) & mod.name!='all1') {
            p = common_legend(p)
            ps = common_legend(ps)
        }
        ## version for the supplimentary
                                        #do.call('grid.arrange', p) ## version for the supplimentary
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT.png'), do.call('grid.arrange', p), width=15, height=10, dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT.pdf'), do.call('grid.arrange', p), width=15, height=10)

        pred.1=stats.abt(mod, fitMethod=1, analysis = analyses[[a]]$groups[[f]])
        
        save(pred.1, file=paste0('data/pred.1_',mod.name,'_',resp,'.RData'))
        
        optim=summarize_values(pred.1$optim, type='optim') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.optim')
        
        rel.inf = summarize_values(pred.1$rel.imp, type='Rel.inf') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.rel.inf')
        
        R2=summarize_values(pred.1$R2.value, 'R2') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.R2')
        
        stats = optim %>% full_join(R2) %>% full_join(rel.inf) %>%
            left_join(var.lookup %>% dplyr::select(Var=Abbreviation, pretty.name))
        
        save(stats, file=paste0('data/stats_',mod.name,'_',resp,'.RData'))
        write.csv(stats %>% as.data.frame, file=paste0('output/data/stats_',mod.name,'_',resp,'.csv'), quote=FALSE, row.names=FALSE)

        ## Version full model with just the relative importance and the substantial panels
        ps1 = arrangeGrob(ps[[1]],ps[[length(ps)]], widths=c(3,1))
        numberOfPlots = length(ps)-2 #minus 1 since one of the items is the legend and minus one for the relative importance plot
        numberOfPlotRows = ceiling(numberOfPlots / 2)
        if (length(ps)<3) {
          g=ps1 ##arrangeGrob(ps1, ps[2])    #do.call('arrangeGrob', c(ps[c(-1, -length(ps))], list(ncol=2))), heights=c(2,numberOfPlotRows))
        } else {
          g=arrangeGrob(ps1, do.call('arrangeGrob', c(ps[c(-1, -length(ps))], list(ncol=2))), heights=c(2,numberOfPlotRows))
        }
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT_short.png'), g, width=10, height=2+(1.7*numberOfPlotRows), dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT_short.pdf'), g, width=10, height=2+(1.7*numberOfPlotRows))

        ## Version with a grid of partials in the lower right corner of the influence figure
        ## The location and size of the figure in figure will be determined by the scale of relative influence
        ## ideally, we want the subfigure to be 3/4 of the width and 3/4 of the height
        ymax=max(pretty(as.numeric(as.character(stats$upper.rel.inf))))
        ymin=ymax*1/4
        xmax=(length(p)-2)*3/4
        if (ymin <= 100/(length(p)-2)) {  # if the left side of the subfigure is too close to the dashed vertical line
            ymax = ymax + 5
            ymin=ymax*1/4
            ps[[1]] = ps[[1]] + scale_y_continuous('Relative Importance', limits=c(0,ymax))
        }
        ps[2:(length(ps)-1)]=lapply(ps[2:(length(ps)-1)], function(f) f+theme(axis.title.y=element_blank())) 
        ps2=do.call('arrangeGrob', c(ps[c(length(ps),2:(length(ps)-1))], list(ncol=2)))
                                        #g= ps[[1]] + annotation_custom(grob=ps2, xmin=1, xmax=15, ymin=10, ymax=Inf)
        g= ps[[1]] + annotation_custom(grob=ps2, xmin=1, xmax=xmax, ymin=ymin, ymax=Inf)
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT_short_in.png'), g, width=10, height=8, dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt.',mod.name,'_',resp,'_ABT_short_in.pdf'), g, width=10, height=8)
        
        ## Version for vertical table of regional relative importance and substantial panels
        if (length(ps)<2) {
          g = ps[[1]]
        } else {
          g=do.call('arrangeGrob', c(ps[1:4], list(ncol=1, heights=c(1.5, rep(1,3)))))
        }
                                        #save(g, file=paste0('data/g.abt_',mod_num,'_',resp,'.RData'))
        
                                        #save(data.all.abt, file=paste0('data/data.all.abt_',mod_num,'_',resp,'.RData'))
        
        
                                        #save(rel.importance, file=paste0('data/rel.importance_',mod_num,'_',resp,'.RData'))
                                        #save(thresholds, file=paste0('data/thresholds_', mod_num,'_',resp,'.RData'))
                                        #write.table(bind_rows(thresholds, .id='Var'), file=paste0('data/thresholds_',mod_num,'_',resp,'.csv'), quote=FALSE, row.names=FALSE)
        if (1==2) { # all this is to attempt to address co-authors comments
            ## Many of the predictors vary over space and time.  As is, the models just indicate whether there is a relationship between
            ## the response and the various predictors.  Although fitting the models separately for the different Regions does partly tease
            ## apart the spatial from temporal, it does not completely do this (as there is still spatiality within Regions.
            ## One way we might be able to address this is to take the 'important predictors' from a given model and then refit with
            ## specifically centred versions of the predictors.  For example, if we centre a predictor within a site, then we effectively remove
            ## the spatial component of this variable.  Similarly, if we centre on time (year), then we remove the temporal element.

            ## 1. start by identifying the important predictors
            wch=rel.inf$Var[which(as.numeric(rel.inf$Mean.rel.inf)> (100/nrow(rel.inf)))]
            wch.n = wch[unlist(lapply(wch, function(x) is.numeric(fish.sub[,x])))]
            wch.c = wch[unlist(lapply(wch, function(x) is.factor(fish.sub[,x])))]

            ## Remove spatial element
            fish.sub1 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.t')) %>%
                dplyr::select(!!c('TFD','YEAR',paste0(wch.n,'.t')))
            fish.sub1.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                summarize_at(vars(wch.n), function(x) mean(x)) 

            ## #fish.sub1 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC, color=SITE)) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC, color=factor(YEAR))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=SITE)) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=factor(YEAR))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=YEAR, color=factor(SITE))) + facet_wrap(~REGION)
            ## fish.sub1 %>% ggplot() + geom_line(aes(y=LHC.t,x=YEAR, color=factor(SITE))) + facet_wrap(~REGION)
            ## #fish.sub %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC)) + facet_wrap(~REGION)
            ## fish.sub %>% ggplot() + geom_line(aes(y=LHC,x=YEAR, color=factor(SITE))) + facet_wrap(~REGION)
            
            ## Remove temporal element
            fish.sub2 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'YEAR')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.s')) %>%
                dplyr::select(!!c('TFD','SITE',paste0(wch.n,'.s')))

            fish.sub3 = fish.sub1 %>% full_join(fish.sub2)
            
            fish.sub2.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'YEAR')))) %>%
                summarize_at(vars(wch.n), function(x) mean(x)) 


            ## fish.sub1 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC.t, color=factor(YEAR))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=factor(SITE))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=factor(YEAR))) + facet_wrap(~REGION)
            ## fish.sub2 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC.s, color=factor(YEAR))) + facet_wrap(~REGION)
            ## fish.sub %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC)) + facet_wrap(~REGION)

            ## fish.sub1 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC.s, color=factor(SITE))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=factor(SITE))) + facet_wrap(~REGION)
            ## #fish.sub1 %>% ggplot() + geom_line(aes(y=log(TFD),x=LHC, color=factor(YEAR))) + facet_wrap(~REGION)
            ## fish.sub2 %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC, color=factor(SITE))) + facet_wrap(~REGION)
            ## fish.sub %>% ggplot() + geom_point(aes(y=log(TFD),x=LHC)) + facet_wrap(~REGION)


            
            ## fish.sub2 %>% ggplot() + geom_line(aes(y=LHC,x=YEAR, color=factor(SITE))) + facet_wrap(~REGION)
            ## fish.sub2 %>% ggplot() + geom_line(aes(y=LHC,x=as.numeric(SITE), color=factor(YEAR))) + facet_wrap(~REGION)

            p.names=colnames(fish.sub3)[!colnames(fish.sub3) %in% c('SITE','TFD','YEAR')]
            ff = update(analyses[[a]]$formulas[[f]], . ~ 1)
            ff=as.formula(c(gsub('1','',deparse(ff)), paste(p.names, collapse='+')))
            MONOTONE = assignMonotone(fish.sub3, ff)
            set.seed(123)
            mod = abt(ff, data=fish.sub3, distribution=analyses[[a]]$family,
                      cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                      var.monotone=as.vector(MONOTONE))
            summary(mod)
            ## pdf(file='TFD.W.test.pdf', width=10, height=30)
            ## par(mfrow=c(5,2))
            ## #for (i in c(2,7,3,8,4,9,5,10,6,11,1)) plot(mod,i, ylim=c(-1.5,1.5))
            ## for (i in c(1,6,2,7,3,8,4,9,5,10)) plot(mod,i, ylim=c(-1,1))
            ## dev.off()

            ## And now including the other predictors as well
            fish.sub1 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.t')) %>%
                dplyr::select(!!c('TFD','YEAR',paste0(wch.n,'.t')))
            fish.sub2 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'YEAR')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.s')) %>%
                dplyr::select(!!c('TFD','SITE',paste0(wch.n,'.s')))

            fish.sub3 = fish.sub1 %>% full_join(fish.sub2) %>% full_join(fish.sub)
            
            ff = deparse(analyses[[a]]$formulas[[f]])
            for (w in wch.n) ff = gsub(w,paste0(w,'.t + ',w,'.s'), ff)
            MONOTONE = assignMonotone(fish.sub3, ff)
            set.seed(123)
            mod = abt(ff, data=fish.sub3, distribution=analyses[[a]]$family,
                      cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                      var.monotone=as.vector(MONOTONE))
            summary(mod)
            ## Capture the list of influential spatial and temporal predictors...


            
            ## pdf(file='TFD.W.test.pdf', width=10, height=30)
            ## pdf(file='TFD.Region.test.pdf', width=10, height=30)
            ## par(mfrow=c(5,2))
            ##                             #for (i in c(2,7,3,8,4,9,5,10,6,11,1)) plot(mod,i, ylim=c(-1.5,1.5))
            ##                             #for (i in c(5,6,7,8,2,3,19,20,11,22)) plot(mod,i, ylim=c(-1,1))
            ## for (i in c(3,4,8,9,20,21,16,1,6,7)) plot(mod,i, ylim=c(-1,1))
            ## dev.off()

        }
        
        rm(mod,pred.1,p,ps,ps1,g,stats,optim,rel.inf,R2,fit)
        gc()
    }
}
#summary(mod)
#plot(mod,c(15,1))

## ----end

##========================Species specific====================================================

## ---- analysis.list.species
load('data/fish.species.RData')
load('data/var.lookup.species.RData')
resp.lookup.species = var.lookup.species %>% filter(Type=='Response') %>% droplevels
pred.lookup.species = var.lookup.species %>% filter(Type=='Predictor') %>% droplevels

groupings.all.species = vector('list', nrow(pred.lookup.species))
names(groupings.all.species) = pred.lookup.species$Abbreviation
groupings.region.species = vector('list', nrow(pred.lookup.species))
names(groupings.region.species) = pred.lookup.species$Abbreviation
for (i in 1:nrow(pred.lookup.species)) {
    pred = as.character(pred.lookup.species[i,'Abbreviation'])
    groupings.all.species[[pred]] = ifelse(pred=='REGION',NA,'REGION')
    groupings.region.species[[pred]] = ifelse(pred=='NTR.Pooled',NA,'NTR.Pooled')
}
groupings.all.species = do.call('c',groupings.all.species)
groupings.region.species = do.call('c',groupings.region.species)

gfun = function(f, form) {
    #print(f)
    #print(form[[f]])
    form = attr(terms(form[[f]]),'term.labels')
    if (f=='all') return(groupings.all.species[form])
    else return(groupings.region.species[form])
}

analyses.species = vector('list',nrow(resp.lookup.species))
names(analyses.species) = resp.lookup.species$Abbreviation
for (i in 1:nrow(resp.lookup.species)) {
    resp = as.character(resp.lookup.species[i,'Abbreviation'])
    fun = as.character(resp.lookup.species[i,'Transform'])
    fam = as.character(resp.lookup.species[i,'Family'])
    
    analyses.species[[resp]] = list('formulas' = lapply(formulas, function(f) update(f, paste0(fun,'(',resp,') ~.'))),
                            'family' = fam)
    #if (resp %in% c('PI','PTD','PLD')) analyses[[resp]][['formulas']] = lapply(analyses[[resp]][['formulas']], function(f) f=update(f, .~.+PREY.DENSITY))
    #if (resp %in% c('PTB','PLB')) analyses[[resp]][['formulas']] = lapply(analyses[[resp]][['formulas']], function(f) f=update(f, .~.+PREY.BIOMASS))
    analyses.species[[resp]][['groups']] = sapply(c('all','Palm','Magnetic','Whitsunday','Keppel'), gfun, form=analyses.species[[resp]][['formulas']], USE.NAMES = TRUE,simplify=FALSE)
}
save(analyses.species, file='data/analyses.species.RData')
## ----end


## ---- ABT.species.HPC
library(ssh)
session = ssh_connect('mlogan@hpc-login', keyfile='~/.ssh/hpc_id_rsa')
scp_upload(session,
           files=c('FISH_functions.R',
                   'FISH_HPC_gbm.R',
                   'FISH_HPC_gbm.batch'),
           to = "~/tmp/fish", verbose = TRUE)
scp_upload(session,
           files=c('data/analyses.species.RData',
                   'data/fish.species.RData',
                   'data/speciesList.RData',
                   'data/var.lookup.species.RData'),
           to = "~/tmp/fish/data", verbose = TRUE)
ssh_exec_wait(session,
              command = "qsub -l mem=1gb -l nodes=1:ppn=20 ~/tmp/fish/FISH_HPC_gbm.batch\n")
ssh_exec_wait(session, command = "qstat -a")
#ssh_exec_wait(session, command = "qdel 71728")
ssh_disconnect(session)
## ----end

## ---- ABT.species.HPC.receive
session = ssh_connect('mlogan@hpc-login', keyfile='~/.ssh/hpc_id_rsa')
scp_download(session,
             files=c('/export/home/l-p/mlogan/tmp/fish/output/figures/*_in.p*'),
             to = "output/figures")
scp_download(session,
             files=c('/export/home/l-p/mlogan/tmp/fish/data/pred.1*.RData'),
             to = "data/")
scp_download(session,
             files=c('/export/home/l-p/mlogan/tmp/fish/data/stats*.RData'),
             to = "data/")
scp_download(session,
             files=c('/export/home/l-p/mlogan/tmp/fish/data/thresholds*.RData'),
             to = "data/")
# Note, for some reason the above returns an error, even through it is successful
ssh_disconnect(session)
## ----end



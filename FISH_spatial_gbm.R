source('FISH_functions.R')

## ---- VarlookupSpatial
load('data/fish.RData')
load('data/var.lookup.RData')

 var.lookup = var.lookup %>%
    rbind(var.lookup %>%
          mutate(Field.name=paste0(Field.name,'.t'),
                 Abbreviation=paste0(Abbreviation,'.t'))) %>%
    rbind(var.lookup %>%
          mutate(Field.name=paste0(Field.name,'.s'),
                 Abbreviation=paste0(Abbreviation,'.s')))
save(var.lookup, file='data/var.lookup_spatial.RData')
## ----end


## ---- IdentifyInfluentialPredictors
load(file='data/analyses.RData')
new_analyses = analyses
#fish.sub=fish.sub %>% as.data.frame
for (a in 1:length(new_analyses)) {
    resp=names(new_analyses)[a]
    print(paste('Response =',resp))
    for (f in 1:length(new_analyses[[a]]$formulas)) {
        #new_analyses = analyses
        mod.name = names(new_analyses[[a]]$formulas)[f]
        print(paste('Model =',mod.name))

        if (mod.name=='all') {fish.sub=fish
        } else {fish.sub = fish %>% filter(REGION==mod.name)}
        if (any(fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))]==0)) {
            val = fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))]
            val=min(val[val>0], na.rm=TRUE)
            fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))] = fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))] + val
        }

        fish.sub = fish.sub %>% mutate_if(is.character,  as.factor)

        load(file=paste0('data/pred.1_',mod.name,'_',resp,'.RData'))
        rel.inf = summarize_values(pred.1$rel.imp, type='Rel.inf') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.rel.inf')

        ## 1. start by identifying the important predictors
        ## Actually, lets generate spatial and temporal versions of all variables
        ##wch=rel.inf$Var[which(as.numeric(rel.inf$Mean.rel.inf)> (100/nrow(rel.inf)))]
        wch=rel.inf$Var
        if(f==1) wch = unique(c(wch, 'REGION'))
        if(f>1) wch = unique(c(wch, 'NTR.Pooled'))
        wch.n = wch[unlist(lapply(wch, function(x) is.numeric(fish.sub[,x])))]
        wch.c = wch[unlist(lapply(wch, function(x) is.factor(fish.sub[,x])))]

        ## ---- fitPurelyInfluential
        if (1==2) {
            ## 2. Remove spatial element (to produce a temporal version)
            fish.sub1 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x, na.rm=TRUE)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.t')) %>%
                dplyr::select(!!c('TFD','YEAR',paste0(wch.n,'.t')))
            fish.sub1.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                summarize_at(vars(wch.n), function(x) mean(x, na.rm=TRUE)) 

            ## 3. Remove temporal element (to produce a spatial version)
            ## But also ensure that it is not grouped by REGION
            fish.sub2 = fish.sub %>% group_by_at(vars(one_of(c(wch.c[wch.c!='REGION'],'YEAR')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.s')) %>%
                dplyr::select(!!c('TFD','SITE',paste0(wch.n,'.s')))
            fish.sub2.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c[wch.c!='REGION'],'YEAR')))) %>%
                summarize_at(vars(wch.n), function(x) mean(x)) 
            
            ## 4. Join spatial and temporal elements together into a single data frame
            fish.sub3 = fish.sub1 %>% full_join(fish.sub2) %>% ungroup
            
            ## 5. Refit BRT with only influential (spatial and temporal) elements
            p.names=colnames(fish.sub3)[!colnames(fish.sub3) %in% c('SITE','TFD','YEAR')]
            ff = update(new_analyses[[a]]$formulas[[f]], . ~ 1)
            ff=as.formula(c(gsub('1','',deparse(ff)), paste(p.names, collapse='+')))
            MONOTONE = assignMonotone(fish.sub3, ff)
            set.seed(123)
            mod = abt(ff, data=fish.sub3, distribution=new_analyses[[a]]$family,
                      cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                      var.monotone=as.vector(MONOTONE))
            summary(mod)
        }

        ## ----end
        ## ---- fitAllPredictors
        ## 6. And now including the other predictors as well
        ## Remove spatial element (to produce a temporal version)
        fish.sub1 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
            mutate_at(vars(wch.n), function(x) x-mean(x, na.rm=TRUE)) %>%
            dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.t')) %>%
            dplyr::select(!!c('TFD','YEAR',paste0(wch.n,'.t')))
        fish.sub1.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                summarize_at(vars(wch.n), function(x) mean(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.t'))
        fish.sub1.mins = fish.sub %>% group_by_at(vars(one_of(c(wch.c)))) %>%
                summarize_at(vars(wch.n), function(x) min(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.s'))
        fish.sub1.maxs = fish.sub %>% group_by_at(vars(one_of(c(wch.c)))) %>%
                summarize_at(vars(wch.n), function(x) max(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.s'))

        ## Remove temporal element (to produce a spatial version)
        ## But also ensure that it is not grouped by REGION
        fish.sub2 = fish.sub %>% group_by_at(vars(one_of(c(wch.c[wch.c!='REGION'],'YEAR')))) %>%
            mutate_at(vars(wch.n), function(x) x-mean(x, na.rm=TRUE)) %>%
            dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.s')) %>%
            dplyr::select(!!c('TFD','SITE',paste0(wch.n,'.s')))
        fish.sub2.means = fish.sub %>% group_by_at(vars(one_of(c(wch.c[wch.c!='REGION'],'YEAR')))) %>%
            summarize_at(vars(wch.n), function(x) mean(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.s'))
        fish.sub2.mins = fish.sub %>% group_by_at(vars(one_of(c(wch.c)))) %>%
                summarize_at(vars(wch.n), function(x) min(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.t'))
        fish.sub2.maxs = fish.sub %>% group_by_at(vars(one_of(c(wch.c)))) %>%
                summarize_at(vars(wch.n), function(x) max(x, na.rm=TRUE)) %>%
            rename_if(is.numeric, function(x) paste0(x,'.t'))

        new_analyses[[a]]$means <- list()
        new_analyses[[a]]$means[['temporal']] <- fish.sub1.means %>% ungroup
        new_analyses[[a]]$means[['spatial']] <- fish.sub2.means %>% ungroup
        new_analyses[[a]]$mins <- list()
        new_analyses[[a]]$mins[['temporal']] <- fish.sub2.mins %>% ungroup
        new_analyses[[a]]$mins[['spatial']] <- fish.sub1.mins %>% ungroup
        new_analyses[[a]]$maxs <- list()
        new_analyses[[a]]$maxs[['temporal']] <- fish.sub2.maxs %>% ungroup
        new_analyses[[a]]$maxs[['spatial']] <- fish.sub1.maxs %>% ungroup

        fish.sub3 = fish.sub1 %>% full_join(fish.sub2) %>% full_join(fish.sub) %>% ungroup
        
        ff = deparse(new_analyses[[a]]$formulas[[f]])
        for (w in wch.n) {
            ff = gsub(paste(paste0(w,' '),'|',paste0(w,'$'), collapse='', sep=''),paste0(w,'.t + ',w,'.s '), ff)
        }
        MONOTONE = assignMonotone(fish.sub3, ff)
        
        fish.sub3 = fish.sub3 %>% mutate_if(is.character,  as.factor)
        set.seed(123)
        mod = abt(ff, data=fish.sub3, distribution=new_analyses[[a]]$family,
                  cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                  var.monotone=as.vector(MONOTONE))
        summary(mod)
        
        ## Modify the analysis groups to incorporate the spatial and temporal versions.
        for (n in names(new_analyses[[a]]$groups[[f]])) {
            if (n %in% wch.n) {
                new_analyses[[a]]$groups[[f]] = c(new_analyses[[a]]$groups[[f]], setNames(new_analyses[[a]]$groups[[f]][n], paste0(n,'.t')))
                new_analyses[[a]]$groups[[f]] = c(new_analyses[[a]]$groups[[f]], setNames(new_analyses[[a]]$groups[[f]][n], paste0(n,'.s')))
                new_analyses[[a]]$groups[[f]] = new_analyses[[a]]$groups[[f]][-which(names(new_analyses[[a]]$groups[[f]])==n)]
            }
        }
        pred.1=stats.abt(mod, fitMethod=1, analysis = new_analyses[[a]]$groups[[f]])
        rel.inf = summarize_values(pred.1$rel.imp, type='Rel.inf') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.rel.inf')
        VARS=rel.inf %>%
            mutate(Mean.rel.inf = as.numeric(Mean.rel.inf),
                   Flag = Mean.rel.inf > 100/n()) %>%
            arrange(-Mean.rel.inf) %>%
            filter(Flag) %>% pull(Var)
        
        optim=summarize_values(pred.1$optim, type='optim') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.optim')
        
        rel.inf = summarize_values(pred.1$rel.imp, type='Rel.inf') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.rel.inf')
        
        R2=summarize_values(pred.1$R2.value, 'R2') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.R2')
        
        stats = optim %>% full_join(R2) %>% full_join(rel.inf) %>%
            left_join(var.lookup %>% dplyr::select(Var=Abbreviation, pretty.name))
        
        save(stats, file=paste0('data/stats_spatial_temporal_',mod.name,'_',resp,'.RData'))
        
        ## ----end
        ## ---- SpatialAnalyses
        ## remove temporal elements
        if (length(grep('\\.t', VARS))>0) {
          VARS = VARS[-grep('\\.t',VARS)]
        }
        ff = update(as.formula(ff), as.formula(paste0('. ~',paste0(VARS, collapse=' + '))))
        if (f==1) {
            ff= update(ff, ~ . + REGION)
            VARS = unique(c(VARS, 'REGION'))
        }
        if (f > 1) {
            ff = update(ff, ~ . + NTR.Pooled)
            VARS = unique(c(VARS, 'NTR.Pooled'))
        }
        MONOTONE = assignMonotone(fish.sub3, ff)
        set.seed(123)
        mod = abt(ff, data=fish.sub3, distribution=new_analyses[[a]]$family,
                  cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                  var.monotone=as.vector(MONOTONE))
        summary(mod)
        
        ## Modify the analysis groups to incorporate just the selected spatial versions.
        new_analyses[[a]]$groups[[f]] = new_analyses[[a]]$groups[[f]][which(names(new_analyses[[a]]$groups[[f]]) %in% VARS)]

        if (f>1) {
          var.lookup1 = var.lookup %>%
            mutate(Field.name=ifelse(Field.name %in% c('PCO1', 'PCO2'), paste0(Field.name, 'r'), Field.name),
                   Abbreviation=ifelse(Abbreviation %in% c('PCO1', 'PCO2'), paste0(Abbreviation, 'r'), Abbreviation))
        } else {
          var.lookup1 = var.lookup
        }

        p=plot.abts(mod, var.lookup1, center=FALSE, type='response', return.grid=TRUE,
                    groupby=na.omit(unique(new_analyses[[a]]$groups[[f]])),pt.size=14,
                    mins=new_analyses[[a]]$mins[['spatial']],
                    maxs=new_analyses[[a]]$maxs[['spatial']])#
        thresholds=p$thresholds
        save(thresholds, file=paste0('data/thresholds_spatial_',mod.name,'_',resp,'.RData'))
        ps = p[['ps']]
        p = p[['p']]
        if (length(levels(fish$REGION))>1 || length(levels(fish$NTR.Pooled))>1) {
            p = common_legend(p)
            ps = common_legend(ps)
        }
        ## version for the supplimentary
                                        #do.call('grid.arrange', p) ## version for the supplimentary
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT.png'), do.call('grid.arrange', p), width=15, height=10, dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT.pdf'), do.call('grid.arrange', p), width=15, height=10)

        pred.1=stats.abt(mod, fitMethod=1, analysis = new_analyses[[a]]$groups[[f]],
                         mins=new_analyses[[a]]$mins[['spatial']],
                         maxs=new_analyses[[a]]$maxs[['spatial']])
        
        save(pred.1, file=paste0('data/pred.1_spatial_',mod.name,'_',resp,'.RData'))
        optim=summarize_values(pred.1$optim, type='optim') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.optim')
        
        rel.inf = summarize_values(pred.1$rel.imp, type='Rel.inf') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.rel.inf')
        
        R2=summarize_values(pred.1$R2.value, 'R2') %>%
            dplyr::rename_at(vars(-contains('Var')), paste0, '.R2')
        
        stats = optim %>% full_join(R2) %>% full_join(rel.inf) %>%
            left_join(var.lookup %>% dplyr::select(Var=Abbreviation, pretty.name))
        
        save(stats, file=paste0('data/stats_spatial_',mod.name,'_',resp,'.RData'))
        write.csv(stats %>% as.data.frame, file=paste0('output/data/stats_spatial_',mod.name,'_',resp,'.csv'), quote=FALSE, row.names=FALSE)

        ## Version full model with just the relative importance and the substantial panels
        ps1 = arrangeGrob(ps[[1]],ps[[length(ps)]], widths=c(3,1))
        numberOfPlots = length(ps)-2 #minus 1 since one of the items is the legend and minus one for the relative importance plot
        numberOfPlotRows = ceiling(numberOfPlots / 2)
        g=arrangeGrob(ps1, do.call('arrangeGrob', c(ps[c(-1, -length(ps))], list(ncol=2))), heights=c(2,numberOfPlotRows))
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT_short.png'), g, width=10, height=2+(1.7*numberOfPlotRows), dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT_short.pdf'), g, width=10, height=2+(1.7*numberOfPlotRows))

        ## Version with a grid of partials in the lower right corner of the influence figure
        ## The location and size of the figure in figure will be determined by the scale of relative influence
        ## ideally, we want the subfigure to be 3/4 of the width and 3/4 of the height
        ymax=max(pretty(as.numeric(as.character(stats$upper.rel.inf))))
        ymin=ymax*1/4
        xmax=(length(p)-2)*3/4
        if (ymin <= 100/(length(p)-2)) {  # if the left side of the subfigure is too close to the dashed vertical line
            ymax = ymax + 5
            ymin=ymax*1/4
            ymin = 100/(length(p)-2) + 1
            ps[[1]] = ps[[1]] + scale_y_continuous('Relative Importance', limits=c(0,ymax))
        }
        ps[2:(length(ps)-1)]=lapply(ps[2:(length(ps)-1)], function(f) f+theme(axis.title.y=element_blank())) 
        ps2=do.call('arrangeGrob', c(ps[c(length(ps),2:(length(ps)-1))], list(ncol=2)))
                                        #g= ps[[1]] + annotation_custom(grob=ps2, xmin=1, xmax=15, ymin=10, ymax=Inf)
        g= ps[[1]] +
            #scale_y_continuous(limits=c(0,max(50,max(pretty(as.numeric(as.character(stats$upper.rel.inf))))))) +
            annotation_custom(grob=ps2, xmin=1, xmax=xmax, ymin=ymin, ymax=Inf)
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT_short_in.png'), g, width=10, height=8, dpi=300)
        ggsave(filename=paste0('output/figures/data.all.abt_spatial.',mod.name,'_',resp,'_ABT_short_in.pdf'), g, width=10, height=8)

        ## Version for vertical table of regional relative importance and substantial panels
        g=do.call('arrangeGrob', c(ps[1:4], list(ncol=1, heights=c(1.5, rep(1,3)))))
        ## ----end

        ## ---- fitAllPredictorsWithoutRegion
        if (1==2) {
            ## 6. And now including the other predictors as well
            ## Remove spatial element (to produce a temporal version)
            fish.sub1 = fish.sub %>% group_by_at(vars(one_of(c(wch.c,'SITE')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x, na.rm=TRUE)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.t')) %>%
                dplyr::select(!!c('TFD','YEAR',paste0(wch.n,'.t')))
            ## Remove temporal element (to produce a spatial version)
            ## But also ensure that it is not grouped by REGION
            fish.sub2 = fish.sub %>% group_by_at(vars(one_of(c(wch.c[wch.c!='REGION'],'YEAR')))) %>%
                mutate_at(vars(wch.n), function(x) x-mean(x, na.rm=TRUE)) %>%
                dplyr::rename_at(vars(wch.n), ~paste0(wch.n,'.s')) %>%
                dplyr::select(!!c('TFD','SITE',paste0(wch.n,'.s')))
            
            fish.sub3 = fish.sub1 %>% full_join(fish.sub2) %>% full_join(fish.sub) %>% ungroup

            ff = deparse(analyses[[a]]$formulas[[f]])
            for (w in wch.n) ff = gsub(w,paste0(w,'.t + ',w,'.s'), ff)
            ff = gsub('REGION +','', ff)
            MONOTONE = assignMonotone(fish.sub3, ff)
            set.seed(123)
            mod = abt(ff, data=fish.sub3, distribution=analyses[[a]]$family,
                      cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                      var.monotone=as.vector(MONOTONE))
            summary(mod)
        }
        ## ----end

        rm(mod,pred.1,p,ps,ps1,g,stats,optim,rel.inf,R2,fit)
        gc()
    }
}
#save(new_analyses, file='data/new_analyses.RData')
save(new_analyses, file='data/new_analyses1.RData')
## ----end



fish.sub3 = fish.sub1 %>% full_join(fish.sub2) %>% full_join(fish.sub) %>% ungroup
fish.sub3 = fish.sub3 %>% mutate(SSTMEAN=SSTMEAN-mean(SSTMEAN, na.rm=TRUE))
MONOTONE = assignMonotone(fish.sub3, ff)
set.seed(123)
ff = formula(log(TFD) ~ REGION+SSTMEAN)
MONOTONE = assignMonotone(fish.sub3, ff)
mod = abt(ff, data=fish.sub3, distribution=analyses[[a]]$family,
          cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
          var.monotone=as.vector(MONOTONE))
summary(mod)






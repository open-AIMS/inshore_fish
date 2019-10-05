source('FISH_functions.R')



## ---- formulas
formulas = list(
    all = ~REGION + NTR.Pooled + LHC + SC + MA + TURF + UC + BR + CMD + SLOPE + RUG + SCI + CHL + KD490 + SSTMEAN + SSTANOM + WAVE + DEPTH + DHW + CYCLONE + EXP,
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
names(groupings.all) = pred.lookup$Abbreviation
groupings.region = vector('list', nrow(pred.lookup))
names(groupings.region) = pred.lookup$Abbreviation
for (i in 1:nrow(pred.lookup)) {
    pred = as.character(pred.lookup[i,'Abbreviation'])
    groupings.all[[pred]] = ifelse(pred=='REGION',NA,'REGION')
    groupings.region[[pred]] = ifelse(pred=='NTR.Pooled',NA,'NTR.Pooled')
}
groupings.all = do.call('c',groupings.all)
groupings.region = do.call('c',groupings.region)

gfun = function(f, form) {
    #print(f)
    #print(form[[f]])
    form = attr(terms(form[[f]]),'term.labels')
    if (f=='all') return(groupings.all[form])
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
    analyses[[resp]][['groups']] = sapply(c('all','Palm','Magnetic','Whitsunday','Keppel'), gfun, form=analyses[[resp]][['formulas']], USE.NAMES = TRUE,simplify=FALSE)
}
## ----end



## ---- ABT
for (a in 1:length(analyses)) {
    resp=names(analyses)[a]
    print(paste('Response =',resp))
    for (f in 1:length(analyses[[a]]$formulas)) {
        mod.name = names(analyses[[a]]$formulas)[f]
        print(paste('Model =',mod.name))
        MONOTONE = assignMonotone(fish, analyses[[a]]$formulas[[f]])
        if (mod.name=='all') {fish.sub=fish
        } else {fish.sub = fish %>% filter(REGION==mod.name)}
        if (any(fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))]==0)) {
            val = fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))]
            val=min(val[val>0], na.rm=TRUE)
            fish.sub[, as.character(get_response(analyses[[a]]$formulas[[f]]))] = fish.sub[,as.character(get_response(analyses[[a]]$formulas[[f]]))] + val
        }
        set.seed(123)
        mod = abt(analyses[[a]]$formulas[[f]], data=fish.sub, distribution=analyses[[a]]$family,
                  cv.folds=10,interaction.depth=10,n.trees=10000, shrinkage=0.001, n.minobsinnode=2,
                  var.monotone=as.vector(MONOTONE))
        p=plot.abts(mod, var.lookup, center=FALSE, type='response', return.grid=TRUE,
                    groupby=na.omit(unique(analyses[[a]]$groups[[f]])))#
        thresholds=p$thresholds
        save(thresholds, file=paste0('data/thresholds_',mod.name,'_',resp,'.RData'))
        ps = p[['ps']]
        p = p[['p']]
        if (length(levels(fish$REGION))>1 || length(levels(fish$NTR.Pooled))>1) {
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
        g=arrangeGrob(ps1, do.call('arrangeGrob', c(ps[c(-1, -length(ps))], list(ncol=2))), heights=c(2,numberOfPlotRows))
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
        g=do.call('arrangeGrob', c(ps[1:4], list(ncol=1, heights=c(1.5, rep(1,3)))))
                                        #save(g, file=paste0('data/g.abt_',mod_num,'_',resp,'.RData'))
        
                                        #save(data.all.abt, file=paste0('data/data.all.abt_',mod_num,'_',resp,'.RData'))
        
                
        #save(rel.importance, file=paste0('data/rel.importance_',mod_num,'_',resp,'.RData'))
        #save(thresholds, file=paste0('data/thresholds_', mod_num,'_',resp,'.RData'))
                                        #write.table(bind_rows(thresholds, .id='Var'), file=paste0('data/thresholds_',mod_num,'_',resp,'.csv'), quote=FALSE, row.names=FALSE)
        rm(mod,pred.1,p,ps,ps1,g,stats,optim,rel.inf,R2,fit)
        gc()
    }
}
#summary(mod)
#plot(mod,c(15,1))

## ----end

source('FISH_functions.R')

## ---- readFish
#fish = read.csv('data/Fish benthic physical sitelevel 2019.csv', strip.white=TRUE)
fish = read.csv('data/Selected fish benthic physical sitelevel 2019.csv', strip.white=TRUE)
fish %>% glimpse
## ----end

## ---- nameLookup
var.lookup = rbind(
    data.frame(pretty.name='Total density', Field.name='Total.fish.density', Abbreviation='TFD', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Total species richness', Field.name='Total.fish.species.richness', Abbreviation='TFSR', Family='gaussian', Type='Response', Transform='I', Groupby=''),
    data.frame(pretty.name='Benthic invertivores', Field.name='BE', Abbreviation='BE', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Grazers', Field.name='GRAZ', Abbreviation='GRAZ', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Grazers2', Field.name='GRAZ2', Abbreviation='GRAZ2', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Parrot', Field.name='Parrot', Abbreviation='PA', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Corallivores', Field.name='COR', Abbreviation='COR', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Omnivores', Field.name='OM', Abbreviation='OM', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Planktivores', Field.name='PL', Abbreviation='PL', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Carnivores', Field.name='CA', Abbreviation='CA', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Piscivores', Field.name='PI', Abbreviation='PI', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Farmers', Field.name='FA', Abbreviation='FA', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Plectropomus total density', Field.name='Plectropomus.total.density', Abbreviation='PTD', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Plectropomus total biomass', Field.name='Plectropomus.total.biomass', Abbreviation='PTB', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Plectropomus legal density', Field.name='Plectropomus.legal.density', Abbreviation='PLD', Family='gaussian', Type='Response', Transform='log', Groupby=''),
    data.frame(pretty.name='Plectropomus legal biomass', Field.name='Plectropomus.legal.biomass', Abbreviation='PLB', Family='gaussian', Type='Response', Transform='log', Groupby=''),

    data.frame(pretty.name='Region', Field.name='REGION', Abbreviation='REGION', Family=NA, Type='Predictor', Transform='I', Groupby=''),
    data.frame(pretty.name='NTR Pooled', Field.name='NTR.Pooled', Abbreviation='NTR.Pooled', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='LHC % (live hard coral cover)', Field.name='LHC_.', Abbreviation='LHC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SC % (soft coral)', Field.name='SC_.', Abbreviation='SC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='MA % (macroalgal cover)', Field.name='MAC_.', Abbreviation='MA', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Turf %', Field.name='Turf_.', Abbreviation='TURF', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Unconsolidated %', Field.name='Unconsolidated_.', Abbreviation='UC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Benthic richness', Field.name='Benthic.richness', Abbreviation='BR', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Coral morphological diversity', Field.name='Coral_Morph.Diversity', Abbreviation='CMD', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Slope', Field.name='slope', Abbreviation='SLOPE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Rugosity', Field.name='rugosity', Abbreviation='RUG', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SCI (structural complexity)', Field.name='SCI', Abbreviation='SCI', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Chla', Field.name='ChlA', Abbreviation='CHL', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Kd490', Field.name='kd490', Abbreviation='KD490', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SST mean', Field.name='SSTmean', Abbreviation='SSTMEAN', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SST anom', Field.name='SSTanom', Abbreviation='SSTANOM', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Wave exposure', Field.name='wave.exposure.index', Abbreviation='WAVE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Corrected depth', Field.name='Corrected.depth', Abbreviation='DEPTH', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Max DHW', Field.name='maxDHW', Abbreviation='DHW', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Cyclone', Field.name='Cyclone', Abbreviation='CYCLONE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Exposure to primary weeks', Field.name='Exposure.to.primary.weeks', Abbreviation='EXP', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Prey density', Field.name='Prey.density', Abbreviation='PREY.DENSITY', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Prey biomass', Field.name='Prey.biomass', Abbreviation='PREY.BIOMASS', Family=NA, Type='Predictor', Transform='I', Groupby='Region')
)
save(var.lookup, file='data/var.lookup.RData') 
names = with(var.lookup, setNames(as.character(Field.name), Abbreviation))
## exclude those whose names equal their values otherwise there will be duplicate fields created in the fish data
names = names[names(names)!=names]
## ----end

## ---- AbbreviatedNames
## Create duplicates of the fields that are not already abbreviated 
fish = fish %>%
    mutate(SSTmean=ifelse(as.character(SSTmean)=='#N/A',NA,as.numeric(as.character(SSTmean)))) %>%
    mutate_at(as.character(var.lookup$Field.name), list(A=~I)) %>%
    rename(!!! gsub('(.*)','\\1_A',names)) %>%
    mutate(REGION=factor(REGION, levels=c('Palm','Magnetic','Whitsunday','Keppel')))
save(fish, file='data/fish.RData')
## ----end


## Exploratory data analyses

## ---- EDA.TFD
EDA_histograms(var='TFD', dat=fish, var.lookup=var.lookup)
EDA_density(var='TFD', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.TFSR
EDA_histograms(var='TFSR', dat=fish, var.lookup=var.lookup)
EDA_density(var='TFSR', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.BE
EDA_histograms(var='BE', dat=fish, var.lookup=var.lookup)
EDA_density(var='BE', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.GRAZ
EDA_histograms(var='GRAZ', dat=fish, var.lookup=var.lookup)
EDA_density(var='GRAZ', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.GRAZ2
EDA_histograms(var='GRAZ2', dat=fish, var.lookup=var.lookup)
EDA_density(var='GRAZ2', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PA
EDA_histograms(var='PA', dat=fish, var.lookup=var.lookup)
EDA_density(var='PA', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.COR
EDA_histograms(var='COR', dat=fish, var.lookup=var.lookup)
EDA_density(var='COR', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.OM
EDA_histograms(var='OM', dat=fish, var.lookup=var.lookup)
EDA_density(var='OM', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PL
EDA_histograms(var='PL', dat=fish, var.lookup=var.lookup)
EDA_density(var='PL', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.CA
EDA_histograms(var='CA', dat=fish, var.lookup=var.lookup)
EDA_density(var='CA', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PI
EDA_histograms(var='PI', dat=fish, var.lookup=var.lookup)
EDA_density(var='PI', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.FA
EDA_histograms(var='FA', dat=fish, var.lookup=var.lookup)
EDA_density(var='FA', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PTD
EDA_histograms(var='PTD', dat=fish, var.lookup=var.lookup)
EDA_density(var='PTD', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PTB
EDA_histograms(var='PTB', dat=fish, var.lookup=var.lookup)
EDA_density(var='PTB', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PLD
EDA_histograms(var='PLD', dat=fish, var.lookup=var.lookup)
EDA_density(var='PLD', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PLB
EDA_histograms(var='PLB', dat=fish, var.lookup=var.lookup)
EDA_density(var='PLB', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end


## ---- readSelections
selections = read.csv('data/Species selection.csv', strip.white=TRUE)
selections %>% glimpse
## ----end


## Predictors

## ---- EDA.REGION
EDA_histograms(var='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.NTR.Pooled
EDA_histograms(var='NTR.Pooled', dat=fish, var.lookup=var.lookup)
EDA_histograms(var='NTR.Pooled', group='REGION',dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.LHC
EDA_histograms(var='LHC', dat=fish, var.lookup=var.lookup)
EDA_density(var='LHC', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.SC
EDA_histograms(var='SC', dat=fish, var.lookup=var.lookup)
EDA_density(var='SC', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.MA
EDA_histograms(var='MA', dat=fish, var.lookup=var.lookup)
EDA_density(var='MA', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.TURF
EDA_histograms(var='TURF', dat=fish, var.lookup=var.lookup)
EDA_density(var='TURF', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.UC
EDA_histograms(var='UC', dat=fish, var.lookup=var.lookup)
EDA_density(var='UC', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.BR
EDA_histograms(var='BR', dat=fish, var.lookup=var.lookup)
EDA_density(var='BR', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.CMD
EDA_histograms(var='CMD', dat=fish, var.lookup=var.lookup)
EDA_density(var='CMD', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.SLOPE
EDA_histograms(var='SLOPE', dat=fish, var.lookup=var.lookup)
EDA_density(var='SLOPE', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.RUG
EDA_histograms(var='RUG', dat=fish, var.lookup=var.lookup)
EDA_density(var='RUG', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.SCI
EDA_histograms(var='SCI', dat=fish, var.lookup=var.lookup)
EDA_density(var='SCI', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.CHL
EDA_histograms(var='CHL', dat=fish, var.lookup=var.lookup)
EDA_density(var='CHL', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.KD490
EDA_histograms(var='KD490', dat=fish, var.lookup=var.lookup)
EDA_density(var='KD490', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.SSTMEAN
EDA_histograms(var='SSTMEAN', dat=fish, var.lookup=var.lookup)
EDA_density(var='SSTMEAN', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.SSTANOM
EDA_histograms(var='SSTANOM', dat=fish, var.lookup=var.lookup)
EDA_density(var='SSTANOM', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.WAVE
EDA_histograms(var='WAVE', dat=fish, var.lookup=var.lookup)
EDA_density(var='WAVE', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.DEPTH
EDA_histograms(var='DEPTH', dat=fish, var.lookup=var.lookup)
EDA_density(var='DEPTH', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.DHW
EDA_histograms(var='DHW', dat=fish, var.lookup=var.lookup)
EDA_density(var='DHW', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.CYCLONE
EDA_histograms(var='CYCLONE', dat=fish, var.lookup=var.lookup)
EDA_density(var='CYCLONE', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.EXP
EDA_histograms(var='EXP', dat=fish, var.lookup=var.lookup)
EDA_density(var='EXP', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PREY.DENSITY
EDA_histograms(var='PREY.DENSITY', dat=fish, var.lookup=var.lookup)
EDA_density(var='PREY.DENSITY', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end

## ---- EDA.PREY.BIOMASS
EDA_histograms(var='PREY.BIOMASS', dat=fish, var.lookup=var.lookup)
EDA_density(var='PREY.BIOMASS', group='REGION', dat=fish, var.lookup=var.lookup)
## ----end


## ---- readSpeciesList
speciesList = read.csv('data/Species selection.csv', strip.white=TRUE)
speciesList %>% glimpse
speciesList = speciesList %>% dplyr::select(-X) %>% dplyr::rename(`all.model`=`All.regions.model`, `Palm.model`=`Palms.model`)
## ----end

## ---- nameLookup.species
speciesNames = levels(unlist(speciesList))
speciesNames = speciesNames[speciesNames!=""]
## The names dont quite line up between the two data sets
## The following lines are to help them match
speciesNames = gsub(' ','.',speciesNames)
speciesNames = speciesNames[speciesNames!="pse.tuka"]


speciesList = sapply(speciesNames, function(x) colSums(speciesList==x[1]), simplify=FALSE,USE.NAMES = TRUE)
save(speciesList, file='data/speciesList.RData')

var.lookup.species = do.call('rbind',
                             lapply(speciesNames, function(x) data.frame(pretty.name=x, Field.name=x,
                                                                              Abbreviation=x, Family='poisson',
                                                                              Type='Response', Transform='I',
                                                                              Groupby='')
                                    )
                             )
var.lookup.species = rbind(
    var.lookup.species,
    data.frame(pretty.name='Region', Field.name='REGION', Abbreviation='REGION', Family=NA, Type='Predictor', Transform='I', Groupby=''),
    data.frame(pretty.name='NTR Pooled', Field.name='NTR.Pooled', Abbreviation='NTR.Pooled', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='LHC % (live hard coral cover)', Field.name='LHC_.', Abbreviation='LHC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SC % (soft coral)', Field.name='SC_.', Abbreviation='SC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='MA % (macroalgal cover)', Field.name='MAC_.', Abbreviation='MA', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Turf %', Field.name='Turf_.', Abbreviation='TURF', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Unconsolidated %', Field.name='Unconsolidated_.', Abbreviation='UC', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Benthic richness', Field.name='Benthic.richness', Abbreviation='BR', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Coral morphological diversity', Field.name='Coral_Morph.Diversity', Abbreviation='CMD', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Slope', Field.name='slope', Abbreviation='SLOPE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Rugosity', Field.name='rugosity', Abbreviation='RUG', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SCI (structural complexity)', Field.name='SCI', Abbreviation='SCI', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Chla', Field.name='ChlA', Abbreviation='CHL', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Kd490', Field.name='kd490', Abbreviation='KD490', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SST mean', Field.name='SSTmean', Abbreviation='SSTMEAN', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='SST anom', Field.name='SSTanom', Abbreviation='SSTANOM', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Wave exposure', Field.name='wave.exposure.index', Abbreviation='WAVE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Corrected depth', Field.name='Corrected.depth', Abbreviation='DEPTH', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Max DHW', Field.name='maxDHW', Abbreviation='DHW', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Cyclone', Field.name='Cyclone', Abbreviation='CYCLONE', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Exposure to primary weeks', Field.name='Exposure.to.primary.weeks', Abbreviation='EXP', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Prey density', Field.name='Prey.density', Abbreviation='PREY.DENSITY', Family=NA, Type='Predictor', Transform='I', Groupby='Region'),
    data.frame(pretty.name='Prey biomass', Field.name='Prey.biomass', Abbreviation='PREY.BIOMASS', Family=NA, Type='Predictor', Transform='I', Groupby='Region')
)
save(var.lookup.species, file='data/var.lookup.species.RData') 
names.species = with(var.lookup.species, setNames(as.character(Field.name), Abbreviation))
## exclude those whose names equal their values otherwise there will be duplicate fields created in the fish data
names.species = names.species[names(names.species)!=names.species]
## ----end

## ---- AbbreviatedNamesSpecies
## Create duplicates of the fields that are not already abbreviated
fish = read.csv('data/Fish benthic physical sitelevel 2019.csv', strip.white=TRUE)

fish.species = fish %>%
    mutate(SSTmean=ifelse(as.character(SSTmean)=='#N/A',NA,as.numeric(as.character(SSTmean)))) %>%
    mutate_at(as.character(var.lookup.species$Field.name), list(A=~I)) %>%
    rename(!!! gsub('(.*)','\\1_A',names.species)) %>%
    mutate(REGION=factor(REGION, levels=c('Palm','Magnetic','Whitsunday','Keppel')))
save(fish.species, file='data/fish.species.RData')
## ----end




## ---- SpeciesEDALoopOld
load('data/var.lookup.species.RData')
resp.lookup = var.lookup.species %>% filter(Type=='Response') %>% droplevels

plots = lapply(speciesNames, function(x) {
    p1=EDA_histograms_grob(var=x, dat=fish, var.lookup=var.lookup.species)
    p2=EDA_density_grob(var=x, group='REGION', dat=fish, var.lookup=var.lookup.species)
    list('Hist'=p1,'Dens'=p2)
})
## Need to make this into a single level list
plots=unlist(plots, recursive = FALSE)
paths = paste0(apply(expand.grid(c('eda_','dens_'), speciesNames), 1, paste0, collapse=""), '.png')
pwalk(list(paths,plots), ggsave, path='fishanalysis_species_files/figure-html', width=8, height=5, dpi=300)
## ----end



## ---- SpeciesEDALoop
load('data/var.lookup.species.RData') 
resp.lookup = var.lookup.species %>% filter(Type=='Response') %>% droplevels

for (s in speciesNames) {
                                        #print(s)
    resp=s
    cat(paste('## ',resp.lookup$pretty.name[resp.lookup$Abbreviation==resp],' {.tabset .tabset-pills} \n\n'))
    EDA_histograms(var=s, dat=fish, var.lookup=var.lookup.species)
    cat('\n\n')
    EDA_density(var=s, group='REGION', dat=fish, var.lookup=var.lookup.species)
    cat('\n\n')
}

## ----end



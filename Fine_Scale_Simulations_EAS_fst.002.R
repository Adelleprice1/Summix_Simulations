
setwd("/storage/math/projects/gnomad-public/gnomADv3.1.2/chr22")

#Create vector of fine-scale ancestries to be simulated
anc_list <- c("dai_cdx", "lahu", "naxi_yizu", "oroqen", "hezhen", "she", "tujia", "xibo", "cambodian", "miaozu", "daur_mongola", "tu", "jpt_japanese", "khv", "chb_han", "chs")
#Create parameter dataset

col_names <- paste(c("A", "S" , "E"), rep(1:length(anc_list), each = 3), sep = "")
a <- c("dai_cdx", "0", "0", "lahu", "0", "0", "naxi_yizu", "0", "0", "oroqen", "0", "0", "hezhen", "0", "0", "she", "0", "0", "tujia", "0", "0", "xibo", "0", "0", "cambodian", "0", "0", "miaozu", "0", "0", "daur_mongola", "0", "0", "tu", "0", "0", "jpt_japanese", ".499", ".501", "khv", "0", "0", "chb_han", ".499", ".501", "chs", "0", "0")
params <- data.frame(matrix(a, nrow = 1, ncol = length(col_names), dimnames = list(c(), col_names)))
b <- c("dai_cdx", ".099", ".101", "lahu", "0", "0", "naxi_yizu", ".099", ".101", "oroqen", "0", "0", "hezhen", "0", "0", "she", "0", "0", "tujia", "0", "0", "xibo", "0", "0", "cambodian", "0", "0", "miaozu", "0", "0", "daur_mongola", ".099", ".101", "tu", "0", "0", "jpt_japanese", ".199", ".201", "khv", ".099", ".101", "chb_han", ".199", ".201", "chs", ".199", ".201")
c <- c("dai_cdx", ".0624", ".0626", "lahu", ".0624", ".0626", "naxi_yizu", ".0624", ".0626", "oroqen", ".0624", ".0626", "hezhen", ".0624", ".0626", "she", ".0624", ".0626", "tujia", ".0624", ".0626", "xibo", ".0624", ".0626", "cambodian", ".0624", ".0626", "miaozu", ".0624", ".0626", "daur_mongola", ".0624", ".0626", "tu", ".0624", ".0626", "jpt_japanese", ".0624", ".0626", "khv", ".0624", ".0626", "chb_han", ".0624", ".0626", "chs", ".0624", ".0626")
params <- as.data.frame(rbind(params, b, c))
params[,grep('S|E', names(params))] <- lapply(params[,grep('S|E', names(params))], as.character)
params[,grep('S|E', names(params))] <- lapply(params[,grep('S|E', names(params))], as.numeric)

#read in merged and AF filtered data
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
library(Summix)
chr22 <- read.table("chr22AFfilteredallNAremoved.txt.gz", header = TRUE)
chr22$AFjpt_japanese <- (chr22$ANjpt * chr22$AFjpt + chr22$ANjapanese * chr22$AFjapanese) / (chr22$ANjpt+ chr22$ANjapanese)
chr22$AFnaxi_yizu <- (chr22$ANnaxi * chr22$AFnaxi + chr22$ANyizu * chr22$AFyizu) / (chr22$ANnaxi+ chr22$ANyizu)
chr22$AFdaur_mongola <- (chr22$ANdaur * chr22$AFdaur + chr22$ANmongola * chr22$AFmongola) / (chr22$ANdaur+ chr22$ANmongola)
chr22$AFchb_han <- (chr22$ANchb * chr22$AFchb + chr22$ANhan * chr22$AFhan) / (chr22$ANchb+ chr22$ANhan)
chr22$AFdai_cdx <- (chr22$ANdai * chr22$AFdai + chr22$ANcdx * chr22$AFcdx) / (chr22$ANdai+ chr22$ANcdx)

chr22$AF_EAS <- (chr22$ANcdx * chr22$AFcdx + chr22$ANchb * chr22$AFchb + chr22$ANchs * chr22$AFchs + chr22$ANjpt * chr22$AFjpt + chr22$ANkhv * chr22$AFkhv + chr22$ANcambodian * chr22$AFcambodian + chr22$ANdai * chr22$AFdai + chr22$ANdaur * chr22$AFdaur + chr22$ANhan * chr22$AFhan + chr22$ANhezhen * chr22$AFhezhen + chr22$ANjapanese * chr22$AFjapanese + chr22$ANlahu * chr22$AFlahu + chr22$ANmiaozu * chr22$AFmiaozu + chr22$ANmongola * chr22$AFmongola + chr22$ANnaxi * chr22$AFnaxi + chr22$ANoroqen * chr22$AForoqen + chr22$ANshe * chr22$AFshe + chr22$ANtu * chr22$AFtu + chr22$ANtujia * chr22$AFtujia + chr22$ANxibo * chr22$AFxibo + chr22$ANyizu * chr22$AFyizu)/ (chr22$ANcdx + chr22$ANchb + chr22$ANchs + chr22$ANjpt + chr22$ANkhv + chr22$ANcambodian + chr22$ANdai + chr22$ANdaur + chr22$ANhan + chr22$ANhezhen + chr22$ANjapanese + chr22$ANlahu  + chr22$ANmiaozu + chr22$ANmongola + chr22$ANnaxi + chr22$ANoroqen + chr22$ANshe + chr22$ANtu + chr22$ANtujia + chr22$ANxibo + chr22$ANyizu)
#filter for MAF<.01 in simulated observed ancestry
chr22flt <- chr22[chr22$AF_EAS>0.01 & chr22$AF_EAS<0.99,]
#remove full chr22 file
rm(chr22)

# Number which calls specific paramaters for simulation, increasing
ivalnum = 1

# Number of simulations to run
testnum = 100
# Number of people to simulate
ntot = 10000


# Load in reference data and paramaters
parameters = params
nodecontrol = ceiling(dim(parameters)[1]/20)
tivec = c(1, round(seq(1, dim(parameters)[1], by = dim(parameters)[1]/nodecontrol))[-1], dim(parameters)[1])

# Pulls reference group names from parameter file
AncFrame = data.frame(parameters[,grep('A', names(parameters))])
# Pulls proportion windows from parameter file
SEframe = data.frame(parameters[,grep('S|E', names(parameters))])
# Calculates observed intervals
testint = c(ifelse(ivalnum == 1, tivec[ivalnum], tivec[ivalnum] + 1), tivec[ivalnum+1])

#Create final simulation data frame output
finalframe = data.frame(matrix(vector(), 0, (6 + 3*length(anc_list) + dim(parameters)[2] + 1),
                               dimnames=list(c(), c('P_Num', 'T_Num', 'Seed', names(parameters), 
                                               paste(anc_list, rep("_parameter", each = length(anc_list)), sep = ""),
                                               'Sum_obj', 'Sum_iterations', 'Sum_time', 'Sum_filt',
                                               paste(anc_list, rep("_estimate", each = length(anc_list)), sep = ""),  
                                               paste(anc_list, rep("_accuracy", each = length(anc_list)), sep = "")))), stringsAsFactors=F)



# Simulations
for (m in testint[1]:testint[2]){
  
  # Pull data from reference
  propvec = SEframe[m,]
  ancvec = AncFrame[m,]
  ancnum = dim(AncFrame)[2]
  outframe = finalframe[FALSE,]
  
  for (j in 1:testnum){
    
    # Set Seed
    seed = as.integer(Sys.time())
    set.seed(seed)


    # Select 100K SNPS from reference data
    refdat = chr22flt %>% 
        sample_n(100000) %>% 
        select(CHROM, POS, REF, ALT, paste(rep("AF", each = length(anc_list)), anc_list, sep="")) 
    names(refdat) <- gsub("AF", "", names(refdat))


    # Determines proportions of population
    popvecprop = function(propvec){
      prop = numeric(length = (dim(propvec)[2]/2))

      starts = propvec %>% 
        select(starts_with('S'))
      ends = propvec %>% 
        select(starts_with('E'))
      starts <- as.numeric(starts)
      ends <- as.numeric(ends)

      for (i in (1:length(prop))){
        prop[i] = runif(1, min = starts[i], max = ends[i])
      }

      return(prop)
    }


        # Generate vector of proportions of ancestry
        prop = popvecprop(propvec)


        popnum = numeric(ancnum)
        for (i in 1:ancnum){
            popnum[i] = floor((ntot * prop[i]))
          }

        # Initialize matrix for simulated population
        pop_matrix = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = 3))
        # Generate allele counts for population using multinom 
        for (i in 1:ancnum){
          popmatrixadd = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, popnum[i], prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
          pop_matrix = pop_matrix + popmatrixadd
        }
        # MAF threshold to filter by
        MAF_thresh = 0.01

        # Calcluate allele frequencies
        master_frame_gen1 <- data.frame(refdat[,c(1:4)], pop_matrix)
        master_frame_gen1$AF <- (2 * master_frame_gen1[,5] + master_frame_gen1[,6]) / (2 * ntot)
        # Filter by MAF
        master_frame_gen2 <- master_frame_gen1[master_frame_gen1$AF > MAF_thresh & master_frame_gen1$AF < (1-MAF_thresh),]


        # Pull reference data
        refdatm = refdat %>% 
          select(POS, REF, ALT, all_of(anc_list))
        # Set Observed data
        obsvecm = master_frame_gen2 %>% 
          select(POS, REF, ALT, AF)

        # Merge reference and observed
        mergeframe = merge(refdatm, obsvecm, by = c("POS","REF","ALT")) %>% 
          select(POS, REF, ALT, all_of(anc_list), AF)


    # Run Summix
        R_sum = summix(mergeframe, reference = all_of(anc_list), observed = "AF")

    # Create dataframe with set parameter values for each fine-scale reference ancestry
    ancs <- anc_list
    propinfo <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        props = prop[i]
        propinfo[,i] = props
    }
    names(propinfo)<-c(paste((anc_list), rep("_prop", each=length(anc_list)), sep="_"))

    # Create dataframe with Summix estimations for each fine-scale reference ancestry
    R_ests <- R_sum[,5:length(R_sum)]
    Sum_est <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        ests = R_ests[i]
        Sum_est[,i] = ests
    }
    names(Sum_est)<-c(paste((anc_list), rep("_est", each=length(anc_list)), sep="_"))

    # Create dataframe with accuracy of Summix estimations (estimations-parameters)
    Sum_acc <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        acc = (Sum_est[,i] - propinfo[,i])
        Sum_acc[,i] = acc
    }
    names(Sum_acc)<-c(paste((anc_list), rep("_acc", each=length(anc_list)), sep="_"))

    # Save test info for each simulation
    testinfo = data.frame(P_Num = m, T_Num = j, Seed = seed)
    # Save parameters
    parinfo = parameters[m,]
    # Convert to character
    ancvec[] <- lapply(ancvec, as.character)

    # Bind Data
    outline = cbind(testinfo, parinfo, propinfo, R_sum[1], R_sum[2], R_sum[3], R_sum[4], Sum_est, Sum_acc)
    outframe[j,] = outline
  }
  
  # Bind Data
  finalframe = rbind(finalframe, outframe)
  }

write.table(finalframe, file="Fine_Scale_Simulations_EAS_fst.002.testrun.text", sep="\t", quote = FALSE, row.names=F,col.names=T)

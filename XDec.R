#################################################################
##                                                             ##
## R Script to run the deconvolution algorithm of RNA-Seq data ##
##                                                             ##
## Author: Oscar Murillo (odmurill@bcm.edu)                    ##
##                                                             ##
## Version 1.0.1 (2019-01-10)                                  ##
## Version 1.0.2 (2019-03-31)                                  ##
##                                                             ##
#################################################################
library(optparse)
library(EDec)
library(gplots)
library(ggplot2)
library(reshape2)

## RNA Sets
OriginalSet.RNA = c("RP11-386I8.6:lincRNA","RP11-267L5.1:lincRNA","hsa_piR_001042","RP11-760N9.1:lincRNA","CTC-297N7.8:lincRNA","hsa-miR-487a-3p","RBFADN:lincRNA","LINC01067:lincRNA","RP11-1036E20.9:lincRNA","SNORA3:snoRNA","SNORD115-31:snoRNA","SNORD116-1:snoRNA","hsa-miR-3150b-3p","CTD-2135D7.3:lincRNA","hsa_piR_004420","hsa-miR-301a-3p","CTD-2017D11.1:lincRNA","Trp_tRNA","hsa-miR-24-3p","hsa-miR-27a-3p","hsa-let-7i-5p","hsa-miR-1301-3p","RNU2-68P:snRNA","hsa-miR-106a-5p","hsa-miR-155-5p","hsa-miR-127-5p","RP11-3B12.3:lincRNA","LINC00511:lincRNA","RNU2-26P:snRNA","RNU2-57P:snRNA","CTD-2044J15.2:lincRNA","hsa-miR-379-5p","hsa-miR-410-3p","CTD-2540B15.12:lincRNA","RNU2-28P:snRNA","RNU6-671P:snRNA","RNU2-16P:snRNA","RNU6-761P:snRNA","RNU6-1329P:snRNA","SCARNA11:snoRNA","RNU6-891P:snRNA","RNU6-21P:snRNA","hsa-miR-1283","hsa_piR_017033","RNU2-59P:snRNA","Y_RNA.320:misc_RNA","Y_RNA.348:misc_RNA","SNORD42A:snoRNA","RNU6-42P:snRNA","RNU6-1:snRNA","RP11-277P12.10:lincRNA","RNU2-37P:snRNA","RNU2-33P:snRNA","RNU2-7P:snRNA","RNU2-48P:snRNA","RP11-256I9.2:lincRNA","LINC00265:lincRNA","RP5-1070A16.1:lincRNA","hsa-miR-199a-3p|hsa-miR-199b-3p","CTD-2501M5.1:lincRNA","RP11-539G18.2:lincRNA","AC092415.1:lincRNA","RP11-346D14.1:lincRNA","RP11-212P7.2:lincRNA","Y_RNA.555:misc_RNA","Y_RNA.95:misc_RNA","Y_RNA.605:misc_RNA","Y_RNA.331:misc_RNA","Y_RNA.389:misc_RNA","hsa-miR-548d-3p","RNU4ATAC16P:snRNA","Y_RNA.328:misc_RNA","RNU6-1058P:snRNA","RNU6-1139P:snRNA","Y_RNA.513:misc_RNA","hsa_piR_014620","Y_RNA.362:misc_RNA","Ser_tRNA","Y_RNA.31:misc_RNA","hsa-miR-744-5p","RNU3P3:snoRNA")
OptiPrep.RNA = c("hsa-miR-181b-5p","hsa-miR-107","hsa-miR-181a-5p","hsa-miR-22-3p","Val_tRNA","Y_RNA.148:misc_RNA","RP11-160N1.9:lincRNA","AC010982.1:lincRNA","MIR181A1HG:lincRNA","RP11-757O6.1:lincRNA","SNORA70.14:snoRNA","SNORA70.16:snoRNA","hsa_piR_012287","RP11-238K6.2:lincRNA","RP11-363J20.2:lincRNA","hsa_piR_011015","RP11-96L7.2:lincRNA","hsa_piR_015471","hsa_piR_014547","RP11-284F21.11:lincRNA","RP11-25O10.2:lincRNA","RP11-366L20.4:lincRNA","RP11-145E17.2:lincRNA","hsa_piR_004325","AC018731.3:lincRNA","RP11-430H10.3:lincRNA","RP11-355F22.1:lincRNA","RP11-473E2.4:lincRNA","AC064834.2:lincRNA","AC019117.2:lincRNA","RP11-292F9.1:lincRNA","TTTY8:lincRNA","TTTY8B:lincRNA","RP11-347K2.2:lincRNA","RP11-760D2.1:lincRNA","RP11-503L19.1:lincRNA","hsa_piR_004137","RP11-268P4.5:lincRNA","RP4-680D5.9:lincRNA","hsa_piR_005828","RP11-91I20.2:lincRNA","RP11-15I11.3:lincRNA","hsa_piR_022275","RP11-136K7.1:lincRNA","snoU13.357:snoRNA","Y_RNA.244:misc_RNA","hsa_piR_002037","hsa_piR_007228","RP11-128B16.3:lincRNA","hsa_piR_017647","hsa_piR_009450","hsa_piR_006005","RP3-463P15.1:lincRNA","RP11-53L24.1:lincRNA","RNU6-1279P:snRNA","RP11-849N15.4:lincRNA","RP11-444A22.1:lincRNA","AC092669.2:lincRNA","RP11-309L24.9:lincRNA","RP11-223C24.1:lincRNA","RP11-300M6.1:lincRNA","RP11-129J12.1:lincRNA","RP11-486P11.1:lincRNA","AC018866.1:lincRNA","hsa_piR_011396","hsa_piR_013166","hsa_piR_020863","hsa_piR_014062","hsa_piR_000029","hsa_piR_018126","SNORD112.18:snoRNA","RP11-370P15.2:lincRNA","RNU1-57P:snRNA","CTA-481E9.4:lincRNA","RP11-554D15.4:lincRNA","RP11-1072N2.3:lincRNA","RP11-109D9.4:lincRNA","RP4-625H18.2:lincRNA","RP11-394A14.2:lincRNA","FAM157B:lincRNA")

## Read in User Parameters
option_list = list(
  make_option(c("-i", "--inputdir"), type = "character", default = "input/", help = "input folder directory [default= %default]", metavar = "character"),
  make_option(c("-b", "--backgrounddir"), type = "character", default = "background/", help = "background data folder directory [default= %default]", metavar = "character"),
  make_option(c("-m", "--miRNA"), type = "character", default = list.files("input/", pattern = "^.*miRNA_ReadsPerMillion.txt$"), help = "input miRNA file name [default= %default]", metavar = "character"),
  make_option(c("-p", "--piRNA"), type = "character", default = list.files("input/", pattern = "^.*piRNA_ReadsPerMillion.txt$"), help = "input piRNA file name [default= %default]", metavar = "character"),
  make_option(c("-t", "--tRNA"), type = "character", default = list.files("input/", pattern = "^.*tRNA_ReadsPerMillion.txt$"), help = "input tRNA file name [default= %default]", metavar = "character"),
  make_option(c("-g", "--gencode"), type = "character", default = list.files("input/", pattern = "^.*gencode_ReadsPerMillion.txt$"), help = "input gencode file name [default= %default]", metavar = "character"),
  make_option(c("-meta", "--metadata"), type = "character", default = "exceRpt_sample_descriptors.txt", help = "metadata file name [default= %default]", metavar = "character"),
  make_option(c("-r", "--RNASet"), type = "character", default = OriginalSet.RNA, help = "informative RNA set [default= %default]", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "output/", help = "output folder directory [default= %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Create folders for results
dir.create(opt$outdir, showWarnings = FALSE, recursive = FALSE)
Stage1 = paste(opt$outdir,"Stage_1/",sep = "")
Stage2 = paste(opt$outdir,"Stage_2/",sep = "")
UpdateAtlas = paste(opt$outdir,"UpdateAtlas/",sep = "")
dir.create(Stage1, showWarnings = FALSE, recursive = FALSE)
dir.create(Stage2, showWarnings = FALSE, recursive = FALSE)
dir.create(UpdateAtlas, showWarnings = FALSE, recursive = FALSE)

## Transformation functions
quantile_normalisation = function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
logistic = function(x){
  a = 1/max(x)
  1 - exp(1)^(-a*x)
}

## Global Parameters
my_palette = colorRampPalette(c("white","lightgreen","darkgreen"))(n = 29)
col_breaks = c(seq(0.0000,0.3330, length = 10), seq(0.3331,0.6660, length = 10), seq(0.6661,1.0000, length = 10))

## All exRNA Atlas profiles
print("Reading in exRNA Atlas")
Atlas.RPM = read.table(file = paste(opt$backgrounddir, "Atlas_Update.txt", sep = ""), quote = "\"", sep = "\t", check.names = FALSE, comment.char = "!", header = TRUE, row.names = 1)

## Read in metadata file
input.meta = read.table(file = paste(opt$inputdir, opt$metadata, sep = ""), header = TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)
input.meta[,colnames(input.meta)] = data.frame(apply(input.meta[colnames(input.meta)], 2, as.factor))
names(input.meta) = gsub(pattern = " ", replacement = "_", x = names(input.meta))

## Deconvolution results - 75 cargo profiles
Atlas.Results = read.table(file = paste(opt$backgrounddir, "Atlas_DeconvolutionResults.txt", sep = ""), quote = "\"", sep = "\t", check.names = FALSE, comment.char = "!", header = TRUE, row.names = 1)
Predictions.CTSubtype = c("CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT1","CT2","CT2","CT2","CT2","CT2","CT2","CT2","CT2","CT2","CT2","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3A","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3B","CT3C","CT3C","CT3C","CT3C","CT3C","CT3C","CT3C","CT3C","CT3C","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4","CT4")
colors.1 = rep("white",length(Predictions.CTSubtype))
colors.1[Predictions.CTSubtype=="CT1"] = "#2FB24B"
colors.1[Predictions.CTSubtype=="CT2"] = "#B49E31"
colors.1[Predictions.CTSubtype=="CT3A"] = "#11BCC2"
colors.1[Predictions.CTSubtype=="CT3B"] = "#353635"
colors.1[Predictions.CTSubtype=="CT3C"] = "#703B96"
colors.1[Predictions.CTSubtype=="CT4"] = "#F3746D"
colors.1 = adjustcolor(colors.1, alpha.f = 0.5)
colnames(Atlas.Results) = c(1:75)

## Read in RNA files, and grep gencode file for additional RNAs
miRNA.RPM = read.table(file = paste(opt$inputdir, opt$miRNA, sep = ""), header = TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)
piRNA.RPM = read.table(file = paste(opt$inputdir, opt$piRNA, sep = ""), header = TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)
rownames(piRNA.RPM) = gsub("\\|gb\\|.*","", rownames(piRNA.RPM))
tRNA.RPM = read.table(file = paste(opt$inputdir, opt$tRNA, sep = ""), header = TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)
rownames(tRNA.RPM) = gsub("$","_tRNA", rownames(tRNA.RPM))
gencode.RPM = read.table(file = paste(opt$inputdir, opt$gencode, sep = ""), header = TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)
YRNA.RPM = gencode.RPM[grepl("*Y_RNA*", rownames(gencode.RPM)), ]
lincRNA.RPM = gencode.RPM[grepl("*lincRNA*", rownames(gencode.RPM)), ]
snoRNA.RPM = gencode.RPM[grepl("*snoRNA*", rownames(gencode.RPM)), ]
snRNA.RPM = gencode.RPM[grepl("*snRNA*", rownames(gencode.RPM)), ]

## Merge all RNA subtypes into single file, remove all zero rows
allncRNA.RPM = rbind(miRNA.RPM, piRNA.RPM, tRNA.RPM, YRNA.RPM, lincRNA.RPM, snoRNA.RPM, snRNA.RPM)
allncRNA.RPM.nonZero = allncRNA.RPM[rowSums(allncRNA.RPM[,]) > 0,]
input.excerpt = allncRNA.RPM.nonZero
input.samples = intersect(rownames(input.meta), colnames(allncRNA.RPM.nonZero))
allncRNA.RPM.nonZero = allncRNA.RPM.nonZero[,input.samples]

## Check if >= to 40 samples are provided by the user
if(length(input.samples) < 40){
  print("Please provide greater than 40 samples")
  quit(status = 1)
}

## Merge to master Atlas profiles matrix
Atlas.input = merge(Atlas.RPM, input.excerpt, by = "row.names", all.x = TRUE)
Atlas.input[is.na(Atlas.input)] <- 0
rownames(Atlas.input) = Atlas.input[,1]
Atlas.input[,1] = NULL
write.table(Atlas.input, file = paste(UpdateAtlas, "Atlas_Update.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)

## Transform merged matrix
print("Transformation...")
input.excerpt.QN = as.data.frame(t(quantile_normalisation(t(Atlas.input))))
input.excerpt.QN.logistic = t(apply(t(input.excerpt.QN), 2, logistic))
input.excerpt.QN.logistic.nonZero = input.excerpt.QN.logistic[rowSums(input.excerpt.QN.logistic[, ]) > 0, ]
input.excerpt.QN.logistic.nonZero.max = (1/max(input.excerpt.QN.logistic.nonZero)) * input.excerpt.QN.logistic.nonZero

## Subset to user defined samples
input.deconv = input.excerpt.QN.logistic.nonZero.max[,colnames(input.excerpt)]

if(mean(input.deconv[OriginalSet.RNA,]) <= 0.05){
  print("Warning: Mean coverage of informative RNA set is low")
} else if(mean(input.deconv[OriginalSet.RNA,]) <= 0.015){
  print("Mean coverage of informative RNA set is too low")
  quit(status = 1)
}

## Run stability function and define k
print("Running Stability Function")
find.k = estimate_stability(meth_bulk_samples = input.deconv, informative_loci = opt$RNASet, possible_num_ct = c(3:6), subset_prop = 0.8, num_subsets = 3, reps_per_subset = 2, max_its = 1000, rss_diff_stop = 1e-8)
nCts.find.k = find.k$most_stable_num_ct

print(paste("Modeling k = ", nCts.find.k, sep = ""))

## Run Stage 1
print("Running Stage 1")
input.deconv.Stage1 = run_edec_stage_1(meth_bulk_samples = input.deconv, informative_loci = opt$RNASet, num_cell_types = nCts.find.k, max_its = 2000, rss_diff_stop = 1e-10)
input.deconv.Stage1.methy = input.deconv.Stage1$methylation
input.deconv.Stage1.methy.chosenProbes = input.deconv.Stage1.methy[opt$RNASet,]
input.deconv.Stage1.prop = as.data.frame(t(input.deconv.Stage1$proportions))
input.deconv.Stage1.Props = input.deconv.Stage1$proportions

## Output tables - transformed RPM values & proportions
colnames(input.deconv.Stage1.methy) = c(1:nCts.find.k)
write.table(input.deconv.Stage1.methy, file = paste(Stage1, "Stage1_Results_Expression.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(input.deconv.Stage1.prop, file = paste(Stage1 ,"Stage1_Results_Proportions.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)

## Initiate Stage 1 PDF
pdfPath = paste(paste(Stage1, "Stage1_Results_Heatmaps.pdf", sep = ""))
pdf(file = pdfPath, width = 30, height = 28)

## Correlation Heatmap
CorMatrix.input.deconv.Stage1 = cor(Atlas.Results, input.deconv.Stage1.methy.chosenProbes)
par(cex.main = 3)
heatmap.2(as.matrix(t(CorMatrix.input.deconv.Stage1)), main = "\n\nStage 1 - Profile Correlations", 
          srtCol = 360, adjCol = c(0.5,1), trace = "row", tracecol = "black", linecol = NULL, col = my_palette, breaks = col_breaks,
          cexCol = 1.5, cexRow = 7, dendrogram = 'none', Colv = FALSE, Rowv = FALSE, ColSideColors = colors.1,
          key = T, key.title = "", key.ylab = "correlation", key.par = list(cex = 1.5), density.info = "none", keysize = 1)
legend("topright", legend = unique(Predictions.CTSubtype), col = unique(colors.1), lty = 1, lwd = 8, cex = 1.8)

## Proportions Heatmap
par(cex.main = 3)
heatmap.2(as.matrix(input.deconv.Stage1.prop), main = "\nStage 1 - Proportions", 
          trace = "none", col = my_palette, breaks = col_breaks, labCol = "samples", srtCol = 360, cexCol = 4, cexRow = 7, 
          key = T, key.title = "", key.xlab = "proportion", key.par = list(cex = 1.5), density.info = "none", keysize = 1)
dev.off()

## Split samples based on metadata
input.deconv.Stage1.prop.sort = input.deconv.Stage1.prop[,order(names(input.deconv.Stage1.prop))]
rownames(input.deconv.Stage1.prop.sort) = c(1:length(rownames(input.deconv.Stage1.prop.sort)))
input.meta.sort = input.meta[order(rownames(input.meta)),]

## Check all samples are included in metadata in the same order
stopifnot(colnames(input.deconv.Stage1.prop.sort) == rownames(input.meta.sort))

## Plot proportions based on disease and biofluid
input.deconv.Stage1.prop.t = as.data.frame(t(input.deconv.Stage1.prop.sort))
input.deconv.Stage1.prop.t.meta = merge(input.deconv.Stage1.prop.t, input.meta.sort, by = 0, all.x = TRUE)
rownames(input.deconv.Stage1.prop.t.meta) = input.deconv.Stage1.prop.t.meta$Row.names
input.deconv.Stage1.prop.t.meta$Row.names = NULL
input.deconv.Stage1.prop.t.meta[is.na(input.deconv.Stage1.prop.t.meta)] = 0 
input.deconv.Stage1.prop.t.Melt = melt(input.deconv.Stage1.prop.t.meta, id.vars = colnames(input.meta))
colnames(input.deconv.Stage1.prop.t.Melt) = c(colnames(input.meta), "profile", "proportion")

theme_set(theme_grey(base_size = 15))
pdfPath = paste(paste(Stage1, "Stage1_Results_Boxplots.pdf", sep = ""))
pdf(file = pdfPath)
for(i in names(input.meta)){
  if(length(unique(input.meta[[i]])) <= 10){
    metadata = input.deconv.Stage1.prop.t.Melt[[i]]
    print(ggplot(data = input.deconv.Stage1.prop.t.Melt, aes(x = profile, y = proportion)) + geom_boxplot(aes(fill = metadata)) + ggtitle(paste("Stage 1 Proportions - ", i, sep = "")) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")))
  }
}
dev.off()

print("Running Stage 2")
input.Stage2 = input.deconv.Stage1.prop.t.meta
props.names = list()
props.matrix = list()
RPM.matrix = list()
Stage2.matrix = list()

for(i in names(input.meta)){
  meta.cols = c(table(input.meta[[i]]))
  for(j in names(meta.cols)){
    if(meta.cols[[j]] >= 20){
      props.names[[j]] = rownames(input.Stage2[which(input.Stage2[[i]] == j),])
      props.matrix[[j]] = as.matrix(input.deconv.Stage1.Props[props.names[[j]],])
      colnames(props.matrix[[j]]) = c(1:length(rownames(input.deconv.Stage1.prop.sort)))
      ##Subset RPM raw data based on metadata
      RPM.matrix[[j]] = as.matrix(allncRNA.RPM.nonZero[grep("miR|let", rownames(allncRNA.RPM.nonZero)), props.names[[j]]])
      ##Run Stage 2 on each subset
      Stage2.matrix[[j]] = run_edec_stage_2(gene_exp_bulk_samples = RPM.matrix[[j]], cell_type_props = props.matrix[[j]])
      Stage2.filename = paste(Stage2, "Stage2_",i,"_",j,"_miRNA_RPM.txt",sep = "")
      Stage2.filename = gsub(" ", "_", x = Stage2.filename)
      write.table(Stage2.matrix[[j]], file = Stage2.filename, sep = "\t", quote = FALSE, row.names = TRUE)
      print(paste("Stage 2: Complete for: ",i,"-",j, sep = ""))
    }
  }
}

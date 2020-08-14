library(ggplot2)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)

################################## Importing data ##################################


chrNum = "chr7"
START = 115000000
END = 117000000



# Importing BED files
# BED: information about peak locations
ImportBed_fn <- function(file) {
        import.bed(file, genome = "hg19")
}

pk_siNT1 <- ImportBed_fn("siNT_ER_E2_r1_SRX176856_peaks.bed")
pk_siNT2 <- ImportBed_fn("siNT_ER_E2_r2_SRX176858_peaks.bed")
pk_siNT3 <- ImportBed_fn("siNT_ER_E2_r3_SRX176860_peaks.bed")


pk_siGATA1 <- ImportBed_fn("siGATA_ER_E2_r1_SRX176857_peaks.bed")
pk_siGATA2 <- ImportBed_fn("siGATA_ER_E2_r2_SRX176859_peaks.bed")
pk_siGATA3 <- ImportBed_fn("siGATA_ER_E2_r3_SRX176861_peaks.bed")



# Importing BAM files
# BAM: alignment between read sequences and the reference genome in a compressed binary format

ImportBam_fn <- function(file, chr) {
        bam <- readGAlignments(file)
        bam <- bam[bam@seqnames == chr]
        return(bam)
}

bm_siNT1 <- ImportBam_fn("siNT_ER_E2_r1_SRX176856_sort.bam", chrNum)
bm_siNT2 <- ImportBam_fn("siNT_ER_E2_r2_SRX176858_sort.bam", chrNum)
bm_siNT3 <- ImportBam_fn("siNT_ER_E2_r3_SRX176860_sort.bam", chrNum)

bm_siGATA1 <- ImportBam_fn("siGATA_ER_E2_r1_SRX176857_sort.bam", chrNum)
bm_siGATA2 <- ImportBam_fn("siGATA_ER_E2_r2_SRX176859_sort.bam", chrNum)
bm_siGATA3 <- ImportBam_fn("siGATA_ER_E2_r3_SRX176861_sort.bam", chrNum)


################################## Peak Inspection ##################################

library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(chromstaR)


# Importing blacklisted region
BLFilePath <- "consensusBlacklist.bed"
BlackList <- import.bed(BLFilePath)
BlackListedOverlap <- findOverlaps(pk_siNT3, BlackList, type = "within")

ideo <- IdeogramTrack(chrNum, "hg19")
ax <- GenomeAxisTrack()



# a function for building a coverage track 
CT_fn <- function(cover, label) {
        DataTrack(cover,
                  window = 10000,
                  name = label, 
                  type = "h")
        
}

# a function for creating GRanges object from a BAM file
ConvertToGRanges_fn <- function(BamFile, chromo) {
        
        # Creating a GRanges object
        GR <-  readBamFileAsGRanges(BamFile,
                                    bamindex = BamFile,
                                    chromosomes = chromo)
        
        # Converting strand data
        GR@strand@values <- rep("*", 
                                times = length(GR@strand@values))
        
        return(GR)
}


# GRanges objects
GRange_NT3 <- ConvertToGRanges_fn("siNT_ER_E2_r3_SRX176860_sort.bam", 
                                  chrNum)

GRange_GATA3 <- ConvertToGRanges_fn("siGATA_ER_E2_r3_SRX176861_sort.bam", 
                                    chrNum)


# coverage tracks
CoverTrack1 <- CT_fn(GRange_NT3, "cov_NT3")
CoverTrack2 <- CT_fn(GRange_GATA3, "cov_GATA3") 



# Peak 
PeakTrack1 <- AnnotationTrack(pk_siNT1,
                              name = "Pk_NT3")
PeakTrack2 <- AnnotationTrack(pk_siGATA1,
                              name = "Pk_GATA3")


# Blacklist 
Black <- AnnotationTrack(BlackList,
                         name = "BL")

# Gene
# transcriptAnnotation = "symbol", "gene", "transcript", "exon", or "feature"
tx <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
                      chromosome = chrNum,
                      start = START,
                      end = END,
                      name = "Genes",
                      transcriptAnnotation = "gene")

plotTracks(list(ideo,
                ax,
                CoverTrack1,
                CoverTrack2,
                PeakTrack1,
                PeakTrack2, 
                Black,
                tx),
           from = START,
           to = END,  
           chromosome = chrNum,
           type = "heatmap",
           background.title = "darkgrey",
           main = "Chr7:115000000-117000000")



################################## Quality Control ##################################

library(ChIPQC)
library(BiocParallel)
register(DoparParam())
registered() 
bpparam("SerialParam")

SampleID = c("siNT1", 
             "siNT2",
             "siNT3",
             "siGATA1",
             "siGATA2",
             "siGATA3")

Treatment <- c(rep("Control", 3),
               rep("Knockdown", 3))

Replicate <- rep(1:3, times = 2)


# bam.bai (index) files have to be in the same folder
bamReads <- c("siNT_ER_E2_r1_SRX176856_sort.bam",
              "siNT_ER_E2_r2_SRX176858_sort.bam",
              "siNT_ER_E2_r3_SRX176860_sort.bam",
              "siGATA_ER_E2_r1_SRX176857_sort.bam",
              "siGATA_ER_E2_r2_SRX176859_sort.bam",
              "siGATA_ER_E2_r3_SRX176861_sort.bam")


Peaks <- c("siNT_ER_E2_r1_SRX176856_peaks.bed",
           "siNT_ER_E2_r2_SRX176858_peaks.bed",
           "siNT_ER_E2_r3_SRX176860_peaks.bed",
           "siGATA_ER_E2_r1_SRX176857_peaks.bed",
           "siGATA_ER_E2_r2_SRX176859_peaks.bed",
           "siGATA_ER_E2_r3_SRX176861_peaks.bed")


sample_meta <- data.frame(
        SampleID = SampleID,
        Tissue = "MCF7",
        Factor = NA, 
        Condition = NA,
        Treatment = Treatment,
        Replicate = Replicate,
        bamReads = bamReads,
        Peaks = Peaks,
        PeakCaller = "macs")


# creating GC report for chromosome 6
sample_exp <- ChIPQC(experiment = sample_meta, 
                     annotation = "hg19",
                     chromosomes = chrNum,
                     blacklist = BLFilePath)

ChIPQCreport(sample_exp, facetBy = "Treatment")

QC_Plot_fn <- function(object, tit) {
        object + theme_bw() +
                ggtitle(tit)
}

CoverageHisto <- QC_Plot_fn(plotCoverageHist(sample_exp, 
                                             facetBy = "Treatment"), 
                            "Coverage Histogram")

CrossCoverage <- QC_Plot_fn(plotCC(sample_exp, 
                                   facetBy = "Treatment"), 
                            "Cross Coverage")

library(gridExtra)
grid.arrange(CoverageHisto, CrossCoverage)

RelativeEnrichment <- QC_Plot_fn(plotRegi(sample_exp, 
                                          facetBy = "Treatment"), 
                                 "Relative Enrichment") + 
        theme(axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.3,
                                         size = 11),
              axis.text.y = element_text(size = 11))

PeakProfiles <- QC_Plot_fn(plotPeakProfile(sample_exp, 
                                           facetBy = "Treatment"), 
                           "Peak Profiles")

ReadsOverlappingPeaks <- QC_Plot_fn(plotRap(sample_exp, 
                                            facetBy = "Treatment"), 
                                    "Reads Overlapping Peaks") 

ReadsOverlappingBlacklist <- QC_Plot_fn(plotFribl(sample_exp, 
                                            facetBy = "Treatment"), 
                                    "Reads Overlapping Blacklist")

CorrelationHeatmap <- plotCorHeatmap(sample_exp, 
                                    attributes = "Treatment")

SamplePCA <- plotPrincomp(sample_exp, 
                    attributes = "Treatment")




################# Removing blacklisted regions



# a function removing blacklisted region
RemoveBlacklist_fn <- function(peaks) {
        
        # Find all overlaps between peaks and blacklisted regions
        BL <- findOverlaps(peaks, BlackList, type="within")
        
        # Remove all blacklisted peaks
        clean_peaks <- peaks[-from(BL)]
        return(clean_peaks)
}



# Updating peaks
pk_siNT1 <- RemoveBlacklist_fn(pk_siNT1)
pk_siNT2 <- RemoveBlacklist_fn(pk_siNT2)
pk_siNT3 <- RemoveBlacklist_fn(pk_siNT3)
pk_siGATA1 <- RemoveBlacklist_fn(pk_siGATA1)
pk_siGATA2 <- RemoveBlacklist_fn(pk_siGATA2)
pk_siGATA3 <- RemoveBlacklist_fn(pk_siGATA3)




################# Removing low quality alignments


RemoveLQ_fn <- function(bam, chr, threshold) {
        
        # Load reads with mapping qualities by requesting the "mapq" entries
        # Requires bam.bai files 
        rd <- readGAlignments(bam, 
                              param = ScanBamParam(what = "mapq"))
        rd <- rd[rd@seqnames == chr]
        
        # Identify good quality alignments
        HQ <- mcols(rd)$mapq >= threshold
        
        # Examine mapping quality distribution for high and low quality alignments
        rd_df <- data.frame(read = rd, high_quality = factor(HQ))
        
        # Remove low quality alignments
        GoodQuality <- subset(rd, HQ)
        return(GoodQuality)
}

# Updating coverage
bm_siNT1 <- RemoveLQ_fn("siNT_ER_E2_r1_SRX176856_sort.bam", chrNum, 20)
bm_siNT2 <- RemoveLQ_fn("siNT_ER_E2_r2_SRX176858_sort.bam", chrNum, 20)
bm_siNT3 <- RemoveLQ_fn("siNT_ER_E2_r3_SRX176860_sort.bam", chrNum, 20)
bm_siGATA1 <- RemoveLQ_fn("siGATA_ER_E2_r1_SRX176857_sort.bam", chrNum, 20)
bm_siGATA2 <- RemoveLQ_fn("siGATA_ER_E2_r2_SRX176859_sort.bam", chrNum, 20)
bm_siGATA3 <- RemoveLQ_fn("siGATA_ER_E2_r3_SRX176861_sort.bam", chrNum, 20)



################# Assessing Enrichment

Extend_fn <- function(bam, QCReport, Sample) {
        
        # Data loading 
        gr <- granges(bam)
        
        # average fragment length 
        len <- fragmentlength(QCReport)[Sample]
        
        # Extending reads and computing coverage
        ext <- resize(gr, width = len)
        cover <- coverage(ext)
        
        return(list(ext, cover))
}

# Extended reads
ext_bm_siNT1 <- Extend_fn(bm_siNT1, sample_exp, "siNT1")[[1]]
ext_bm_siNT2 <- Extend_fn(bm_siNT2, sample_exp, "siNT2")[[1]]
ext_bm_siNT3 <- Extend_fn(bm_siNT3, sample_exp, "siNT3")[[1]]
ext_bm_siGATA1 <- Extend_fn(bm_siGATA1, sample_exp, "siGATA1")[[1]]
ext_bm_siGATA2 <- Extend_fn(bm_siGATA2, sample_exp, "siGATA2")[[1]]
ext_bm_siGATA3 <- Extend_fn(bm_siGATA3, sample_exp, "siGATA3")[[1]]



# Creating 200bp bins
TotalBins <- tileGenome(seqinfo(ext_bm_siGATA3), 
                        tilewidth = 200, 
                        cut.last.tile.in.chrom = TRUE)

# a function finding bins 
CountBins_fn <- function(bins, bam, peak, chromo) {
        
        # Bins overlapping peaks 
        OverlapBins <- findOverlaps(bins, peak)
        OverlapBins <- bins[from(OverlapBins), ]
        OverlapBins$score <- countOverlaps(OverlapBins, bam)
        OverlapBins <- OverlapBins[OverlapBins@seqnames == chromo]
        return(OverlapBins)
}



# a function calculating coverage
PkBlBkg_coverage_fn <- function(bins, bam, peak, black, chromo) {
        
        # Getting Peak bins
        PeakBins <- CountBins_fn(bins, bam, peak, chromo)
        
        # Getting Blacklist bins 
        BlacklistBins <- CountBins_fn(bins, bam, black, chromo)
        
        # Removing Peak and Blacklist bins 
        Bkg <- subset(bins,
                      !bins %in% PeakBins & !bins %in% BlacklistBins)
        
        # Getting Background bins
        BkgBins <- CountBins_fn(bins, bam, Bkg, chromo)
        
        # Creating a coverage score table
        pk <- data.frame(category = "Peaks", score = PeakBins$score)
        bl <- data.frame(category = "blacklist", score = BlacklistBins$score)
        bkg <- data.frame(category = "background", score = BkgBins$score)
        
        TotalScores <- rbind(pk, bl, bkg)
        
        return(TotalScores)
}

samples <- c("siNT1",
             "siNT2",
             "siNT3",
             "siGATA1",
             "siGATA2",
             "siGATA3")


# coverage data for chromosome 6 by sample
cover_siNI1 <- PkBlBkg_coverage_fn(TotalBins, 
                                   ext_bm_siNT1, 
                                   pk_siNT1, 
                                   BlackList, 
                                   chrNum) %>%
        mutate(sampleID = samples[1])

cover_siNT2 <- PkBlBkg_coverage_fn(TotalBins, 
                                   ext_bm_siNT2, 
                                   pk_siNT2, 
                                   BlackList, 
                                   chrNum) %>%
        mutate(sampleID = samples[2])

cover_siNT3 <- PkBlBkg_coverage_fn(TotalBins, 
                                   ext_bm_siNT3, 
                                   pk_siNT3, 
                                   BlackList,
                                   chrNum) %>%
        mutate(sampleID = samples[3])

cover_siGATA1 <- PkBlBkg_coverage_fn(TotalBins, 
                                     ext_bm_siGATA1, 
                                     pk_siGATA1, 
                                     BlackList, 
                                     chrNum) %>%
        mutate(sampleID = samples[4])

cover_siGATA2 <- PkBlBkg_coverage_fn(TotalBins, 
                                     ext_bm_siGATA2, 
                                     pk_siGATA2, 
                                     BlackList, 
                                     chrNum) %>%
        mutate(sampleID = samples[5])

cover_siGATA3 <- PkBlBkg_coverage_fn(TotalBins, 
                                     ext_bm_siGATA3, 
                                     pk_siGATA3, 
                                     BlackList,
                                     chrNum) %>%
        mutate(sampleID = samples[6])


# coverage table across the samples
total_coverage <- rbind(cover_siNI1,
                        cover_siNT2,
                        cover_siNT3,
                        cover_siGATA1,
                        cover_siGATA2,
                        cover_siGATA3) %>% 
        mutate(sampleID = factor(sampleID, 
                                 levels = samples))

# creating a boxplot for coverage 
coverage_plot <- 
        ggplot(total_coverage, 
               aes(x = category, y = score)) + 
        geom_boxplot(fill = "grey") + 
        facet_grid(. ~ sampleID, scale = "free_y") + 
        theme_bw() +
        ylab("Coverage Score (Number of Overlaping Fragments)") + 
        xlab("Category") + 
        ggtitle("Coverage Distribution") 

# a boxplot with larger mag
coverage_plot_LargerMag <- coverage_plot + 
        ylim(0, 100)


grid.arrange(coverage_plot,
             coverage_plot_LargerMag,
             ncol = 1)

################################## Comparing differential binding ##################################  


library(DiffBind)

#################### Analysis

# Creating a dba object (counting reads in a peak set)
PeakCounts <- dba.count(sample_exp, summit = 250)

# Establishing a contrast 
PeakCounts <- dba.contrast(PeakCounts,
                           categories = DBA_TREATMENT)

# Running analysis
DBind <- dba.analyze(PeakCounts)


#################### Inspection

# PCA: investigating all peaks
dba.plotPCA(DBind, DBA_TREATMENT)

# PCA: investigating differentially bound peaks
DBind_PCA <- dba.plotPCA(DBind, 
                         DBA_TREATMENT, 
                         contrast = 1)

# Heatmap: investigating all peaks
DBind_Heatmap_All <- dba.plotHeatmap(DBind, DBA_TREATMENT, correlation = FALSE)

# Heatmap: investigating differentially bound peaks
DBind_Heatmap_Diff <- dba.plotHeatmap(DBind, 
                                      DBA_TREATMENT, 
                                      correlation = FALSE, 
                                      contrast = 1,
                                      main = "Differentially Bound Peaks")




################################## Visualization of differential binding ##################################


# MA plot
DBind_MA <- dba.plotMA(DBind)

# volcano plot
Dbind_volcano <- dba.plotVolcano(DBind) 

# Boxplot: distribution of affinity (= peak intensity)
DBind_affinity <- dba.plotBox(DBind, notch = FALSE)



################################## Interpretation ##################################

####################### Annotating peaks

# Prepping annotation (as a GRanges object)
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- genes[genes@seqnames == chrNum]

# Extracting peaks from ChIPQCexperiment object
PCalls <- peaks(sample_exp)

library(ChIPpeakAnno)

# Annotating peaks 
AnnotatedPeaks <- list()

for (i in 1:6) {
        
        AnnotatedPeaks[[i]] <- annoPeaks(PCalls[[i]], 
                                         genes, 
                                         bindingType = "startSite",
                                         bindingRegion = c(-5000, 5000))  
}


####################### Investigating Overlapping peaks


# Find overlapping peaks 
# findOverlapsOfPeaks function only accepts 4 GRanges objects maximum
OverlappingPeaks1 <- findOverlapsOfPeaks(PCalls[[1]],
                                         PCalls[[2]],
                                         PCalls[[4]],
                                         PCalls[[5]],
                                         maxgap = 50)

OverlappingPeaks2 <- findOverlapsOfPeaks(PCalls[[2]],
                                         PCalls[[3]],
                                         PCalls[[5]],
                                         PCalls[[6]],
                                         maxgap = 50)

# Getting merged peaks
OverlappingPeaks1_merged <- OverlappingPeaks1$mergedPeaks
OverlappingPeaks2_merged <- OverlappingPeaks2$mergedPeaks

# annotating overlapping peaks
OverlappingPeaks1_anno <- 
        annoPeaks(OverlappingPeaks1_merged, 
                  genes, 
                  bindingType = "startSite",
                  bindingRegion = c(-5000, 5000))  

OverlappingPeaks2_anno <- 
        annoPeaks(OverlappingPeaks2_merged, 
                  genes, 
                  bindingType = "startSite",
                  bindingRegion = c(-5000, 5000))  

# location of the peak relative to the gene
MakeTable_fn <- function(feature, replic) {
        df <- as.data.frame(prop.table(table(feature))) %>%
                mutate(Replicates = replic)
        names(df) <- c("Features", "Proportion", "Replicates")
        return(df)
}

PeakLocation <- rbind(MakeTable_fn(OverlappingPeaks1_anno$insideFeature, 
                                   "1/2/4/5"),
                      MakeTable_fn(OverlappingPeaks2_anno$insideFeature, 
                                   "2/3/5/6")) %>%
        mutate(Proportion = round(Proportion, digits = 2))


# statistical illustration of peak location relative to the gene
PeakLocationPlot <- 
        ggplot(PeakLocation, aes(x = Features,
                                 y = Proportion,
                                 fill = Replicates)) + 
        geom_bar(stat = "identity", 
                 position = "dodge") + 
        theme_bw() + 
        ggtitle("Peak Location Relative to Each Gene") + 
        xlab("Peak Location") + 
        theme(axis.text.x = element_text(size = 11),
              legend.text = element_text(size = 11))




####################### Visualizing Overlapping peaks


# Binding Site Overlaps (Venn diagram)
BindingSiteOVerlaps_Venn1 <- dba.plotVenn(DBind, 
                                         mask = c(1, 2, 4, 5),
                                         main = "Binding Site Overlaps in 1/2/4/5")

BindingSiteOVerlaps_Venn2 <- dba.plotVenn(DBind, 
                                          mask = c(2, 3, 5, 6),
                                          main = "Binding Site Overlaps in 2/3/5/6")



# Binding Site Overlaps (Upset plot) 
library(UpSetR)

CalledPeaks <- as.data.frame(DBind@DBA$called)

BindingSiteOVerlaps_upset <-
        upset(CalledPeaks, 
              sets = colnames(CalledPeaks),
              keep.order = TRUE,
              order.by = "freq",
              main.bar.color = "darkblue") 



################################## Interpreting gene lists ##################################

library(chipenrich)



################### Data Cleaning

# I need a data table for comparing peaks across the samples!!! 

#       seqnames     start       end score_NT1 score_NT2 score_NT3 score_GATA1 score_GATA2 score_GATA3
# 1:     chr6    321001    321200        13         4        27          17           7          55
# 2:     chr6    321201    321400        24         6        51          22          25          85
# 3:     chr6    321401    321600        10         3        12           3           6          11
# 4:     chr6    379001    379200        18         7        33          27          13          87
# 5:     chr6    379201    379400        24        40       130          68          62         278


# Getting binned peaks
BindingPeaks_NT1 <- CountBins_fn(TotalBins, ext_bm_siNT1, pk_siNT1, chrNum)
BindingPeaks_NT2 <- CountBins_fn(TotalBins, ext_bm_siNT2, pk_siNT1, chrNum)
BindingPeaks_NT3 <- CountBins_fn(TotalBins, ext_bm_siNT3, pk_siNT1, chrNum)

BindingPeaks_GATA1 <- CountBins_fn(TotalBins, ext_bm_siGATA1, pk_siGATA1, chrNum)
BindingPeaks_GATA2 <-CountBins_fn(TotalBins, ext_bm_siGATA2, pk_siGATA2, chrNum)
BindingPeaks_GATA3 <-CountBins_fn(TotalBins, ext_bm_siGATA3, pk_siGATA3, chrNum)


# functions for data cleaning
cleanDF_fn <- function(obj) {
        
        as.data.table(obj)[, -c("width", "strand")]
}

InnerJoin_fn <- function(p1, p2, s1, s2) {
        
        inner_join(cleanDF_fn(p1),
                   cleanDF_fn(p2),
                   by = c("seqnames", "start", "end"),
                   suffix = c(s1, s2)) 

}

InnerJoin_simpler_fn <- function(p1, p2, s1, s2) {
        
        inner_join(p1,
                   p2,
                   by = c("seqnames", "start", "end")) 
        
}



# Inner-joining iteratively 
BindingPeaks_12 <- InnerJoin_fn(BindingPeaks_NT1, 
                                BindingPeaks_NT2,
                                "_NT1", 
                                "_NT2")

BindingPeaks_34 <- InnerJoin_fn(BindingPeaks_NT3, 
                                BindingPeaks_GATA1,
                                "_NT3", 
                                "_GATA1")

BindingPeaks_56 <- InnerJoin_fn(BindingPeaks_GATA2, 
                                BindingPeaks_GATA3,
                                "_GATA2", 
                                "_GATA3")

BindingPeaks_1234 <- InnerJoin_simpler_fn(BindingPeaks_12,
                                          BindingPeaks_34)

BindingPeaks_123456 <- InnerJoin_simpler_fn(BindingPeaks_1234,
                                            BindingPeaks_56)




# Normalizing peak size (compared to read depth)
BindingPeaks_123456_norm <- BindingPeaks_123456[, 
                            c("score_NT1",
                              "score_NT2",
                              "score_NT3",
                              "score_GATA1",
                              "score_GATA2",
                              "score_GATA3") := 
                                    .(score_NT1 / sum(score_NT1) * 100000,
                                      score_NT2 / sum(score_NT2) * 100000,
                                      score_NT3 / sum(score_NT3) * 100000,
                                      score_GATA1 / sum(score_GATA1) * 100000,
                                      score_GATA2 / sum(score_GATA2) * 100000,
                                      score_GATA3 / sum(score_GATA3) * 100000)]

# confirming normalized peak counts
colSums(BindingPeaks_123456[, 4:9])


# Finding decreased or increased peaks
Peaks_NTs <- rowSums(BindingPeaks_123456_norm[, c("score_NT1", 
                                                  "score_NT2", 
                                                  "score_NT3")])
Peaks_GATAs <- rowSums(BindingPeaks_123456_norm[, c("score_GATA1", 
                                                    "score_GATA2", 
                                                    "score_GATA3")])

ChangeName_fn <- function(DataTable) {
        
        names(DataTable) <- c("chrom", "start", "end")
        return(as.data.frame(DataTable))
        
}


GATADominant_Peaks <- ChangeName_fn(BindingPeaks_123456_norm[Peaks_NTs < Peaks_GATAs, 
                                                             1:3])

NTDominant_Peaks <- ChangeName_fn(BindingPeaks_123456_norm[Peaks_NTs > Peaks_GATAs, 
                                                           1:3])



################### Enrichment Test (Inspection)

# I need a peak format shown below which has been creased by the data cleaning step
#
#     chrom    start      end
# 1    chr6   321001   321200
# 2    chr6   321201   321400
# 3    chr6   321401   321600
# 4    chr6   379001   379200



# Distribution of distances between peaks and transcription start sites
DistToTss_GATADominant_Peaks <- plot_dist_to_tss(GATADominant_Peaks, 
                                             genome = "hg19")

DistToTss_NTDominant_Peaks <- plot_dist_to_tss(NTDominant_Peaks, 
                                             genome = "hg19")


# Relationship btw gene length and presence of peaks
plot_chipenrich_spline(GATADominant_Peaks, 
                       genome = "hg19", 
                       locusdef = 'nearest_tss')

plot_chipenrich_spline(NTDominant_Peaks, 
                       genome = "hg19", 
                       locusdef = 'nearest_tss')


################### Enrichment Test (GSEA)


# GSEA 
EnrichPath_fn <- function(peak_object, geneset) {
        
        # supported_genomes(): available genomes 
        # supported_locusdefs()
        # supported_genesets(): available genesets
        # supported_methods() -> "chipenrich", "fet", "broadenrich" (check out vignette)
        chipenrich(peak_object, 
                   genome="hg19", 
                   genesets = geneset, 
                   out_name = NULL, 
                   locusdef = "nearest_tss", 
                   method = "chipenrich",
                   qc_plots = FALSE)
        
}

# Extracting result tables
GATADominant_Peaks_Enrich <- EnrichPath_fn(GATADominant_Peaks, 
                                           "GOBP")
NTDominant_Peaks_Enrich <- EnrichPath_fn(NTDominant_Peaks, 
                                         "GOBP")

# filtering significant pathways
GATADominant_Peaks_Enrich_SigResult <- GATADominant_Peaks_Enrich$results[1:5, ] %>% 
        mutate(Dominant_Peaks_in = "siGATA")

NTDominant_Peaks_Enrich_SigResult <- NTDominant_Peaks_Enrich$results[1:5, ] %>%
        mutate(Dominant_Peaks_in = "siNT")


Combined_EnrichPath <- rbind(NTDominant_Peaks_Enrich_SigResult,
                             GATADominant_Peaks_Enrich_SigResult) %>% 
        as.data.table()

Combined_EnrichPath_Table <- Combined_EnrichPath[, c("Geneset.ID",
                                                     "Description",
                                                     "P.value",
                                                     "FDR",
                                                     "Status",
                                                     "Geneset.Peak.Genes",
                                                     "Dominant_Peaks_in")][, FDR := round(FDR, 4)]
                                                                           

library(formattable)

Combined_EnrichPath_Table_disp <- 
        formattable(Combined_EnrichPath_Table)

library(pathview)

pathview(gene.data  = GATADominant_Peaks_Enrich$results,
         pathway.id = "hsa00480", 
         species    = "hsa",
         limit      = list(gene = 1, cpd = 1))

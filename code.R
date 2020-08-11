library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)

################################## Importing data ##################################


chrNum = "chr6"
START = 132000000
END = 133000000



# Importing BED files

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
library(GenomicRanges)

# Importing blacklisted region
BlackList <- import.bed("dukeExcludeRegions.bed")
BlackListedOverlap <- findOverlaps(pk_siNT1, BlackList, type = "within")

ideo <- IdeogramTrack(chrNum, "hg19")
ax <- GenomeAxisTrack()



# Coverage

CT_fn <- function(peak, label) {
        DataTrack(peak,
                  window = 10000,
                  name = "label", 
                  type = "h")
        
}
CoverTrack1 <- CT_fn(pk_siNT1, "coverage1")
CoverTrack2 <- CT_fn(pk_siGATA1, "coverage2") 


# Peak 
PeakTrack1 <- AnnotationTrack(pk_siNT1,
                              name = "Peaks1")
PeakTrack2 <- AnnotationTrack(pk_siGATA1,
                              name = "Peaks2")


# Blacklist 
Black <- AnnotationTrack(BlackList,
                         name = "Blacklist")

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
           background.title = "darkgrey")



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

sample_exp <- ChIPQC(experiment = sample_meta, 
                     annotation = "hg19",
                     chromosomes = chrNum,
                     blacklist = "dukeExcludeRegions.bed")

ChIPQCreport(sample_exp, facetBy = "Treatment")


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
TotalBins <- tileGenome(seqinfo(bm_siNT1), tilewidth = 200, cut.last.tile.in.chrom = TRUE)


# a function finding bins 
CountBins_fn <- function(bins, bam, peak) {
        
        # Bins overlapping peaks 
        OverlapBins <- findOverlaps(bins, peak)
        OverlapBins <- bins[from(OverlapBins), ]
        OverlapBins$score <- countOverlaps(OverlapBins, bam)
        return(OverlapBins)
}

# a function calculating coverage
PkBlBkg_coverage_fn <- function(bins, bam, peak, black) {
        
        # Getting Peak bins
        PeakBins <- CountBins_fn(bins, bam, peak)
        
        # Getting Blacklist bins 
        BlacklistBins <- CountBins_fn(bins, bam, black)
        
        # Removing Peak and Blacklist bins 
        Bkg <- subset(bins,
                      !bins %in% PeakBins & !bins %in% BlacklistBins)
        
        # Getting Background bins
        BkgBins <- CountBins_fn(bins, bam, Bkg)
        
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

cover_siNI1 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siNT1, pk_siNT1, BlackList) %>%
        mutate(sampleID = samples[1])
cover_siNT2 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siNT2, pk_siNT2, BlackList) %>%
        mutate(sampleID = samples[2])
cover_siNT3 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siNT3, pk_siNT3, BlackList) %>%
        mutate(sampleID = samples[3])
cover_siGATA1 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siGATA1, pk_siGATA1, BlackList) %>%
        mutate(sampleID = samples[4])
cover_siGATA2 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siGATA2, pk_siGATA2, BlackList) %>%
        mutate(sampleID = samples[5])
cover_siGATA3 <- PkBlBkg_coverage_fn(TotalBins, ext_bm_siGATA3, pk_siGATA3, BlackList) %>%
        mutate(sampleID = samples[6])

total_coverage <- rbind(cover_siNI1,
                        cover_siNT2,
                        cover_siNT3,
                        cover_siGATA1,
                        cover_siGATA2,
                        cover_siGATA3) %>% 
        mutate(sampleID = factor(sampleID, 
                                 levels = samples))


coverage_plot <- 
        ggplot(total_coverage, 
               aes(x = category, y = score)) + 
        geom_boxplot() + 
        facet_grid(. ~ sampleID, scale = "free_y") + 
        theme_bw() +
        ylab("Coverage Score (Number of Overlaping Fragments)") + 
        xlab("Category") + 
        ggtitle("Coverage Distribution") 

coverage_plot_LargerMag <- coverage_plot + 
        ylim(c(0, 300))



################################## Comparing differential binding ##################################  


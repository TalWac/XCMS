
mzMLs <- dir('C:/Users/BCDDNB6/Documents/Projects/STD MIX/convert2format', full.names = TRUE,
            recursive = TRUE)

pd <- data.frame(sample_name = sub(basename(mzMLs), pattern = ".mzML",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("SA", 5), rep("SB", 5)),
                 stringsAsFactors = FALSE) 


raw_data <- readMSData(files = mzMLs, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk") 
# filter the data set to retention time range
raw_data<- filterRt(raw_data, c(1800, 2500))

head(rtime(raw_data))
head(mz(raw_data))

mzs <- mz(raw_data)

## Split the list by file
mzs_by_file <- split(mzs, f = fromFile(raw_data))

length(mzs_by_file) 

## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(raw_data, aggregationFun = "max")
## Define colors for the two groups
group_colors <- brewer.pal(3, "Set1")[1:2]
names(group_colors) <- c("SA", "SB")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data$sample_group])

bpi_1 <- bpis[1, 1]
## Get the total ion current by file
tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current") 

# https://www.bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html

## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))



## Define the rt and m/z range of the peak area
rtr <- c(1800, 2400)
mzr <- c(334.9, 335.1)
## extract the chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group]) 

raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")

cwp <- CentWaveParam(peakwidth = c(30, 80), noise = 1000)
xdata <- findChromPeaks(raw_data, param = cwp) 

head(chromPeaks(xdata)) 


summary_fun <- function(z) {
  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
}
T <- lapply(split.data.frame(chromPeaks(xdata),
                             f = chromPeaks(xdata)[, "sample"]),
            FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(xdata))
pandoc.table(T,
             caption = paste0("Summary statistics on identified chromatographic",
                              " peaks. Shown are number of identified peaks per",
                              " sample and widths/duration of chromatographic ",
                              "peaks.")) 

plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata, border = group_colors[chr_raw$sample_group],
                    lty = 3, rt = rtr, mz = mzr) 

pander(chromPeaks(xdata, mz = mzr, rt = rtr),
       caption = paste("Identified chromatographic peaks in a selected ",
                       "m/z and retention time range."))

## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities distribution per sample")
grid(nx = NA, ny = NULL) 

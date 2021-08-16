#### Initialize ----

zap()
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(maplet)
library(dplyr)

#### Global Parameters ----

# pathway of interest
pw <- "Cysteine, methionine, SAM, taurine metabolism" 
# maximum missingness allowed for normalization
max_miss_norm <- 0.2
# maximum missingness allowed per grade
max_miss <- 0.85
# significance level
alpha <- 0.05

##### Define Analysis Functions ----

diff_analysis_tau <- function(D, outvar, name, alpha) {
  D %>%
    # Kendall's correlation for ordinal variables
    mt_stats_univ_cor(in_col = outvar,
                      method = "kendall",
                      stat_name = name,
                      samp_filter =(!is.na(!!sym(outvar)))) %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = name, method = "BH") %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name     = name,
                     x = statistic,
                     feat_filter =  p.adj < alpha,
                     colour       = p.adj < alpha) %>%
    # boxplots of significant results
    mt_plots_box_scatter(plot_type = "box",
                         stat_name = name,
                         x    = !!sym(outvar),
                         fill  = !!sym(outvar),
                         feat_filter = p.adj < 0.05,
                         feat_sort = p.value,
                         restrict_to_used_samples = T,
                         jitter = "jitter",
                         annotation   = "{sprintf('P-value: %.1e', p.value)}\nP.adj: {sprintf('%.1e', p.adj)}")
}

#### Load Data ----

# data file
file <- "data/Glioma_Metabolon_MTfriendly.xlsx"

D <- mt_reporting_heading(heading = "Load Data", lvl = 1) %>%
  # load data sheet
  mt_load_xls(file=file, sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>% 
  # load sample annotation sheet
  mt_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_id_col="SAMPLE_NAME") %>% 
  # load metabolite annotation sheet
  mt_anno_xls(file=file, sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col="name") %>%
  # convert grade to numeric
  mt_anno_mutate(anno_type = "samples", col_name = "GRADE.tau", term = case_when(GRADE=="II"~2,
                                                                                 GRADE=="III"~3,
                                                                                 GRADE=="IV"~4)) %>%
  # convert run day variable to factor for batch correction
  mt_anno_mutate(anno_type = "samples",col_name = "RUN.DAY", term = as.factor(RUN.DAY))

#### Preprocessing ----

D <- D %>%
  mt_reporting_heading(heading="Preprocessing", lvl=1) %>%
  # sample boxplots before batch correction
  mt_plots_sample_boxplot(color=RUN.DAY, title='before batch correction - RUN DAY',plot_logged=T) %>%
  mt_plots_sample_boxplot(color=Group, title='before batch correction - Group',plot_logged=T) %>% 
  # batch correction
  mt_pre_batch_median(batch_col = "RUN.DAY") %>%
  # sample boxplots after batch correction
  mt_plots_sample_boxplot(color=RUN.DAY, title='after batch correction - RUN DAY',plot_logged=T) %>%
  mt_plots_sample_boxplot(color=Group, title='after batch correction - Group',plot_logged=T) %>%
  # normalization
  mt_pre_norm_quot(feat_max = max_miss_norm) %>% 
  mt_plots_dilution_factor(in_col = "Group") %>%
  # sample boxplots after normalization
  mt_plots_sample_boxplot(color=RUN.DAY, title='after normalization - RUN DAY',plot_logged=T) %>%
  mt_plots_sample_boxplot(color=Group, title='after normalization - Group', plot_logged=T) %>%
  # log transformation
  mt_pre_trans_log()

#### Filtering ----

D <- D %>%
  mt_reporting_heading(heading = "Select Cysteine Pathway", lvl=1) %>%
  # keep only metabolites in cysteine pathway
  mt_modify_filter_features(filter = SUB.PATHWAY %in% pw) %>%
  mt_reporting_heading(heading = "Filtering", lvl=1) %>%
  # filter out metabolites with too many missing values (if any)
  mt_plots_missingness(feat_max=max_miss) %>%
  mt_pre_filter_missingness(feat_max=max_miss, group_col = "Group") %>%
  mt_plots_missingness(feat_max=max_miss)

#### Perform Differential Analysis ----

D <- D %>%
  mt_reporting_heading(heading = "Differential Analysis - metabolites", lvl=1) %>%
  # perform differential analysis
  diff_analysis_tau(outvar="GRADE.tau", name="Grade II vs. III vs. IV", alpha=alpha)

#### Write to File ----

D %>%
  # write report to html
  mt_reporting_html(file="Glioma_GradeDifferentialAnalysis_CysteinePathway.html", 
                    title="Analysis Pipeline")

D %>%  
  # write statistical results to file
  mt_write_stats(file = "DifferentialAnalysisResults.xlsx")


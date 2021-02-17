#' Non-Small Cell Lung Cancer MAF Data
#'
#' @description A pre-loaded mutation dataset from Campbell et. al (2016), downloaded from The Cancer Genome Atlas.
#'
#' @format An annotated mutation dataframe with 6 columns and 299855 rows:
#'   \describe{
#'     \item{Tumor_Sample_Barcode}{A sample id for each mutation.}
#'     \item{Hugo_Symbol}{The name of the gene location for each mutation.}
#'     \item{Variant_Classification}{The mutation type for each mutation.}
#'     \item{Chromosome}{Chromosome on which the mutation occurred.}
#'     \item{Start_Position}{Start nucleotide location for mutation.}
#'     \item{End_Position}{End nucleotide location for mutation.}
#'    }
#' @source \url{https://www.cbioportal.org/study/summary?id=nsclc_tcga_broad_2016}
#'
"nsclc_maf"

#' Non-Small Cell Lung Cancer Survival and Clinical Data
#'
#' @description A pre-loaded clinical dataset containing survival and clinical data from Campbell et. al (2016), downloaded from The Cancer Genome Atlas.
#'
#' @format An annotated mutation dataframe with 23 columns and 1144 rows. Each row
#' corresponds to a sample, and details clinical and survival information about the
#' patient from whom the sample was derived. Its columns are as follows:
#'   \describe{
#'     \item{CASE_ID}{}
#'     \item{AGE}{}
#'     \item{AGE_AT_SURGERY}{}
#'     \item{CANCER_TYPE}{}
#'     \item{CANCER_TYPE_DETAILED}{}
#'     \item{DAYS_TO_DEATH}{}
#'     \item{DAYS_TO_LAST_FOLLOWUP}{}
#'     \item{FRACTION_GENOME_ALTERED}{}
#'     \item{HISTORY_NEOADJUVANT_TRTYN}{}
#'     \item{HISTORY_OTHER_MALIGNANCY}{}
#'     \item{MUTATION_COUNT}{}
#'     \item{M_STAGE}{}
#'     \item{N_STAGE}{}
#'     \item{ONCOTREE_CODE}{}
#'     \item{OS_MONTHS}{}
#'     \item{OS_STATUS}{}
#'     \item{SAMPLE_COUNT}{}
#'     \item{SEX}{}
#'     \item{SMOKING_HISTORY}{}
#'     \item{SMOKING_PACK_YEARS}{}
#'     \item{SOMATIC_STATUS}{}
#'     \item{STAGE}{}
#'     \item{T_STAGE}{}
#'    }
#' @source \url{https://www.cbioportal.org/study/clinicalData?id=nsclc_tcga_broad_2016}
#'
"nsclc_survival"

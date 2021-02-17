#' Group and Filter Mutation Types
#'
#' @description A function to create a mutation dictionary to group and filter mutation types: this can be useful for computational practicality. It is often not practical to model
#' each distinct mutation type together, so for practicality one may group multiple classes together (e.g. all indel mutations, all nonsynonymous mutations).
#' Additionally, some mutation types may be excluded from modelling (for example, one may wish not to use synonymous mutations in the model fitting process).
#'
#' @param for_biomarker (string)
#' Specify some standard groupings of mutation types, corresponding the the coarsest groupings
#' of nonsynonymous mutations required to evaluate the biomarkers TMB and TIB. If "TMB", groups all nonsynonymous
#' mutations together, if "TIB" groups indel mutations together and all other mutations together.
#' @param include_synonymous (logical)
#' Determine whether synonymous mutations should be included in the dictionary.
#' @param maf (dataframe)
#' An annotated mutation table containing the column 'Variant_Classification', only used to check if the dictionary specified does not contain all the variant types
#' in your dataset.
#' @param dictionary (character)
#' Directly specify the dictionary, in the form of a vector of grouping values. The names of the vector should correspond to the set
#' of variant classifications of interest in the mutation annotated file (MAF).
#'
#' @return
#' A vector of characters, with values corresponding to the grouping labels for mutation types, and with names corresponding to the mutation types as they
#' will be referred to in a mutation annotated file (MAF). See examples.
#' @export
#'
#' @examples
#' # To understand the dictionary format, note that the following code
#' dictionary <- get_mutation_dictionary(for_biomarker = "TMB")
#' # is equivalent to
#' dictionary <- c(rep("NS",9), rep("S", 8))
#' names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
#' 'Splice_Site', 'Translation_Start_Site',
#' 'Nonstop_Mutation', 'In_Frame_Ins',
#' 'In_Frame_Del', 'Frame_Shift_Del',
#' 'Frame_Shift_Ins', 'Silent',
#' 'Splice_Region', '3\'Flank', '5\'Flank',
#' 'Intron', 'RNA', '3\'UTR', '5\'UTR')
#' # where the grouping levels are chosen to be "NS" and "S" for
#' # nonsynonymous and synonymous mutations respectively.
#' # the code
#' dictionary <- get_mutation_dictionary(for_biomarker = "TIB", include_synonymous = FALSE)
#' # is equivalent to
#' dictionary <- dictionary <- c(rep("NS",7), rep("I", 2))
#' names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
#'                        'Splice_Site', 'Translation_Start_Site',
#'                       'Nonstop_Mutation', 'In_Frame_Ins',
#'                       'In_Frame_Del', 'Frame_Shift_Del',
#'                       'Frame_Shift_Ins')
#' # where now "I" is used as a label to refer to indel mutations,
#' # and synonymous mutations are filtered out.

get_mutation_dictionary <- function(for_biomarker = "TIB", include_synonymous = TRUE, maf = NULL, dictionary = NULL) {

  if (is.null(dictionary)) {
    if(for_biomarker == "TMB") {
      if (include_synonymous) {
        dictionary <- c(rep("NS",9), rep("S", 8))
        names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
                               'Splice_Site', 'Translation_Start_Site',
                               'Nonstop_Mutation', 'In_Frame_Ins',
                               'In_Frame_Del', 'Frame_Shift_Del',
                               'Frame_Shift_Ins', 'Silent',
                               'Splice_Region', '3\'Flank', '5\'Flank',
                               'Intron', 'RNA', '3\'UTR', '5\'UTR')
      }

      else {
        dictionary <- c(rep("NS",9))
        names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
                               'Splice_Site', 'Translation_Start_Site',
                               'Nonstop_Mutation', 'In_Frame_Ins',
                               'In_Frame_Del', 'Frame_Shift_Del',
                               'Frame_Shift_Ins')
      }
    }

    else if (for_biomarker == "TIB") {
      if (include_synonymous) {
        dictionary <- c(rep("NS",7), rep("I", 2), rep("S", 8))
        names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
                               'Splice_Site', 'Translation_Start_Site',
                               'Nonstop_Mutation', 'In_Frame_Ins',
                               'In_Frame_Del', 'Frame_Shift_Del',
                               'Frame_Shift_Ins', 'Silent',
                               'Splice_Region', '3\'Flank', '5\'Flank',
                               'Intron', 'RNA', '3\'UTR', '5\'UTR')
      }

      else {
        dictionary <- c(rep("NS",7), rep("I", 2))
        names(dictionary) <- c('Missense_Mutation', 'Nonsense_Mutation',
                               'Splice_Site', 'Translation_Start_Site',
                               'Nonstop_Mutation', 'In_Frame_Ins',
                               'In_Frame_Del', 'Frame_Shift_Del',
                               'Frame_Shift_Ins')

      }
    }

    else {
      stop("Not a recognised biomarker")
    }
  }

  if(!is.null(maf)){
    maf_mut_types <- unique(maf$Variant_Classification)

    if (length(setdiff(maf_mut_types, names(dictionary))) > 0) {
      if(include_synonymous) {
        stop(paste0(c("The following mutation types in your MAF file are not in your dictionary: ", paste0(setdiff(maf_mut_types, names(dictionary)), collapse = ", "),
                    "."), collapse = ""))
      }

      else {
        warning(paste0(c("The following mutation types in your MAF file are not in your dictionary: ", paste0(setdiff(maf_mut_types, names(dictionary)), collapse = ", "),
                       ". This is likely because they are synonymous mutation types. If any of them aren't, check your mutation type dictionary."), collapse = ""))
      }
    }
  }

  return(dictionary)

}

#' Produce a Mutation Matrix from a MAF
#'
#' @description A function to, given a mutation annotation dataset with columns for sample barcode, gene name and mutation type, to reformulate this as a mutation matrix,
#' with rows denoting samples, columns denoting gene/mutation type combinations, and the individual entries giving the number of mutations observed. This
#' will likely be very sparse, so we save it as a sparse matrix for efficiency.
#'
#' @param maf (dataframe)
#' A table of annotated mutations containing the columns 'Tumor_Sample_Barcode', 'Hugo_Symbol', and 'Variant_Classification'.
#' @param sample_list (character)
#' Optional parameter specifying the set of samples to include in the mutation matrix.
#' @param gene_list (character)
#' Optional parameter specifying the set of genes to include in the mutation matrix.
#' @param acceptable_genes (character)
#' Optional parameter specifying a set of acceptable genes, for example those which are in an ensembl databse.
#' @param for_biomarker (character)
#' Used for defining a dictionary of mutations. See the function get_mutation_dictionary() for details.
#' @param include_synonymous (logical)
#' Optional parameter specifying whether to include synonymous mutations in the mutation matrix.
#' @param dictionary (character)
#' Optional parameter directly specifying the mutation dictionary to use. See the function get_mutation_dictionary() for details.
#' @return
#' A list with the following entries:
#' * matrix: A mutation matrix, a sparse matrix showing the number of mutations present in each sample, gene and mutation type.
#' * sample_list: A vector of characters specifying the samples included in the matrix: the rows of the mutation matrix correspond to each of these.
#' * gene_list: A vector of characters specifying the the genes included in the matrix.
#' * mut_types_list: A vector of characters specifying the mutation types (as grouped into an appropriate dictionary) to be included in the matrix.
#' * col_names: A vector of characters identifying the columns of the mutation matrix. Each entry will be comprised of two parts separated by the character
#'   '_', the first identifying the gene in question and the second identifying the mutation type. E.g. 'GENE1_NS" where 'GENE1' is an element of gene_list,
#'   and 'NS' is an element of the dictionary vector.
#'
#' @export
#'
#' @examples
#' # We use the preloaded maf file example_maf_data
#' # Now we make a mutation matrix
#' table <- get_table_from_maf(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))
#'
#' print(names(table))
#' print(table$matrix[1:10,1:10])
#' print(table$col_names[1:10])

get_table_from_maf <- function(maf, sample_list = NULL, gene_list = NULL, acceptable_genes = NULL, for_biomarker = "TIB", include_synonymous = TRUE, dictionary = NULL) {

  if (is.null(sample_list)) {
    sample_list <- unique(maf$Tumor_Sample_Barcode)
  }

  n_samples <- length(sample_list)
  sample_ids <- 1:n_samples
  names(sample_ids) <- sample_list


  if (is.null(gene_list)) {
    gene_list <- unique(maf$Hugo_Symbol)
  }


  if (!is.null(acceptable_genes)) {
    wrong_genes <- setdiff(gene_list, acceptable_genes)
    if (length(wrong_genes) > 0) {
      warning(paste("Eliminating the following", length(wrong_genes), "genes not in acceptable_genes: ", paste0(wrong_genes, collapse = ", ")))
      gene_list <- intersect(gene_list, acceptable_genes)
    }
  }


  if (!is.null(acceptable_genes)) {
    gene_list <- intersect(gene_list, acceptable_genes)
  }

  n_genes <- length(gene_list)
  gene_ids <- 1:n_genes
  names(gene_ids) <- gene_list

  dictionary <- get_mutation_dictionary(for_biomarker = for_biomarker, include_synonymous = include_synonymous, maf = maf, dictionary = dictionary)

  n_mut_types <- length(unique(dictionary))
  mut_type_ids <- 1:n_mut_types
  names(mut_type_ids) <- unique(dictionary)

  abridged_maf <- dplyr::filter(maf, maf$Tumor_Sample_Barcode %in% sample_list,
                              maf$Hugo_Symbol %in% gene_list,
                              maf$Variant_Classification %in% names(dictionary))

  maf_ids <- data.frame(sample_id = sample_ids[abridged_maf$Tumor_Sample_Barcode],
                        gene_id = gene_ids[abridged_maf$Hugo_Symbol],
                        mut_type_id = mut_type_ids[dictionary[abridged_maf$Variant_Classification]])

  mutation_vector_locations <- (maf_ids$gene_id - 1)*n_samples*n_mut_types + (maf_ids$mut_type_id - 1)*n_samples + maf_ids$sample_id
  mutation_matrix <- Matrix::sparseMatrix(i = mutation_vector_locations,
                                          j = rep(1, length(mutation_vector_locations)),
                                          x = 1,
                                          dims = c(n_samples*n_genes*n_mut_types, 1))
  dim(mutation_matrix) <- c(n_samples, n_genes*n_mut_types)
  col_names <- paste0(rep(gene_list, each = n_mut_types), "_", rep(unique(dictionary), times = n_genes))

  return(list(matrix = mutation_matrix, sample_list = sample_list, gene_list = gene_list, mut_types_list = unique(dictionary), col_names = col_names))

}

#' Produce Training, Validation and Test Matrices
#'
#' @description This function allows for i) separation of a mutation dataset into training, validation and testing components, and ii) conversion from annotated mutation format
#' to sparse mutation matrices, as described in the function get_table_from_maf().
#'
#' @param maf (dataframe)
#' A table of annotated mutations containing the columns 'Tumor_Sample_Barcode', 'Hugo_Symbol', and 'Variant_Classification'.
#' @param split (double)
#' A vector of three positive values with names 'train', 'val' and 'test'. Specifies the proportions into which to split the dataset.
#' @param sample_list sample_list (character)
#' Optional parameter specifying the set of samples to include in the mutation matrices.
#' @param gene_list (character)
#' Optional parameter specifying the set of genes to include in the mutation matrices.
#' @param acceptable_genes (character)
#' Optional parameter specifying a set of acceptable genes, for example those which are in an ensembl databse.
#' @param for_biomarker (character)
#' Used for defining a dictionary of mutations. See the function get_mutation_dictionary() for details.
#' @param include_synonymous (logical)
#' Optional parameter specifying whether to include synonymous mutations in the mutation matrices.
#' @param dictionary (character)
#' Optional parameter directly specifying the mutation dictionary to use. See the function get_mutation_dictionary() for details.
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#'
#' @return
#' A list of three items with names 'train', 'val' and 'test'. Each element will contain a sparse mutation matrix for the samples in that branch, alongside other information
#' as described as the output of the function get_table_from_maf().
#' @export
#'
#' @examples
#' tables <- get_mutation_tables(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))
#'
#' print(names(tables))
#' print(names(tables$train))

get_mutation_tables <- function(maf, split = c(train = 0.7, val = 0.15, test = 0.15), sample_list = NULL, gene_list = NULL, acceptable_genes = NULL,
                                for_biomarker = "TIB", include_synonymous = TRUE, dictionary = NULL, seed_id = 1234) {

  set.seed(seed_id)
  if (is.null(sample_list)) {
    sample_list <- unique(maf$Tumor_Sample_Barcode)
  }

  if (is.null(gene_list)) {
    gene_list <- unique(maf$Hugo_Symbol)
  }
  split <- split / sum(split)
  train_samples <- sample(sample_list, size = floor(split['train']*length(sample_list)), replace = FALSE)
  val_samples <- sample(setdiff(sample_list, train_samples), floor(split['val']*length(sample_list)), replace = FALSE)
  test_samples <- setdiff(sample_list, c(train_samples, val_samples))

  split_matrices <- purrr::map(list(train = train_samples, val = val_samples, test = test_samples), ~ get_table_from_maf(maf = maf, sample_list = ., gene_list = gene_list,
                                                                                                       acceptable_genes = acceptable_genes, for_biomarker = for_biomarker,
                                                                                                       include_synonymous = include_synonymous, dictionary = dictionary))

  return(split_matrices)


  }

#' Produce a Table of Biomarker Values from a MAF
#'
#' @description A function to recover true biomarker values from a mutation annotation file.
#'
#' @param maf (dataframe)
#' A table of annotated mutations containing the columns 'Tumor_Sample_Barcode', 'Hugo_Symbol', and 'Variant_Classification'.
#' @param biomarker (character)
#' Which biomarker needs calculating? If "TMB" or "TIB", then appropriate mutation types will be selected. Otherwise, will be interpreted as a vector of characters
#' denoting mutation types to include.
#' @param sample_list (character)
#' Vector of characters giving a list of values of Tumor_Sample_Barcode to include.
#' @param gene_list (character)
#' Vector of characters giving a list of genes to include in calculation of biomarker.
#' @param biomarker_name (character)
#' Name of biomarker. Only needed if biomarker is not "TMB" or "TIB"
#'
#' @return
#' A dataframe with two columns, 'Tumor_Sample_Barcode' and values of the biomarker specified.
#'
#' @export
#'
#' @examples print(head(get_biomarker_from_maf(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))))

get_biomarker_from_maf <- function(maf, biomarker = "TIB", sample_list = NULL, gene_list = NULL, biomarker_name = NULL) {
  if (is.null(sample_list)) {
    sample_list <- unique(maf$Tumor_Sample_Barcode)
  }

  if (is.null(gene_list)) {
    gene_list <- unique(maf$Hugo_Symbol)
  }

  if (biomarker == "TMB") {
    mut_types <- c('Missense_Mutation', 'Nonsense_Mutation',
                       'Splice_Site', 'Translation_Start_Site',
                       'Nonstop_Mutation', 'In_Frame_Ins',
                       'In_Frame_Del', 'Frame_Shift_Del',
                       'Frame_Shift_Ins')

    if (is.null(biomarker_name)) {
      biomarker_name <- "TMB"
    }
  }

  else if (biomarker == "TIB") {
    mut_types <- c('Frame_Shift_Ins', 'Frame_Shift_Del')

    if (is.null(biomarker_name)) {
      biomarker_name <- "TIB"
    }
  }

  else {
    mut_types <- biomarker
    if (is.null(biomarker_name)) {
      biomarker_name <- "biomarker"
    }
  }

  abridged_maf <- dplyr::filter(maf, maf$Tumor_Sample_Barcode %in% sample_list,
                                     maf$Hugo_Symbol %in% gene_list,
                                     maf$Variant_Classification %in% mut_types)

  table <- dplyr::ungroup(dplyr::count(dplyr::group_by(abridged_maf, abridged_maf$Tumor_Sample_Barcode)))
  colnames(table) <- c("Tumor_Sample_Barcode", biomarker_name)

  zero_samples <- setdiff(sample_list, table$Tumor_Sample_Barcode)

  if (length(zero_samples) > 0) {
    zeros_table <- data.frame(Tumor_Sample_Barcode = zero_samples, biomarker = rep(0, length(zero_samples)), stringsAsFactors = FALSE)
    names(zeros_table) <- c("Tumor_Sample_Barcode", biomarker_name)

    table <- dplyr::bind_rows(table, zeros_table)
  }

  table <- as.data.frame(table)
  rownames(table) <- table$Tumor_Sample_Barcode
  table <- table[sample_list,]

  return(table)
}

#' Get True Biomarker Values on Training, Validation and Test Sets
#'
#' A function, similar to get_mutation_tables(), but returning the true biomarker values for a training, validation and test sets.
#'
#' @param maf (dataframe)
#' A table of annotated mutations containing the columns 'Tumor_Sample_Barcode', 'Hugo_Symbol', and 'Variant_Classification'.
#' @param biomarker (character)
#' Which biomarker needs calculating? If "TMB" or "TIB", then appropriate mutation types will be selected. Otherwise, will be interpreted as a vector of characters
#' denoting mutation types to include.
#' @param sample_list (character)
#' Vector of characters giving a list of values of Tumor_Sample_Barcode to include.
#' @param gene_list (character)
#' Vector of characters giving a list of genes to include in calculation of biomarker.
#' @param biomarker_name (character)
#' Name of biomarker. Only needed if biomarker is not "TMB" or "TIB"
#' @param tables (list)
#' Optional parameter, the output of a call to get_mutation_tables(), which already has a train/val/test split.
#' @param split (numeric)
#' Optional parameter directly specifying the proportions of a train/test/val split.
#' @param seed_id (numeric)
#' Input value for the function set.seed().
#'
#' @return A list of three objects: 'train', 'val' and 'test. Each comprises a dataframe with two columns, denoting sample ID and biomarker value.
#' @export
#'
#' @examples print(head(get_biomarker_tables(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))))

get_biomarker_tables <- function(maf, biomarker = "TIB", sample_list = NULL, gene_list = NULL, biomarker_name = NULL,
                                 tables = NULL, split = c(train = 0.7, val = 0.15, test = 0.15), seed_id = 1234) {
  set.seed(seed_id)
  if (is.null(tables)) {
    if (is.null(sample_list)) {
      sample_list <- unique(maf$Tumor_Sample_Barcode)
    }

    if (is.null(gene_list)) {
      gene_list <- unique(maf$Hugo_Symbol)
    }
    split <- split / sum(split)
    train_samples <- sample(sample_list, size = floor(split['train']*length(sample_list)), replace = FALSE)
    val_samples <- sample(setdiff(sample_list, train_samples), floor(split['val']*length(sample_list)), replace = FALSE)
    test_samples <- setdiff(sample_list, c(train_samples, val_samples))
  }
  else {
    train_samples <- tables$train$sample_list
    val_samples <- tables$val$sample_list
    test_samples <- tables$test$sample_list
  }

  split_biomarker_tables <- purrr::map(list(train = train_samples, val = val_samples, test = test_samples),
                                       ~ get_biomarker_from_maf(maf = maf, biomarker = biomarker, sample_list = ., gene_list = gene_list, biomarker_name = biomarker_name))

  return(split_biomarker_tables)
}

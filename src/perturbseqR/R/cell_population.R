CellPopulation <- R6::R6Class(
  "CellPopulation",
  public = list(
    matrix = NULL,            # expression matrix: rows = cells, cols = genes
    normalized_matrix = NULL, # optional normalized matrix (same dim/rownames/colnames)
    genes = NULL,             # gene metadata (data.table or data.frame) indexed by gene id
    cells = NULL,             # cell metadata (data.table or data.frame) indexed by barcode
    source = NULL,
    guides = NULL,

    initialize = function(matrix,
                          cell_list,
                          gene_list,
                          source = "arrays",
                          normalized_matrix = NULL,
                          calculate_statistics = TRUE) {

      # matrix: matrix, dgCMatrix, or data.frame/matrix-like (rows=cells, cols=genes)
      # cell_list & gene_list: data.frame/data.table; they will be coerced to data.table
      self$matrix <- matrix
      self$normalized_matrix <- normalized_matrix
      self$source <- source

      # coerce metadata to data.table and ensure rownames available
      if (!is.data.table(cell_list)) cell_list <- as.data.table(cell_list, keep.rownames = FALSE)
      if (!is.data.table(gene_list)) gene_list <- as.data.table(gene_list, keep.rownames = FALSE)

      # ensure gene_list rows are named by gene id (use rownames or a column gene_id)
      if (!is.null(rownames(gene_list)) && any(nzchar(rownames(gene_list)))) {
        rownames(gene_list) <- rownames(gene_list)
      } else if ("gene_id" %in% colnames(gene_list)) {
        # set rownames from gene_id
        rownames(gene_list) <- as.character(gene_list$gene_id)
      } else {
        # leave as-is (user might index by position)
      }

      # compute statistics
      if (calculate_statistics) {
        message("Generating summary statistics...")
        # handle sparse/dense automatically: use Matrix::colMeans / colSds if sparse
        if (inherits(self$matrix, "dgCMatrix")) {
          gene_means <- Matrix::colMeans(self$matrix)
          # colSds is not in Matrix; approximate using apply on dense conversion if small, else compute var via crossprod
          gene_vars <- Matrix::colSums((self$matrix)^2)/nrow(self$matrix) - gene_means^2
          gene_sd <- sqrt(pmax(0, gene_vars))
        } else {
          gene_means <- colMeans(self$matrix, na.rm = TRUE)
          gene_sd <- apply(self$matrix, 2, sd, na.rm = TRUE)
        }
        gene_list$mean <- gene_means[rownames(gene_list)]
        gene_list$std <- gene_sd[rownames(gene_list)]
        gene_list$cv <- gene_list$std / gene_list$mean
        gene_list$fano <- (gene_list$std^2) / gene_list$mean
        gene_list$in_matrix <- TRUE
      }

      # assign metadata (as data.tables keyed by rownames if possible)
      self$genes <- as.data.table(gene_list, keep.rownames = "gene_id")
      if ("gene_id" %in% colnames(self$genes)) {
        setkey(self$genes, gene_id)
      }

      # cell metadata
      # if cell_list has rownames of barcodes, bring in as column cell_barcode
      if (!is.null(rownames(cell_list)) && any(nzchar(rownames(cell_list)))) {
        cell_list <- as.data.table(cell_list, keep.rownames = "cell_barcode")
      } else if ("cell_barcode" %in% colnames(cell_list)) {
        cell_list <- as.data.table(cell_list)
      } else {
        cell_list <- as.data.table(cell_list)
      }
      if ("cell_barcode" %in% colnames(cell_list)) setkey(cell_list, cell_barcode)
      self$cells <- cell_list

      if (calculate_statistics) {
        # gem_group by parsing barcode suffix after '-'
        if ("cell_barcode" %in% colnames(self$cells)) {
          self$cells[, gem_group := as.integer(sub(".*-([0-9]+)$", "\\1", cell_barcode))]
        }
        # guides: unique guide_identity among single_cell cells
        if ("single_cell" %in% colnames(self$cells) && "guide_identity" %in% colnames(self$cells)) {
          self$guides <- sort(unique(self$cells[single_cell == TRUE, guide_identity]))
        } else {
          self$guides <- character(0)
        }
      }

      message("Done.")
    },

    # HDF I/O using hdf5r
    from_hdf = function(filename) {
      # class method analog: returns a new CellPopulation
      # This is provided as an instance method to follow the asker's pattern.
      # Use: pop <- CellPopulation$new(...)$from_hdf("file.h5")
      t0 <- proc.time()[3]
      fh <- H5File$new(filename, mode = "r")
      on.exit(fh$close_all())

      message("Loading matrix...")
      # Expect datasets: /matrix (dense) or /matrix/dgCMatrix components; here we expect dense matrix
      matrix <- tryCatch({
        fh[["matrix"]][,]
      }, error = function(e) {
        # maybe stored as sparse components (i, j, x)
        if (all(c("i", "j", "x", "Dim") %in% names(fh))) {
          # reconstruct sparse
          di <- fh[["Dim"]][]
          i <- fh[["i"]][]
          j <- fh[["j"]][]
          x <- fh[["x"]][]
          # note: hdf5 storage conventions vary; this is a best-effort
          sparseMatrix(i = i + 1, j = j + 1, x = x, dims = di)
        } else stop("matrix dataset not found or unsupported format in HDF5 file")
      })

      normalized_matrix <- NULL
      if ("normalized_matrix" %in% names(fh)) {
        message("Loading normalized matrix...")
        normalized_matrix <- fh[["normalized_matrix"]][,]
      }

      message("Loading metadata...")
      gene_list <- as.data.frame(fh[["gene_list"]][,])
      cell_list <- as.data.frame(fh[["cell_list"]][,])
      fh$close_all()
      message(sprintf("Done in %.2fs.", proc.time()[3] - t0))

      CellPopulation$new(matrix, cell_list, gene_list, source = filename, normalized_matrix = normalized_matrix, calculate_statistics = FALSE)
    },

    to_hdf = function(filename, store_normalized_matrix = FALSE) {
      t0 <- proc.time()[3]
      fh <- H5File$new(filename, mode = "w")
      on.exit(fh$close_all())
      message("Writing matrix...")

      # store dense matrix directly if dense, else store sparse components
      if (inherits(self$matrix, "dgCMatrix")) {
        # store sparse representation
        sp <- self$matrix
        fh[["i"]] <- as.integer(sp@i)
        fh[["j"]] <- as.integer(rep(seq_len(ncol(sp)) - 1, diff(sp@p)))
        fh[["x"]] <- as.numeric(sp@x)
        fh[["Dim"]] <- as.integer(dim(sp))
      } else {
        fh[["matrix"]] <- as.matrix(self$matrix)
      }

      if (store_normalized_matrix && !is.null(self$normalized_matrix)) {
        message("Writing normalized matrix...")
        fh[["normalized_matrix"]] <- as.matrix(self$normalized_matrix)
      }

      message("Writing metadata...")
      fh[["cell_list"]] <- as.data.frame(.strip_cat_cols(as.data.frame(self$cells)))
      fh[["gene_list"]] <- as.data.frame(.strip_cat_cols(as.data.frame(self$genes)))

      fh$close_all()
      message(sprintf("Done in %.2fs.", proc.time()[3] - t0))
    },

    # Load from 10x CellRanger Matrix Market format
    from_file = function(directory, genome = NULL, filtered = TRUE, raw_umi_threshold = 2000) {
      # This is a convenience factory to return a new CellPopulation.
      # Use: pop <- CellPopulation$public_methods$from_file("path")
      if (is.null(genome)) genome <- "GRCh38"

      # Candidate dirs
      cand1 <- file.path(directory, "outs", "filtered_gene_bc_matrices", genome)
      cand2 <- file.path(directory, "outs", "raw_gene_bc_matrices", genome)
      cand3 <- file.path(directory, "outs", "filtered_gene_bc_matrices_mex", genome)
      cand4 <- file.path(directory, "outs", "raw_gene_bc_matrices_mex", genome)

      if (dir.exists(cand1)) {
        matrix_directory <- if (filtered) cand1 else cand2
      } else {
        matrix_directory <- if (filtered) cand3 else cand4
      }

      message(sprintf("Loading digital expression data: %s", file.path(matrix_directory, "matrix.mtx")))
      genes_path <- file.path(matrix_directory, "genes.tsv")
      gene_list <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
      colnames(gene_list) <- c("gene_id", "gene_name")

      barcodes_path <- file.path(matrix_directory, "barcodes.tsv")
      cell_barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
      colnames(cell_barcodes) <- "cell_barcode"

      # read matrix market file
      mm_path <- file.path(matrix_directory, "matrix.mtx")
      raw_mm <- Matrix::readMM(mm_path) # usually genes x barcodes or barcodes x genes depending on 10x version

      message("Densifying matrix...")
      if (filtered) {
        # CellRanger filtered typically genes x barcodes -> transpose to cells x genes
        mat <- t(raw_mm)
        # ensure column names are gene ids
        colnames(mat) <- gene_list$gene_id
        rownames(mat) <- cell_barcodes$cell_barcode
        mat <- as(mat, "dgCMatrix") # keep sparse
      } else {
        # raw matrices might be huge; filter barcodes by total UMIs
        umi_per_barcode <- as.integer(colSums(raw_mm))
        keep_idx <- which(umi_per_barcode >= raw_umi_threshold)
        message(sprintf("Filtering cell barcodes with fewer than %d UMIs...", raw_umi_threshold))
        mat <- raw_mm[, keep_idx]
        mat <- t(mat)
        colnames(mat) <- gene_list$gene_id
        rownames(mat) <- cell_barcodes$cell_barcode[keep_idx]
        mat <- as(mat, "dgCMatrix")
      }

      # gene_list: set rownames
      rownames(gene_list) <- gene_list$gene_id

      # guide identities file
      identity_filename <- if (filtered) "cell_identities.csv" else "raw_cell_identities.csv"
      message(sprintf("Loading guide identities: %s", file.path(directory, "outs", identity_filename)))
      guide_path <- file.path(directory, "outs", identity_filename)
      guide_identities <- read.csv(guide_path, stringsAsFactors = FALSE, check.names = FALSE)
      # normalize column names similar to Python code
      rename_map <- c(
        "cell BC" = "cell_barcode",
        "read count" = "guide_read_count",
        "UMI count" = "guide_UMI_count",
        "coverage" = "guide_coverage",
        "cell_BC" = "cell_barcode",
        "read_count" = "guide_read_count",
        "UMI_count" = "guide_UMI_count"
      )
      for (old in names(rename_map)) {
        if (old %in% colnames(guide_identities)) colnames(guide_identities)[colnames(guide_identities) == old] <- rename_map[[old]]
      }
      # set index
      if (!"cell_barcode" %in% colnames(guide_identities)) stop("guide_identities lacks cell_barcode column")
      rownames(guide_identities) <- guide_identities$cell_barcode

      # merge cell_barcodes and guide identities
      cb_df <- data.frame(cell_barcode = cell_barcodes$cell_barcode, stringsAsFactors = FALSE)
      if (filtered) {
        cell_list <- merge(cb_df, guide_identities, by = "cell_barcode", all.x = TRUE, sort = FALSE)
      } else {
        # keep only those in keep_idx
        cell_list <- merge(cb_df[keep_idx, , drop = FALSE], guide_identities, by = "cell_barcode", all.x = TRUE, sort = FALSE)
      }
      rownames(cell_list) <- cell_list$cell_barcode
      # compute guide_target
      if ("guide_identity" %in% colnames(cell_list)) {
        cell_list$guide_target <- vapply(cell_list$guide_identity, function(x) strsplit(as.character(x), "_")[[1]][1], FUN.VALUE = character(1))
      } else {
        cell_list$guide_target <- NA_character_
      }

      # single_cell, UMI_count
      cell_list$single_cell <- (cell_list$number_of_cells == 1) & (cell_list$good_coverage == TRUE) & (cell_list$guide_identity != "*")
      cell_list$UMI_count <- Matrix::rowSums(mat)

      # return new object
      CellPopulation$new(mat, cell_list, gene_list, source = directory)
    },

    # where
    where = function(cells = NULL,
                     genes = NULL,
                     normalized = FALSE,
                     gene_names = FALSE,
                     densify = TRUE,
                     return_query_str = FALSE,
                     dropna = FALSE,
                     ...) {

      # choose matrix
      mat <- if (normalized) self$normalized_matrix else self$matrix
      if (is.null(mat)) stop("Requested matrix not available (normalized = ", normalized, ")")

      # densify if requested and sparse
      if (densify && inherits(mat, "dgCMatrix")) {
        if (!normalized) message("(where) Densifying matrix...")
        else message("(where) Densifying normalized matrix...")
        mat <- as.matrix(mat)
      }

      # quick return if no queries:
      if (is.null(cells) && is.null(genes)) {
        if (gene_names) {
          out_mat <- mat
          colnames(out_mat) <- self$genes$gene_name[match(colnames(out_mat), self$genes$gene_id)]
          if (return_query_str) return(list(out_mat, c("", "")) ) else return(out_mat)
        } else {
          if (return_query_str) return(list(mat, c("", ""))) else return(mat)
        }
      }

      # figure out whether complex (string) queries or simple vectors
      complex_cell_indexing <- is.character(cells) && length(cells) == 1
      complex_gene_indexing <- is.character(genes) && length(genes) == 1
      cell_query_str <- ""
      gene_query_str <- ""

      # cell_index
      if (complex_cell_indexing) {
        idx <- .eval_query_index(as.data.frame(self$cells), cells, extra = list(...))
        cell_index <- rownames(self$cells)[idx]
        cell_query_str <- paste0("| ", cells, " ")
      } else if (!is.null(cells)) {
        # assume vector of barcodes
        cell_index <- as.character(cells)
        cell_query_str <- "| cells in list "
      } else {
        cell_index <- NULL
      }

      # gene_index
      if (complex_gene_indexing) {
        idxg <- .eval_query_index(as.data.frame(self$genes), genes, extra = list(...))
        gene_index <- rownames(self$genes)[idxg]
        gene_query_str <- paste0("| ", genes, " ")
      } else if (!is.null(genes)) {
        # if provided gene names (not ENSG), map to ids
        if (!is.character(genes)) stop("genes must be a character vector of gene IDs or a query string")
        test_gene <- genes[[1]]
        if (nchar(test_gene) >= 4 && substring(test_gene, 1, 4) != "ENSG") {
          # assume gene symbols; map using gene_name column
          if (!"gene_name" %in% colnames(self$genes)) stop("Mapping from gene symbols to ids requires self$genes$gene_name column")
          gene_ids <- self$genes$gene_id[match(genes, self$genes$gene_name)]
          gene_index <- gene_ids[!is.na(gene_ids)]
        } else {
          gene_index <- genes
        }
        gene_query_str <- "| genes in list "
      } else {
        gene_index <- NULL
      }

      # subset cases
      if (is.null(gene_index) && !is.null(cell_index)) {
        out_mat <- mat[cell_index, , drop = FALSE]
      } else if (!is.null(gene_index) && is.null(cell_index)) {
        out_mat <- mat[, gene_index, drop = FALSE]
      } else if (!is.null(gene_index) && !is.null(cell_index)) {
        out_mat <- mat[cell_index, gene_index, drop = FALSE]
      } else {
        out_mat <- mat
      }

      # optionally convert column names to gene symbols
      if (gene_names) {
        colnames(out_mat) <- self$genes$gene_name[match(colnames(out_mat), self$genes$gene_id)]
      }

      # drop columns with NaN/Inf if requested
      if (dropna) {
        bad <- apply(out_mat, 2, function(col) any(is.infinite(col) | is.nan(col)))
        out_mat <- out_mat[, !bad, drop = FALSE]
      }

      if (return_query_str) {
        return(list(out_mat, c(cell_query_str, gene_query_str)))
      } else {
        return(out_mat)
      }
    },

    # subpopulation
    subpopulation = function(cells = NULL, genes = NULL, normalized_matrix = NULL, ...) {
      res <- self$where(cells = cells, genes = genes, return_query_str = TRUE, ...)
      new_matrix <- res[[1]]
      query_strs <- res[[2]]
      cell_query_str <- query_strs[1]
      gene_query_str <- query_strs[2]

      # subset metadata
      new_cells <- as.data.frame(self$cells)[rownames(new_matrix), , drop = FALSE]
      new_genes <- as.data.frame(self$genes)[colnames(new_matrix), , drop = FALSE]

      new_pop <- CellPopulation$new(new_matrix, new_cells, new_genes, source = paste0(self$source, "  ||", gene_query_str, cell_query_str, " || "))

      if (identical(normalized_matrix, "inherit")) {
        message("Inheriting from parent normalized matrix...")
        if (is.null(self$normalized_matrix)) {
          new_pop$normalized_matrix <- NULL
        } else {
          new_pop$normalized_matrix <- self$normalized_matrix[rownames(new_matrix), colnames(new_matrix), drop = FALSE]
        }
      } else {
        new_pop$normalized_matrix <- normalized_matrix
      }

      # if some genes were removed (in_matrix), add them back
      if ("in_matrix" %in% colnames(self$genes) && !all(self$genes$in_matrix)) {
        missing <- as.data.frame(self$genes)[!self$genes$in_matrix, , drop = FALSE]
        combined <- rbind(missing, as.data.frame(new_pop$genes))
        # reorder to original gene order if possible
        if ("gene_id" %in% colnames(self$genes)) {
          new_pop$genes <- as.data.table(combined)[match(self$genes$gene_id, combined$gene_id), ]
        } else {
          new_pop$genes <- as.data.table(combined)
        }
      }

      return(new_pop)
    },

    # fit
    fit = function(transformer, y = NULL, cells = NULL, genes = NULL, normalized = FALSE, ...) {
      mat <- self$where(cells = cells, genes = genes, normalized = normalized, ...)
      # if transformer is R6-like with $fit
      if (is.environment(transformer) || is.list(transformer) || inherits(transformer, "R6")) {
        if (!is.null(transformer$fit)) {
          transformer$fit(mat, y)
          return(invisible(transformer))
        }
      }
      # if transformer is a function: call as transformer(mat, y)
      if (is.function(transformer)) {
        transformer(mat, y)
        return(invisible(transformer))
      }
      stop("Unsupported transformer: must have $fit or be a function")
    },

    fit_transform = function(transformer, y = NULL, cells = NULL, genes = NULL, normalized = FALSE, prefix = NULL, return_dataframe = TRUE, ...) {
      mat <- self$where(cells = cells, genes = genes, normalized = normalized, ...)
      Z <- NULL

      # R6-like transformer with $fit_transform
      if (is.environment(transformer) || is.list(transformer) || inherits(transformer, "R6")) {
        if (!is.null(transformer$fit_transform)) {
          Z <- transformer$fit_transform(mat, y)
        } else if (!is.null(transformer$fit) && !is.null(transformer$transform)) {
          transformer$fit(mat, y)
          Z <- transformer$transform(mat)
        }
      } else if (is.function(transformer)) {
        # assume function performs fit_transform
        Z <- transformer(mat, y)
      } else {
        stop("Unsupported transformer for fit_transform")
      }

      # if Z is data.frame-like, return as-is or its matrix
      if (is.data.frame(Z) || is.matrix(Z)) {
        if (return_dataframe) {
          if (is.matrix(Z)) {
            Z <- as.data.frame(Z)
            rownames(Z) <- rownames(mat)
          }
          return(Z)
        } else {
          if (is.data.frame(Z)) return(as.matrix(Z)) else return(Z)
        }
      } else {
        # e.g. numeric matrix result
        if (return_dataframe) {
          Zdf <- as.data.frame(Z, stringsAsFactors = FALSE)
          rownames(Zdf) <- rownames(mat)
          if (!is.null(prefix)) colnames(Zdf) <- paste0(prefix, seq_len(ncol(Zdf)))
          return(Zdf)
        } else {
          return(Z)
        }
      }
    },


    groupby = function(key_name, densify = TRUE, ...) {
      # returns a list of named matrices: list(key1 = mat1, key2 = mat2, ...)
      args <- list(...)
      data_normalized <- if ("normalized" %in% names(args)) args$normalized else FALSE

      # if sparse and densify requested, convert once
      sparsify <- FALSE
      if (!data_normalized && inherits(self$matrix, "dgCMatrix") && densify) {
        sparsify <- TRUE
        orig_mat <- self$matrix
        self$matrix <- as.matrix(self$matrix)
        on.exit({ if (sparsify) self$matrix <- orig_mat }, add = TRUE)
      } else if (data_normalized && inherits(self$normalized_matrix, "dgCMatrix") && densify) {
        sparsify <- TRUE
        orig_nmat <- self$normalized_matrix
        self$normalized_matrix <- as.matrix(self$normalized_matrix)
        on.exit({ if (sparsify) self$normalized_matrix <- orig_nmat }, add = TRUE)
      }

      if (!("cells" %in% names(args))) {
        keys <- sort(unique(self$cells[[key_name]]))
      } else {
        cells_query <- args$cells
        data <- self$where(cells = cells_query, densify = densify, ...)
        keys <- sort(unique(self$cells[rownames(data), get(key_name)]))
      }

      out <- list()
      for (k in keys) {
        key_barcodes <- rownames(self$cells)[self$cells[[key_name]] == k]
        # cell_str expression to be used by where: 'index in @key_barcodes' style mapping not native here
        data_k <- self$where(cells = key_barcodes, densify = densify, ...)
        out[[as.character(k)]] <- data_k
      }
      return(out)
    },

    groupby_values = function(key_name, ...) {
      gb <- self$groupby(key_name, ...)
      unname(gb)
    }
  )
)
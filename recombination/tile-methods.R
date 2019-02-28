# tile-methods.R

## TODO: validation methods?

# physical tiles for a single chromosome
.physicalTiles <- function(x, chrom, tilewidth) {
  tiles <- unlist(tileGenome(seqlengths(x@ranges)[chrom], tilewidth=tilewidth))
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  # convert these tiles to a list of overlapping row loci indices
  ranges_chrom <- x@ranges[this_chrom]
  lapply(split(tiles, start(tiles)), function(tile_i) {
         hits <- findOverlaps(ranges_chrom, tile_i)
         loci <- queryHits(hits)
         loci
  })
}

# SNP tiles for a single chromosome
.snpTiles <- function(x, chrom, nsnps) {
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  ranges_chrom <- x@ranges[this_chrom]
  nloci <- length(ranges_chrom)
  m <- floor(nloci/nsnps)
  group <- rep(seq_len(m), each=nsnps)
  # append small end group
  group <- c(group, rep(m+1, nloci - length(group)))
  stopifnot(length(group) == length(ranges_chrom))
  split(seq_len(nloci), group)
}

#' Load genetic map of four columns: marker, chrom, position, interval, cumulative
loadGeneticMap <- function(mapfile) {
  cols <- c(marker="factor", seqnames="factor", position="integer",
            interval="numeric", cumulative="numeric")
  read.delim(mapfile, colClasses=cols, col.names=names(cols))
}

#' This uses a monotonic polynomial regression (9th order by default) to smooth
#' over the genentic map markers
approxGeneticMap <- function(genmap, degree=9) {
  # check if multicore available
  ncores <- getOption("mc.cores")
  lpfun <- if (!is.null(ncores) && ncores > 1) mclapply else lapply
  gm <- genmap
  chrs <- split(gm, gm$seqnames)
  if (!all(!sapply(chrs, function(x) is.unsorted(x$position) & is.unsorted(x$cumulative))))
      stop("Genetic map unsorted; either position or cumulative distance not sorted.")
  ## create some points at origin to force it through (TODO)
  # this is hack because damn MonoPoly doesn't obey y ~ x + 0
  n <- 13
  origin <- data.frame(NA, NA, rep(0, n), NA, rep(0, n))
  colnames(origin) <- colnames(gm)
  approx <- lpfun(chrs, function(d) {
      fit <- tryCatch(MonoPoly::monpol(cumulative ~ position, algorithm="Hawkins",
                                       data=rbind(origin, d), degree=degree),
                      warning=function(x) { warning(x); NULL })
  })
  if (any(sapply(approx, is.null)))
   stop("error: fits from approxGeneticMap() not converged!")
  names(approx) <- unique(gm$seqnames)
  list(fits=approx, map=genmap)
}

#' Use Monotonic Polynomial Regression Fits to Predict Genetic Distance Between
#' Arbitary Markers
predictGeneticMap <- function(approx, ranges, as_df=TRUE) {
  # make data dataframe from range object
  data <- data.frame(seqnames=seqnames(ranges), position=start(ranges))
  tmp <- unique(data$seqnames)
  keep_chrs <- tmp %in% unique(approx$map$seqnames)
  if (!all(keep_chrs)) {
      msg <- "%d chromosomes in ProgenyArray object do not matching chromosomes in genetic map: %s"
      #warning(sprintf(msg, sum(keep_chrs), paste(tmp[!keep_chrs], sep=", ")))
  }
  # drop chroms not in gen map
  data <- data[data$seqnames %in% unique(approx$map$seqnames), ]
  data$seqnames <- droplevels(data$seqnames)
#  rang <- approx$map %>% group_by(seqnames) %>%
#              summarise(min=min(position), max=max(position))

  pred <- unlist(lapply(split(data, data$seqnames), function(d) {
      chr <- as.character(d$seqnames[1])
      fit <- approx$fits[[chr]]
      res <- unname(predict(fit, d)[, 1])
      # we need to handle edge cases, since silly MonoPoly won't
      # allow regression through the origin
      ## minval <- approx$map$cumulative[approx$map$seqnames == chr][1]
      ## res <- ifelse(d$position < rang$min[rang$seqnames == chr], minval, res)

      # now we buffered fit with (0,0), so just force any slightly neg points 
      # positive. Yes this is an ugly hack, but MonoPoly doesn't do regression
      # through origin
      res <- ifelse(res < 0, 0, res)
      res
  }))
  if (!as_df) return(pred)
  data$smoothed_cumulative <- pred
  data
}

#' Function to inspect the Montonic Polynomial Regression Smoothing of Genetic Map
#'
#' @param tiles a PhasingTiles object
#' @export 
plotGeneticMap <- function(tiles) {
  genmap <- tiles@info$genetic_map
  approx <- tiles@info$smoothed_genetic_map
  p <- ggplot(approx) 
  p <- p + facet_wrap(~seqnames, scales="free")
  p <- p + geom_point(data=genmap, aes(x=position, y=cumulative), size=1)
  p <- p + geom_line(aes(x=position, y=smoothed_cumulative), color="blue")
  if (!is.null(tiles)) {
      dt <- tiles@info$physical_tiles
      dt$y <- rowMeans(dt[, c("start", "end")])
      p <- p + geom_segment(data=dt, color="red",
                            aes(x=phys_start, xend=phys_end, y=y, yend=y))
  }
  p
}

#' Create a PhasingTiles object for tiles based on genetic distances
#'
#' @param x a ProgenyArray object
#' @param mapfile a tab-delimited file (no header) of marker, chromosome, position, cumulative genetic distance
#' @param length length of each tile in centiMorgans (actual windows may be a bit larger than this)
#'
#' @export
geneticTiles <- function(x, mapfile, length) {
    gm <- loadGeneticMap(mapfile)
    message("Approximating genetic map with monotonic polynomial regression...")
    ap <- approxGeneticMap(gm)
    message("Smoothing genetic map...")
    pgm <- predictGeneticMap(ap, x@ranges)
    chrs <- split(pgm, pgm$seqnames)
    # create tiles for chroms
    rngs <- x@ranges
    message("Creating tiles...")
    tmp <- lapply(names(chrs), function(nam) {
        chr <- chrs[[nam]]
        # total gen dist
        dist <- max(chr$smoothed_cumulative)
        tiles <- cut(chr$smoothed_cumulative, breaks=floor(dist/length), dig.lab=10)
        chr$tiles <- tiles
        # turn these tiles groups into loci groups
        # by generating sequence of length of num SNPs
        tile_ints <- split(seq_along(rngs[seqnames(rngs) == nam]), as.integer(tiles))
        # parse the cut labels to get the start/end positions
        bins <- extractCutLabels(tiles)
        bins$tiles <- tiles
        bins$seqnames <- nam
        # find min and max physical positions for these genetic tiles
#        pos <- chr %>% group_by(tiles) %>%
#                   summarise(start=min(position), end=max(position))
        # match everything up, to add as columns to bins
        i <- match(bins$tiles, pos$tiles)
        bins$phys_start <- pos$start[i]
        bins$phys_end <- pos$end[i]
        list(bins, tile_ints)
    })
    tiles <- lapply(tmp, `[[`, i=2)
    names(tiles) <- names(chrs)
    bins <- do.call(rbind, lapply(tmp, `[[`, i=1))
    # calculate actual lengths
    mn_width <- mean(apply(bins[, 1:2], 1,
                           function(x) abs(diff(x))))
    info <- list(physical_tiles=bins, mapfile=mapfile, approx=ap,
                 genetic_map=gm, smoothed_genetic_map=pgm, width=mn_width)
    new("PhasingTiles", type="genetic", width=as.integer(length), tiles=tiles, info=info)
}

#' Create a PhasingTiles object for tiles based on physical distances
#' @param x a ProgenyArray object
#' @param width width in bases of each tile (except last)
#'
#' @export
physicalTiles <- function(x, width) {
  nchroms <- nlevels(seqnames(x@ranges))
  tiles <- lapply(seqlevels(x@ranges),
                  function(chrom) .physicalTiles(x, chrom, width))
  names(tiles) <- seqlevels(x@ranges)
  new("PhasingTiles", type="physical", width=as.integer(width), tiles=tiles)
}

#' Create a PhasingTiles object for tiles based on fixed number of SNPs per tile
#' @param x a ProgenyArray object
#' @param nsnps number of SNPs per tile
#'
#' @export
snpTiles <- function(x, nsnps) {
  nchroms <- nlevels(seqnames(x@ranges))
  tiles <- lapply(seqlevels(x@ranges),
                  function(chrom) .snpTiles(x, chrom, nsnps))
  names(tiles) <- seqlevels(x@ranges)
  new("PhasingTiles", type="snp", width=as.integer(nsnps), tiles=tiles)
}

#' Pretty-print a PhasingTiles object
#'
#' @param object a PhasingTiles object
#'
#' @export
#setMethod("show",
#          c(object="PhasingTiles"),
#          function(object) {
#              cat(sprintf("PhasingTiles (%s)\n", object@type))
#              cat(sprintf(" number of chromosomes: %d\n", length(object@tiles)))
#              if (object@type != "genetic") {
#                  cat(sprintf(" width: %d\n", object@width))
#              } else {
#                  cat(sprintf(" width: %d cM (specified), %0.3f cM (actual) \n", object@width, object@info$width))
#              }
#          })

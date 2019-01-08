#!/usr/bin/env Rscript

# maskaln.r
#
# Uses the Biostrings package to mask an alignment, both column- and row-wise.
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

SCRIPT_VERSION = "0.9.2"

# Arguments for testing:
# opt = list(options = list(min_blockwidth = 1, min_gap_fraction = 0.9, max_prop_gaps = 0.5), args = c('maskaln.00.alnfaa', 'maskaln.00.out'))
# Get arguments
option_list = list(
  make_option(
    '--max_prop_gaps', type = 'double', default = 0.5, 
    help = 'Maximum proportion of gap characters in a sequence to delete it, default: %default. The lower, the stricter.'
  ),
  make_option(
    '--min_blockwidth', type = 'integer', default = 1,
    help = 'Minimum widht of a block gaps in a column to mask the position, default: %default.'
  ),
  make_option(
    '--min_gap_fraction', type = 'double', default = 0.9, 
    help = 'Minimum proportion of gaps in a column to mask the position, default: %default. Set to 1.0 to not mask anything.'
  ),
  make_option(
    '--protect_taxa', type = 'character', default = '',
    help = 'Name of file containing taxa to protect, i.e. not filter out. One taxon per line.'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] input.alnfaa output.alnfaa\n\n\tThe input alignment will be masked and output as fasta.", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
  quit('no')
}

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

if ( opt$options$protect_taxa != '' ) {
  protect <- read_tsv(opt$options$protect_taxa, col_names = c('taxon'), col_types = cols(taxon = col_character()))
} else {
  protect <- tibble(taxon = character())
}

logmsg(sprintf("Reading alignment in %s", opt$args[1]))
seqs <- readAAMultipleAlignment(opt$args[1]) %>%
  maskGaps(
    min.fraction = opt$options$min_gap_fraction, 
    min.block.width = opt$options$min_blockwidth
  )

# Mask rows of the alignment with mostly gaps, i.e. allow only sequences with at most 20% gaps.
st <- data.frame(seqname = rownames(seqs), seq = as.character(seqs), stringsAsFactors = FALSE) %>%
  filter(
    seqname %in% protect$taxon |
    str_count(seq, fixed('-'))/str_length(seq) < opt$options$max_prop_gaps
  ) %>%
  tibble::column_to_rownames('seqname') %>%
  as.matrix()

if ( dim(st)[1] > 0 ) {
  ss <- as(st, 'AAStringSet')
  names(ss) <- rownames(st)

  logmsg(sprintf("Writing %d sequences of %d length to %s", nrow(st), str_length(st[1,1]), opt$args[2]))
  writeXStringSet(ss, opt$args[2])
} else {
  write("Nothing left after masking", stderr())
}

logmsg("Done")

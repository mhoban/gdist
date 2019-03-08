#!/usr/bin/env littler

options(warn=-1)

lt <- function(mat)
{
  return(mat[lower.tri(mat)])
}

jc <- function(mat,lower=F)
{
  if (lower) {
    mat <- lt(mat)
  }
  return((-3/4)*log(1-(4*mat/3)))
}

suppressMessages(library(docopt))
# check the options and such before we load the packages that take time
doc <- 
"Usage:
  gdist.r (--field1=<cat1> --field2=<cat2>) [--sequence-id=<seq> --field-name=<field>] <alignment> <datafile>

Options:
  -1 --field1=<cat1>  First group category
  -2 --field2=<cat2>  Second group category
  -s --sequence-id=<seq>  Column index/name for sequence ids [default: 1]
  -f --field-name=<field>  Column index/name for group categories [default: 2]"

opt <- docopt(doc)

# opt <- list("sequence-id"="1","field-name"="2",field1="lineopunctatus",field2="quagga",alignment="trees/cirripectes_id.fasta",datafile="trees/clade.tsv")

# make sure alignment and data files exist
if (!file.exists(opt[['alignment']])) {
  stop("Alignment file doesn't exist")
}

if (!file.exists(opt[['datafile']])) {
  stop("Data file doesn't exist")
}

# make sure our field and categories are treated correctly (integer or string)
seq.id <- ifelse(!is.na(as.integer(opt[['sequence-id']])),as.integer(opt[['sequence-id']]),opt[['sequence-id']])
field <- ifelse(!is.na(as.integer(opt[['field-name']])),as.integer(opt[['field-name']]),opt[['field-name']])
cat1 <- opt[['field1']]
cat2 <- opt[['field2']]

cat("Loading packages (may take a minute)...\n")
suppressMessages(library(tidyverse))
suppressMessages(library(pegas))

# load the alignment
alignment <- read.dna(opt[['alignment']],format="fasta")
# load the data table
samples <- suppressMessages(read_tsv(opt[['datafile']]))

# get column types
col.types <- samples %>% summarize_all(class)

# make sure the field exists in the data file and is the correct type
if (is.character(field) & !(field %in% names(samples))) {
 stop(paste(field,"not found in data file"))
} 
if (col.types[[field]] != "factor" & col.types[[field]] != "character") {
  stop(paste("Field",field,"must be character or factor"))
}

samples <- samples %>% dplyr::filter(.[[field]] %in% c(cat1,cat2))
samples <- samples %>% dplyr::filter(.[[seq.id]] %in% labels(alignment))
alignment <- alignment[labels(alignment) %in% dplyr::pull(samples,seq.id),]

groupone <- samples %>% dplyr::filter(.[[field]] == cat1) %>% pull(seq.id)
grouptwo <- samples %>% dplyr::filter(.[[field]] == cat2) %>% pull(seq.id)
rawdist <- dist.dna(alignment,model="raw",variance=T)

raw.between <- as.matrix(rawdist)[groupone,grouptwo]
raw.within.one <- as.matrix(rawdist)[groupone,groupone]
raw.within.two <- as.matrix(rawdist)[grouptwo,grouptwo]

pi.k <- mean(lt(raw.within.one))
pi.n <- mean(lt(raw.within.two))
pi.k.jc <- mean(jc(raw.within.one,lower=T))
pi.n.jc <- mean(jc(raw.within.two,lower=T))

dxy <- mean(raw.between)
dxy_jc <- mean(jc(raw.between))

da <- dxy - mean(pi.k,pi.n)
da_jc <- dxy_jc - mean(pi.k.jc,pi.n.jc)

cat("Genetic distances between",cat1,"and",cat2,"\n")
cat("Dxy:",dxy,"\t","Dxy (JC):",dxy_jc,"\n")
cat("Da:",da,"\t","Da (JC)",da_jc,"\n")



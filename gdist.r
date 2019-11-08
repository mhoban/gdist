#!/usr/bin/env littler

# the methods of doing this stuff I think I swiped from the program called satDNA

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

# save(opt,file="opt.data")
# stop()
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
suppressPackageStartupMessages(
  {
    library(dplyr)
    library(readit)
    library(readr)
    library(ape)
  }
)

# load the alignment
alignment <- read.dna(opt[['alignment']],format="fasta")


# try to load the data table using various methods
tryCatch(
  {
     samples <- suppressMessages(readit(opt[['datafile']]))
  },
  error = function(problem) {
    tryCatch(
      {
         samples <<- suppressMessages(read_tsv(opt[['datafile']]))
        if (ncol(samples) == 1) stop("single-column tab-sep")   
      },
      error = function(problem2) {
        tryCatch(
          {
             samples <<- suppressMessages(read_csv(opt[['datafile']]))   
            if (ncol(samples) == 1) stop("single-column comma-sep")
          },
          error = function(problem3) {
            #final failure condition
             samples <<- NA
          }
        )
      }
    )
  },
  finally = {
    if (is.na(samples)) {
      stop("Couldn't load data file")
    }
  }
)

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

rawdist <- dist.dna(alignment,model="raw",variance=T,as.matrix=T)
jcdist  <- dist.dna(alignment,model="JC69",variance=T,as.matrix=T)

dxy <- mean(rawdist[groupone,grouptwo])
dxy_jc <- mean(jcdist[groupone,grouptwo])

da <- dxy - mean(mean(lt(rawdist[groupone,groupone])),mean(lt(rawdist[grouptwo,grouptwo])))
da_jc <- dxy_jc - mean(mean(lt(jcdist[groupone,groupone])),mean(lt(jcdist[grouptwo,grouptwo])))

cat("Genetic distances between",cat1,"and",cat2,"\n")
cat("Dxy:",dxy,"\t","Dxy (JC):",dxy_jc,"\n")
cat("Da:",da,"\t","Da (JC):",da_jc,"\n")




# MEEM suite and R
## MEME Suite Installation on WSL2 (e.g. Ubuntu)

<https://meme-suite.org>

### Verify WSL2 configuration
```
wsl -l -v
# In case of old WSL1 - convert your Ubuntu distribution to WSL2
wsl --set-version Ubuntu-xx.yy 2
# start WSL2 instance
wsl -d Ubuntu-xx.yy
````
### Install required system dependencies
```
sudo apt update
sudo apt install -y \
    build-essential wget perl python3 \
    zlib1g-dev libxml2-dev libxslt1-dev \
    libgd-dev libexpat1-dev libpcre3-dev \
    libssl-dev imagemagick ghostscript \
    libxml-parser-perl
```

### Download and extract MEME Suite (e.g. 5.5.9
```)
wget https://meme-suite.org/meme/meme-software/5.5.9/meme-5.5.9.tar.gz --no-check-certificate

tar -xzf meme-5.5.9.tar.gz
cd meme-5.5.9
```
### Configure the build
```
# Install everything into /usr/local/meme instead of the system default
# Build MEME with its own internal copy of libxml2 instead of relying on the system version, with its own internal libxslt library (for XML/XSLT transformations), with GD so it can generate motif logos and images
./configure --prefix=/usr/local/meme \
    --enable-build-libxml2 \
    --enable-build-libxslt \
    --with-gd
```
### Compile, test, and install
```
make -j4 # run up to 4 compilation jobs at the same time
make test
make install
````
### Check that MEME Suite is installed correctly
```
wsl bash -lc 'which meme'
wsl bash -lc 'meme -version'
wsl bash -lc 'echo $PATH'
# you should see /usr/local/meme/bin in your PATH.
```
### Confirm MEME tools are available
```
which meme fimo tomtom streme dreme
```
### Optional - copy binaries to the current directory
```
cp /usr/local/meme/bin/* .
chmod +x *
```

## in R
- WSL can read your Windows files via /mnt/c/...
- R can call WSL directly

### Confirm that R sees MEME
```
wsl = function(cmd) {
  system(sprintf('wsl --cd ~ bash -lc "%s"', cmd))
}
wsl("meme -version")
```
### Normalise names
```
to_wsl <- function(path) {
  path = normalizePath(path, winslash = "/", mustWork = TRUE)
  sub("^([A-Za-z]):", "/mnt/\\L\\1", path, perl = TRUE)
}

pos_file = file.path(
  out.dir,
  paste0(gsub('\\..*', '', my.file.1),
         '_positives_', upstream.width, 'nt_upstream-regions.fasta')
)

neg_file = file.path(
  out.dir,
  paste0(gsub('\\..*', '', my.file.1),
         '_negatives_', upstream.width, 'nt_upstream-regions.fasta')
)
pos_wsl = to_wsl(pos_file)
neg_wsl = to_wsl(neg_file)
out_wsl = to_wsl(paste0(out.dir, '/MEME'))
```
### Run MEME
```
run_meme_wsl <- function(pos, neg, outdir,
                         minw = 6, maxw = 24, nmotifs = 3, mod = "zoops") {

  inner = sprintf(
    'meme %s -neg %s -objfun de -dna -mod %s -minw %d -maxw %d -nmotifs %d -oc %s',
    shQuote(pos),
    shQuote(neg),
    mod,
    minw,
    maxw,
    nmotifs,
    shQuote(outdir)
  )

  cmd = sprintf('"%s"', inner)

  cat("CMD SENT TO WSL:\n", cmd, "\n\n")

  system2("wsl", c("bash", "-lc", cmd))
}
run_meme_wsl(pos_wsl, neg_wsl, out_wsl)
```
### Run STREME
```
out_wsl = to_wsl(file.path(out.dir, 'STREME'))

run_streme_wsl <- function(pos, neg, outdir,
                           minw = 6, maxw = 24) {

  inner = sprintf(
    'streme --p %s --n %s --oc %s --minw %d --maxw %d --dna',
    shQuote(pos),
    shQuote(neg),
    shQuote(outdir),
    minw,
    maxw
  )

  cmd = sprintf('"%s"', inner)

  cat("CMD SENT TO WSL:\n", cmd, "\n\n")

  system2("wsl", c("bash", "-lc", cmd))
}

run_streme_wsl(
  pos = pos_wsl,
  neg = neg_wsl,
  outdir = out_wsl,
  minw = 6,
  maxw = 24
)
```

# Motif scan with TFBSTools and Jaspar in R
### packages
```
library(TFBSTools)
library(JASPAR2022)
library(universalmotif)
library(Biostrings)
library(pbapply)
library(ggseqlogo)
```
### Ath TF Motifs (JASPAR2022)
```
# Load Arabidopsis PFMs
opts = list(species = 3702)
pfms = TFBSTools::getMatrixSet(JASPAR2022, opts)
```
### TF-name lookup
```
tf_lookup <- sapply(names(pfms), function(id) {
  nm <- pfms[[id]]@tags$remap_tf_name
  if (is.null(nm) || nm == "") id else nm
})
```
### Promoter positive set
```
pos_file = file.path(
  out.dir,
  paste0(gsub('\\..*', '', my.file.1),
         '_positives_', upstream.width, 'nt_upstream-regions.fasta')
)
pos = readDNAStringSet(pos_file)
```
### Filter out low-info/low-qual motifs
```
IC = sapply(pfms, function(pfm) {
  mat = pfm@profileMatrix
  mat = sweep(mat, 2, colSums(mat), "/")  # normalize to probabilities

  colIC = apply(mat, 2, function(col) {
    2 - (-sum(col * log2(col + 1e-6)))    # 2 bits max for DNA
  })

  sum(colIC)
})

pfms_filtered = pfms[IC > 6] # change as needed
```
### Build Log-Odds PWMs for scanning
```
pwms_filtered = lapply(pfms_filtered, function(pfm) {
  u = convert_motifs(pfm, class = "universalmotif")
  u_pwm = convert_type(u, type = "PWM", pseudocount = 0.1)
  PWMatrix(
    ID = u_pwm@name,
    name = u_pwm@name,
    profileMatrix = u_pwm@motif,  # log-odds matrix
    tags = list()
  )
})

names(pwms_filtered) = names(pfms_filtered)
```
### Scan promoters with high-qual motifs
```
# adjust min.score as needed
all_hits = pblapply(pwms_filtered, function(pwm) {
  searchSeq(pwm, pos, min.score = "95%", strand = "*")
})

names(all_hits) = names(pwms_filtered)
```
### Summarise hit counts per motif
```
hit_counts = sapply(all_hits, function(h) sum(lengths(h)))
hit_counts_sorted = sort(hit_counts, decreasing = TRUE)
hit_counts_sorted[1:10]
# readjust min.score as needed 
```
### Plot top motifs logos
```
top = names(hit_counts_sorted)[1:20]

for (m in top) {
  
  print(m)
  print(as.data.frame(t(unlist(pfms[[m]]@tags))))
  cat('\n')

  
  mat = pfms_filtered[[m]]@profileMatrix
  mat = sweep(mat, 2, colSums(mat), "/")  # convert to probabilities

  print(
    ggseqlogo(mat, method = "prob") +
      ggtitle(tf_lookup[m])   # pretty name
  )
}
```
### Extract all hits
```
extract_hits = function(all_hits, seq_names) {
  nested = pblapply(seq_along(all_hits), function(m) {
    motif = names(all_hits)[m]
    hits_motif = all_hits[[m]]

    lapply(seq_along(hits_motif), function(i) {
      if (length(hits_motif[[i]]) > 0) {
        df = as.data.frame(hits_motif[[i]])
        df$motif_id = motif
        df$motif = tf_lookup[motif]   # pretty name
        df$sequence = seq_names[i]
        return(df)
      } else {
        return(NULL)
      }
    })
  })
  do.call(rbind, unlist(nested, recursive = FALSE))
}

hit_table = extract_hits(all_hits, names(pos))

# Collect all possible tag names across all motifs
all_tag_names = unique(unlist(lapply(pfms, function(x) names(x@tags))))
all_tag_names
# Build metadata table
motif_metadata = do.call(
  rbind,
  lapply(names(pfms), function(id) {
    tags = pfms[[id]]@tags
    # Create row with all tag names
    row = setNames(as.list(rep(NA, length(all_tag_names))), all_tag_names)
    # Fill in existing tags, collapsing multi-value fields
    for (tag in names(tags)) {
      value = tags[[tag]]
      # collapse vectors into a single string
      if (length(value) > 1) {
        value = paste(value, collapse = "; ")
      }
      row[[tag]] = value
    }
    # To data.frame
    df = as.data.frame(row, stringsAsFactors = FALSE)
    df$motif_id = id
    df
  })
)

rownames(motif_metadata) = NULL
hit_table = merge(hit_table, motif_metadata, by = "motif_id", all.x = TRUE, all.y = FALSE)
```
### Write hits to file
```
file.output.name = paste0(gsub('\\..*', '', my.file.1), '_subset_promoter_', upstream.width, '_motifs-in-R_TFBSTools_JASPAR2022.xlsx')
openxlsx::write.xlsx(hit_table, 
                     file = file.path(out.dir, file.output.name),
                     asTable = TRUE)
```

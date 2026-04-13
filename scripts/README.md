
# MEEM suite and R
## MEME Suite Installation on WSL2 (e.g. Ubuntu)

<https://meme-suite.org>

### Verify WSL2 configuration
```
wsl -l -v
# In case of of WSL1 convert your Ubuntu distribution to WSL2
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
./configure --prefix=/usr/local/meme \
    --enable-build-libxml2 \
    --enable-build-libxslt \
    --with-gd
```
### Compile, test, and install
```
make -j4
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

### confirm that R sees MEME
```
wsl = function(cmd) {
  system(sprintf('wsl --cd ~ bash -lc "%s"', cmd))
}
wsl("meme -version")
```
### normalise names
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

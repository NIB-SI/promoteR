
# MEME Suite Installation on WSL2 (e.g. Ubuntu 20.04)

<https://meme-suite.org>

## Verify WSL2 configuration
```
wsl -l -v
# In case of of WSL1 convert your Ubuntu distribution to WSL2
wsl --set-version Ubuntu-20.04 2
# start WSL2 instance
wsl -d Ubuntu-20.04
````
## Install required system dependencies
```
sudo apt update
sudo apt install -y \
    build-essential wget perl python3 \
    zlib1g-dev libxml2-dev libxslt1-dev \
    libgd-dev libexpat1-dev libpcre3-dev \
    libssl-dev imagemagick ghostscript \
    libxml-parser-perl
```

## Download and extract MEME Suite (e.g. 5.5.9
```)
wget https://meme-suite.org/meme/meme-software/5.5.9/meme-5.5.9.tar.gz --no-check-certificate

tar -xzf meme-5.5.9.tar.gz
cd meme-5.5.9
```
## Configure the build
```
./configure --prefix=/usr/local/meme \
    --enable-build-libxml2 \
    --enable-build-libxslt \
    --with-gd
```
## Compile, test, and install
```
make -j4
make test
make install
````
## Check that MEME Suite is installed correctly
```
wsl bash -lc 'which meme'
wsl bash -lc 'meme -version'
wsl bash -lc 'echo $PATH'
# you should see /usr/local/meme/bin in your PATH.
```
## Confirm MEME tools are available
```
which meme fimo tomtom streme dreme
```
## Optional - copy binaries to the current directory
```
cp /usr/local/meme/bin/* .
chmod +x *
```


#!/bin/bash

atRoot=$False

if [ $(basename "`pwd`") != "scripts" ]; then
    cd scripts
    atRoot=$True  
fi

cd ..

if [ ! -f finalproject.tar.gz ]; then
    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1bXAHETXs5_UlzRMJEQr3QhdYygPe-wmv' -O finalproject.tar.gz 
fi

if [ -d data ]; then
    rm -rf data
fi

tar -xvf finalproject.tar.gz

cd data

if [ -d htseq.counts ];then
    rm -rf htseq.counts
fi

mkdir htseq.counts

cd htseq.counts

tar xvf ../gdc_download_20211110_070038.465477.tar

# for annotations in $(ls **/annotations.txt)
# do
#     echo "\n" >> $annotations
# done

cd ..

if [ -d clinical ];then
    rm -rf clinical
fi

mkdir clinical

cd clinical

tar -xvf ../clinical.cart.2021-11-09.tar.gz 

cd ..

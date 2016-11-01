#!/bin/bash -e


for file in $(find html -type f)
do
    curl --ftp-create-dirs -T $file -u "$FTP_USER:$FTP_PASSWORD" $FTP_TARGET/$file
done; # file


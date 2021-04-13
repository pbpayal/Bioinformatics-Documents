# will give each sig
find ./* -type f|xargs -n1 md5  >/tmp/output.sigs
# gives number of files
find ./* -type f |wc
# give lines words chars in a directory
ls -R |wc 



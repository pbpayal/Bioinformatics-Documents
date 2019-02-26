## This script creates a index of your genome for bwa
## after -p chr14 denotes the name of your outfiles. e.g. chr14.ann and chr.bwt
## chr14.fa is your genome
## The & symbol sets the task to the background

bwa_dir=/u/local/apps/bwa/0.6.2/
$bwa_dir/bwa index -a bwtsw -p chr14.fa chr14.fa &


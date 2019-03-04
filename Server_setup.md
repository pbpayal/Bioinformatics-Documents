### Centos
First check you have root priviliges

* Then install wget, but in most cases use **sudo yum install**
* Check python version, mostly comes with python2.7
* To install python3 - https://linuxize.com/post/how-to-install-python-3-on-centos-7/
  * sudo yum install centos-release-scl
  * sudo yum install rh-python36
    * After typing python --version , if it outputs the latest 3.7 its fine or else type **scl enable rh-python36 bash**
      

**Install Anaconda** - https://linuxize.com/post/how-to-install-anaconda-on-centos-7/
* conda config --add channels bioconda
* conda config --add channels conda-forge
* Add at least 2 environmnets python2 and python3

**Install R** - https://linuxize.com/post/how-to-install-r-on-centos-7/
* sudo yum install epel-release
* sudo yum install R
* R --version

** Setting up WGS/WES pipeline requiremnets
* Conda Install the following tools-
  * Bwa-mem
  * GATK
  * Picard tools
  * Samtools
  ** Got an error "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory" 
  **conda install -c bioconda r-rcurl openssl=1.0**
  * Platypus
  * Snpeff and Snpsift
  * Gemini - flexible framework for exploring genetic variation, popultaion data <https://gemini.readthedocs.io/en/latest/index.html#>
   For this however craete a sperate environmnet, like **conda create -n gemeini**, because its downgrading a lot of other dpeendant packages, which might cause error while running other packages!
   
   

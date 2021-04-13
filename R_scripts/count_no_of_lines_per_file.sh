cd /Users/pbanerjee/Documents/CBU/CBU_Projects/Neuroscience/Jiaqi/Original

for file in $(ls *txt)
do
	wc -l $file
done


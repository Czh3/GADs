### mouse

for i in `ls ../Juicer/merge/*.hic`
do
	python Gene_Interaction.py $i > ${i/.hic/.gene_IA} &
done



### human


for i in `ls /lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/*.hic/*allValidPairs.hic`
do
	python Gene_Interaction_hg19.py $i > ${i/.allValidPairs.hic/.gene_IA} &
done

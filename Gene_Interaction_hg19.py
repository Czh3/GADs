import sys
import os
#import multiprocessing as mp

def strawHic(pos1, pos2, hic):
		cmd = "./straw KR %s %s %s BP 10000" % (hic, pos1, pos2)
		#cmd = "/lustre/user/liclab/zhangc/tools/juicer-master/CPU/common/juicer_tools dump oe KR %s %s %s  BP 10000" % (hic, pos1, pos2)
		response = os.popen(cmd).readlines()
		#result=sum([int(i.strip().split("\t")[2]) for i in response])
		n=0
		result=0
		#if len(response) < 30:
		#	return(-1)
		for i in response:
			j = i.strip().split("\t")
			if j[0] == j[1]:
				continue
			if j[2] == 'nan':continue
			if int(j[1]) - int(j[0]) <= 30000:
				continue
			#if(float(j[2])<1):continue
			result = result+float(j[2])
			n = n+1
		if n==0:
			return(0)
		return(result/n)

#print strawHic("10:64320000:64320001", "10:64620000:64720001", "../Juicer/merge/ST.hic")


window = 1000000
for i in open("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/gene_body.hg19.bed"):
	j = i.strip().split("\t")
	if (j[3]!="protein_coding"):
		continue
	gene_length = int(j[2]) - int(j[1])
	if gene_length < 100000:
		continue
	pos1 = j[0].strip("chr")+":"+j[1]+":"+j[2]
	#pos2 = j[0].strip("chr")+":"+ str(int(j[1])-window) +":"+str(int(j[1])+window)
	res = strawHic(pos1, pos1, sys.argv[1])
	#print j[-1]+"\t"+str(res)
	
	pos_left = j[0].strip("chr")+":"+ str(int(j[1])-gene_length) +":"+j[1]
	pos_right = j[0].strip("chr")+":"+ j[2] +":"+ str(int(j[2])+gene_length)
	control_left = strawHic(pos_left, pos_left, sys.argv[1])
	control_right = strawHic(pos_right, pos_right, sys.argv[1])
	#print j[-1]+"\t"+str(2*res/(control_left+control_right))
	print("\t".join([j[-1], str(res), str(control_left), str(control_right)]))

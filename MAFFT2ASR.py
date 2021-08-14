class MAFFTTOASR():
	def __init__(self, inputfilename,outputfilename):
		mafft_output = []
		mafft_form = []
		mafft_formchange = []
		inputphylo = []
		outputphylo = []
		mafft_output = "mafft-align.out"
		mafft_form = "fasta"
		mafft_formchange = "mafft-formchange.out"
		inputphylo = mafft_formchange+"_phyml_tree.txt"
		outputphylo = "paml-"+mafft_formchange+"-tree.txt"
	
		os.system("mafft %(inputfilename)s > %(mafft_output)s"%vars())
		self.ALIGNCONV_MAFFT2PAML(mafft_output,mafft_form,mafft_formchange)
		os.system("phyml -i %(mafft_formchange)s -b 100 -m JTT -d aa"%vars())
		self.ALIGNCONV(mafft_output,mafft_form, outputfilename)
		self.PHYLOCONV(inputphylo,outputphylo)
		
		self.CODEMLGENERATOR(outputfilename,outputphylo)
		os.system("codeml codeml_param.txt")

	#Subroutine "ALIGNCONV" generates nexus format of alignment data which can apply to run PAML.
	#The generated file was named as defined in "-OUTPUTFILE".
	def ALIGNCONV(self, inputfilename,inputformat,outputfilename):
		outputformat = []
		outputformat = "nexus"
		align = AlignIO.read(open(inputfilename,"rU"),inputformat,alphabet = IUPAC.extended_protein)
		merged_data = ''.join(align.format(outputformat).split("\t"))
		num_of_seqs,residue_num,sequences = self.DATAPICK(merged_data)
		self.OUTPUTDATA(num_of_seqs,residue_num,sequences,outputfilename)

	def DATAPICK(self,merged_data):
		tmp_data = [];num_of_seqs = int(0);residue_num = int(0);seq_flag = int(0);sequences = []
		tmp_data = merged_data.split("\n")
		for i in range(len(tmp_data)):
			tmp_data1 = []
			tmp_data1 = tmp_data[i]
			if re.search("^dimension.*",tmp_data1):
				tmp1 = []
				tmpa = []
				tmpa = ''.join(re.split(";",tmp_data1))
				tmp1 = re.split("\s+",tmpa)
				#print tmp1
				for j in range(len(tmp1)):
					tmp_param = []
					tmp_param = tmp1[j]
					if re.search("ntax.*",tmp_param):
						#print re.split("=",tmp_param)
						num_of_seqs = int(re.split("=",tmp_param)[1])
					elif re.search("nchar.*",tmp_param):
						residue_num = int(re.split("=",tmp_param)[1])
			
			if re.search("matrix",tmp_data1):
				seq_flag = int(1)
				continue
	
			if re.search(";",tmp_data1):
				seq_flag = int(0)
			
			if seq_flag == int(1):
				tmp_seq = []
				tmp_seq = re.split("\s+",tmp_data1)
				#print tmp_seq
				for k in range(len(tmp_seq)):
					tmp_seq1 = []
					tmp_seq1 = tmp_seq[k]
					sequences.append(tmp_seq1)
					sequences.append("\n")
				sequences.append("\n")
		
		return num_of_seqs,residue_num,sequences	

	def OUTPUTDATA(self,num_of_seqs,residue_num,sequences,outputfilename):
		output_data = open(outputfilename,"w")
		output_data.write("  %(num_of_seqs)i    %(residue_num)i\n\n\n"%vars())
		for i in range(len(sequences)):
			tmp_output = []
			tmp_output = sequences[i]
			output_data.write(tmp_output)
		output_data.close()

	def PHYLOCONV(self,inputphylo,outputphylo):
		phyloformat = [];phylo_flag = int(0);phylo_data = [];seq_num = int(0)
		phyloformat = "newick"
		Phylo.convert(inputphylo,phyloformat,"temporary.nex","nexus")
		tmp_genphylo = ''.join(open("temporary.nex","rU").read().split("\t"))
		tmp_genphylo = tmp_genphylo.split("\n")
		#print tmp_genphylo
		for i in range(len(tmp_genphylo)):
			tmp_data = []
			tmp_data = tmp_genphylo[i]
			if re.search(".*NTax.*",tmp_data):
				tmp_num = []
				tmp_num = re.split("=",tmp_data)[1]
				seq_num = int(''.join(re.split(";",tmp_num)))
			if re.search(".*Tree tree.*=.*",tmp_data):
				tmp_data1 = []
				tmp_data1 = re.split("=",tmp_data)[1]
				tmp_data1 = ''.join(re.split(";",tmp_data1))
				phylo_data.append(tmp_data1)
		phylo_data.append(";")
		phylo_data = ''.join(phylo_data)	
	

		output_data = open(outputphylo,"w")
		output_data.write("%(seq_num)i   1\n"%vars())
		output_data.write("%(phylo_data)s\n"%vars())
	
		output_data.close()

	def CODEMLGENERATOR(self,outputfilename,outputphylo):
		parameter_filename = "codeml_param.txt"
		ancestor_file = "ancestral_seq.out"
		output_data = open(parameter_filename,"w")
		output_data.write("seqfile = %(outputfilename)s\n"%vars())
		output_data.write("treefile = %(outputphylo)s\n"%vars())
		output_data.write("outfile = %(ancestor_file)s\n\n"%vars())
		output_data.write("noisy = 9   * 0,1,2,3,9: how much rubbish on the screen\n\
	      verbose = 1   * 1: detailed output, 0: concise output\n\
	      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic * 3: StepwiseAddition; (4,5):PerturbationNNI \n\
	      seqtype = 2   * 1:codons; 2:AAs; 3:codons-->AAs\n\
	    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n\
	        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree\n\
	        model = 2   * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n\
		aaDist = 0\n\
		aaRatefile = /Users/snakano/bioinformatics/paml4.8/dat/jones.dat \n\
	      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive\n\
	        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below\n\
	    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated\n\
	        kappa = 2   * initial or fixed kappa\n\
	    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate\n omega = .4   * initial or fixed omega, for codons or codon-transltd AAs\n\
	    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha\n \
	        alpha = 0.  * initial or fixed alpha, 0:infinity (constant rate)\n \
	       Malpha = 0   * different alphas for genes\n \
	        ncatG = 3   * # of categories in the dG or AdG models of rates\n \
		fix_rho = 1 \n\
		rho = 0.\n \
	        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n \
	 	RateAncestor = 1   * (1/0): rates (alpha>0) or ancestral states (alpha=0)\n \
	       Small_Diff = .5e-6\n \
		*  fix_blength = 1  * 0: ignore, -1: random, 1: initial, 2: fixed\n \
	       method = 0   * 0: simultaneous; 1: one branch at a time)")
		
		output_data.close()

	def ALIGNCONV_MAFFT2PAML(self,inputfilename,inputformat,outputfilename):
		outputformat = []
		outputformat = "phylip-relaxed"
		align = AlignIO.read(open(inputfilename,"rU"),inputformat,alphabet = IUPAC.extended_protein)
		print(align.format(outputformat))
		merged_data = ''.join(align.format(outputformat))
	
		output = open(outputfilename,"w")
		for i in range(len(merged_data)):
			tmp = []
			tmp = merged_data[i]
			output.write("%(tmp)s"%vars())
		
		output.close()
			
		#num_of_seqs,residue_num,sequences = DATAPICK(merged_data)
		#OUTPUTDATA(num_of_seqs,residue_num,sequences,outputfilename)
	
class PAMLLOGANALYSIS():
	def __init__(self, inputfilename, inputlibrary):
		initnode = int(0)
		lastnode = int(0)

		initnode, lastnode = self.NODENUMBER(inputfilename)
		ancseqs = self.ANCSEQS(inputfilename)
		nodelist, ppscores = self.LOGPARAM(initnode,lastnode,inputfilename)
		ancnode,nodes,maxident,minident,commonanc = self.COMMONANC(ancseqs,inputlibrary)
		self.OUTPUTLOG(nodelist, ppscores, ancseqs, ancnode,nodes,maxident,minident,commonanc)

	def ANGLEROUTCONVERT(self,inputlibrary):
		data = [];data2 = []
		data = open(inputlibrary,'r').read().split("\n")
		data = self.NULLELIM(data)
		del data[0]
		#print data
		#data2 = '\n'.join(data2)
		newlibname = []
		newlibname = "mod-"+inputlibrary
		output = open(newlibname,"w")
		for i in range(len(data)):
			tmp_data2 = []
			tmp_data2 = data[i]
			output.write("%(tmp_data2)s\n"%vars())
		output.close()
		return newlibname

	def NULLELIM(self, cvs):
		aft_data = []
		for i in range(len(cvs)):
			temp = []
			temp = cvs[i]
			if re.search('.*\w+.*',temp):
				aft_data.append(temp)
			else:
				continue
		return aft_data

	def NODENUMBER(self, inputfilename):
		alldata = [];initnum = int(0);lastnum = int(0)
		alldata = open(inputfilename,'r').read().split("\n")
		for i in range(len(alldata)):
			tmp = []
			tmp = alldata[i]
			if re.search("^Nodes.*ancestral$",tmp):
				#print tmp
				initnum = int(re.split("\s+",tmp)[1])
				lastnum = int(re.split("\s+",tmp)[3])
				break
		return initnum,lastnum

	def CONVFASTA(self,nexseqs):
		data = []
		data = re.split("\s+",nexseqs)
		data.insert(2,"\n")
		data.insert(0,">")
		data = ''.join(data)
		#print data
		return data
	

	def ANCSEQS(self,inputfilename):
		alldata = [];flag = int(0);sequences = []
		alldata = open(inputfilename,'r').read().split("\n")
		for i in range(len(alldata)):
			tmp = []
			tmp = alldata[i]
			if re.search("Overall accuracy of the 1 ancestral sequences.*",tmp):
				break
			if re.search("List of extant and reconstructed sequences",tmp):
				flag += int(1)
				continue
			if flag == int(1):
				sequences.append(tmp)
			
		del sequences[0:2]
		sequences = '\n'.join(sequences)
	
		seqtofasta = [];ancseqs = []
		seqtofasta = sequences.split("\n")
		for j in range(len(seqtofasta)):
			tmp2 = []
			tmp2 = seqtofasta[j]
			if re.search("^node.*",tmp2):
				tmptofasta = []
				tmptofasta = self.CONVFASTA(tmp2)
				ancseqs.append(tmptofasta)
	
		return ancseqs

	def LOGPARAM(self,initnode,lastnode,inputfilename):
		nodes = []
		nodelist = range(initnode,lastnode+1)
		for i in range(len(nodelist)):
			tmp = []
			tmp = "node#"+str(nodelist[i])
			nodes.append(tmp)
	
		alldata = [];flag = int(0);pplist = []
		alldata = open(inputfilename,'r').read().split("\n")
		#print alldata
	
		for j in range(len(alldata)):
			tmp = []
			tmp = alldata[j]
			if re.search("Prob of best state at each node, listed by site",tmp):
				flag = int(1)
				continue
				
			if re.search("Summary of changes along branches.",tmp):
				flag = int(0)			
				break
	
			if flag == int(1):
				pplist.append(tmp)

		print(pplist)
	
		del pplist[0:3]
		del pplist[-1]
		
		pplist = ','.join(pplist)
		
		pplist2 = []
		pplist2 = re.split(",",pplist)
	
		scorearray = np.zeros([len(nodes),len(pplist2)],float)
	
		#print pplist2
	
		for k in range(len(pplist2)):
			tmp2 = [];tmpdata2 = []
			tmp2 = pplist2[k]
			tmpdata = re.split(":",tmp2)[1]
			tmpdata2 = re.split("\s+",tmpdata)
			tmpdata2 = self.NULLELIM(tmpdata2)
		
			for l in range(len(tmpdata2)):
				tmp_param = [];tmp_param2 = float(0);tmp3 = []
				tmp3 = tmpdata2[l]
				tmp_param = re.split("\(",tmp3)[1]
				tmp_param2 = float(re.split("\)",tmp_param)[0])
				scorearray[l,k] = tmp_param2
		
		return nodes,scorearray 

	def CALCPARAM(self, tmp_ppscores):
		numof090 = int(0); numof080 = int(0); numof070 = int(0); numof050 = int(0); averagepp = float(0); allresnum = int(0)
		allresnum = len(tmp_ppscores)
		averagepp = np.average(tmp_ppscores)
		for i in range(len(tmp_ppscores)):
			tmp_ppscore_column = float(0)
			tmp_ppscore_column = tmp_ppscores[i]
	
			if tmp_ppscore_column >= float(0.50):
				numof050 += int(1.0)
				if tmp_ppscore_column >= float(0.70):
					numof070 += int(1.0)
					if tmp_ppscore_column >= float(0.80):
						numof080 += int(1.0)
						if tmp_ppscore_column >= float(0.90):
							numof090 += int(1.0)
	
		return averagepp,allresnum,numof090,numof080,numof070,numof050	

	def OUTPUTLOG(self,nodelist, ppscores, ancseqs, ancnode, nodes, maxident, minident, commonanc):
	
		if os.path.isdir("result-PAMLlog"):
			shutil.rmtree("result-PAMLlog")
		os.mkdir("result-PAMLlog")
	
		avepp_list = np.zeros(len(nodelist),float)	
			
		for i in range(len(nodelist)):
			resnum = int(0)
			tmp_nodename = []
			tmp_nodecheck = []
			tmp_nodecheck = nodelist[i]
	
			for k in range(len(ancseqs)):
				tmpseq = [];tmp_seq2 = []
				tmpseq = ancseqs[k].split("\n")[0]
				if re.search(tmp_nodecheck,tmpseq):
					tmp_seq2 = ancseqs[k]
					break
	
			tmp_nodename = ''.join(re.split("#",nodelist[i]))
			tmp_outname = "result-PAMLlog/"+tmp_nodename+"-ancestral-data.log"
			tmp_ppscores = ppscores[i]
			averagepp,allresnum,numof090,numof080,numof070,numof050 = self.CALCPARAM(tmp_ppscores)
			
			avepp_list[i] = averagepp
	
			tmp_outputfile = open(tmp_outname,"w")
			tmp_outputfile.write("#The ancestral sequence name: %(tmp_nodename)s\n"%vars())
			tmp_outputfile.write("#The average posterior probability values was: %(averagepp)7.3f\n"%vars())
			tmp_outputfile.write("#Total number of residues was: %(allresnum)i\n"%vars())
			tmp_outputfile.write("#Number of residues which bear >0.90 of PP value: %(numof090)i\n"%vars())
			tmp_outputfile.write("#Number of residues which bear >0.80 of PP value: %(numof080)i\n"%vars())
			tmp_outputfile.write("#Number of residues which bear >0.70 of PP value: %(numof070)i\n"%vars())
			tmp_outputfile.write("#Number of residues which bear >0.50 of PP value: %(numof050)i\n"%vars())
			tmp_outputfile.write("#The residue number, posterior probability\n")
			for j in range(len(tmp_ppscores)):
				resnum += int(1)
				tmp_ppval = float(0)
				tmp_ppval = tmp_ppscores[j]
				tmp_outputfile.write("  %(resnum)i, %(tmp_ppval)7.3f\n"%vars())
			tmp_outputfile.write("\n")
			tmp_outputfile.write("%(tmp_seq2)s\n"%vars())
			tmp_outputfile.close()

		print(nodelist)	
			
		caindex = int(0);comanc_pp = float(0)
		comname = commonanc.split("\n")[0]
		comname2 = re.split(">",comname)[1]
		caindex = int(nodelist.index(comname2))
		comanc_pp = avepp_list[caindex]
		commonanc_tmp = commonanc.split("\n")
		commonanc_tmp.insert(1,";"+str(comanc_pp))
		commonanc_tmp.insert(2,"\n")
		commonanc = ''.join(commonanc_tmp)
		output_commonanc = "result-PAMLlog/commonanc.fasta"
		outputCA = open(output_commonanc,"w")
		outputCA.write("%(commonanc)s"%vars())
		outputCA.close()
	
		summary_output = "result-PAMLlog/summary.out"
		outputsummary = open(summary_output,"w")
		outputsummary.write("#Summary of analysis data\n")
		outputsummary.write("#Ancestral sequences which are near to common ancestor in the input library are as following order:")
	
		for l in range(len(ancnode)):
			tmp_node = []
			tmp_node = ancnode[l]
			outputsummary.write(" %(tmp_node)s,"%vars())
	
		outputsummary.write("\n")
		outputsummary.write("#Sequence identity which bears the highest value when we compared between an ancestral sequence and a library sequence: ")
		for m in range(len(maxident)):
			tmp_maxident = float(0);tmp_node1 = []
			tmp_maxident = maxident[m]
			tmp_node1 = re.split(">",nodes[m])[1]
			outputsummary.write(" %(tmp_maxident)6.3f;%(tmp_node1)s,"%vars())
		outputsummary.write("\n")
	
		outputsummary.write("#Sequence identity which bears the lowest value when we compared between an ancestral sequence and a library sequence: ")
		for n in range(len(minident)):
			tmp_minident = float(0);tmp_node2 = []
			tmp_minident = minident[n]
			tmp_node2 = re.split(">",nodes[n])[1]
			outputsummary.write(" %(tmp_minident)6.3f;%(tmp_node2)s,"%vars())
		outputsummary.write("\n")
		outputsummary.close()

	def IDENTITYCALC(self,tmp_outname):
		mafftout = [];seq1 = [];seq2 = [];flag = int(0)
		mafftout = open(tmp_outname,'r').read().split("\n")
		#print mafftout
		for i in range(len(mafftout)):
			tmp = []
			tmp = mafftout[i]
			
			if re.search("^>.*",tmp):
				if flag == int(0):
					flag = int(1)
					continue
			if re.search("^>.*",tmp):
				#print "test"
				if flag == int(1):
					flag = int(2)
					continue
			if flag == int(1):
				seq1.append(tmp)
			
			if flag == int(2):
				seq2.append(tmp)
		seq1 = ''.join(seq1)
		seq2 = ''.join(seq2)
		
		seq1 = list(seq1)
		seq2 = list(seq2)
	
		total_resnum = int(0)
		total_resnum = int(len(seq1))
		identresnum = int(0)		
		
		#print len(seq1)
		#print len(seq2)
	
		for j in range(len(seq1)):
			tmp_seq1 = seq1[j]
			tmp_seq2 = seq2[j]
			if re.search(tmp_seq1,tmp_seq2):
				identresnum += 1
		identity = float(0)
		identity = float(identresnum)/float(total_resnum)*float(100.0)
		
		return identity				

	def COMMONANC(self,ancseqs,inputlibrary):
	
		libraryseqs = [];flag = int(0);alignedseqs = []
		libraryseqs = open(inputlibrary,'r').read().split("\n")
		for i in range(len(libraryseqs)):
			tmp = []
			tmp = libraryseqs[i]
			if flag == int(1):
				if re.search("^>.*",tmp):
					tmpseq = ''.join(tmpseq)
					alignedseqs.append(tmpseq)
	
			if re.search("^>.*",tmp):
				tmpseq = []
				flag = int(1)
				tmpseq.append(tmp)
				tmpseq.append("\n")
				continue
	
			if flag == int(1):
				tmpseq.append(tmp)
				continue
		#In "alignedseqs", each sequences are contained in the divided lists.
		
		maxident = np.zeros(len(ancseqs),float)
		minident = np.zeros(len(ancseqs),float)
		diffident = np.zeros(len(ancseqs),float)
		nodes = []
	
		for j in range(len(ancseqs)):
			tmpanc = [];tmp_node = [];tmp_identarray = np.zeros(len(alignedseqs),float)
			tmpanc = ancseqs[j]
			tmp_node = tmpanc.split("\n")[0]
			nodes.append(tmp_node)
	
			for l in range(len(alignedseqs)):
				tmpalignedseq = [];tmpidentity = float(0)
				tmpalignedseq = alignedseqs[l]
	
				if os.path.isfile("temporary.txt"):os.remove("temporary.txt")
				if os.path.isfile("temporary.out"):os.remove("temporary.out")
	
				tmpoutputfile = open("temporary.txt","w")
				tmpoutputfile.write("%(tmpanc)s\n%(tmpalignedseq)s\n"%vars())
				tmpoutputfile.close()
				os.system("mafft temporary.txt > temporary.out"%vars())
				
				tmp_outname = []
				tmp_outname = "temporary.out"		
				tmpidentity = self.IDENTITYCALC(tmp_outname)
				tmp_identarray[l] = tmpidentity
	
			tmpmax = float(0);tmpmin = float(0)
			tmpmax = np.max(tmp_identarray)
			tmpmin = np.min(tmp_identarray)
			maxident[j] = tmpmax
			minident[j] = tmpmin
		
		diffident = maxident-minident
		order = [];commonanc = []
		for m in range(len(diffident)):
			tmp = []
			tmp = str(np.argmin(diffident))
	
			if m == int(0):
				#commonanc.append(nodes[int(tmp)])
				commonanc.append(ancseqs[int(tmp)])
				commonanc = '\n'.join(commonanc)
	
			diffident[int(tmp)] = float(100.0)
			order.append(tmp)
	
		ancnode = []
		for n in range(len(order)):
			tmp1 = []
			tmp1 = order[n]
			ancnode.append(nodes[int(tmp1)])
		
		#print ancnode
		
		return ancnode,nodes,maxident,minident,commonanc

def ANGLEROUTCONVERT(inputlibrary):
	data = [];data2 = []
	data = open(inputlibrary,'r').read().split("\n")
	data = NULLELIM(data)
	del data[0]
	newlibname = []
	newlibname = "mod-"+inputlibrary
	output = open(newlibname,"w")
	for i in range(len(data)):
		tmp_data2 = []
		tmp_data2 = data[i]
		output.write("%(tmp_data2)s\n"%vars())
	output.close()
	return newlibname

def NULLELIM(cvs):
	aft_data = []
	for i in range(len(cvs)):
		temp = []
		temp = cvs[i]
		if re.search('.*\w+.*',temp):
			aft_data.append(temp)
		else:
			continue
	return aft_data

from Bio import AlignIO,SeqIO
from Bio import Phylo
from Bio.Alphabet import IUPAC
import os,sys,re,shutil
import numpy as np

anglerflag = int(0)

while len(sys.argv)>1:
	option = sys.argv[1]
	del sys.argv[1]
	if option == "-INPUTFILE":
		inputfilename = sys.argv[1]
		del sys.argv[1]
	elif option == "-OUTPUTFILE":
		outputfilename = sys.argv[1]
		del sys.argv[1]
	elif option == "-IANGLER":
		anglerflag = int(1)

inputlogname = []
inputlogname = "rst"

if anglerflag == int(1):
	inputfilename = ANGLEROUTCONVERT(inputfilename)

MAFFTTOASR(inputfilename,outputfilename)
PAMLLOGANALYSIS(inputlogname,inputfilename)


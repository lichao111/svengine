#!/usr/bin/env python
import csv,time,argparse 
class CommentedFile:
	def __init__(self, f, commentstring="#"):
		self.f = csv.reader(f,delimiter='\t')
		self.commentstring = commentstring
	def next(self):
		line = self.f.next()
		while line[0].startswith(self.commentstring):
			line = self.f.next()
		return line
	def __iter__(self):
		return self

def run():
	parser = argparse.ArgumentParser(description="spike-in mutations to fasta files, \
                                                we insert mutations into a provided reference \
                                                this script is haplotype insensitive and ligation less")
	parser.add_argument('-a', '--aligner', dest='aligner', default='bwamem',
			                        help="which aligner to use")
	parser.add_argument('-p', '--prefix', dest='prefix', default='test',
			                        help="in/out prefix")
	parser.add_argument('-r', '--reffile', dest='reffile', default='example.fna',
			                        help="reffile")
	parser.add_argument('-m', '--mutfile', dest='mutfile', default='adeno.fna',
			                        help="mutfile")
	parser.add_argument('-s', '--simfile', dest='simfile', default='test.txt',
			                        help="simfile")
	args = parser.parse_args()
	main(args)

def echocmd(cmd):
	return "echo [sim]"+cmd+";"+cmd

def main(args):
	i=1
	aligner = args.aligner
	prefix = args.prefix; rf = args.reffile; mf = args.mutfile; sf = args.simfile
	for line in CommentedFile(open(sf,'rU')):
		freq1 = str(int(100*float(line[0]))); cvg = line[1]; rpn1 = line[2]; rpn0 = line[3]
		fp = '%s.c%s.f%s' % (prefix,cvg,freq1)
		timecmd1='SECONDS=0'
		wgsimcmdm='wgsim -d 300 -s 30 -N %s -1 100 -2 100 %s %s.m.fq1 %s.m.fq2 1>/dev/null' % (rpn1,mf,fp,fp)
		wgsimcmdr='wgsim -d 300 -s 30 -N %s -1 100 -2 100 %s %s.r.fq1 %s.r.fq2 1>/dev/null' % (rpn0,rf,fp,fp)

		if aligner == 'bwamem':
			aligncmdm='bwa mem -t 16 %s %s.m.fq1 %s.m.fq2 2>/dev/null | samtools view -T %s -bS -o %s.m -' % (rf,fp,fp,rf,fp)
			aligncmdr='bwa mem -t 16 %s %s.r.fq1 %s.r.fq2 2>/dev/null | samtools view -T %s -bS -o %s.r -' % (rf,fp,fp,rf,fp)
		else:
			raise "unimplemented aligner method"

		rmcmd1='rm -f %s.m.fq1 %s.m.fq2 %s.r.fq1 %s.r.fq2' % (fp,fp,fp,fp)
		if int(rpn0) == 0: #ref is nothing
			aligncmdr = 'sleep 0' 
			mergecmd = 'mv %s.m %s.bam' % (fp,fp)
		elif int(rpn1) == 0: #mut is nothing
			aligncmdm = 'sleep 0' 
			mergecmd = 'mv %s.r %s.bam' % (fp,fp)
		else:
			mergecmd='samtools merge -f %s.bam %s.m %s.r' % (fp,fp,fp)

		rmcmd2='rm -f %s.m %s.r' % (fp,fp)
		timecmd2='echo $SECONDS'
		print "\"(%s)\"" % ';'.join([timecmd1, wgsimcmdm, wgsimcmdr, aligncmdm, aligncmdr, rmcmd1, mergecmd, rmcmd2, timecmd2])
		#print "echo \"(%s)\" | qsub -d $PWD -N j%d" % (';'.join([timecmd1, wgsimcmdm, wgsimcmdr, aligncmdm, aligncmdr, rmcmd1, mergecmd, rmcmd2, timecmd2]), i)
		i+=1

if __name__ == '__main__':
	starttime=time.time()
	run()
	#print "total runtime: ", time.time()-starttime, " seconds"

#novoindex -k 14 -s 1 celegans.ndx elegans.dna.fa
#novoalign -f sim_l.fastq sim_r.fastq -d chrX

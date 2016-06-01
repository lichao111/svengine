#!/usr/bin/env python
#DEBUG: import svengine.mf.mutforge as mf
#NOTE: currently mergemax for mergefq is disabled

import copy, time, re, os, shutil, sys, random, csv, subprocess, traceback, argparse, tempfile, itertools, operator, collections, ConfigParser, json, imp, string
from itertools import izip
from collections import Counter
from multiprocessing import Pool
from functools import partial
import pysam, pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import shelve

starttime=time.time()
devnull=open(os.devnull, 'w')
io_quote = "this is for end user"
code_quote = 'this is for program'
maxrefnum = 24 			#for debugging and human, should set to a big number afterwards
#burnin = 100000			#burnin region for each chromosome, to avoid N rich regions
varfile_names = ['#VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)']
varbed_names = ['#CHROM','START','END','VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)']

def myTemporaryFile(prefix): #this provides a unique temporary file name only
	return os.path.join(tempfile.tempdir,prefix+next(tempfile._get_candidate_names())+"tmp") 
	#add suffix here to avoid gzip error: already has _z suffix -- unchanged

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

class Bunch(object):
	def __init__(self, adict):
		self.__dict__.update(adict)

#paradigm to make dictionary from local variables
#for i in ('apple', 'banana', 'carrot'):
#		fruitdict[i] = locals()[i]
#paradigm to generate random sequence
#''.join([random.choice('AGTC') for x in range(1200)])
#paradigm to get fasta file size
#pysam.Fastafile('example.fna').lengths
#fas2bed(fasfile)

def merge_file(inlist,output):
	with open(output, 'w') as outfile:
		for infile in inlist:
			shutil.copyfileobj(open(infile), outfile)

def deovlpvar(dic,maxdeovlp=2):
	nowdeovlp=0; 
	varlist=copy.deepcopy(dic['varlist'])
	#print varlist
	vartab=var2bed(Bunch(dic)) #vartab is a BedTool
	#print vartab
	while chkvar(vartab) != True and nowdeovlp<maxdeovlp: 
		varlist=copy.deepcopy(dic['varlist'])
		ovlptab=vartab.intersect(vartab,c=True,f=0).to_dataframe(names=varbed_names+["ovlap"]) 
		#NOTE: f=0 is important so that simple touching is considerred an overlap, which make chkvar work
		print ovlptab,len(ovlptab)
		#print varlist.keys()
		print >>sys.stderr, "total vars before deovlapping", len(ovlptab); assert len(ovlptab)==len(varlist)
		ovlpkey=[ ovlptab.iloc[j,3] for j in range(0,len(ovlptab)) if int(ovlptab.iloc[j,len(varbed_names)]) > 1 ]
		print >>sys.stderr, "Warn:", " overlapped events removed: ", ovlpkey 
		print >>sys.stderr, "Warn:", "complex events such as translocaitons maybe broken due to deovlapping"	
		for key in ovlpkey:
			del varlist[key] #exclude overlapping entries
		dic['varlist']=varlist #this is needed for var2bed
		dic['varcnt']=len(varlist) #this is needed for continuously counting
		vartab=var2bed(Bunch(dic))
		print >>sys.stderr, "total vars after deovlapping", len(vartab); assert len(varlist)==len(vartab)
		nowdeovlp=nowdeovlp+1
	return varlist
	
def runcmd(cmd):
	#print >>sys.stderr, cmd
	run=subprocess.Popen(cmd.split(),env=os.environ,cwd=os.getcwd(),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	run_info=run.communicate() 
	return run_info

def pipecmd(cmd1,cmd2):
	#print >>sys.stderr, cmd1, "|", cmd2
	#run=subprocess.Popen(cmd.split(),env=os.environ,cwd=os.getcwd(),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	#run_info=run.communicate() 
	run_info=os.system(" ".join([cmd1,"|",cmd2]))
	return run_info

def syscmd(cmd):
	run_info=os.system(cmd)
	return run_info

def run2cmd(cmd1,cmd2):
	#print >>sys.stderr, cmd1, "|", cmd2
	run2=subprocess.Popen(cmd2.split(),env=os.environ,cwd=os.getcwd(),stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	run1=subprocess.Popen(cmd1.split(),env=os.environ,cwd=os.getcwd(),stdout=run2.stdin)
	run2_info=run2.communicate()
	run1_info=run1.wait()
	return run1_info,run2_info

def CountUnique(mylist):
	h = {}
	for x in mylist: h[x] = h.get(x, 0) + 1
	return h

def isint(s):
	try: int(s); return True;
	except ValueError: return False;

def locatemin(a):
	smallest = min(a)
	return smallest, [index for index, element in enumerate(a) if smallest == element]

def demultiplex(m):
	return [i for i, x in enumerate([int(x) for x in bin(m)[2:][::-1]]) if x == 1]

#def CountUnique(mylist):
#	h = {}
#	for x in mylist: h[x] = h.get(x, 0) + 1
#	return h

#def isint(s):
#	try: int(s); return True;
#	except ValueError: return False;

#def locatemin(a):
#	smallest = min(a)
#	return smallest, [index for index, element in enumerate(a) if smallest == element]

#def demultiplex(m):
#	return [i for i, x in enumerate([int(x) for x in bin(m)[2:][::-1]]) if x == 1]

def test_map(a,b):
	print a,b

def main(args):

	#NOTE: varlist var location is varlist[var][1]+nligation, not converted
	#NOTE: varbed var location is varbed[var][1], converted
	#NOTE: sequence position is always converted

	nprocs = args.nprocs
	layout = args.layout
	outbam = args.outbam
	burnin = args.burnin
	edgein = args.edgein
	mergemax = args.mergemax
	tempfile.tempdir = args.tmpdir
	try:
		os.mkdir(tempfile.tempdir)
	except:
		print >>sys.stderr, tempfile.tempdir, "error!, either it exists or cann't be created! please change path"
		quit()
	print >>sys.stderr, "tmpdir=", tempfile.tempdir
	pybedtools.set_tempdir(tempfile.tempdir)

	if args.metafile != None and args.varfile != None:
		print >>sys.stderr, "either a .meta or a .var file can be specified but not both"
		quit()

	runmode=0 #1 fasta+meta, 2 bam+meta, 3 fasta+var, 4 bam+var
	if args.inbamfilenames!=None and args.metafile!=None:
		runmode=2;
	elif args.inbamfilenames!=None and args.varfile!=None:
		runmode=4;
	elif args.hapfiles!=None and args.metafile!=None:
		runmode=1;
	elif args.hapfiles!=None and args.varfile!=None:
		runmode=3;
	else:
		raise Exception("must provide meta/var file and bam/hap files!")

	if runmode==0:
		sys.stderr.write("inccorrect input specification, please use on of following: \n \
													(1) fasta+meta, (2) bam+meta, (3) fasta+var, (4) bam+var")

	hapfiles = args.hapfiles.split(':') 		#this is required input
	gapfile = args.gapfile 								#this is required input
	parfile = args.parfile								#parfile controls WGS
	reffile = args.reffile 								#this is required input
	nploid = args.nploid 									#this input has default
	plansize = args.plansize							#this input has default
	nligation = args.nligation						#this input has default
	parlist = ConfigParser.ConfigParser()
	parlist.readfp(parfile)
	metafile = args.metafile
	varfile = args.varfile
	#if varfile != None: print varfile.name.rstrip("var").rstrip(".") #weird behavior of rstrip(".meta") -> exampl
	#if metafile != None: print metafile.name.rstrip("meta").rstrip(".")
	oname = metafile.name.rstrip("meta").rstrip(".") if metafile != None else varfile.name.rstrip(".var").rstrip(".")
	oprefix = oname if not args.oprefix else args.oprefix

	if len(hapfiles) == 1: #reuse hapfiles for haplotype
		hapseq=[pysam.Fastafile(hapfiles[0])]*nploid
	elif len(hapfiles)==nploid:
		hapseq=[pysam.Fastafile(hapfiles[i]) for i in range(0,nploid)]
		for ref in freseq:
			assert ref.lengths == hapseq[0].lengths, "hapfile seqs must be of the same length cross haplos"
	else:
		raise Exception("incorrect number of hapfiles provided")

	try:
		sizefile = fas2size(hapseq[0].filename)		 #require hap file to be the same size and gap applies to all
		#FIXME: this sizefile has also to be hapseq dependent 
		gaptab=pybedtools.BedTool(gapfile.read(),from_string=True)
		gaptab=gaptab.flank(g=sizefile, b=edgein)							#need to merge
	except pybedtools.cbedtools.BedToolsFileError:
		raise Exception("incorrect gapfile provided, must in BED format")

	dic={} #common dictionary interface for passing parameters
	for i in ('hapseq','gaptab','reffile','plansize','nligation','nploid','oprefix','nprocs','edgein','burnin','mergemax'):
		dic[i] = locals()[i]

	print >>sys.stderr, "entering runmode", runmode
	#print dic

	if runmode==1 or runmode==2: #have meta input, no need of deovlapping
		metatab = CommentedFile(metafile)
		dic['varlist']=collections.OrderedDict()
		dic['varcnt']=0
		varlist,varcnt = makevar(Bunch(dic),metatab)
		dic['varlist']=varlist
		dic['varcnt']=varcnt
		
		var2file(Bunch(dic))
		print >>sys.stderr, "done var2file"

	if runmode==3 or runmode==4: #have var input, assume deovlapping is done externally
		dic['varfile'] = varfile
		varlist = file2var(Bunch(dic))
		dic['varlist'] = varlist
		dic['varcnt'] = len(varlist)

	#NOTE: we cannot allow deovlpvar here bc. varlist may well contain overlaps due to hetero zygosity
	#			 this deovlapping is skiiped for mutforge for now
	#varlist=deovlpvar(dic); dic['varcnt']=len(varlist); dic['varlist']=varlist; vartab=var2bed(Bunch(dic))
	#assert chkvar(vartab) == True, "==Warn: varfile still contains self overlaps, which is not allowed, bail out"
	vartab=var2bed(Bunch(dic))
	var2bed(Bunch(dic))
	print >>sys.stderr, "done var2bed"

	if layout:
		print >>sys.stderr, "done layout only, now quit"
		quit()
	
	if runmode==1 or runmode==3:
		dic['parlist'] = parlist
		ffq = var2fq(Bunch(dic))
		#print ffq
		print >>sys.stderr, "done var2fq", time.time()-starttime, " seconds"

		if outbam:
			fbam = fq2bam(ffq,Bunch(dic))
			print >>sys.stderr, "done fq2bam", time.time()-starttime, " seconds"

	if runmode==2 or runmode==4:
		inbamfilenames = args.inbamfilenames.split(':')
		var2bs(Bunch(dic)) # use bam surgeon, not implemented yet

	shutil.rmtree(tempfile.tempdir) #clean temp files
	print >>sys.stderr, "done mutforge"

	return None #all done

def fq2bam(ffq,bun):
	reffile=bun.reffile.name
	nlibrary=bun.parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(bun.parlist.get('xwgsim', 'libnames'))
	nprocs=bun.nprocs
	mergemax=bun.mergemax
	oprefix=bun.oprefix
	libbams = [None] * nlibrary
	for lib in range(0,nlibrary):
		fq1 = ffq[lib][0]; fq2 =	ffq[lib][1]
		#sbam = '.'.join([oprefix,libnames[lib],"bam"])
		tbam = tempfile.NamedTemporaryFile('w',delete=False).name
		#ssam = tempfile.NamedTemporaryFile('w',delete=False).name
		#bwasw_cmd1="bwa mem -t %s %s %s %s >%s" % (nprocs,reffile,fq1,fq2,ssam) 
		#bwasw_cmd2="samtools view -bS -o %s %s" % (sbam,ssam) #it auto gives bam if not provided
		bwasw_cmd1="bwa mem -t %s %s %s %s" % (nprocs,reffile,fq1,fq2) 
		bwasw_cmd2="samtools view -bS -o %s -" % (tbam) #it auto gives bam if not provided
		#bwasw_cmd3="rm -f %s" % ssam
		#print >>sys.stderr, bwasw_cmd1, "\n", bwasw_cmd2
		pipecmd(bwasw_cmd1,bwasw_cmd2)
		#bwasw_cmd="bwa mem -t %s %s %s %s | samtools view -bS -o %s - " % (nprocs,reffile,fq1,fq2,sbam) 
		#print >>sys.stderr, bwasw_cmd #
		#runcmd(bwasw_cmd)
		libbams[lib] = tbam #post processing: librarywise merging bams

	fbam = [None] * nlibrary
	for lib in range(0,nlibrary):
		sbam = '.'.join([oprefix,libnames[lib]])
		tbam = tempfile.NamedTemporaryFile('w',delete=False).name
		mbam = mergebam([libbams[lib]],tbam,mergemax)
		sort_cmd = "samtools sort -@ %s %s %s" % (nprocs, mbam, sbam)
		#print >>sys.stderr, sort_cmd #
		syscmd(sort_cmd)
		fbam[lib] = sbam+".bam" #post processing: librarywise merging bams

	return fbam

def var2bs(bun):
	#print >>sys.stderr, "this part will use bamsurgeon, it has not been implemented yet"; quit()
	numinbamfiles = len(inbamfilenames)
	oprefix = args.oprefix.split(':') if args.oprefix!=None else args.oprefix.split(':')
	return None

def fas2bed(fasfile):
	fasseq = pysam.Fastafile(fasfile)
	bed = [ fasseq.references[i]+"\t0\t"+str(fasseq.lengths[i]) for i in range(0,min(fasseq.nreferences,maxrefnum)) ]
	return pybedtools.BedTool('\n'.join(bed),from_string=True)

def fas2size(fasfile):
	sizefile = tempfile.NamedTemporaryFile('w',delete=False,prefix="fasize_")
	fasseq = pysam.Fastafile(fasfile)
	size = [ fasseq.references[i]+'\t'+str(fasseq.lengths[i]) for i in range(0,min(fasseq.nreferences,maxrefnum)) ]
	sizefile.write('\n'.join(size))
	sizefile.close()
	return sizefile.name

def chkvar(varbed):
	#print varbed.sort().merge()
	#print varbed.sort()
	if varbed.sort().merge().count() == varbed.sort().count(): return True;
	else: return False

def ppbam(bun,bamfiles):
	return None

def var2wgs(bun):
	hapseq=bun.hapseq
	nprocs=bun.nprocs
	mergemax=bun.mergemax
	oprefix=bun.oprefix
	parlist=bun.parlist
	nploid=bun.nploid
	nligation=bun.nligation
	nlibrary=bun.parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	varlist=bun.varlist
	reffile = bun.reffile.name
	parlist = bun.parlist
	#parfile in python ConfigParser format
	#LibraryNumber = 3
	#Have to use very complex merging

	#0, if no variant to be inserted, bail out. If some control is needed, try use fasforge for now
	if len(varlist)==0: #control case, this has to be implemented
		print >>sys.stderr, "no variant specified, use fasforge/wgsim directly, bail out..."
		quit()
	
	#1, construct continuous contig sets with ligations (nploid)*(variant/non-variant)
	#hapvar{0: bed_features}
	hapvar = collections.OrderedDict(); hapbed = collections.OrderedDict();
	#print varlist.values()
	for var in varlist.values():
		#print var
		#print [var[4],var[5],var[5]+var[7],var[0],var[1],var[2],var[6],var[8]]
		#CHR ST	ED	VID	MID	FREQ	DEL	INS; adding one base size to ensure flanking regions are separate
		hapvar[var[3]] = hapvar.get(var[3],[]) + [[var[4],var[5],var[5]+var[7]+2*nligation+1,var[0],var[1],var[2],var[7],formatins(var[9])]]
	for hap in hapvar.keys():
		hapvar[hap] = pybedtools.BedTool('\n'.join([ '\t'.join([str(c) for c in x]) for x in hapvar[hap] ]), from_string=True)
		#assert chkvar(hapvar[hap]) == True, "varfile contains self overlaps, which is not allowed"
	for hap in hapvar.keys():
		sizefile = fas2size(hapseq[hap].filename)		 #create genome size file
		hapbed[hap] = hapvar[hap].complement(g=sizefile) #create complementary blocks
		hapbed[hap] = hapbed[hap].slop(g=sizefile,b=nligation)
		hapbed[hap] = pybedtools.BedTool('\n'.join([ f[0]+'\t'+str(f[1])+'\t'+str(f[2])+'\tNone\tNone\t0\t0\tNone' for f in hapbed[hap]]),from_string=True)
		hapvar[hap] = hapvar[hap].cat(hapbed[hap],postmerge=False)

	hapiter=[]
	#print len(hapvar[0])
	#print hapvar[0][0]

	idx = 1; total = 0; 
	for hap in hapvar.keys():
		total += len(hapvar[hap])
	for hap in hapvar.keys():
		hapfas=bun.hapseq[hap].filename
		for feat in hapvar[hap]:
			hapiter.append([hap,hapfas,nligation,parlist,reffile,nploid,idx,total,feat,mergemax])
			idx += 1
	#for it in hapiter:
	#	print it

	#NOTE: debug code for series execution
	bedbams = []
	if nprocs<2:
		for hapit in hapiter:
			bedbams.append(fakebam(hapit))
	else:
		pool = Pool(processes=int(nprocs))
		bedbams = pool.map(fakebam,hapiter,1)
	#print bedbams
	print >>sys.stderr, "done fakebam", 	time.time()-starttime, " seconds"	

	libbams = {} #bamfiles to be merged: [ [hap1lib1, hap1lib2, ...], [hap2lib1, hap2lib2, ...] ]
	#libcnt0 = 0; libcnt1 = 0
	for i in xrange(0,len(bedbams)):
		for lib in xrange(0,nlibrary):
			#print bedbams[i][1][lib];
			#if lib == 0: 
			#	libcnt0 += 1
			#else: 
			#	libcnt1 += 1;
			if libbams.get(lib,None) == None:
				libbams[lib]=[bedbams[i][1][lib]]
			else:
				libbams[lib].append(bedbams[i][1][lib]) 

	#print "len bedbams:", len(bedbams)
	#for lib in range(0,nlibrary):
	#	print "len libbams:", len(libbams[lib])
	#print	libcnt0,libcnt1
	#print bedbams
	#print libbams
	#quit()

	fbam = [None] * nlibrary
	for lib in range(0,nlibrary):
		sbam = '.'.join([oprefix,libnames[lib]])
		tbam = tempfile.NamedTemporaryFile('w',delete=False).name
		mbam = mergebam([libbams[lib]],tbam,mergemax)
		sort_cmd = "samtools sort -@ %s %s %s" % (nprocs, mbam, sbam)
		#print >>sys.stderr, sort_cmd #
		syscmd(sort_cmd)
		fbam[lib] = sbam+".bam" #post processing: librarywise merging bams

	print >>sys.stderr, "done merge to final bam"

	return fbam

def var2fq(bun):
	hapseq=bun.hapseq
	nprocs=bun.nprocs
	oprefix=bun.oprefix
	parlist=bun.parlist
	nploid=bun.nploid
	nligation=bun.nligation
	nlibrary=bun.parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	varlist=bun.varlist
	reffile = bun.reffile.name
	parlist = bun.parlist
	mergemax = bun.mergemax
	#parfile in python ConfigParser format
	#LibraryNumber = 3
	#Have to use very complex merging

	#0, if no variant to be inserted, bail out. If some control is needed, try use fasforge for now
	if len(varlist)==0: #control case, this has to be implemented
		print >>sys.stderr, "no variant specified, use fasforge/wgsim directly, bail out..."
		quit()
	
	#1, construct continuous contig sets with ligations (nploid)*(variant/non-variant)
	#hapvar{0: bed_features}
	hapvar = collections.OrderedDict(); hapbed = collections.OrderedDict();
	#print varlist.values()
	for var in varlist.values():
		#print var
		#print [var[4],var[5],var[5]+var[7],var[0],var[1],var[2],var[6],var[8]]
		#CHR ST	ED	VID	MID	FREQ	DEL	INS; adding one base size to ensure flanking regions are separate
		hapvar[var[3]] = hapvar.get(var[3],[]) + [[var[4],var[5],var[5]+var[7]+2*nligation+1,var[0],var[1],var[2],var[7],formatins(var[9])]]
	for hap in hapvar.keys():
		hapvar[hap] = pybedtools.BedTool('\n'.join([ '\t'.join([str(c) for c in x]) for x in hapvar[hap] ]), from_string=True)
		#assert chkvar(hapvar[hap]) == True, "varfile contains self overlaps, which is not allowed"
	for hap in hapvar.keys():
		sizefile = fas2size(hapseq[hap].filename)		 #create genome size file
		hapbed[hap] = hapvar[hap].complement(g=sizefile) #create complementary blocks
		hapbed[hap] = hapbed[hap].slop(g=sizefile,b=nligation)
		hapbed[hap] = pybedtools.BedTool('\n'.join([ f[0]+'\t'+str(f[1])+'\t'+str(f[2])+'\tNone\tNone\t0\t0\tNone' for f in hapbed[hap]]),from_string=True)
		hapvar[hap] = hapvar[hap].cat(hapbed[hap],postmerge=False)

	hapiter=[]
	#print len(hapvar[0])
	#print hapvar[0][0]

	idx = 1; total = 0; 
	for hap in hapvar.keys():
		total += len(hapvar[hap])
	for hap in hapvar.keys():
		hapfas=bun.hapseq[hap].filename
		for feat in hapvar[hap]:
			hapiter.append([hap,hapfas,nligation,parlist,reffile,nploid,idx,total,feat,mergemax])
			idx += 1
	#for it in hapiter:
	#	print it

	#NOTE: debug code for series execution
	bedfqs = []
	if nprocs<2:
		for hapit in hapiter:
			bedfqs.append(fakefq(hapit))
	else:
		pool = Pool(processes=int(nprocs))
		bedfqs = pool.map(fakefq,hapiter,1)
	#print bedfqs
	print >>sys.stderr, "done fakefq", 	time.time()-starttime, " seconds"	

	libfqs = {} #fqfiles to be merged: [ [hap1lib1, hap1lib2, ...], [hap2lib1, hap2lib2, ...] ]
	#hap1lib1=[hap1lib1fq1,hap1lib1fq2]
	#libcnt0 = 0; libcnt1 = 0
	for i in xrange(0,len(bedfqs)):
		for lib in xrange(0,nlibrary):
			if libfqs.get(lib,None) == None:
				libfqs[lib]=[bedfqs[i][1][lib]]
			else:
				libfqs[lib].append(bedfqs[i][1][lib]) 

	ffq = [None,None] * nlibrary
	for lib in range(0,nlibrary): #each libfqs[lib] has everything to be merged to ffq[lib]
		sfq = '.'.join([oprefix,libnames[lib]])
		ffq[lib] = mergefq(libfqs[lib],sfq,mergemax)

	print >>sys.stderr, "done merge to final fq"

	return ffq


def fakefq(hapit): # makebam from a pair of seqs of representing the original and sved haplotype
	hap=hapit[0]; hapfas=hapit[1]; nligation=hapit[2]; parlist=hapit[3]; reffile=hapit[4]; nploid=hapit[5]; bed=hapit[8]; mergemax=hapit[9]
	#print >>sys.stderr, "processing %d out of %d regions" % (hapit[6],hapit[7])
	#print >>sys.stderr, "bed=", bed
	tmpfas=pysam.Fastafile(hapfas)
	seq0=SeqRecord(Seq(tmpfas.fetch(str(bed[0]),int(bed[1])+1,int(bed[2])), generic_dna),id=bed[0], name="", description="")
	tmpfas.close()
	seq1=seq0
	#CHR ST	ED	VID	MID	FREQ	DEL	INS
	if int(bed[6])!=0: #work on deletion
		seq1.seq=seq1.seq[:nligation]+seq1.seq[nligation+int(bed[6]):];
	if bed[7]!="None": #work on insertion
		info = readins(bed[7])
		tmpfas=pysam.Fastafile(info[0])
		seqi = SeqRecord(Seq(pysam.Fastafile(info[0]).fetch(str(info[1]),info[2],info[2]+info[3]), generic_dna),id=info[0], name="", description="")
		tmpfas.close() #reduce open handles
		#['example.fasta', '22', 25202000, 1000, 2, False]
		if info[5]:
			seqi.seq = seqi.seq.reverse_complement();
		seqt=seqi
		for i in range(1,info[4]):
			seqi.seq += seqt.seq
		seq1.seq=seq1.seq[:nligation]+seqi.seq+seq1.seq[nligation:]
	nlibrary=parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	rl1=json.loads(parlist.get('xwgsim', 'read1'))
	rl2=json.loads(parlist.get('xwgsim', 'read2'))
	freq=float(bed[5])

	fqs = [None] * nlibrary
	for lib in range(0,nlibrary):
		libn=libnames[lib]; rl=min(rl1[libn],rl2[libn])
		freq0 = 1-freq; empty0=False; empty1=False;
		# a successful region would be no N in middle + not all N in ligation; goodseq holds if Ns are consecutive
		def goodseq(seq, nl, rl):
			gd = nl - seq[:nl].count('N') > rl and nl - seq[(len(seq)-nl):len(seq)].count('N') > rl and len(seq)-2*nl-seq[nl:len(seq)-nl].count('N') > rl
			return gd
		#print seq0.seq
		if not goodseq(seq0.seq, nligation, rl): empty0=True
		#print goodseq(seq0.seq, nligation, rl)
		if empty0:
			wgs0 = [None, None]
		else:
			wgs0 = fqwgs(parlist,nligation,reffile,nploid,[seq0],lib,freq0,mergemax)
		#print seq1.seq
		if not goodseq(seq1.seq, nligation, rl): empty1=True #FIXME, these N regions has to be masked from input
		#print goodseq(seq1.seq, nligation, rl)
		if empty1:
			wgs1 = [None, None]
		else:
			wgs1 = fqwgs(parlist,nligation,reffile,nploid,[seq1],lib,freq,mergemax)
		fqopre = myTemporaryFile(prefix="fakefq_")
		#only need the temporary filename here
		#print "wgs=",wgs0,wgs1
		fqs[lib] = mergefq([wgs0,wgs1],fqopre,mergemax)
		#print "fqs[lib]=",fqs[lib]
		#assert all([ os.path.isfile(x) for x in fqs[lib] if x != None ])
		#assert all([ os.path.isfile(x[0]) for x in fqs[lib] if x[0] != None ])
		#assert all([ os.path.isfile(x[1]) for x in fqs[lib] if x[1] != None ])
		
	print >>sys.stderr, "processed %d out of %d regions" % (hapit[6],hapit[7])
	return [hap,fqs]

def fqwgs(parlist,nligation,reffile,nploid,seq,lib,freq,mergemax): #create bam file for 1 library given seq and freq
	#passing both filename and contents to save space
	tmpfa=False #fa is temporary created or constant
	if(type(seq)==type({})):
		fa=seq['filename']; seqs=seq['seq'] #for fasforge calls
	else:
		fa=myTemporaryFile(prefix="fawgs_") #this has to be removed
		tmpfa=True
		seqs=seq
	seql = sum( [ len(s) for s in seqs ] ) #len of all seqs
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	coverage=json.loads(parlist.get('xwgsim', 'coverage'))
	isize=json.loads(parlist.get('xwgsim', 'isize'))
	rl1=json.loads(parlist.get('xwgsim', 'read1'))
	rl2=json.loads(parlist.get('xwgsim', 'read2'))
	par=' '.join(['%s %s' % (key, value) for (key, value) in json.loads(parlist.get('xwgsim', 'par')).items()])
	libn=libnames[lib]
	nrdp = seql*coverage[libn]/nploid*freq/(rl1[libn]+rl2[libn])
	nmode = len(isize[libn])
	mode = [ None ] * nmode
	if freq == 0 or nrdp == 0:
		return [ None, None ]
	else:
		for m in range(0,nmode):
			mden=isize[libn][m][0]; mnread = int(nrdp * mden)
			if mnread>0:
				SeqIO.write(seqs,fa,'fasta') #we need this temp fa file
				misize=isize[libn][m][1]; misizesd=isize[libn][m][2]
				fq1=myTemporaryFile(prefix="fq1wgs_")
				fq2=myTemporaryFile(prefix="fq2wgs_")
				wgs_cmd="xwgsim -h -d %d -s %d -N %d -1 %d -2 %d -l %d %s %s %s %s " % (misize,misizesd,mnread,rl1[libn],rl2[libn],nligation,par,fa,fq1,fq2)
				#print >>sys.stderr, wgs_cmd
				runcmd(wgs_cmd) #fastq file
				#we opt for speed instead of space
				#gzfq1_cmd="gzip %s" % fq1 #	fq1.gz, fq1 auto deleted
				#gzfq2_cmd="gzip %s" % fq2 #	fq2.gz, fq2 auto deleted
				#print >>sys.stderr, gzfq1_cmd
				#print >>sys.stderr, gzfq2_cmd
				#runcmd(gzfq1_cmd) #fastq file
				#runcmd(gzfq2_cmd) #fastq file
				#clean up
				if tmpfa: os.remove(fa)
				#fq3,fq4=filterfq(fq1,fq2,nligation,len(seq)-nligation) #half ligation region reads
				#print (fq1,fq2,fq3,fq4,nligation,len(seq)-nligation)
				#bwasw_cmd="bwa mem -t 4 %s %s %s | samtools view -bS -o %s - " % (reffile,fq1,fq2,mbam+".bam") #it auto gives bam if not provided
				#print bwasw_cmd
				#runcmd(bwasw_cmd)
				#TODO: for tmpfile in (fq1,fq2,fq3,fq4): os.remove(tmpfile)
				#mode[m]=[fq1+".gz",fq2+".gz"]
				mode[m]=[fq1,fq2]
			else:
				mode[m]=[None,None]
	#assert all([ os.path.isfile(x[0]) for x in mode if x[0] != None ])
	#assert all([ os.path.isfile(x[1]) for x in mode if x[1] != None ])
	libfq = mergefq(mode,myTemporaryFile(prefix="libfq_"),mergemax)
	#print "libfq=", libfq
	#assert all([ os.path.isfile(x) for x in libfq if x != None ])
	#quit()
	#assert all([ os.path.isfile(x[0]) for x in libfq if x[0] != None ])
	#assert all([ os.path.isfile(x[1]) for x in libfq if x[1] != None ])
	#TODO: for tmpfile in mode + [fa]: os.remove(tmpfile)
	return libfq

def mergefq(fqs,oprefix,mergemax):
	fq1s = [x[0] for x in fqs if x[0] is not None] #remove empty components
	fq2s = [x[1] for x in fqs if x[1] is not None] #remove empty components
	assert len(fq1s)==len(fq2s)	 #all reads have to be properly paired
	#check all fq1s are good
	nfqs = len(fq1s); ffq1=oprefix+".fq1"; ffq2=oprefix+".fq2";
	#print "nfqs,fq1s,fq2s=",nfqs,fq1s,fq2s
	#assert all([ os.path.isfile(x) for x in fq1s ])
	#assert all([ os.path.isfile(x) for x in fq2s ])
	if nfqs==0:
		return [None, None]
	elif nfqs==1:
		#print >>sys.stderr, "mv %s %s" % (fq1s[0],ffq1)
		#print >>sys.stderr, "mv %s %s" % (fq2s[0],ffq2)
		shutil.move(fq1s[0],ffq1) #if merge only one file, just move it to new name
		shutil.move(fq2s[0],ffq2) #if merge only one file, just move it to new name
	else:
		merge_file(fq1s,ffq1)
		merge_file(fq2s,ffq2)
		for f in fq1s: os.remove(f)
		for f in fq2s: os.remove(f)
	#else: #trunk merge nfqs 
	#	for i in range(0,(nfqs/mergemax)+1):
	#		if i == 0: #first iteration
	#			cfq1 = ""
	#			cfq2 = ""
	#			mfq1 = myTemporaryFile(prefix="mfq1_")
	#			mfq2 = myTemporaryFile(prefix="mfq2_")
	#		if i == nfqs/mergemax: #last iteration
	#			if cfq1 != "": cfq1 = mfq1
	#			if cfq2 != "": cfq2 = mfq2
	#			mfq1 = ffq1
	#			mfq2 = ffq2
	#		else:
	#			cfq1 = mfq1
	#			cfq2 = mfq2
	#			mfq1 = myTemporaryFile(prefix="mfq1_")
	#			mfq2 = myTemporaryFile(prefix="mfq2_")
	#		merge1_cmd = "cat %s %s >%s" % ( cfq1, " ".join([fq1 for fq1 in fq1s[i*mergemax:min((i+1)*mergemax,nfqs)]]), mfq1 )
	#		merge2_cmd = "cat %s %s >%s" % ( cfq2, " ".join([fq2 for fq2 in fq2s[i*mergemax:min((i+1)*mergemax,nfqs)]]), mfq2 )
	#		#print >>sys.stderr, merge1_cmd #unsorted
	#		#print >>sys.stderr, merge2_cmd #unsorted
	#		runcmd(merge1_cmd)
	#		runcmd(merge2_cmd)
	#		#clean temp files
	#		#for f in fq1s[i*mergemax:min((i+1)*mergemax,nfqs)]: os.remove(f)
	#		#for f in fq2s[i*mergemax:min((i+1)*mergemax,nfqs)]: os.remove(f)
	#		if cfq1 != "": os.remove(cfq1)
	#		if cfq2 != "": os.remove(cfq2)
	#FIXME: who do we clean up this mess of temporary files here?
	return [ffq1,ffq2]

def fakebam(hapit): # makebam from a pair of seqs of representing the original and sved haplotype
	hap=hapit[0]; hapfas=hapit[1]; nligation=hapit[2]; parlist=hapit[3]; reffile=hapit[4]; nploid=hapit[5]; bed=hapit[8]; mergemax=hapit[9]; #print mergemax
	print >>sys.stderr, "processing %d out of %d regions" % (hapit[6],hapit[7])
	seq0=SeqRecord(Seq(pysam.Fastafile(hapfas).fetch(str(bed[0]),int(bed[1])+1,int(bed[2])), generic_dna),id=bed[0], name="", description="")
	seq1=seq0
	#CHR ST	ED	VID	MID	FREQ	DEL	INS
	if int(bed[6])!=0: #work on deletion
		seq1.seq=seq1.seq[:nligation]+seq1.seq[nligation+int(bed[6]):];
	if bed[7]!="None": #work on insertion
		#print bed[7]
		info = readins(bed[7])
		#print info
		seqi = SeqRecord(Seq(pysam.Fastafile(info[0]).fetch(str(info[1]),info[2],info[2]+info[3]), generic_dna),id=info[0], name="", description="")
		#['example.fasta', '22', 25202000, 1000, 2, False]
		if info[5]:
			seqi.seq = seqi.seq.reverse_complement();
		seqt=seqi
		for i in range(1,info[4]):
			seqi.seq += seqt.seq
		seq1.seq=seq1.seq[:nligation]+seqi.seq+seq1.seq[nligation:]
	nlibrary=parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	freq=float(bed[5])

	bams = [None] * nlibrary
	for lib in range(0,nlibrary):
		freq0 = 1-freq
		bam0 = runwgs(parlist,nligation,reffile,nploid,[seq0],lib,freq0)
		bam1 = runwgs(parlist,nligation,reffile,nploid,[seq1],lib,freq)
		bams[lib] = mergebam([bam0,bam1],tempfile.NamedTemporaryFile('w',delete=False).name,mergemax)
	print >>sys.stderr, "processed %d out of %d regions" % (hapit[6],hapit[7])
	return [hap,bams]

def runwgs(parlist,nligation,reffile,nploid,seq,lib,freq): #create bam file for 1 library given seq and freq
	fa=tempfile.NamedTemporaryFile('w',delete=False).name
	SeqIO.write(seq,fa,'fasta')
	seql = sum( [ len(s) for s in seq ] ) #len of all seqs
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	coverage=json.loads(parlist.get('xwgsim', 'coverage'))
	isize=json.loads(parlist.get('xwgsim', 'isize'))
	rl1=json.loads(parlist.get('xwgsim', 'read1'))
	rl2=json.loads(parlist.get('xwgsim', 'read2'))
	par=' '.join(['%s %s' % (key, value) for (key, value) in json.loads(parlist.get('xwgsim', 'par')).items()])
	libn=libnames[lib]
	nrdp = seql*coverage[libn]/nploid*freq/(rl1[libn]+rl2[libn])
	nmode = len(isize[libn])
	mode = [ None ] * nmode
	#print "run wgs for:", libn, fa, coverage[libn], freq, nrdp
	if freq == 0 or nrdp == 0:
		return None
	else:
		for m in range(0,nmode):
			mden=isize[libn][m][0]; mnread = int(nrdp * mden)
			if mnread>0:
				misize=isize[libn][m][1]; misizesd=isize[libn][m][2]
				fq1=tempfile.NamedTemporaryFile('w',delete=False).name
				fq2=tempfile.NamedTemporaryFile('w',delete=False).name
				mbam=tempfile.NamedTemporaryFile('w',delete=False).name
				#print (misize,misizesd,mnread,rl1[libn],rl2[libn],fa,fq1,fq2)
				wgs_cmd="xwgsim -h -d %d -s %d -N %d -1 %d -2 %d -l %d %s %s %s %s " % (misize,misizesd,mnread,rl1[libn],rl2[libn],nligation,par,fa,fq1,fq2)
				#print >>sys.stderr, wgs_cmd
				runcmd(wgs_cmd) #fastq file
				#fq3,fq4=filterfq(fq1,fq2,nligation,len(seq)-nligation) #half ligation region reads
				#print (fq1,fq2,fq3,fq4,nligation,len(seq)-nligation)
				bwasw_cmd1="bwa mem -t 4 %s %s %s" % (reffile,fq1,fq2) 
				bwasw_cmd2="samtools view -bS -o %s -" % mbam+".bam" #it auto gives bam if not provided
				#print >>sys.stderr, bwasw_cmd1, " | ", bwasw_cmd2
				pipecmd(bwasw_cmd1,bwasw_cmd2)
				#TODO: for tmpfile in (fq1,fq2,fq3,fq4): os.remove(tmpfile)
				mode[m]=mbam+".bam"
			else:
				mode[m]=None
	libbam = mergebam(mode,tempfile.NamedTemporaryFile('w',delete=False).name,mergemax)
	#TODO: for tmpfile in mode + [fa]: os.remove(tmpfile)
	return libbam

def filterfq(fq1,fq2,leftcut,rightcut): #half ligation region reads, this is no longer needed for xwgsim
	fq3seq = []; fq4seq=[]; thincnt=0
	for rfq1,rfq2 in zip(SeqIO.parse(open(fq1,'rU'),'fastq'),SeqIO.parse(open(fq2,'rU'),'fastq')):
		assert rfq1.id[:-1] == rfq2.id[:-1], "unmatching fastq1 and fastq2, error in read simulation, quit"
		#@fins.seq1_461_964_1:0:0_1:0:0_0/1
		rdpl = int(rfq1.id.split('_')[1]); rdpr = int(rfq1.id.split('_')[2])
		if (rdpl < leftcut and rdpr < leftcut) or (rdpl > rightcut and rdpr > rightcut):
			#TODO: currently thin when both in ligation region, future may need check not thin if at chrom start/end
			thincnt+=1
			if thincnt % 2 == 1: #odd counts not thinned, thin exactly 1/2
				fq3seq.append(rfq1); fq4seq.append(rfq2)
		else:
			fq3seq.append(rfq1); fq4seq.append(rfq2)
	fq3=tempfile.NamedTemporaryFile('w',delete=False).name
	fq4=tempfile.NamedTemporaryFile('w',delete=False).name
	SeqIO.write(fq3seq,fq3,'fastq')
	SeqIO.write(fq4seq,fq4,'fastq')
	return fq3,fq4 #NOTE: these are freeed in upper level

def mergebam(bams,oprefix,mergemax):
	bams = [x for x in bams if x is not None] #remove empty components
	nbams = len(bams); fbam=oprefix+".bam";
	#print nbams
	#samtools merge can take up to 4092 files in one merge
	if nbams==1:
		shutil.move(bams[0],fbam) #if merge only one file, just move it to new name
	else:
		for i in range(0,nbams/mergemax+1):
			if i == 0: #first iteration
				cbam = ""
				mbam = tempfile.NamedTemporaryFile('w',delete=False).name
			if i == nbams/mergemax: #last iteration
				if cbam != "": cbam = mbam
				mbam = fbam
			else:
				cbam = mbam
				mbam = tempfile.NamedTemporaryFile('w',delete=False).name
			merge_cmd = "samtools merge -f %s %s %s " % (mbam, cbam, " ".join([bam for bam in bams[i*mergemax:min((i+1)*mergemax,nbams)]]))
			#print >>sys.stderr, merge_cmd #unsorted
			runcmd(merge_cmd) #merging is possible with unsorted bams
	return fbam

def lumpy2var(bp,ref,nligation):
	idx=1; delidx=1; insidx=1; traidx=1; invidx=1; dupidx=1; varlist=collections.OrderedDict()
	for line in bp:
		if line[10] == "TYPE:DELETION":
			mid="DEL_"+str(delidx); var="VAR_"+str(idx); idx+=1; delidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[1])-nligation #VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
			varlist[var][8] = False	#VAR INS
			varlist[var][9] = None	#VAR INS SIZE
		if line[10] == "TYPE:INVERSION":
			mid="INV_"+str(invidx); var="VAR_"+str(idx); idx+=1; invidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[1])-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	line[0], int(line[1]), int(line[4])-int(line[1])+1, 1, True)	#VAR INS SIZE
		if line[10] == "TYPE:DUPLICATION": #something weird with lumpy's output here, start is bigger than end
			mid="DUP_"+str(dupidx); var="VAR_"+str(idx); idx+=1; dupidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[4])-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[1])-int(line[4])+1 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	line[0], int(line[4]), int(line[1])-int(line[4])+1, 2, False)	#VAR INS SIZE
		if line[10] == "TYPE:TRANSLOCATION":
			line2 = bp.next() #this can be ingnored
			if int(line[2])-int(line[1])==1:#where insertion happens
				delst=int(line[4]); deled=int(line[5]); delchr=line[3]; insst=int(line[1]); insted=int(line[2]); inschr=line[0]
			else:
				delst=int(line[1]); deled=int(line[2]); delchr=line[0]; insst=int(line[4]); insted=int(line[5]); inschr=line[3]
			#do insertion
			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; 
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = inschr		#VAR CHR
			varlist[var][5] = insst-nligation	#VAR POS
			varlist[var][6] = False	#VAR DEL
			varlist[var][7] = 0 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	delchr, delst, deled-delst+1, 1, False)	#VAR INS SIZE
			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; traidx+=1
			#do deletion
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = delchr		#VAR CHR
			varlist[var][5] = delst-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = deled-delst+1 #VAR DEL SIZE
			varlist[var][8] = False	#VAR INS
			varlist[var][9] = None	#VAR INS SIZE
	return varlist

def file2var(bun):
	ifile=CommentedFile(bun.varfile); hapseq=bun.hapseq; nligation=bun.nligation
	#NOTE: CommentedFile can only be iterated once. Need seek to reiterate
	varlist=collections.OrderedDict()
	for row in ifile: 
		varlist[row[0]]=row;
		#print row ,varlist
	for var in varlist.keys():
		varlist[var][2] = float(varlist[var][2])
		varlist[var][3] = int(varlist[var][3])
		varlist[var][5] = int(varlist[var][5])-nligation
		varlist[var][6] = varlist[var][6] in ['True']
		varlist[var][7] = int(varlist[var][7])
		varlist[var][8] = varlist[var][8] in ['True']
		varlist[var][9] = readins(varlist[var][9])
	return varlist

def readins(cell):
	if cell == 'None': return None;
	tmp = cell.split(','); chr = tmp[1].split(':')[0]; copy = int(tmp[2])
	start = int(tmp[1].split(':')[1].split('-')[0])-1; span = int(tmp[1].split(':')[1].split('-')[1])-start
	return [ tmp[0], chr, start, span, copy, True if tmp[3] == 'r' else False ]

def formatins(cell):
	#cell = ['example.fasta', '22', 25202000, 1000, 2, False] #0-based index cell[1] and 1-based cell[2]
	#print cell
	if not cell: return 'None';
	revcomp = 'r' if cell[5] else 'f'
	coord = cell[1]+':'+str(cell[2]+1)+'-'+str(cell[2]+cell[3])
	return ','.join([cell[0],coord,str(cell[4]),revcomp])

def var2file(bun):
	#IMPORTANT: OrderedDict and List are all mutable
	varcopy=copy.deepcopy(bun.varlist); nligation=bun.nligation; ofile=csv.writer(open(bun.oprefix+'.out.var','w'),delimiter='\t')
	ofile.writerow(varfile_names)
	for var in varcopy.keys():
		varcopy[var][9]=formatins(varcopy[var][9])
	for var in varcopy.keys():
		varcopy[var][5]=varcopy[var][5]+nligation #actual site
	ofile.writerows(varcopy.values())
	return None

def var2bed(bun):
	#IMPORTANT: OrderedDict and List are all mutable
	if len(bun.varlist) == 0: return None
	varcopy=copy.deepcopy(bun.varlist); nligation=bun.nligation; 
	for var in varcopy.keys():
		varcopy[var][9]=formatins(varcopy[var][9])
		varcopy[var][5]=varcopy[var][5]+nligation #actual site
	#print varcopy
	varbed=pybedtools.BedTool("\n".join([ str(varcopy[var][4])+"\t"+str(varcopy[var][5])+"\t"+str(varcopy[var][5]+varcopy[var][7])+"\t"+"\t".join([str(x) for x in varcopy[var]]) for var in varcopy.keys() ]), from_string=True)
	#print varbed
	varbed=varbed.sort()
	varbed.saveas(bun.oprefix+'.out.bed',trackline='\t'.join(varbed_names))
	#varbed.saveas(bun.oprefix+'.bed',trackline='\t'.join(['#CHROM','START','END','VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)']))
	#ofile=csv.writer(open(bun.oprefix+'.bed','w'),delimiter='\t')
	#ofile.writerow(['#CHROM','START','END','VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)'])
	#for var in varcopy.keys():
	#	ofile.writerow([varcopy[var][4],varcopy[var][5],varcopy[var][5]+varcopy[var][7]+1]+varcopy[var])
	return varbed

def makevar(bun,metatab,fas=False):
	#metatab = bun.metatab
	#FIXME: this has to be modified to disallow streteching Ns and random selection in similar size
	gaptab = bun.gaptab
	hapseq = bun.hapseq
	nploid = bun.nploid
	plansize = bun.plansize
	nligation = bun.nligation
	burnin = bun.burnin
	varlist = bun.varlist
	varcnt = bun.varcnt
	#print >>sys.stderr, "input varcnt:", varcnt												#total number of input vars
	trunks = maketrunks(bun)	 													#trunks is in bed format
	trunks,stretches = blocktrunks(trunks,gaptab)				#block trunks that not usable
	trunks=trunks.to_dataframe()
	# seqfile,sizestring #random seq from seqfile using sizestring
	# seqfile #use seqfile exactly
	metacnt={}; delvec=set(); finsvec=set(); dinsvec=set(); dupvec=set(); travec=set(); idupvec=set(); invvec=set(); itravec=set(); idinsvec=set();
	for var in varlist:
		meta=varlist[var][1].split("_")[0]
		if meta == "DEL":
			delvec.add(varlist[var][1])
		elif meta == "FINS":
			finsvec.add(varlist[var][1])
		elif meta == "DINS":
			dinsvec.add(varlist[var][1])
		elif meta == "IDINS":
			idinsvec.add(varlist[var][1])
		elif meta == "DUP":
			dupvec.add(varlist[var][1])
		elif meta == "TRA":
			travec.add(varlist[var][1])
		elif meta == "IDUP":
			idupvec.add(varlist[var][1])
		elif meta == "INV":
			invvec.add(varlist[var][1])
		elif meta == "ITRA":
			itravec.add(varlist[var][1])
	metacnt['DEL']=len(delvec);
	metacnt['FINS']=len(finsvec);
	metacnt['DINS']=len(dinsvec);
	metacnt['DUP']=len(dupvec);
	metacnt['TRA']=len(travec);
	metacnt['IDUP']=len(idupvec);
	metacnt['INV']=len(invvec);
	metacnt['IDINS']=len(idinsvec);
	metacnt['ITRA']=len(itravec);
	
	rows=[]
	for row in metatab:
		if fas: row[4:8]=['fix_1']*4	#override by 'fix_1' for fasforge
		rows.append(row)
		
	rowix=range(0,len(rows))
	#print rowix, rows
	#random.shuffle(rowix) #randomized between metatypes
	
	for r in rowix:
		sizevec=[]; copyvec=[]; donorvec=[]; recptvec=[]; donorfreq=[]; recptfreq=[]
		row = rows[r]; metatype = row[0]; cnt = int(row[1])
		
		if metatype != 'FINS':
			sizedist = expanddist(row[2],cnt)
		else:
			sizeopt = row[2].split(',')
			foreignseq = pysam.Fastafile(sizeopt[0])
			if len(sizeopt)==1: #use exact size
				randseq = False
				if(foreignseq.nreferences!=cnt):
					raise Exception("number of foreign sequences provided not equal to specified in FINS!")
				else:
					sizedist = CountUnique(foreignseq.lengths)
			else:
				randseq = True
				sizedist = expanddist(sizeopt[1],cnt)

		#distributions are marginal/independent and mixed/combinatorics
		copydist = expanddist(row[3],cnt)
		donordist = expanddist(row[4],cnt)
		recptdist = expanddist(row[5],cnt)
		donorfdist = expanddist(row[6],cnt)
		recptfdist = expanddist(row[7],cnt)
		for key in sizedist.keys():
			sizevec.extend([key]*sizedist[key])
		for key in copydist.keys():
			copyvec.extend([key]*copydist[key])
		for key in donordist.keys():
			donorvec.extend([key]*donordist[key])
		for key in recptdist.keys():
			recptvec.extend([key]*recptdist[key])
		for key in donorfdist.keys():
			donorfreq.extend([key]*donorfdist[key])
		for key in recptfdist.keys():
			recptfreq.extend([key]*recptfdist[key])
		### randomized within metatypes ###
		random.shuffle(copyvec)
		random.shuffle(donorvec)
		random.shuffle(recptvec)
		random.shuffle(donorfreq)
		random.shuffle(recptfreq)
		assert(len(sizevec)==cnt)
		assert(len(copyvec)==cnt)
		assert(len(donorvec)==cnt)
		assert(len(recptvec)==cnt)
		assert(len(donorfreq)==cnt)
		assert(len(recptfreq)==cnt)
			
		for i in xrange(0,cnt): #this could be parallelized? yes if we preallocate trunks
			metacnt[metatype] = metacnt.get(metatype,0) + 1
			mid = metatype + "_" + str(metacnt[metatype])
			donor = demultiplex(donorvec[i]) #this is a list of donor haplos, 0-indexed, empty if none
			recpt = demultiplex(recptvec[i]) #this is a list of recpt haplos, 0-indexed, empty if none
			dpick = random.sample(donor,1)[0] #randomly choose a donor sequence
			vartag = 'VAR'
			
			try:
				if metatype == 'DEL':		#this should be only operations to donors
					dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					for d in donor:
						varcnt += 1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,donorfreq[i],d,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1],True,sizevec[i],False,None]
				elif metatype == 'TRA' or metatype == 'ITRA':	#one deletion at donor followed by a insertion at recepient
					rc = True if metatype == 'ITRA' else False
					dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					for d in donor:
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,donorfreq[i],d,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1],True,sizevec[i],False,None]
					rtid,stretches=maketid(stretches,1)
					for r in recpt:		#inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],False,0,True,\
							[hapseq[dpick].filename,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1]+nligation,sizevec[i],copyvec[i],rc]]
				elif metatype == 'INV' or metatype == 'IDUP':	#one deletion at donor followed by insertion of reverse
					dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					for d in donor:
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,donorfreq[i],d,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1],True,sizevec[i],True,
							[hapseq[dpick].filename,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1]+nligation,sizevec[i],copyvec[i],True]]
				elif metatype == 'DINS':
					dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					rtid,stretches=maketid(stretches,1)
					for r in recpt:		#inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],False,0,True,\
							[hapseq[dpick].filename,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1]+nligation,sizevec[i],copyvec[i],False]]
				elif metatype == 'IDINS':
					dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					rtid,stretches=maketid(stretches,1)
					for r in recpt:		#inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],False,0,True,\
							[hapseq[dpick].filename,str(trunks.iloc[dtid,0]),trunks.iloc[dtid,1]+nligation,sizevec[i],copyvec[i],True]]
				elif metatype == 'DUP':
					rtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
					for r in recpt:		#inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],True,sizevec[i],True,\
							[hapseq[r].filename,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1]+nligation,sizevec[i],copyvec[i],False]]
				elif metatype == 'FINS':
					#print donor, recpt, dpick
					rtid,stretches=maketid(stretches,1)
					#print foreignseq.filename
					#print foreignseq.references[i]
					for r in recpt:		#inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
						varcnt+=1
						varid = vartag+"_"+str(varcnt)
						#print foreignseq[i]
						if randseq == False:
							varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],False,0,True,\
								[foreignseq.filename,foreignseq.references[i],0,foreignseq.lengths[i],copyvec[i],False]]
						else:
							rvec = [ j for j in xrange(0,len(foreignseq.references)) if foreignseq.lengths[j] > sizevec[i] ]
							ri = random.sample( rvec, 1 )[0]
							si = random.sample( xrange(0,foreignseq.lengths[ri]-sizevec[i]), 1 )[0]
							rc = random.sample( [True, False], 1 )[0]
							varlist[varid]=[varid,mid,recptfreq[i],r,str(trunks.iloc[rtid,0]),trunks.iloc[rtid,1],False,0,True,\
								[foreignseq.filename,foreignseq.references[ri],si,sizevec[i],copyvec[i],rc]]
				else:
					raise Exception("unknown metatype %s" % metatype)
			except IndexError:
				raise Exception("index out of boundry, best bet: insufficient haplo sequences provided than specified in .meta/.var")

	print >>sys.stderr, "output varcnt", varcnt
	return varlist,varcnt #synonym to varlist

def blocktrunks(trunks,gaptab):
#	for gap in gaptab:
#		print gap
#	for trunk in trunks:
#		print trunk
	remains = trunks.intersect(gaptab,v=True) #only keeping non overlapping trunks
#	for remain in remains:
#		print remain
	#FIXME: sometimes bedtool error will rise if temporary files were deleted, how to fix in safeway?
	nr = len(remains)
	cnt=1; stretches = [0] * nr; j=0;
	for feat in remains:
		if j == 0:
			j += 1; pfeat=feat; continue
		if feat.chrom == pfeat.chrom and feat.start == pfeat.end:
			cnt += 1 #continue counting from j-1 to j
		else:		#stop counting, i.e. j-1 in old, j is new
			stretches[j-cnt] = cnt # e.g. s[0]=2, s[1]=0, s[2]=1
			for k in xrange(j-cnt+1,j): stretches[k] = 0 #marking a stretch by its first
			cnt=1	#start new counting
		if j == nr-1 : #this is forced ending, we need to trace back (cnt-1)
			stretches[j-cnt+1]=cnt
			for k in xrange(j-cnt+2,j): stretches[k] = 0
		j += 1; pfeat=feat
	return remains,stretches
	#remains: trunks that remains, BED formatted
	#stretches: if a trunk remain, how far it stretches;
	#					 counting from 1

def maketrunks(bun):
	#one trunk will be [ sv_pos-buffer, sv_pos+maxsv+buffer ]
	#trunks are stored by [ {1:(chr1,sv_pos),2:(chr2,sv_pos),...}, {1:(chr1,sv_pos),2:(chr2,sv_pos),...}, ... ]
	#all converted to bed format
	hapseq=bun.hapseq
	plansize=bun.plansize
	nligation=bun.nligation
	burnin=bun.burnin
	trunks = [] # trunk[i].start = sv[i].pos - nligation; trunk[i+1] = trunk[i].start + plansize + 2*ligation
	for i in xrange(0,min(hapseq[0].nreferences,maxrefnum)):
		pos=xrange(burnin,hapseq[0].lengths[i]-plansize-2*nligation-burnin,plansize+2*nligation)
		trunks.extend([ '%s\t%s\t%s' % (hapseq[0].references[i], str(p), str(p+plansize+2*nligation)) for p in pos ])
	trunks = pybedtools.BedTool("\n".join([ x for x in trunks ]), from_string=True) #convert to BED
	return trunks

def maketid(stretches,varstretch):
	#print stretches
	#print varstretch
	#print len(stretches)
	tids = [t for t, x in enumerate(stretches) if x >= varstretch] #always random pick among min qualifying tids
	if tids == []: raise Exception("cannot allocate genome region for metafile, please reduce plansize!");
	if varstretch >= 1: #could from any place that has a stretch >= varstretch
		mtids = [varstretch,tids]
	else:
		mtids = locatemin([stretches[i] for i in tids]) #this is not currently used
	#print tids
	#print mtids
	tid = random.sample(mtids[1],1)[0]	 #locate the tid
	vtid = stretches[tid]
	#print tid
	for i in xrange(tid,tid+varstretch): stretches[i]=0;						#mark head of stretch used
	if vtid-varstretch>0: stretches[tid+varstretch]=vtid-varstretch; #mark remaining stretch by its head
	return tid, stretches

#maketid tests
#>>> ts=[2,0,1,1,2,0,3,0,0,1]
#>>> s=copy.deepcopy(ts)
#>>> maketid(s,1)
#(9, [2, 0, 1, 1, 2, 0, 3, 0, 0, 0])
#>>> maketid(s,1)
#(3, [2, 0, 1, 0, 2, 0, 3, 0, 0, 0])
#>>> maketid(s,1)
#(6, [2, 0, 1, 0, 2, 0, 0, 2, 0, 0])
#>>> maketid(s,1)
#(4, [2, 0, 1, 0, 0, 1, 0, 2, 0, 0])
#>>> maketid(s,1)
#(0, [0, 1, 1, 0, 0, 1, 0, 2, 0, 0])
#>>> maketid(s,1)
#(1, [0, 0, 1, 0, 0, 1, 0, 2, 0, 0])
#>>> maketid(s,1)
#(2, [0, 0, 0, 0, 0, 1, 0, 2, 0, 0])
#>>> maketid(s,1)
#(5, [0, 0, 0, 0, 0, 0, 0, 2, 0, 0])
#>>> maketid(s,1)
#(7, [0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
#>>> maketid(s,1)
#(8, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
#>>> s=copy.deepcopy(ts)
#>>> maketid(s,2)
#(4, [2, 0, 1, 1, 0, 0, 3, 0, 0, 1])
#>>> maketid(s,2)
#(6, [2, 0, 1, 1, 0, 0, 0, 0, 1, 1])
#>>> maketid(s,2)
#(0, [0, 0, 1, 1, 0, 0, 0, 0, 1, 1])
#>>> s=copy.deepcopy(ts)
#>>> maketid(s,3)
#(6, [2, 0, 1, 1, 2, 0, 0, 0, 0, 1])

def expanddist(diststring, cnt):
	dist=diststring.split('_')
	uniques = {}
	new_cnt = 0
	if dist[0] == 'fix':
		values = [ int(x) if isint(x) else float(x) for x in dist[1:len(dist)] ]
		uniques = CountUnique(values)
		sorted_uniques = sorted(uniques.items(), key=operator.itemgetter(1)) 
		#increasing order sorting of values
		nvalue = len(values)
		nunique = len(uniques)
		amp = cnt/nvalue
		rem = cnt-amp*nvalue
		for item in sorted_uniques:
			key = item[0]
			if rem>0:
				delta = uniques[key] if rem-uniques[key]>=0 else rem
				uniques[key] = uniques[key] * amp + delta
				rem = rem - delta
			else:
				delta = 0
				uniques[key] = uniques[key] * amp
			new_cnt = new_cnt + uniques[key]
		assert new_cnt == cnt

	elif dist[0] == 'nori':
		#create cnt normal integer; 
		mean = int(dist[1]) 
		sd = int(dist[2])
		values = [ random.gauss(mean, sd) for i in xrange(0,cnt) ]
		values = [ mean if v <= 0 else int(v) for v in values ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	elif dist[0] == 'norf':
		#create cnt normal floats; 
		mean = float(dist[1]); sd = float(dist[2])
		values = [ random.gauss(mean, sd) for i in xrange(0,cnt) ]
		values = [ mean if v<=0 else v for v in values ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	elif dist[0] == 'unii':
		#create cnt uniform integer; 
		low = int(dist[1]); up = int(dist[2])
		values = [ random.uniform(low, up) for i in xrange(0,cnt) ]
		values = [ int(v) for v in values ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	elif dist[0] == 'unif':
		#create cnt uniform floats; 
		low = float(dist[1]); up = float(dist[2])
		values = [ random.uniform(low, up) for i in xrange(0,cnt) ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	elif dist[0] == 'expi':
		#create cnt uniform integer; 
		lam = float(dist[1]);
		values = [ random.expovariate(lam) for i in xrange(0,cnt) ]
		values = [ int(1/lam) if v<=0 else int(v) for v in values ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	elif dist[0] == 'expf':
		#create cnt uniform float; 
		lam = float(dist[1]);
		values = [ random.expovariate(lam) for i in xrange(0,cnt) ]
		values = [ 1/lam if v<=0 else v for v in values ]
		uniques = CountUnique(values)
		assert new_cnt == cnt

	else:
		raise Exception("can't expand specified distribution!")

	return uniques

def run():
	parser = argparse.ArgumentParser(description="spike-in mutations to fasta/bam files, \
																									for bam files first assemble haplo-types \
																									for fasta files we segment fasta files \
																									to confirm if the configuration is viable")
	parser.add_argument('hapfiles', metavar='hapfiles', #example.fasta
												help="inputfile of hap file(s), required")
	parser.add_argument('gapfile', metavar='gapfile', type=argparse.FileType('rU'), #example.bed
												help="inputfile of gap or other banned regions, required")
	parser.add_argument('parfile', metavar='parfile', type=argparse.FileType('rU'), #example.par
												help="inputfile of library-wise parameters, required"),
	parser.add_argument('reffile', metavar='reffile', type=argparse.FileType('rU'), #example.fasta
												help="inputfile of ref file(s), required"),
	parser.add_argument('-m', '--metafile', dest='metafile', type=argparse.FileType('rU'), #example.meta
												help="inputfile of meta specification")
	parser.add_argument('-v', '--varfile', dest='varfile', type=argparse.FileType('rU'),	#example.var/op
												help="inputfile of variant specification")
	parser.add_argument('-b', '--burnin', dest='burnin', type=int, default=1000000,
												help="burnin is the chromosome tip skip size for genome planning to avoid tip padding Ns, default: %(default)s bp")
	parser.add_argument('-e', '--edgein', dest='edgein', type=int, default=100000,
												help="edgein is the gap neibourgh region skip size for genome planning, default: %(default)s bp")
	parser.add_argument('-i', '--bamfile', dest='inbamfilenames', #in1.bam,in2.bam,in3.bam
												help="inputfile of library-wise bam file(s)")
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #output1,output2,output3,...
												help="prefix for output files")
	parser.add_argument('-y', '--nploid', dest='nploid', default=2, type=int,
												help="genome ploidity default= %(default)s ploid")
	parser.add_argument('-f', '--plansize', dest='plansize', default=100000, type=int,
												help="palnsize + ligation*2 is the grain size for genome planning, default: %(default)s bp)")
	parser.add_argument('-g', '--nligation', dest='nligation', default=2000, type=int,
												help="genome trunk ligation (default: %(default)s bp)")
	parser.add_argument('-n', '--nprocs', dest='nprocs', default=1, type=int,
												help="split into multiple processes (default: %(default)s )")
	parser.add_argument('-r', '--mergemax', dest='mergemax', default=4000, type=int,
												help="max number of bamfiles in one iteration of samtool merging (default: %(default)s )")
	parser.add_argument('-d', '--tmpdir', dest='tmpdir', default=os.path.join(os.environ['HOME'],'svetmp_'+''.join([random.choice(string.ascii_letters) for i in range(4)])),
												help="root dir for keeping temporary files, default (last 4 digits random): %(default)s")
	parser.add_argument('--layout', action='store_true', dest='layout', default=False, help="dry run to layout")
	parser.add_argument('--outbam', action='store_true', dest='outbam', default=False, help="wet run to bam")
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
	print >>sys.stderr, "total runtime: ", time.time()-starttime, " seconds"

#def var2vcf(bun):
#	varlist=bun.varlist; ofile=[ csv.writer(open(bun.oprefix+'.'+str('.vcf','w'),delimiter='\t')
#	vcf_meta = call_meta + sample_meta + contig_meta
#	ofile.writerow(vcf_meta)
#	ofile.writerow(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"])
#	metalist={}
#	for var in varlist.keys():
#		if not metalist.get(varlist[var][1],False):
#			metalist[varlist[var][1]].append(varlist[var])
#		else:
#			metalist[varlist[var][1]]=varlist[var]
#	for mid in metalist.keys():
#		sv_type=''.join([i for i in mid if not i.isdigit()])
#		if len(metalist[mid])
	#NOTE: only when reffile == hapfiles vcf file will make sense
	#NOTE: need to go over by MID and report correspondingly

#meta_keys=["fileformat","fileDate","source","reference","phasing"]
#meta_values=["None","None","None","None","partial"]
#info_values=[
#'<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
#'<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
#'<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
#'<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
#'<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
#'<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
#'<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
#'<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">',
#'<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
#'<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
#'<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
#'<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
#'<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
#'<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
#'<ID=SOMATIC,Number=0,Type=Flag,Description="SOMATIC structural variation">',
#'<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">',
#'<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
#'<ID=METRANS,Number=.,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">',
#'<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
#'<ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">',
#'<ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVAR">',
#'<ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">',
#'<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">',
#'<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">',
#'<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
#'<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the length of the inserted material between breakends">',
#'<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">',
#'<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">',
#'<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of adjacency">',
#'<ID=CICN,Number=2,Type=Integer,Description="Confidence interval around copy number for the segment">',
#'<ID=CICNADJ,Number=.,Type=Integer,Description="Confidence interval around copy number for the adjacency">'
#]
#format_values=[
#'<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',
#'<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">',
#'<ID=CNL,Number=.,Type=Float,Description="Copy number genotype likelihood for imprecise events">',
#'<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel with respect to the genome ancestor">',
#'<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">',
#'<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">',
#'<ID=GT,Number=1,Type=String,Description="Genotype">',
#'<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
#'<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
#'<ID=AD,Number=1,Type=Integer,Description="Read Supporting ALT">',#TCGA
#'<ID=BQ,Number=1,Type=Integer,Description="Base Quality ofr Read Supporting ALT">',#TCGA
#'<ID=SS,Number=1,Type=Integer,Description="Somatic Status 0=wt,1=gm,3=so,4=LOH,5=PTM/unknown">',#TCGA
#'<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
#'<ID=GL,Number=G,Type=Integer,Description="Genotype Likelihood">',
#'<ID=lW,Number=1,Type=Integer,Description="Max Coverage Likelihood Ratio">',
#'<ID=lWc,Number=1,Type=Integer,Description="Cumulative Max Coverage Likelihood Ratio">',
#'<ID=lCd,Number=1,Type=Integer,Description="Deletion Insert Size Likelihood Ratio">',
#'<ID=lCi,Number=1,Type=Integer,Description="Insertion Insert Size Likelihood Ratio">',
#'<ID=lDl,Number=1,Type=Integer,Description="Left/Plus Anchored Hang Read Likelihood Ratio">',
#'<ID=lDr,Number=1,Type=Integer,Description="Right/Minus Anchored Hang Read Likelihood Ratio">',
#'<ID=RGl,Number=1,Type=Integer,Description="Left Limit of Voted Region">',
#'<ID=RGr,Number=1,Type=Integer,Description="Right Limit of Voted Region">',
#'<ID=BPl,Number=1,Type=Integer,Description="Softclipped Reads Left Break Point">',
#'<ID=BPr,Number=1,Type=Integer,Description="Softclipped Reads Right Break Point">'
#]
#alt_values=[
#'<ID=DEL,Description="Deletion">',
#'<ID=DEL:ME:ALU,Description="Deletion of ALU element">',
#'<ID=DEL:ME:L1,Description="Deletion of L1 element">',
#'<ID=DUP,Description="Duplication">',
#'<ID=DUP:TANDEM,Description="Tandem Duplication">',
#'<ID=INS,Description="Insertion of novel sequence">',
#'<ID=INS:ME:ALU,Description="Insertion of ALU element">',
#'<ID=INS:ME:L1,Description="Insertion of L1 element">',
#'<ID=INV,Description="Inversion">',
#'<ID=CNV,Description="Copy number variable region">'
#]
#filter_values=[
#'<ID=q10,Description="Quality below 10">',
#'<ID=s50,Description="Less than 50% of samples have data">'
#]
#info_keys= ["INFO"] * len(info_values)
#format_keys= ["FORMAT"] * len(format_values)
#alt_keys= ["ALT"] * len(alt_values)
#filter_keys= ["FILTER"] * len(filter_values)
#sv_keys= info_keys + format_keys + alt_keys + filter_keys
#sv_values= info_values + format_values + alt_values + filter_values
#
#def sample_meta(id, genomes, mixture, desc):
#	template=['##SAMPLE=<ID=%s,Genomes=%s,Mixture=%s,Description="%s">'] * len(id)
#	return "\n".join( [ template[i] % (id[i], genomes[i], mixture[i], desc[i]) for i in range(0,len(id)) ] )
#
#def contig_meta(id, length, md5, species):
#	template=['##contig=<ID=%s,length=%s,md5=%s,species="%s">'] * len(id)
#	return "\n".join( [ template[i] % (id[i], length[i], md5[i], species[i]) for i in range(0,len(id)) ] )
#
#def call_meta(keys, values, program="mutforge"):
#	keys.append(["fileformat","fileDate","source"])
#	values.append(["VCFv4.1",time.ctime(),program)])
#	return "\n".join([ '##'+key+'='+value for key,value in izip(keys,values) ])

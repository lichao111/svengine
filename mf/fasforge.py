#!/usr/bin/env python
#import svengine.mf.mutforge as mf

import copy, time, re, os, shutil, sys, random, csv, subprocess, traceback, argparse, tempfile, itertools, operator, collections, ConfigParser, json, imp, StringIO, string
from multiprocessing import Pool
import pysam
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import pybedtools.featurefuncs
import shelve
#properly import muforge
try:
	import mf.mutforge as mf #in debugging
except ImportError:
	import mutforge as mf #in production
starttime=time.time()

#def lumpy2var(bp,ref,nligation):
#	idx=1; delidx=1; insidx=1; traidx=1; invidx=1; dupidx=1; varlist=collections.OrderedDict()
#	for line in bp:
#		if line[10] == "TYPE:DELETION":
#			mid="DEL_"+str(delidx); var="VAR_"+str(idx); idx+=1; delidx+=1
#			varlist[var]=[None] * 10
#			varlist[var][0] = var	#VID
#			varlist[var][1] = mid	#MID
#			varlist[var][2] = 1.0 #VAR FREQ
#			varlist[var][3] = 0		#VAR HAP
#			varlist[var][4] = line[0]		#VAR CHR
#			varlist[var][5] = int(line[1])-nligation #VAR POS
#			varlist[var][6] = True	#VAR DEL
#			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
#			varlist[var][8] = False	#VAR INS
#			varlist[var][9] = None	#VAR INS SIZE
#		if line[10] == "TYPE:INVERSION":
#			mid="INV_"+str(invidx); var="VAR_"+str(idx); idx+=1; invidx+=1
#			varlist[var]=[None] * 10
#			varlist[var][0] = var	#VID
#			varlist[var][1] = mid	#MID
#			varlist[var][2] = 1.0 #VAR FREQ
#			varlist[var][3] = 0		#VAR HAP
#			varlist[var][4] = line[0]		#VAR CHR
#			varlist[var][5] = int(line[1])-nligation	#VAR POS
#			varlist[var][6] = True	#VAR DEL
#			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
#			varlist[var][8] = True	#VAR INS
#			varlist[var][9] = (ref.name,	line[0], int(line[1]), int(line[4])-int(line[1])+1, 1, True)	#VAR INS SIZE
#		if line[10] == "TYPE:DUPLICATION": #something weird with lumpy's output here, start is bigger than end
#			mid="DUP_"+str(dupidx); var="VAR_"+str(idx); idx+=1; dupidx+=1
#			varlist[var]=[None] * 10
#			varlist[var][0] = var	#VID
#			varlist[var][1] = mid	#MID
#			varlist[var][2] = 1.0 #VAR FREQ
#			varlist[var][3] = 0		#VAR HAP
#			varlist[var][4] = line[0]		#VAR CHR
#			varlist[var][5] = int(line[4])-nligation	#VAR POS
#			varlist[var][6] = True	#VAR DEL
#			varlist[var][7] = int(line[1])-int(line[4])+1 #VAR DEL SIZE
#			varlist[var][8] = True	#VAR INS
#			varlist[var][9] = (ref.name,	line[0], int(line[4]), int(line[1])-int(line[4])+1, 2, False)	#VAR INS SIZE
#		if line[10] == "TYPE:TRANSLOCATION":
#			line2 = bp.next() #this can be ingnored
#			if int(line[2])-int(line[1])==1:#where insertion happens
#				delst=int(line[4]); deled=int(line[5]); delchr=line[3]; insst=int(line[1]); insted=int(line[2]); inschr=line[0]
#			else:
#				delst=int(line[1]); deled=int(line[2]); delchr=line[0]; insst=int(line[4]); insted=int(line[5]); inschr=line[3]
#			#do insertion
#			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; 
#			varlist[var]=[None] * 10
#			varlist[var][0] = var	#VID
#			varlist[var][1] = mid	#MID
#			varlist[var][2] = 1.0 #VAR FREQ
#			varlist[var][3] = 0		#VAR HAP
#			varlist[var][4] = inschr		#VAR CHR
#			varlist[var][5] = insst-nligation	#VAR POS
#			varlist[var][6] = False	#VAR DEL
#			varlist[var][7] = 0 #VAR DEL SIZE
#			varlist[var][8] = True	#VAR INS
#			varlist[var][9] = (ref.name,	delchr, delst, deled-delst+1, 1, False)	#VAR INS SIZE
#			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; traidx+=1
#			#do deletion
#			varlist[var]=[None] * 10
#			varlist[var][0] = var	#VID
#			varlist[var][1] = mid	#MID
#			varlist[var][2] = 1.0 #VAR FREQ
#			varlist[var][3] = 0		#VAR HAP
#			varlist[var][4] = delchr		#VAR CHR
#			varlist[var][5] = delst-nligation	#VAR POS
#			varlist[var][6] = True	#VAR DEL
#			varlist[var][7] = deled-delst+1 #VAR DEL SIZE
#			varlist[var][8] = False	#VAR INS
#			varlist[var][9] = None	#VAR INS SIZE
#	return varlist

def var2seq(brpt,reffile,nligation,oprefix): #assume varbed is checked for nonoverlapping
	#print brpt,reffile,nligation,oprefix
	ref=pysam.Fastafile(reffile.name); seq=[]
	for i in xrange(0,ref.nreferences): #iterate through
		rname=ref.references[i]; rlen=ref.lengths[i]; segst=0;
		rseq=Seq(ref.fetch(rname,0,rlen), generic_dna)
		if brpt.get(rname,None) == None: #no var at all
			seg = [ rseq ]
		else: #1+ vars
			brs = collections.OrderedDict(sorted(brpt[rname].items(), key=lambda t: t[0]))
			seg = []; nbrs=len(brs) #member of seg are Seqs
			#if i==1: print nbrs
			for j in xrange(0,nbrs):
				#VID MID VARFRQ VARHAP VARCHR VARPOS VARDEL VARDELSZ VARINS VARINSSZ
				try:
					var=brs.values()[j]
					if j == (nbrs-1): #last var, seged is 
						seged=rlen
					else:
						seged=brs.values()[j+1][5]
					if j == 0: #the first seg is appended without var
						seg.append(rseq[segst:var[5]]); segst=var[5] 
					seq1 = rseq[segst:seged] #seg if connected is the new seq 
					if var[6]: #work on deletion
						seq1 = seq1[:nligation]+seq1[nligation+var[7]:]
						segst = seged
					if var[8]: #work on insertion
						#['example.fasta', '22', 25202000, 1000, 2, False]
						info = var[9]
						seqi = Seq(pysam.Fastafile(info[0]).fetch(info[1],info[2],info[2]+info[3]), generic_dna) #st=0-based, ed=0-based+1
						if info[5]:
							seqi = seqi.reverse_complement();
						seqi = Seq(str(seqi)*info[4], generic_dna)
						seq1 = seq1[:nligation]+seqi+seq1[nligation:]
						segst = seged
					seg.append(seq1)
				except TypeError:
					print var
					print >>sys.stderr,"most likely formatting wrong in provided .var file, ignore entry!"
		#print "old:",len(rseq)
		#print "new:",len(SeqRecord( Seq(''.join([str(s) for s in seg]), generic_dna), id=rname, name=rname, description=oprefix+"."+rname ) )
		seq.append( SeqRecord( Seq(''.join([str(s) for s in seg]), generic_dna), id=rname, name=rname, description=oprefix+"."+rname ) )
	SeqIO.write(seq,oprefix+".fnv",'fasta')
	return { 'filename':oprefix+".fnv", 'seq':seq }

def main(args):
	burnin = args.burnin
	edgein = args.edgein
	outbam = args.outbam
	nprocs = args.nprocs
	varfile = args.varfile
	metafile = args.metafile
	reffile = args.reffile
	gapfile = args.gapfile
	parfile = args.parfile
	parlist = ConfigParser.ConfigParser()
	parlist.readfp(parfile)
	plansize = args.plansize
	mergemax = args.mergemax
	layout = args.layout
	tempfile.tempdir = args.tmpdir
	targetfile=args.targetfile
	tar_len = 0
	try:
		os.mkdir(tempfile.tempdir)
	except:
		print >>sys.stderr, tempfile.tempdir, "error!, either it exists or cann't be created! please change path"
		quit()
	print >>sys.stderr, "tmpdir=", tempfile.tempdir
	pybedtools.set_tempdir(tempfile.tempdir)

	runmode=0 #1 meta, 2 var, 3 var+meta
	if varfile!=None and metafile!=None:
		runmode=3;
		oname=varfile.name.rsplit(".")[0]	
		#don't allow mode 3
		print >>sys.stderr, "either a .meta or a .var file can be specified but not both"
		quit()
	elif varfile!=None and metafile==None:
		runmode=2;
		oname=varfile.name.rsplit(".")[0]	
	elif metafile!=None and varfile==None:
		runmode=1;
		oname=metafile.name.rsplit(".")[0]	
	else:
		raise Exception("must provide meta and/or var files!")
	oprefix = oname if not args.oprefix else args.oprefix
	nligation = 0; nploid = 1; #in fasforge, niligation is always 0 and nploid always is 1 
	freq1 = [ float(f) for f in args.mixfreq.split(",") ]
	for f in freq1:
		assert f>=0 and f<=1, "mixfreq must be a float between 0 and 1" 
	try:
		sizefile = mf.fas2size(reffile.name) #create genome size file
		gaptab=pybedtools.BedTool(gapfile.read(),from_string=True)
	#------------------add for targetfile---------------		
		if targetfile != None:
			targfile=pybedtools.BedTool(targetfile.read(),from_string=True)
			targtab=targfile.complement(g=sizefile)
			gaptab=targtab.cat(gaptab.each(pybedtools.featurefuncs.extend_fields,n=targtab.field_count()).saveas())
		gaptab=gaptab.slop(g=sizefile, b=edgein)
		print gaptab.count()
		for i in xrange(gaptab.count()):
		#	print gaptab[i].start
			tar_len += gaptab[i].end-gaptab[i].start
		print tar_len
	except pybedtools.cbedtools.BedToolsFileError:
		raise Exception("incorrect gapfile provided, must in BED format")

	#create a string meta file:
	#FINS	2500	adeno.fna,fix_...	fix_1	fix_1	fix_1	fix_1	fix_1	#FINS=COPY*INS, foreign insertion
	#DINS	2500	fix_...	fix_1	fix_1	fix_1	fix_1	fix_1	#DINS=COPY*INS, domestic insertion
	#inssize="fix_100_200_300_400_500_600_700_800_900_1000_1000_2000_3000_4000_5000_6000_7000_8000_9000_10000"
	#metastr='\n'.join([ "\t".join(['FINS',str(2),'adeno.fna,'+inssize,'fix_1','fix_1','fix_1','fix_1','fix_1' ]),
	#					"\t".join(['DINS',str(2),'fix_100_200','fix_1','fix_1','fix_1','fix_1','fix_1' ]) ])
	varlist = collections.OrderedDict()
	dic={'nligation':nligation,'varlist':varlist,'oprefix':oprefix,'reffile':reffile,'nploid':nploid,'plansize':plansize,'varcnt':0,'edgein':edgein,'gaptab':gaptab,'varfile':varfile,'parfile':parfile,'parlist':parlist,'nprocs':nprocs,'burnin':burnin,'mergemax':mergemax}
	print >>sys.stderr, "runmode=", runmode
	dic['hapseq']=[pysam.Fastafile(reffile.name)];

	if runmode > 1: #have var input
		varintype = varfile.name.split(".")[-1]
		if varintype == 'lumpy':
			bp=mf.CommentedFile(varfile) #blue print
			#bp=csv.reader(varfile,delimiter="\t")
			#somehow we need to set ligation setting 
			#assume that every true var pos out of program is true position
			#every true var pos inside the program is position+nligation
			varlist=mf.lumpy2var(bp,reffile,nligation)
		elif varintype == 'var':
			varlist=mf.file2var(mf.Bunch(dic))
		else:
			raise Exception("varintype %s hasn't been implemented" % varintype)
		dic['varcnt']=len(varlist) #this is needed for continuously counting
		dic['varlist']=varlist #this is needed for continuously counting
		#d.iloc[0:1,0:1]; deoverlapping
		vartab=mf.var2bed(mf.Bunch(dic))
		#print "br1"
		if(len(varlist)>1): #only needed if input is more than 2
				varlist=mf.deovlpvar(dic) 
				dic['varcnt']=len(varlist) 
				dic['varlist']=varlist 
				vartab=mf.var2bed(mf.Bunch(dic))
				assert mf.chkvar(vartab) == True, "==Warn: varfile still contains self overlaps, which is not allowed, bail out"
		#print "br2"
		gaptab=vartab.cat(gaptab.each(pybedtools.featurefuncs.extend_fields,n=vartab.field_count()).saveas())
		dic['gaptab']=gaptab #considerring var tab as pre-existing mask that works for meta-variant 
		#print "br3"

	if runmode%2 == 1: #have meta input, merge var with meta
		metatab = mf.CommentedFile(metafile)
		#for row in metatab: #force input of fix_1 columns 5-8 
		#	row[4:8]=['fix_1']*4
		varlist,varcnt = mf.makevar(mf.Bunch(dic),metatab,fas=True)
		dic['varlist'] = varlist
		dic['varcnt'] = varcnt
		vartab=mf.var2bed(mf.Bunch(dic))
		if(len(varlist)>1): #only needed if input is more than 2
				varlist=mf.deovlpvar(dic) 
				dic['varcnt']=len(varlist) 
				dic['varlist']=varlist 
				vartab=mf.var2bed(mf.Bunch(dic))
				assert mf.chkvar(vartab) == True, "==Warn: varfile still contains self overlaps, which is not allowed, bail out"

	mf.var2file(mf.Bunch(dic))
	#read in bedpe to generate varlist
	#now we need to check whether the	sv has been correctly inserted
	#we also need to check whether the position is properly selected 
	brpt=collections.OrderedDict()
	#print "br4"
	for key in varlist.keys():
		var = varlist[key]
		if brpt.get(var[4],None) == None:
			brpt[var[4]] = collections.OrderedDict({var[5]:var})
		else:
			brpt[var[4]][var[5]]=var	#brpt[chr][pos]
	#print "br5"
	ref = pysam.Fastafile(reffile.name)
#---------------add to dele the target region-------------
	#ref_del = seq2targ(ref,gaptab)
	#reffile_del = file("tar.fnv")
	#ref_del = pysam.Fastafile("tar.fnv")


	seq0 = { 'filename':reffile.name, 'seq':[ SeqRecord(Seq(ref.fetch(ref.references[i],0,ref.lengths[i]) ),id=ref.references[i],name=ref.references[i],description=ref.references[i]) for i in xrange(0,ref.nreferences) ] } #input sequence
	#print "br6"
	seq1 = var2seq(brpt,reffile,nligation,oprefix)
	print >>sys.stderr, "done var2seq", time.time()-starttime, " seconds" 

	nlibrary=parlist.getint('xwgsim', 'nlibrary')
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	#bams = [None] * nlibrary

	if layout:
		print >>sys.stderr, "done layout and fas file generation"
		quit()
	
	pariter=[]
	for lib in range(0,nlibrary):
		for fi in range(0,len(freq1)):
			pariter.append([parlist,nligation,reffile.name,nploid,seq0,lib,freq1[fi],seq0,seq1,lib,oprefix,mergemax,tar_len])

	#NOTE: debug code
	ffq=[]
	if nprocs<2:
		for parit in pariter:
			ffq.append(makefq(parit))
	else:
		pool = Pool(processes=int(nprocs))
		ffq=pool.map(makefq,pariter,1)
	print >>sys.stderr, "done var2fq", time.time()-starttime, " seconds"

	if outbam:
		libbams = mf.fq2bam(ffq,mf.Bunch(dic))
		
		print >>sys.stderr, "done fq2bam", time.time()-starttime, " seconds"

		fbam = [None] * nlibrary
		for lib in range(0,nlibrary): #there will be only one bam in libbams[lib]
			sbam = '.'.join([oprefix,libnames[lib]]) # oprefix.lib.st
			tbam = tempfile.NamedTemporaryFile('w',delete=False).name
			mbam = mf.mergebam([libbams[lib]],tbam,mergemax)
			sort_cmd = "samtools sort -@ %s %s %s" % (nprocs, mbam, sbam)
			#print sort_cmd #
			mf.syscmd(sort_cmd)
			fbam[lib] = sbam+".bam" #post processing: librarywise merging bams;

		print >>sys.stderr, "done merge to final bam:", str(fbam), time.time()-starttime, " seconds"

	shutil.rmtree(tempfile.tempdir) #clean temp files
	print >>sys.stderr, "done fasforge"
	
	return None

#OBSOLETE
def makebam(parit):		
	print "processing lib %s freq %s regions" % (str(parit[5]),str(parit[6]))
	parlist=parit[0]; nligation=parit[1]; rname=parit[2]; nploid=parit[3]; seq0=parit[4]; lib=parit[5]; freq1=parit[6]; freq0=1-freq1; seq0=parit[7]; seq1=parit[8]; lib=parit[9]; oprefix=parit[10]
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	bam0 = mf.runwgs(parlist,nligation,rname,nploid,seq0,lib,freq0)
	bam1 = mf.runwgs(parlist,nligation,rname,nploid,seq1,lib,freq1)
	#print [bam0,bam1]
	bam = mf.mergebam([bam0,bam1],oprefix+"."+libnames[lib])
	print "done lib %s freq %s regions" % (str(parit[5]),str(parit[6]))
	return bam

def makefq(parit):		
	print "processing lib %s freq %s regions" % (str(parit[5]),str(parit[6]))
	parlist=parit[0]; nligation=parit[1]; rname=parit[2]; nploid=parit[3]; lib=parit[5]; freq1=parit[6]; freq0=1-freq1; seq0=parit[7]; seq1=parit[8]; lib=parit[9]; oprefix=parit[10]; mergemax=parit[11];tar_len=parit[12] #seq0=parit[4]; 
	libnames = json.loads(parlist.get('xwgsim', 'libnames'))
	wgs0 = mf.fqwgs(parlist,nligation,rname,nploid,seq0,lib,freq0,mergemax,tar_len)
	wgs1 = mf.fqwgs(parlist,nligation,rname,nploid,seq1,lib,freq1,mergemax,tar_len)
	#bam0 = mf.runwgs(parlist,nligation,rname,nploid,seq0,lib,freq0)
	#bam1 = mf.runwgs(parlist,nligation,rname,nploid,seq1,lib,freq1)
	#print [bam0,bam1]
	ffq = mf.mergefq([wgs0,wgs1],oprefix+"."+libnames[lib],mergemax)
	print "done lib %s freq %s" % (str(parit[5]+1),str(parit[6])), time.time()-starttime, " seconds"
	return ffq

def seq2targ(reffile,gaptab):
	seq=[]
	ref=reffile
	tar_name=[]
	for j in xrange(gaptab.count()):
		while gaptab[j][0]  not in tar_name:
			tar_name.append(gaptab[j][0])
	print tar_name
	for i in xrange(0 , ref.nreferences):
		
		rname=ref.references[i]; rlen=ref.lengths[i]; segst=0;
#		print rname
		rseq=Seq(ref.fetch(rname,0,rlen),generic_dna)
		if rname not in tar_name:
			seg=[rseq]
		else:
			for k in xrange(gaptab.count()):
				print gaptab[k],gaptab[k][0]
				if gaptab[k][0] == rname:
					seg =rseq[:gaptab[k].start]+rseq[gaptab[k].end:]
		seq.append( SeqRecord( Seq(''.join([str(s) for s in seg]), generic_dna), id=rname, name=rname, description='tar'+"."+rname ) )
	SeqIO.write(seq,'tar.fnv','fasta') 
	return { 'filename':"tar.fnv", 'seq':seq }

def run():
	parser = argparse.ArgumentParser(description="spike-in mutations to fasta files, \
																									we insert mutations into a provided reference \
																									this script is haplotype insensitive and ligation less")
	parser.add_argument('gapfile', metavar='gapfile', type=argparse.FileType('rU'), #example.bed
												help="inputfile of gap or other banned regions, required")
	parser.add_argument('parfile', metavar='parfile', type=argparse.FileType('rU'), #example.par
												help="inputfile of library-wise parameters, required"),
	parser.add_argument('reffile', metavar='reffile', type=argparse.FileType('rU'), #example.fasta
												help="inputfile of ref file(s), required"),
	parser.add_argument('-m', '--metafile', dest='metafile', type=argparse.FileType('rU'), #example.meta
												help="inputfile of meta specification")	# variants to be randomly added
	parser.add_argument('-v', '--varfile', dest='varfile', type=argparse.FileType('rU'),	#example.var/op
												help="inputfile of variant specification") # prespecified variants, bedpe format or others
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #output1,output2,output3,...
												help="prefix for output files")
	parser.add_argument('-f', '--plansize', dest='plansize', type=int, default=100000,
												help="palnsize is the grain size for genome planning, exluding ligation buffer (default= %(default)s bp)")
	parser.add_argument('-b', '--burnin', dest='burnin', type=int, default=1000000,
												help="burnin is the chromosome tip skip size for genome planning to avoid tip padding Ns, default: %(default)s bp")
	parser.add_argument('-e', '--edgein', dest='edgein', type=int, default=100000,
												help="edgein is the gap neibourgh region skip size for genome planning, default: %(default)s bp")
	parser.add_argument('-q', '--mixfreq', dest='mixfreq', default="1",
												help="generated spikein-in mixing frequency, homozygous=1")
	parser.add_argument('-n', '--nprocs', dest='nprocs', type=int, default=1,
												help="split into multiple processes (default= %(default)s )")
	parser.add_argument('-r', '--mergemax', dest='mergemax', default=2000, type=int,
												help="max number of files in one iteration of merging (default: %(default)s )") 
	parser.add_argument('-d', '--tmpdir', dest='tmpdir', default=os.path.join(os.environ['HOME'],'svetmp_'+''.join([random.choice(string.ascii_letters) for i in range(4)])),
												help="root dir for keeping temporary files, default (last 4 digits random): %(default)s")
	parser.add_argument('--layout', dest='layout', action='store_true', default=False,
												help="only dry layout")
	parser.add_argument('--outbam', action='store_true', dest='outbam', default=False, 
												help="wet run to bam")
	parser.add_argument('-t','--targetfile',dest='targetfile',type=argparse.FileType('rU'),default=None, help='inputfile of target region, for example Exome region, in bed format')
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
	print "total runtime: ", time.time()-starttime, " seconds"

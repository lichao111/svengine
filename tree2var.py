#!user/bin/env python
import argparse,time,subprocess,csv,tempfile

#try:
#	import mf.mutforge as mf
#except ImportError:
#	import mutforge as mf
start_time=time.time()

class Node:
	def __init__(self,i,j,v,a):#i:int,j:int,v:dict,a:range[0,1]
		self.i=int(i)
		self.j=int(j)
		self.v=v
		self.a=float(a)

#class tree:
#	def __init__(self, f, commentstring="/ "or "\ "):        
#		self.f = csv.reader(f,delimiter='\t')
#		self.commentstring = commentstring
#	def next(self):
#		line = self.f.next()
#		while commentstring in line
#			line = self.f.next()
#		return line
#	def __iter__(self):
#		return self


#def tree2treelist():

def C(Node):
	i=Node.i+1
	j1=2*(Node.j)
	j2=j1+1
	return {Node(i,j1),Node(i,j2)}

def P(Node):	#find Node's parent
	i=Node.i-1
	j=Node.j/2
	return globals()['N'+str(i)+str(j)]

def anchor(Node,n):	#find Node's anchor
	for i in range(n):
		Node=P(Node)
	return Node

def depth(j,generation):	#the depth of j brunch, generation minu the number of None in this brunch
	assert j in range(2**generation),	'%d  brunch in all' %(2**generation)	
#	print Nonelist	
	i=generation
	depth=generation
#	for l in Nonelist:
#		if  not in l[1]:
#			return generation/2
#		else:			
	while [i,j] in Nonelist:
#		print [i,j]
		i=i-1
		j=int(j/2)
		depth=depth-1
	return depth	
	
def get_idname (d,Node):	#d=locals(),but not here
	for key in d:
		if d[key]==Node:
			return key
			
def runcmd(cmd):   
	#print >>sys.stderr, cmd
	run=subprocess.Popen(cmd.split(),env=os.environ,cwd=os.getcwd(),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	run_info=run.communicate() 
	return run_info
	






def main(args):
	treefile=args.treefile
	generation=args.generation
	oprefix=args.oprefix
			
		
	Nodelist=[]	#store all Node class
	id_Nodelist=[]	#store all Node id, for example N00,N10...
	globals()['Nonelist']=[]	#store the 'None' Node
	VAR=[]		#store all kinds of var name
#	with open('/home/lichao/Documents/github/svengine/exam_tree') as treefile:
#	treefile=file

#	with open('/home/lichao/Documents/github/svengine/tree1','w+') as tree1:
#tree1=open('/home/lichao/Documents/github/svengine/tree1','w+')
#tree1=tempfile.TemporaryFile()
#tree1=open('tree1','w+')
	for g,line in enumerate(treefile):
			if '/' not in line and g<=2*generation:
				#tree1.write(line)
				line=line.expandtabs().strip().split(' ')
				while '' in line:
					line.remove('')		
				l=len(line)
				
				if g==0:# the root node
					print line	
					assert  line==['N'],	"the tree should use N as root Node"
					globals()['N'+'0'+'0']=Node(0,0,{'Normal':None},1)
					Nodelist.append(Node(0,0,{'Normal':None},1))
				else:
					for j in range(l) :		
						
						if line[j]=='None':	#check the 'None' condition
							#print line[j],type(line[j])							
							Nonelist.append([g/2,j])
#							print Nonelist							
							globals()['N'+str(g/2)+str(j)]=Node(g/2,j,globals()['N'+str(g/2-1)+str(int(j/2))].v,0.5)
							Nodelist.append(Node(g/2,j,globals()['N'+str(g/2-1)+str(int(j/2))].v,0.5))
						else:
						
							
						
							

							line[j]=line[j].split(';')
							
							var=line[j][0].split(',')
							freq=line[j][1]
							#print	VAR.append			
						
							#print var.split(':')[0],var.split(':')[1],var,type(var)
						
							if var[0][0]=='N' and var[0][1].isdigit() :
								assert len(var)==1,	"Normal Node only have one type 'N'"	
								globals()['N'+str(g/2)+str(j)]=Node(g/2,j,{'Normal':None},freq)	
								Nodelist.append(Node(g/2,j,{'Normal':None},freq))
							else:
								for i in range(len(var)):
									if var[i] not in VAR:
										VAR.append(var[i])#f  var not in VAR:
						
	#						line[j]=line[j].split(';')
	#						b=line[j][0]
								varlist={}
							#freq=line[j][1]
							#print varlist
								
								map(lambda x:varlist.setdefault(x.split(':')[0], int(x.split(':')[1])), var)
						
								globals()['N'+str(g/2)+str(j)]=Node(g/2,j,varlist,freq)
								Nodelist.append(Node(g/2,j,varlist,freq))
	assert len(Nodelist)==2**(generation+1)-1,	'the number of Node should equal to 2**(generation+1)-1, the treefile may include error.'	
	print VAR_frq(VAR,generation)	


def freq(vari,generation):		#calculte the allel freq
	cell_all=0	#the number of all cell
	cell_vari=0	#the number of cell inculde vari: if the vari is heterozygote,the number be half
	count_brunch=[]	#the number of cell each brunch
	brunch_vari=[]	#the brunch id which include vari 
	for i in range(2**generation):
		a=1
		for j in range (generation):
			
			a=a*(anchor(globals()['N'+str(generation)+str(i)],j).a)
#			print a
		
		count_brunch.append(2**(depth(i,generation))*a)
		cell_all =round(sum(count_brunch),2)
	haplist=[]
	for i in range(2**generation):
		if vari in globals()['N'+str(generation)+str(i)].v.keys():
			haplist.append(globals()['N'+str(generation)+str(i)].v[vari])
			brunch_vari.append(i)
	if haplist.count(haplist[0])!=len(haplist):	
		print "Warning! The haplotype of each var must keep unchange,you better check your %s  "	%vari
	hap=haplist[0]
	assert hap in [0,1],	"the haplotype must 1(homozygote) or 0(heterozygote)"
	for j in brunch_vari:
		cell_vari=round((cell_vari+count_brunch[j])*(float(hap+1)/2),2)	#0-->0.5,1-->1
							
								
	return cell_vari,cell_all,round(cell_vari/cell_all,2)				
		
#print freq('var1')#,freq('var2'),freq('var3'),freq('var4'),freq('var5')
def VAR_frq(VAR,generation):
	VAR_frq={}
	for var in VAR:
		VAR_frq[var.split(':')[0]]=freq(var.split(':')[0],generation)	
	return VAR_frq










def run():
	parser = argparse.ArgumentParser(description="translate a binary tree to varfile can be accepted by fasforge and muteforge.")
	parser.add_argument("treefile", metavar="treefile", type=argparse.FileType('rU'), #example_tree	
											help="inputfile of the svolution tree, required")
	parser.add_argument('-g','--generation',dest='generation',type=int,help="the numbre of generation plan to simulation,required")
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #output1,output2,output3,...
												help="prefix for output files")
	args = parser.parse_args()
	main(args)


if __name__ == '__main__':
	run()
	print "tree2var have done! Total runtime: ", time.time()-start_time, " seconds"




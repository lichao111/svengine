#!user/bin/env python
import argparse,time,subprocess,csv,tempfile

try:
	import mf.mutforge as mf
except ImportError:
	import mutforge as mf
start_time=time.time()

class Node:
	def __init__(self,i,j,v,a)
	self.i=i
	self.j=j
	self.v=v
	self.a=a

class tree:
	def __init__(self, f, commentstring="/ "or "\ "):        
		self.f = csv.reader(f,delimiter='\t')
		self.commentstring = commentstring
	def next(self):
		line = self.f.next()
		while commentstring in line
			line = self.f.next()
		return line
	def __iter__(self):
		return self


def tree2treelist():

def C(Node):
	i=Node.i+1
	j1=2*(Node.j)
	j2=j1+1
	return {Node(i,j1),Node(i,j2)}

def p(Node):
	i=Node.i-1
	j=int((Node.j)/2)
	return Node(i,j)


def runcmd(cmd):   
	#print >>sys.stderr, cmd
	run=subprocess.Popen(cmd.split(),env=os.environ,cwd=os.getcwd(),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	run_info=run.communicate() 
	return run_info
	






def main():
	treefile=args.treefile
	generation=args.ganeration
	
	
	










def run():
	parser = argparse.ArgumentParser(description="Simulating allele-specific frequency cancer  evolution data")
	parser.add_argument("treefile", metavar="treefile", type=argparse.FileType('rU'), #example_tree	
											help="inputfile of the svolution tree, required")
	parser,add_argument('-g','--generation',dest='generation',type=int,help="the numbre of generation plan to simulation")
	args = parser.parse_args()
	main(args)


if __name__ == '__main__':
	run()
	print "total runtime: ", time.time()-starttime, " seconds"




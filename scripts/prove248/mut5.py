import sys,importlib,ast; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
mod=importlib.import_module(sys.argv[1])
N=5; fc=0; tot=0
for line in open('../prove247b/z5.txt'):
    c=int(line.split()[3]); zs=set(ast.literal_eval(line[line.index('['):line.rindex(']')+1]))
    for lam,mu,nu,a,b in mod.allA(N):
        if (a,b,lam,mu,nu) in zs or (a,b,mu,lam,nu) in zs:
            tot+=1; fc+=mod.start('A',N,lam,mu,nu,a,b,('eq',c))
print(sys.argv[1],'A N=5 zeros',tot,'falsecert',fc,'PSTAT',mod.PSTAT,flush=True)

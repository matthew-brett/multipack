#!/usr/bin/env python
"""
dpythonrun --- A distiller of the pythonrun.sty documents.

Copyright 1999 Pearu Peterson all rights reserved,
Pearu Peterson <pearu@ioc.ee>          
Permission to use, modify, and distribute this software is given under the
terms of the LGPL.  See http://www.fsf.org

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Pearu Peterson

Usage:
    dpythonrun <infile> <outfile> # Note that `rm -f <outfilename>.*'
                                  # will be executed
    dpythonrun <infile>           # Recommended usage
    dpythonrun                    # The distillation will be not 
                                  # completed (see the code at the EOF) 
    
    Default <outfile> is stdout. Default <infile> is stdin.
Limitations:
    1) A line must contain at most one pyhtonrun.sty's environment or
    command
    2) All pythonrun.sty environments must start at the first position
    in the line
    3) All comments (lines starting with `%') are removed after
    `\begin{document}'
    4) <outfile> must end with `.tex'
"""

__version__ = "0.1"

stdoutflag=0
import sys,os,string,fileinput,re

try: fn=sys.argv[2]
except:
    try: fn='distil_'+sys.argv[1]
    except:
        stdoutflag=1
        fn='tmp_pythonrun_distil_stdout.tex'
try: fi=sys.argv[1]
except: fi=()
sys.stdout=open(fn,'w')
fh=open('tmp_pythonrun_distil_head.tex','w')

nonverb=r'[\w\s\\&=\^\*\.\{\(\)\[\?\+\$/]*(?!\\verb.)'
nonverb2=r'[\w\s&=\^\*\.\{\(\)\[\?\+\$/]*(?!\\verb.)'
usepythonrun=re.compile(nonverb+r'\\usepackage.*\{pythonrun}')
beginverbatim=re.compile(r'^\\begin\s*\{verbatim\*?}')
endverbatim=re.compile(r'^\\end\s*\{verbatim\*?}')
begindocument=re.compile(r'^\\begin\s*\{document}')
enddocument=re.compile(r'^\\end\s*\{document}')
beginpython=re.compile(r'^\\begin\s*\{python(run(withsave|)|call|)\*?}')
endpython=re.compile(r'^\\end\s*\{python(run(withsave|)|call|)\*?}')
pythonreset=re.compile(nonverb2+r'\\python(reset|clean)')
pythoninput=re.compile(nonverb+r'\\python(input|save)\*?\s*\{.*}')
pythoninputshort=re.compile(nonverb+r'\\python(input|save)')
pythonline=re.compile(nonverb+r'\\pythonline\*?\s*\{.*}')
pythonline2=re.compile(nonverb+r'\\pythonline\*?\s*\{.*')
comment=re.compile(r'[^%]*%')
slash=re.compile(r'[^%]*\\')
bcurly=re.compile(r'[^%]*}')
docflag=0
verbatimflag=0
lineflag=0
def addnl(x):
    if string.strip(x)=='': return ''
    return x+'\n'
def raddnl(x):
    if string.strip(x)=='': return ''
    return '\n'+x
for l in fileinput.input(fi):
    l=l[:-1]
    l1=''
    if comment.match(l):
        m=comment.match(l)
        l1=l[m.end()-1:]
        l=l[:m.end()-1]
    if docflag:
        if verbatimflag:
            m=endverbatim.match(l)
            print l
            if m:
                verbatimflag=0
        else:
            if lineflag:
                if bcurly.match(l):
                    lineflag=0
                    m=bcurly.match(l)
                    print addnl(l[:m.end()])+r'\begin{pythonlatex}'+ \
                          raddnl(l[m.end():])
                continue
            if beginverbatim.match(l):
                verbatimflag=1
            elif enddocument.match(l):
                print r'\end{pythonlatex}'
            elif beginpython.match(l):
                print r'\end{pythonlatex}'
            if pythonreset.match(l):
                #sys.stderr.write(l+'\n')
                m=pythonreset.match(l)
                s=slash.match(l[:m.end()])
                print addnl(l[:s.end()-1])+r'\end{pythonlatex}'+'\n'+ \
                      addnl(l[s.end()-1:m.end()])+r'\begin{pythonlatex}'+\
                      raddnl(l[m.end():])
            elif pythoninput.match(l):
                ms=pythoninputshort.match(l)
                m=pythoninput.match(l)
                s=slash.match(l[:ms.end()])
                print addnl(l[:s.end()-1])+r'\end{pythonlatex}' +'\n'+\
                      addnl(l[s.end()-1:m.end()])+r'\begin{pythonlatex}'+\
                      raddnl(l[m.end():])
            elif pythonline.match(l):
                #sys.stderr.write(l+'\n')
                m=pythonline.match(l)
                s=slash.match(l[:m.end()])
                print addnl(l[:s.end()-1])+r'\end{pythonlatex}'+'\n'+\
                      addnl(l[s.end()-1:m.end()])+r'\begin{pythonlatex}'+\
                      raddnl(l[m.end():])
            elif pythonline2.match(l):
                lineflag=1
                m=pythonline2.match(l)
                s=slash.match(l[:m.end()])
                print addnl(l[:s.end()-1])+r'\end{pythonlatex}'+\
                      raddnl(l[s.end()-1:m.end()])
            else:
                print l
            if endpython.match(l):
                print r'\begin{pythonlatex}'
    else:
        if usepythonrun.match(l):
            print r'%'+l
            print r'\usepackage[distil]{pythonrun}'
            fh.write(r'%'+l+'\n')
            fh.write(r'%\usepackage[distil]{pythonrun}'+'\n')
        else:
            print l+l1
            fh.write(l+l1+'\n')
            if begindocument.match(l):
                docflag=1
                print r'\begin{pythonlatex}'
fh.close()

sys.stdout.close()
os.system('latex --shell-escape '+fn+' 1>&2')
os.system('rm -f '+fn[:-3]+'* 1>&2')
os.system('cat tmp_pythonrun_distil_head.tex tmp_pythonrun_distil1.tex > '+fn)
os.system('echo "\\end{document}" >> '+fn)
if stdoutflag:
    os.system('cat '+fn)
os.system('rm -f tmp_pythonrun_distil*.tex 1>&2')

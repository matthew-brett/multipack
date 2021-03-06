#!/usr/bin/env python
"""
setmodules --- Helper program for setting up multipack modules.

Copyright 1999 Pearu Peterson all rights reserved,
Pearu Peterson <pearu@ioc.ee>          
Permission to use, modify, and distribute this software is given under the
terms of the LGPL.  See http://www.fsf.org

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Usage:
    setmodules.py <modul1> <modul2> .. : <ulib1> <ulib2> ..
The following flags can be used:
    -q  --- run in quite mode
    -a  --- create file m_multipackmodule.c (obsolete)
    -m  --- create files m_*module.c
*/
Pearu Peterson
"""
__version__ = "0.4"

import sys,fileinput,re,string

infile='_multipackmodule.c'
outfile='m_multipackmodule.c'
mlist=[]    # user required modules
ulibs=[]    # user provided libraries
cflag=0
q_flag=a_flag=m_flag=0
for m in sys.argv[1:]:
    m=string.strip(m)
    if m==':': cflag=1;continue
    if m=='-q':q_flag=1;continue
    if m=='-a':a_flag=1;continue
    if m=='-m':m_flag=1;continue
    if cflag: ulibs.append(m)
    else: mlist.append(m)
ulibs=map(string.strip,ulibs)
files=[]    # __*.c to be included to 'outfile'
methods=[]  # module methods to be inserted to 'outfile'
mods=[]     # existing submodules
libs=[]     # libraries required by submodules
llibs=[]    # final list of libraries to be linked with multipack module
python=[]   # python wrappers of C-submodules: to be imported in Multipack.py
insmeth=re.compile(r'.*PyMethodDef\s*multipack_module_method.*{')
#modmeth=re.compile(r'.*/\*.*multipack_module_methods\s*:')
modmeth=re.compile(r'.*/\*.*module_methods\s*:')
linklibs=re.compile(r'.*/\*.*link\s*libraries\s*:')
pyfiles=re.compile(r'.*/\*.*python\s*files\s*:')
revision=re.compile(r'.*\$Revision.*\$')
endcomment=re.compile(r'.*\*/')
for m in mlist:
    if not q_flag:
        sys.stderr.write('Scanning module "%s" for methods\n'%m)
    try:
        f=open('__%s.c'%m,'r')
        f.close()
        files.append('__%s.c'%m)
        mods.append(m)
    except:
        sys.stderr.write('Warning: Could not open file __%s.c. Skipping.\n'%m)
        continue
    cflag,eflag=0,[]
    version=''
    if m_flag:
        f=open('m_%smodule.c'%m,'w')
        f.write("""/*
    Multipack project.
    This file is generated by setmodules.py. Do not modify it.
 */
#include \"../multipack.h\"
static PyObject *%s_error;
#include \"../%s\"
static struct PyMethodDef %s_module_methods[] = {
"""%(m,files[-1],m))
    for l in fileinput.input(files[-1]):
        l=l[:-1]
        if cflag:
            mc=endcomment.match(l)
            if mc:
                if cflag not in eflag: eflag.append(cflag)
                l=l[:mc.end()-2]
            l=string.strip(l)
            if l:
                if cflag==1:
                    methods.append(l)
                    if m_flag:
                        f.write(l+'\n')
                elif cflag==2:
                    llibs.append(l)
                    if not ((l in ulibs) or ('!'+l in ulibs)):
                        libs.append(l)
                    else:
                        if not q_flag:
                            sys.stderr.write('Skipping library %s. It is provided by user.\n'%`l`)
                elif cflag==3:
                    try:
                        fp=open(l,'r');fp.close()
                        python.append(l)
                    except:
                        if not q_flag:
                            sys.stderr.write('Warning: Could not open python file %s. Skipping.\n'%`l`)
            if mc: cflag=0
        elif modmeth.match(l): cflag=1
        elif linklibs.match(l): cflag=2
        elif pyfiles.match(l): cflag=3
        elif revision.match(l):
            mc=revision.match(l)
            i=mc.end()-2
            while (not l[i]=='$') and i: i=i-1
            version=l[i+10:mc.end()-1]
        if eflag==[1,2,3]: fileinput.close() # one can comment this out in order to
                                             # force scanning the hole files
                                             # (it takes then 4-5 times longer)
    if m_flag:
        f.write("""{NULL,		NULL, 0, NULL}
};
void init_%s() {
  PyObject *m, *d, *s;
  m = Py_InitModule(\"_%s\", %s_module_methods);
  import_array();
  d = PyModule_GetDict(m);

  s = PyString_FromString(\"%s\");
  PyDict_SetItemString(d, \"__version__\", s);
  %s_error = PyErr_NewException (\"%s.error\", NULL, NULL);
  Py_DECREF(s);
  if (PyErr_Occurred())
    Py_FatalError(\"can't initialize module %s\");
}
        """%(m,m,m,version,m,m,m))
        f.close()
    if eflag==[] and (not q_flag):
        sys.stderr.write('Warning: No module methods or libraries found in "%s".\n'%files[-1])
mflag=0
if a_flag:
    f=open(outfile,'w')
    for l in fileinput.input(infile):
        if insmeth.match(l):
            for i in files:
                f.write('#include "%s"\n'%i)
            f.write(l)
            for m in methods:
                f.write(m+'\n')
        else:
            f.write(l)
    f.close()

llibs=map(string.strip,llibs)
libs=map(string.strip,libs)
def f(ls,i=0):
    ll=[]
    for l in ls:
        if not l[0]=='!': ll.append(l)
        else:
            if i: break
    return ll
i=0
for l in ulibs:
    if l[0]=='!':
        lt=[]
        for ll in llibs:
            if ll==l[1:]:
                lt=lt+f(ulibs[i:])
            else: lt.append(ll)
        llibs=lt
    i=i+1
libs.reverse()
libstmp=[]
for l in libs:
    if l not in libstmp: libstmp.append(l)
libs=libstmp
libs.reverse()
llibs.reverse()
libstmp=[]
for l in llibs:
    if l not in libstmp: libstmp.append(l)
llibs=libstmp
llibs.reverse()
for l in f(ulibs):
    if l not in llibs: llibs.append(l)
print string.join(map(lambda t:'m.'+t,mods),' '),\
      string.join(map(lambda t:'l.'+t,libs),' '),\
      string.join(map(lambda t:'ll.'+t,llibs),' '),\
      string.join(map(lambda t:'p.'+t,python),' ')
if not q_flag:
    sys.stderr.write('setmodules.py found %d methods, %d libraries, %d python files in %d files\n'%(len(methods),len(libs),len(python),len(files)))



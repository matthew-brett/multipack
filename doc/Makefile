# Makefile for compiling Multipack documentation (dvi, ps, html)
# Pearu Peterson <pearu@ioc.ee>

TOP = usersguide
LATEXSRC = main.tex commands.tex fmini.tex nleqns.tex odepack.tex quadpack.tex fitpack.tex 

PYTHONRUN=pythonrun-1.1
LATEX = TEXINPUTS=.:./$(PYTHONRUN): latex
COLLECTINPUT = $(PYTHONRUN)/collectinput.py
DPYTHONRUN = $(PYTHONRUN)/dpythonrun.py

TTH = tth
TTHFILTER = sed -e "s/\vec/\b&/g" | $(TTH) -Ltmp_main -i
TTHMISSING = "\
***************************************************************\n\
Warning:        Could not find tth (a TeX to HTML translator)  \n\
                or error arised by tth\n\
You can download tth from http://hutchinson.belmont.ma.us/tth/  \n\
or\n\
use your favorite LaTeX to HTML translator on file tmp_main.tex\n\
***************************************************************\
"

all: dvi ps html clean
main.dvi: $(LATEXSRC)
	$(LATEX) --shell-escape main.tex
	$(LATEX) --shell-escape main.tex
	$(LATEX) --shell-escape main.tex
$(TOP).dvi: main.dvi
	cp main.dvi $(TOP).dvi
$(TOP).ps: $(TOP).dvi
	dvips $(TOP).dvi -o
$(TOP).html: $(LATEXSRC)
	$(COLLECTINPUT) < main.tex  | $(DPYTHONRUN) > tmp_main.tex
	$(LATEX) tmp_main.tex
	@test `which $(TTH)` && cat tmp_main.tex | $(TTHFILTER) > $(TOP).html\
	 || echo -e $(TTHMISSING)
dvi: $(TOP).dvi
ps: $(TOP).ps
html: $(TOP).html
clean:
	rm -f tmp_main.* main.{aux,dvi,log,toc}
distclean:
	rm -f tmp_main.* main.{aux,dvi,log,toc}
	rm -f $(TOP).{ps,dvi,html}
	rm -f *~

DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: maskaln

clean:
	rm -f *.out

maskaln: maskaln.00 maskaln.01 maskaln.02

maskaln.00:
	../R/maskaln.r $@.alnfaa $@.out
	@$(CHECK)

maskaln.01:
	../R/maskaln.r --verbose $@.alnfaa $@.out

maskaln.02:
	../R/maskaln.r --protect_taxa=$@.protect $@.alnfaa $@.out
	@$(CHECK)

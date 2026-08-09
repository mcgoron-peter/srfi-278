.POSIX:
VERSION=0.2.0

chicken/${VERSION}-6.tar.gz:
	mkdir -p chicken
	tar -czvf chicken/6-${VERSION}.tar.gz srfi-278.egg lib/ tests/ LICENSES

clean:
	rm -f *.so *.o *.link *.install.sh *.import.scm *.build.sh *.log tests/*.log

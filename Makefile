# SPDX-FileCopyrightText: 2026 Peter McGoron
# SPDX-License-Identifier: MIT

.POSIX:
VERSION=0.2.1

chicken/${VERSION}-6.tar.gz:
	mkdir -p chicken
	tar -czvf chicken/${VERSION}-6.tar.gz srfi-278.egg lib/ tests/ LICENSES

clean:
	rm -f *.so *.o *.link *.install.sh *.import.scm *.build.sh *.log tests/*.log

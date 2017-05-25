PATH = /Users/matteo/Documents/MassTodon/MassTodonPy/

install:
	pip install -e $(PATH)

reinstall:
	pip uninstall -y MassTodonPy
	pip install -e $(PATH)


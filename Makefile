PATH = /Users/matteo/Documents/MassTodon/MassTodonPy/

install:
	pip install -e $(PATH)

reinstall:
	pip uninstall -y MassTodonPy
	pip install -e $(PATH)

example_call:
	python2 $(PATH)/Tests/calls/example_call.py
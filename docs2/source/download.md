Download
========

Prerequisites
-------------

MassTodonPy can be used with any UNIX operating system.
It was tested on 
* macOS HighSierra
* Ubuntu 16.04
* Ubuntu 17.10
* Gentoo

It requires a python interpreter.
Both Python 2 and 3 are supported.
However, CVXOPT optimization behaves more stable on Python 2.7.

Installation
-------------

##### PIP

To install MassTodonPy globally, run:
```bash
pip install MassTodonPy
```

To check the installation, type:

```bash
masstodon_example_call
```

Ideally, you will run an example session of MassTodon and see output in your browser.


##### Virtual Environment
If you cannot install dependencies on your machine, but you have access to virtualenv, or you just like virtualenv, because it's neet, then run:
```bash
virtualenv masstodon_ve
source masstodon_ve/bin/activate
pip install MassTodonPy
masstodon_example_call
deactivate
```
[(More on Virtual Environments)](https://python-guide-cn.readthedocs.io/en/latest/dev/virtualenvs.html).

##### GitHub

The package is hosted on [GitHub](https://github.com/MatteoLacki/MassTodonPy/tree/devel): check out the latest devel version!

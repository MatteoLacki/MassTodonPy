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

To install the MassTodon with Virtual Environment, run:

```bash
virtualenv masstodon_ve
source masstodon_ve/bin/activate
pip install MassTodonPy
masstodon_example_call
deactivate
```

This way you might bypass lack of super user privileges.



##### GitHub

The package is hosted on [GitHub](https://github.com/MatteoLacki/MassTodonPy/tree/devel).
To download from github, you will need git.
To get to the *devel* version:

```bash
git clone https://github.com/MatteoLacki/MassTodonPy.git
cd MassTodonPy
git checkout devel
pip install .
```

If you want to install a version that you would locally pimp yourself on branch "my_pimped_masstodon":
```bash
git clone https://github.com/MatteoLacki/MassTodonPy.git
cd MassTodonPy
git checkout -b my_pimped_masstodon
pip install -e .
```

This will only soft-link the MassTodonPy folder with the default folder PIP installs into.

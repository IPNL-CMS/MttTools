# Presentation

This is the analysis tools for mtt study.

## Get the code

Get the code by executing this command:

Using HTTPS:

```bash
myWorkingDir> git clone https://github.com/IPNL-CMS/MttTools.git
```

Using SSH keys (see how to generate a SSH key and add it to GitHub [here](https://help.github.com/articles/generating-ssh-keys)):

```bash
myWorkingDir> git clone git@github.com:IPNL-CMS/MttTools.git
```

You can see that the directory "plotIt" is empty. Retrieve the code from the corresponding sub-repository as shown below:

```bash
myWorkingDir> cd MttTools
MttTools> git submodule init
MttTools> git submodule update
```

## Setup the environment

Every time you want to work with MttTools, you have to setup the environment: 

```bash
MttTools> source setup_lyoserv_env.sh
```

You also have to dowload and compile the external libraries:

```bash
MttTools> cd external/
external> ./build-external.sh
external> cd..
```

## New: make LHAPDF library available

```bash
MttTools> cd external/
external> export PATH=$PWD/bin:$PATH
external> export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
external> export PYTHONPATH=$PWD/lib64/python2.6/site-packages:$PYTHONPATH
```

Check that everything is ok:

```bash
external> lhapdf-config --help
external> cd..
```


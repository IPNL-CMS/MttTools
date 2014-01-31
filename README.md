# Presentation

This is the analysis tools for mtt study.

## Get the code

Get the code by executing this command:

```bash
myWorkingDir> git clone https://github.com/IPNL-CMS/MttTools.git
```

You can see that the directory "plotIt" is empty. Retrieve the code from the corresponding sub-repository as showed below:

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

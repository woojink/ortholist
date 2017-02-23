# Ortholist 2.0
## Overview
This project updates Ortholist ([1.0 here](http://www.greenwaldlab.org/ortholist/)) with results from both additional ortholog prediction programs and updated releases of previously included programs.

## Databases used
* Ensembl Compara (Release 87, 2016-12)
* Homologene (Build 68, 2014-05-06)
* InParanoid (8.0, 2013-12)
* OMA (May 2016)
* OrthoInspector
* OrthoMCL (Version 5, 2015-07-23)

## Instructions
### Pre-requisites
This project makes use of Python 3 and MySQL. If you're on OS X / macOS, you can get these by first installing [Homebrew](https://brew.sh/) by putting the following on your Terminal:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then running
```
brew install python3
brew install mysql
```

Once MySQL installation is complete, get the service started with `mysql.server start`. You can connect to the service using `mysql -uroot` anytime to check the databases within.

### Code
The majority of the code is located under `src/run.py`, which processes the both raw and pre-processed data from ortholog prediction programs stored under `results`.

## Results
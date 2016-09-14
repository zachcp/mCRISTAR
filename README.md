# mCRISTAR

Code for the [mCRISTAR](http://www.mcristar.net/) webapplication.


## App Overview

The basid idea of the mCRISTAR webapp is to facilatate the mCRISTAR cluster refactoring process.
Users upload a genmbank file, edit the promoter sites they would like to change and are given the
CRISPR cassettes and primers needed to generate bridging constructs.

mCristar runs a [Flask](http://flask.pocoo.org/) app which handles the sequence parsing and generation
of primers. The frontend is written in Clojurescript using the [re-frame](https://github.com/Day8/re-frame) system. 


## To Run Locally

```[shell]
pip install -r requirements.txt
python main.py
```

## To Upload
```[shell]
make upload
````

## To Develop
```
#in one terminal window
cd cljs
lein figwheel
#in another window
#python 2
dev_appserver.py.

#python3 osx
/usr/bin/python2.7  /usr/local/bin/dev_appserver.py .
```
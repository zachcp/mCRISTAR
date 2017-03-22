python = /usr/bin/python2.7 
appcfg = /usr/local/bin/appcfg.py #from homebrew: brew install app-engine-python
devserver = /usr/local/bin/dev_appserver.py #from homebrew: brew install app-engine-python


dev:
	#in another window:cd cljs lein figwheel &
	$(python) $(devserver) .

pipinstall:
	pip install -r requirements.txt -t lib

upload:
	cd cljs;lein clean;lein cljsbuild once min
	$(python) $(appcfg) -A mcristar-1125 update app.yaml

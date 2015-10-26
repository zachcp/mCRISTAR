pipinstall:
	pip install -r requirements.txt -t lib

upload:
	cd cljs;lein clean;lein cljsbuild once min
	appcfg.py -A mcristar-1125 update app.yaml

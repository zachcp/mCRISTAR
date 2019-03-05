# from homebrew: 
# brew tap caskroom/cask
# brew cask install google-cloud-sdk


python = /usr/local/bin/python3 

dev:
	#in another window:cd cljs leinfigwheel &
	$(python) main.py

pipinstall:
	pip3 install -r requirements.txt -t lib
	pip3 install -r requirements_bio.txt -t lib --no-deps

compilecljs:
	cd cljs && lein clean && lein cljsbuild once min && cd ..

deploy: compilecljs
	gcloud app deploy

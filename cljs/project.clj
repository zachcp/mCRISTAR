(defproject mCRISTAR "0.1.0-SNAPSHOT"
  :description "FIXME: write this!"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}

  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojure/clojurescript "1.9.227"]
                 [org.clojure/core.async "0.2.374"]
                 [reagent "0.6.1"]
                 [re-frame "0.9.2"]
                 [cljs-ajax "0.5.8"]
                 [cljsjs/jquery "2.2.4-0"]
                 [cljsjs/jquery-ui "1.11.4-0"]]

  :plugins [[lein-cljsbuild "1.1.5"]
            [lein-figwheel "0.5.9"]]

  :source-paths ["src"]

  :clean-targets ^{:protect false} ["../static/js/compiled" "target"]

  :cljsbuild {:builds
              [{:id "dev"
                :source-paths ["src"]

                :figwheel {:on-jsload "mCRISTAR.core/on-js-reload"}

                :compiler {:main mCRISTAR.core
                           :asset-path "../static/js/compiled/out"
                           :output-to "../static/js/compiled/mCRISTAR.js"
                           :optimizations :none
                           :output-dir "../static/js/compiled/out"
                           :foreign-libs [{:file "src/JS/Scribl/Scribl.1.1.5.min.js" :provides ["scribl"]}
                                          {:file "src/JS/utils/dragscrollable.js" :provides ["draggable"]}]
                           :source-map-timestamp true}}
               ;; This next build is an compressed minified build for
               ;; production. You can build this with:
               ;; lein cljsbuild once min
               {:id "min"
                :source-paths ["src"]
                :compiler {:output-to "../static/js/compiled/mCRISTAR.js"
                           :main mCRISTAR.core
                           :optimizations :advanced
                           :foreign-libs [{:file "src/JS/Scribl/Scribl.1.1.5.min.js" :provides ["scribl"]}
                                          {:file "src/JS/utils/dragscrollable.js" :provides ["draggable"]}]
                           :pretty-print false}}]}

  :figwheel {});; :http-server-root "public" ;; default and assumes "resources"
             ;; :server-port 3449 ;; default
             ;; :server-ip "127.0.0.1"

             ;;:css-dirs ["../static/css"] ;; watch and update CSS

             ;; Start an nREPL server into the running figwheel process
             ;; :nrepl-port 7888

             ;; Server Ring Handler (optional)
             ;; if you want to embed a ring handler into the figwheel http-kit
             ;; server, this is for simple ring servers, if this
             ;; doesn't work for you just run your own server :)
             ;; :ring-handler hello_world.server/handler

             ;; To be able to open files in your editor from the heads up display
             ;; you will need to put a script on your path.
             ;; that script will have to take a file path and a line number
             ;; ie. in  ~/bin/myfile-opener
             ;; #! /bin/sh
             ;; emacsclient -n +$2 $1
             ;;
             ;; :open-file-command "myfile-opener"

             ;; if you want to disable the REPL
             ;; :repl false

             ;; to configure a different figwheel logfile path
             ;; :server-logfile "tmp/logs/figwheel-logfile.log"


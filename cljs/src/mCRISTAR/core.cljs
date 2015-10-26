(ns mCRISTAR.core
  (:require [reagent.core :as reagent :refer [atom]]
            [re-frame.core :as re-frame]
            [mCRISTAR.handlers]
            [mCRISTAR.subs]
            [mCRISTAR.views :as views]))

(enable-console-print!)

(println "Edits to this text should show up in your developer console.")

(defn mount-root []
  (reagent/render [views/main-panel]
                  (.getElementById js/document "app")))

(let []
  (re-frame/dispatch-sync [:initialize-db])
  (mount-root))

(defn on-js-reload []
  ;; optionally touch your app-state to force rerendering depending on
  ;; your application
  ;; (swap! app-state update-in [:__figwheel_counter] inc)
)

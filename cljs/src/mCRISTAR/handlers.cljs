(ns mCRISTAR.handlers
  (:require [re-frame.core :as re-frame]
            [mCRISTAR.db :as db]
            [ajax.core :refer [GET POST raw-response-format]]))

(re-frame/register-handler
  :bad-response
  (fn
    [_ [_ resp]]
    (js/alert "notgood!")
    (js/alert (str  resp))
    ((.log js/console resp))))   ;; pure handlers must return a db (unchanged in this case)

(re-frame/register-handler
  :update-gaps
  (fn  [db [_ resp]]
    (assoc db :fill (str resp))))

(re-frame/register-handler
  :initialize-db
  (fn  [_ _]
    db/default-db))

(re-frame/register-handler
  :make-cassettes
  (fn
    [db [_ _]]
    (let [activegaps (re-frame/subscribe [:active-gaps])
          activecount (count @activegaps)]
      (cond
        (< activecount 1) (js/alert "At least one site must be chosen")
        (> activecount 7) (js/alert "At this time a maximum of 7 sites can be chosen.")
        :else
          (POST
            "/makecassettes"
            {:params          {:gaps @activegaps}
             :format :json
             :response-format :json
             :keywords?       true
             :handler         #(re-frame/dispatch [:process-cassettes %1])
             :error-handler   #(re-frame/dispatch [:bad-response %1])}))
      db)))    ;; pure handlers must return a db (unchanged in this case)

(re-frame/register-handler
  :process-cassettes
  (fn
    [db [_ resp]]
    (let [data (js->clj resp)
          cassettes (:cassettes data)
          primers (:primers data)]
      (-> db
          (assoc :cassettes cassettes)
          (assoc :primers primers)
          (assoc :app-cycle 2)))))

(re-frame/register-handler
  :process-gaps
  (fn
    [db [_ resp]]
    (let [data (js->clj resp)
          truefalseconvert (fn [gene]
                             (let [tf (:selected gene)
                                   truefalse (cond
                                               (= "true"  tf) true
                                               (= "false" tf) false)]
                               (assoc gene :selected truefalse)))
          genes  (->> (:genes data)
                      (map truefalseconvert))
          gaps (->> (:gaps data)
                    (map truefalseconvert)
                    (vec))
          clustersize (:end (last genes))
          scale (/ 1000 clustersize)
          scale-min 0
          scale-max (* scale (/ clustersize 1000))
          scale-step (/ (- scale-max scale-min) 100)]
      (js/console.log gaps)
      (-> db
          (assoc :sequence-end clustersize)
          (assoc :svg-scale scale)
          (assoc :svg-scale/min scale-min)
          (assoc :svg-scale/max scale-max)
          (assoc :svg-scale/step scale-step)
          (assoc :genes genes)
          (assoc :gaps gaps)
          (assoc :app-cycle 1)))))

(re-frame/register-handler
  :toggle-gap
  (fn  [appdb [_ id]]
    (js/console.log (str id "  toggled!"))
    (update-in appdb [:gaps id :selected] not)))

(re-frame/register-handler
  :update
  (fn  [appdb [_ key value]]
    (js/console.log (str key " " value))
    (assoc appdb key value)))

(re-frame/register-handler
  :upload-file
  (fn
    [db [_ fd]]
    (POST
      "/upload"
      {:body fd
         ;:response-format (raw-response-format)
       :response-format :json
       :keywords? true
       :handler       #(re-frame/dispatch [:process-gaps %1])
       :error-handler #(re-frame/dispatch [:bad-response %1])})
    db))    ;; pure handlers must return a db (unchanged in this case)


(ns mCRISTAR.subs
  (:require-macros [reagent.ratom :refer [reaction]])
  (:require [re-frame.core :as re-frame]))

(re-frame/register-sub
  :active-gaps
  (fn [db]
    (let [allgaps (reaction (:gaps @db))]
      (reaction (filter #(true? (:selected %)) @allgaps)))))

(re-frame/register-sub
  :app-cycle
  (fn [db]
    (reaction (:app-cycle @db))))

(re-frame/register-sub
  :cassettes
  (fn [db]
    (reaction (:cassettes @db))))

(re-frame/register-sub
  :gaps
  (fn [db]
    (reaction (:gaps @db))))

(re-frame/register-sub
  :genes
  (fn [db]
    (reaction (:genes @db))))

(re-frame/register-sub
  :primers
  (fn [db]
    (reaction (:primers @db))))


(re-frame/register-sub
  :sequence-end
  (fn [db]
    (reaction (:sequence-end @db))))

(re-frame/register-sub
  :sequence-start
  (fn [db]
    (reaction (:sequence-start @db))))

(re-frame/register-sub
  :svg-scale
  (fn [db]
    (reaction (:svg-scale @db))))

(re-frame/register-sub
  :svg-scale/min
  (fn [db]
    (reaction (:svg-scale/min @db))))

(re-frame/register-sub
  :svg-scale/max
  (fn [db]
    (reaction (:svg-scale/max @db))))

(re-frame/register-sub
  :svg-scale/step
  (fn [db]
    (reaction (:svg-scale/step @db))))

(re-frame/register-sub
  :svg-translate
  (fn [db]
    (reaction (:svg-translate @db))))

(re-frame/register-sub
  :svg-width
  (fn [db]
    (reaction (:svg-width @db))))

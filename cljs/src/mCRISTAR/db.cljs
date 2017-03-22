(ns mCRISTAR.db)

(def default-db
  {:app-cycle      0    ;; determine the overall state of the app
                        ;; 0 is for pre-upload
                        ;; 1 is after uplaod with gaps
                        ;; 2 is the table with the cassette information
   :svg-scale      1
   :svg-scale/min  0
   :svg-scale/max  1
   :svg-scale/step 0.1
   :svg-translate  1
   :svg-width      933
   :sequence-start 0
   :sequence-end   1000
   :genes          []
   :gaps           []
   :cassettes      []})


(def test-db
  {:app-cycle      0    ;; determine the overall state of the app
   ;; 0 is for pre-upload
   ;; 1 is after uplaod with gaps
   ;; 2 is the table with the cassette information
   :svg-scale      1
   :svg-scale/min  0
   :svg-scale/max  1
   :svg-scale/step 0.1
   :svg-translate  1
   :svg-width      933
   :sequence-start 0
   :sequence-end   1000
   :genes          [{:id \A :start 0 :end 100 :strand 1}
                    {:id \B :start 150 :end 300 :strand -1}
                    {:id \C :start 325 :end 600 :strand 1}
                    {:id \D :start 625 :end 800 :strand 1}]
   :gaps           [{:id 0
                     :protein1 {:start 10 :end 50 :sequence "ACACAC" :strand 1}
                     :gap      {:start 51 :end 100 :sequence "ACACACACCCCCC" :strand 1}
                     :protein2 {:start 101 :end 140 :sequence "ACACACACCCCCC" :strand 1}
                     :selected true}
                    {:id 1
                     :protein1 {:start 150 :end 180 :sequence "ACACAC" :strand 1}
                     :gap      {:start 190 :end 220 :sequence "ACACACACCCCCC" :strand 1}
                     :protein2 {:start 230 :end 240 :sequence "ACACACACCCCCC" :strand 1}
                     :selected true}
                    :cassettes       [ [{:start 0
                                         :end 10
                                         :sequence "ACACAC"
                                         :type "crisprsite"}]
                                      [{:start 11
                                        :end 20
                                        :sequence "ACACAC"
                                        :type "BSAI site"}]
                                      [{:start 21
                                        :end 30
                                        :sequence "ACACAC"
                                        :type "GoldenGate site"}]]]

   :primers        [{:promoterid "SNP11"
                     :selection  "LEU"
                     :forward "ACACACA"
                     :reverse "ACACAC"
                     :strandorientation -1}
                    {:promoterid "SNP11"
                     :selection  "LEU"
                     :forward "ACACACA"
                     :reverse "ACACAC"
                     :strandorientation -1}]})
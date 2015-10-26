(ns mCRISTAR.views
  (:require [re-frame.core :as re-frame]))


(defn drawgene
  " Draw an SVG Gene by calculating five points.
  Pos: ------------------
       |                 \\
       |                 //
       ------------------
  Neg:   ------------------
       //                  |
       \\                  |
         -----------------
  "
  [gene scalefactor translate]
  (let [start (* scalefactor (- (:start gene) translate))
        end (* scalefactor (- (:end gene) translate))
        height 20
        trackheight 10
        strand (:strand gene)
        arrow-length 10]
    (if (> strand 0)
      ;; Positive Strand
      ;;  --->
      (let [arrowloc (if (<= (- end start) arrow-length )
                      start
                      (- end arrow-length))
            topleft [start, (+ trackheight(/ height 2))]
            bottomleft [start, (- trackheight(/ height 2))]
            topright [arrowloc, (+ trackheight(/ height 2))]
            bottomright [arrowloc, (- trackheight(/ height 2))]
            point [end, trackheight]
            points [topleft topright point bottomright bottomleft]
            pointstr (map (fn [[a b]] (str a "," b)) points)]
        [:polyline {:points (apply str (interpose " " pointstr))
                    :id (:id gene)}])
      ;; Negative Strand
      ;; <-----
      (let [arrowloc (if (<= (- end start) arrow-length )
                       end
                       (+ start arrow-length))
            point [start, trackheight]
            topleft [arrowloc, (+ trackheight(/ height 2))]
            topright [end, (+ trackheight(/ height 2))]
            bottomright [end, (- trackheight(/ height 2))]
            bottomleft [arrowloc, (- trackheight(/ height 2))]
            points [point topleft topright bottomright bottomleft]
            pointstr (map (fn [[a b]] (str a "," b)) points)]
        [:polyline {:points (apply str (interpose " " pointstr))
                    :id (:id gene)}]))))


(defn slider
  [param value min max step]
  (let [width (re-frame/subscribe [:svg-width])]
    [:input {:key param
             :type "range"
             :value @value
             :min min
             :max max
             :step step
             :style {:width @width}
             :on-change #(re-frame/dispatch [:update param (-> % .-target .-value)])}]))


(defn drawgap [gap scalefactor translate]
  "Draw the Gaps in sucha way that they can be toggled on and off."
  (let [id (:id gap)
        klass (if (= (:selected gap) false) "inactive" "active")]
    [:rect {:class (str "gap" " " klass)
            :id id
            :x (* scalefactor (- (:start (:gap gap)) translate))
            :y 30
            :height 20
            :width (* scalefactor (- (:end (:gap gap)) (:start (:gap gap))))
            :fill "red"
            :on-click #(re-frame/dispatch [:toggle-gap id])}]))


(defn form-view []
  [:div {:class "row" }
    [:h3  "Step 1: Upload Your Genbank File"]
    [:hr]
    [:form {:id "gbkform" :method "post" :enc-type "multipart/form-data"}
     [:input {:type "file" :name "file"}]
     [:br]
     [:div {:class "button"
               :on-click (fn []
                           (let [form-data
                                 (js/FormData. (js/document.querySelector "form"))]
                             (re-frame/dispatch [:upload-file form-data])))}
        "Upload"]]])

(defn gap-view [gap]
  "keep track of gaps as a list"
  (let [id (:id gap)
        start (:start (:gap gap))
        end (:end (:gap gap))
        size (- end start)]
    [:div {:class "panel panel-default"}
     (str "Gap ID=" id " Size: " size)

   ;{:id 1
   ; :protein1 {:start 150 :end 180 :sequence "ACACAC" :strand 1}
   ; :gap      {:start 190 :end 220 :sequence "ACACACACCCCCC" :strand 1}
   ; :protein2 {:start 230 :end 240 :sequence "ACACACACCCCCC" :strand 1}
   ]))

(defn minigapssvg [gaps scalefactor]
  (let [widths (map #(* scalefactor (- (:end (:gap %)) (:start (:gap %)))) gaps)
        startpositions (reduce (fn [a b] (conj a  (+ (last a) b))) [] widths)
        status (map #(if (= (:selected %) false) "inactive" "active") gaps)]
    [:svg {:width 933}
      (doall (for [i (range (count gaps))]
               (let [width (get widths i)
                     klass (get status i)
                     start (get startpositions i)]
                 [:rect {:class klass
                         :x start
                         :y 10
                         :width 10
                         :fill "red"}])))]))

(defn svg-main []
  (let [active-gaps (re-frame/subscribe [:active-gaps])
        genes (re-frame/subscribe [:genes])
        gaps (re-frame/subscribe [:gaps])
        translate (re-frame/subscribe [:svg-translate])
        scale (re-frame/subscribe [:svg-scale])
        scale-min (re-frame/subscribe [:svg-scale/min])
        scale-max (re-frame/subscribe [:svg-scale/max])
        scale-step (re-frame/subscribe [:svg-scale/step])
        sequence-end (re-frame/subscribe [:sequence-end])
        width (re-frame/subscribe [:svg-width])]

    [:div
     [:h3  "Step 2: Choose CRISPR Sites"]
     "Gaps are represented below by red blocks. You may choose a maximum of seven
     gaps to design crispr sites for. To inactivate a gap, simple click on it and it will
     change color to blue. When you have selected the gaps of interest, click the
     'Process Gaps' button."
     [:hr]
     "Scale"
     [slider :svg-scale scale @scale-min @scale-max @scale-step]
     "Translate"
     [slider :svg-translate translate (- @sequence-end)  @sequence-end 1]
     [:svg {:width @width}
      (doall (for [gene @genes] [drawgene gene @scale @translate]))
      (doall (for [gap @gaps] [drawgap gap @scale @translate]))]

     ;list of active gaps
     [:h4 "Active Gap Count: " (count @active-gaps)]

     [:div {:class "button"
            :on-click #(re-frame/dispatch [:make-cassettes])} "Process Gaps"]]))


(defn cassettes-main []
  (let [cassettes (re-frame/subscribe [:cassettes])]
    [:div {:class "row"}
     [:h3 {:style {:text-align "center"}} "Crispr Constructs"]
     [:table {:class "u-full-width"}
      [:thead  [:tr [:th "CrisprArray"] [:th "Sequence"]]]
      [:tbody
       [:tr [:th] [:th [:img {:src "/static/img/cassettes-01.png"
                              :style {:text-align "center" :padding "15px"}}]]]

       (for [[idx cass] (keep-indexed vector @cassettes)]
         [:tr
          [:td (str "Cassette " idx)]
          [:td
            (for [feat cass]
              [:span {:class (:type feat)} (:sequence feat)])]])]]]))

(defn primers-main []
  (let [primers (re-frame/subscribe [:primers])]
    [:div
     [:h3 {:style {:text-align "center"}} "Bridge Primers"]
     (for [primer @primers]
       [:div
         [:div {:class "row"}
          [:div {:class "four columns"} "Promoter Cassete: " (:promoterid primer)]
          [:div {:class "six columns"} (:selection primer) ]]
         [:table {:class "u-full-width"}
          [:tr [:td "Forward Primer"] [:td (:forward primer)]]
          [:tr [:td "Reverse Primer"] [:td (:reverse primer)]]]])]))

(defn main-panel []
  (let [app-cycle (re-frame/subscribe [:app-cycle])]
    [:div
     (when (= 0 @app-cycle) [form-view])
     (when (= @app-cycle 1) [svg-main])
     (when (> @app-cycle 1)
       [:div
        [:h3  "Step 3: Obtain Data"]
        "There are two basic outputs to mCRISTAR. The first is the crispr array that will be
        cloned into the pCRCT vector using golden gate assembly. Most likely you will order these
        as gene blocks from a vendor. Although only the sequence is required, the arrays are
        color coded as per the cartoon below. The second output are the bridge primers needed to generate
        the bridge fragments which will supply the new promoters for your CRISPR-cut vector. These
        will be cotransformed into your yeast strain harboring your pCRCT-CRISPRARRAY plasmid."
        [:hr]
        [cassettes-main]
        [primers-main]])]))


var ss1 = null, ss2 = null;
    
$( document ).ready(function() {
        
  var seq = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA";

  var features = [
        {
            category : 'secondary structure',
            type : 'helix',
            start : 34,
            end : 64
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 71,
            end : 100
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 106,
            end : 140
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 150,
            end : 172
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 200,
            end : 230
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 241,
            end : 276
        }, {
            category : 'secondary structure',
            type : 'helix',
            start : 286,
            end : 309
        }];

      var seqEntry = new pviz.SeqEntry({
          sequence : seq
      });

      var view = new pviz.SeqEntryAnnotInteractiveView({
          model : seqEntry,
          el : '#pviz'
      })
      view.render();

      /*
       * add a mouseover/mouseout call back based on the feature type
       */
      pviz.FeatureDisplayer.addMouseoverCallback(['helix', 'turn', 'beta_strand'], function(ft) {
          console.log("mouseover: "+JSON.stringify(ft))
      }).addClickCallback(['helix', 'turn', 'beta_strand'], function(ft) {
          console.log("click: "+JSON.stringify(ft));
          var tm_index;
          $.each(features, function(i, f){
            if(f.start == ft.start && f.end == ft.end && f.category == ft.category && f.type == ft.type){
                tm_index = i+1
            }
          });
          if(ss1 === null || ss2 != null){
            ss1 = tm_index;
            ss2 = null;
            create_circle_flow('#chart','TM'+ss1+'_TM'+ss1+'_interactions.json');
          }
          else{
            ss2 = tm_index;
            create_circle_flow('#chart','TM'+ss1+'_TM'+ss2+'_interactions.json');
          }
          });

      seqEntry.addFeatures(features);
    });
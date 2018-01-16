var viewer;
var structure;
var geom;
var sheets = true;
var strands = false;

function color_face(sheet_num, sheet_residues, sheet_chain, color, recenter, showLines){
  var residues = [].concat.apply([], sheet_residues);
  var sheet = structure.select({chain: sheet_chain, rindices:residues});
  if(showLines){
    viewer.lines("res_lines_"+sheet_num, sheet);
  }
  geom.colorBy(pv.color.uniform(color), sheet);
  if(recenter){
    viewer.centerOn(sheet);
  }
}

function color_faces(face1_residues, face1_chain, face2_residues, face2_chain) {
  if(!sheets && strands){
    viewer.rm("res_lines_*");
  }
  geom.colorBy(pv.color.uniform("grey"), structure);
  color_face(1, face1_residues, face1_chain, "#5cb85c", true, (!sheets && strands));
  color_face(2, face2_residues, face2_chain, "#428bca", false, (!sheets && strands));
  viewer.centerOn(structure);
  viewer.requestRedraw();
}

function reset_colors(){
  viewer.rm("res_lines_*");
  geom.colorBy(pv.color.uniform("grey"), structure);
}

function show_error(msg){
  $("#interaction_chart").html('<div class="alert alert-danger" role="alert">'+msg+'</div>');
}

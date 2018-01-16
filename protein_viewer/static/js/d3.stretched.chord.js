////////////////////////////////////////////////////////////
/////////////// Custom Chord Function //////////////////////
//////// Pulls the chords pullOutSize pixels apart /////////
////////////////// along the x axis ////////////////////////
////////////////////////////////////////////////////////////
///////////// Created by Nadieh Bremer /////////////////////
//////////////// VisualCinnamon.com ////////////////////////
////////////////////////////////////////////////////////////
//// Adjusted from the original d3.svg.chord() function ////
///////////////// from the d3.js library ///////////////////
//////////////// Created by Mike Bostock ///////////////////
////////////////////////////////////////////////////////////

stretchedChord = function() {
    var source = d3_source, 
    target = d3_target, 
    radius = d3_svg_chordRadius, 
    startAngle = d3_svg_arcStartAngle, 
    endAngle = d3_svg_arcEndAngle,
    pullOutSize = 0;
    
  var π = Math.PI,
    halfπ = π / 2;

    function subgroup(self, f, d, i) {
    var subgroup = f.call(self, d, i), 
      r = radius.call(self, subgroup, i), 
      a0 = startAngle.call(self, subgroup, i) - halfπ, 
      a1 = endAngle.call(self, subgroup, i) - halfπ;
      return {
        r: r,
        a0: [a0],
        a1: [a1],
        p0: [ r * Math.cos(a0), r * Math.sin(a0)],
        p1: [ r * Math.cos(a1), r * Math.sin(a1)]
      };
    }

    function arc(r, p, a) {
    var sign = (p[0] >= 0 ? 1 : -1);
    return "A" + r + "," + r + " 0 " + +(a > π) + ",1 " + (p[0] + sign*pullOutSize) + "," + p[1];
    }


    function curve(p1) {
    var sign = (p1[0] >= 0 ? 1 : -1);
    return "Q 0,0 " + (p1[0] + sign*pullOutSize) + "," + p1[1];
    }
  
  /*
  M = moveto
  M x,y
  Q = quadratic Bézier curve
  Q control-point-x,control-point-y end-point-x, end-point-y
  A = elliptical Arc
  A rx, ry x-axis-rotation large-arc-flag, sweep-flag  end-point-x, end-point-y
  Z = closepath

  M251.5579641956022,87.98204731514328
  A266.5,266.5 0 0,1 244.49937503334525,106.02973926358392
  Q 0,0 -177.8355222451483,198.48621369706098
  A266.5,266.5 0 0,1 -191.78901944612068,185.0384338992728
  Q 0,0 251.5579641956022,87.98204731514328
  Z
  */  
    function chord(d, i) {
    var s = subgroup(this, source, d, i), 
      t = subgroup(this, target, d, i);
          
  return "M" + (s.p0[0] + pullOutSize) + "," + s.p0[1] + 
      arc(s.r, s.p1, s.a1 - s.a0) + 
      curve(t.p0) + 
      arc(t.r, t.p1, t.a1 - t.a0) + 
      curve(s.p0) + 
      "Z";
   }//chord

    chord.radius = function(v) {
      if (!arguments.length) return radius;
      radius = d3_functor(v);
      return chord;
    };
    chord.pullOutSize = function(v) {
      if (!arguments.length) return pullOutSize;
      pullOutSize = v;
      return chord;
    };
    chord.source = function(v) {
      if (!arguments.length) return source;
      source = d3_functor(v);
      return chord;
    };
    chord.target = function(v) {
      if (!arguments.length) return target;
      target = d3_functor(v);
      return chord;
    };
    chord.startAngle = function(v) {
      if (!arguments.length) return startAngle;
      startAngle = d3_functor(v);
      return chord;
    };
    chord.endAngle = function(v) {
      if (!arguments.length) return endAngle;
      endAngle = d3_functor(v);
      return chord;
    };


  function d3_svg_chordRadius(d) {
      return d.radius;
  }

  function d3_source(d) {
    return d.source;
  }
    
  function d3_target(d) {
      return d.target;
  }

  function d3_svg_arcStartAngle(d) {
      return d.startAngle;
  }
    
  function d3_svg_arcEndAngle(d) {
      return d.endAngle;
  }

  function d3_functor(v) {
    return typeof v === "function" ? v : function() {
      return v;
    };
  }

  return chord;

}//stretchedChord
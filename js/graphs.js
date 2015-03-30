
/* global d3, _, MathJax, Matrix */

function updateTex() {
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

$(function() {
  var width, height;

  var svg = d3.select('#canvas');

  width = $('#canvas').width();
  height = $('#canvas').height();

  // TODO: add option to freeze vertices
  var force = d3.layout.force()
          .size([width, height])
          .linkDistance(100)
          .linkStrength(0.1)
          .charge(-300)
          .on('tick', function() {
            linkObjs
                    .attr('x1', function(d) {
                      return d.source.x;
                    })
                    .attr('y1', function(d) {
                      return d.source.y;
                    })
                    .attr('x2', function(d) {
                      return d.target.x;
                    })
                    .attr('y2', function(d) {
                      return d.target.y;
                    });
            nodeObjs
                    .attr('cx', function(d) {
                      return d.x;
                    })
                    .attr('cy', function(d) {
                      return d.y;
                    });
          });

  var nodes = force.nodes();
  var links = force.links();

  var drag = force.drag();

  function selectNode(node) {
    node.classed('node-selected', true)
            .transition()
            .duration(300)
            .attr('r', 13);
  }

  function deselectNode(node) {
    node.classed('node-selected', false)
            .transition()
            .duration(100)
            .attr('r', 8);
  }

  function clear() {
    nodes.splice(0, nodes.length);
    links.splice(0, links.length);
    deselectNode(svg.selectAll('.node-selected'));
  }

  var canvasGroup = svg.append('g');

  canvasGroup.append('rect')
          .attr('class', 'bg')
          .attr('width', width)
          .attr('height', height)
          .on('click', function() {
            var selected = svg.selectAll('.node-selected');
            if (selected.empty()) {
              var m = d3.mouse(this);
              nodes.push({ x: m[0], y: m[1], vertNum: nodes.length });
              update();
            } else {
              deselectNode(selected);
            }
          });

  var nodeObjs = canvasGroup.selectAll('.node');
  var linkObjs = canvasGroup.selectAll('.link');

  function update() {
    linkObjs = linkObjs.data(links);
    linkObjs.enter()
            .insert('line', '.node')
            .attr('class', 'link')
            .attr('opacity', 0)
            .transition()
            .duration(200)
            .attr('opacity', 1);
    linkObjs.exit()
            .transition()
            .duration(200)
            .attr('opacity', 0)
            .remove();

    nodeObjs = nodeObjs.data(nodes);
    nodeObjs.enter()
            .append('circle')
            .attr('class', 'node')
            .attr('r', 0)
            .call(drag)
            .on('click', function() {
              if (d3.event.defaultPrevented)
                return;
              console.log('clicked');
              var selected = svg.selectAll('.node-selected');
              var a = d3.select(this).datum(), b = selected.empty() ? null : selected.datum();
              if (d3.event.ctrlKey) {
                console.log('remove');
                if (!selected.empty() && a === b)
                  deselectNode(selected);
                nodes.splice(nodes.indexOf(a), 1);
                nodes.forEach(function(node, i) {
                  node.vertNum = i;
                });
                links = _.filter(links, function(link) {
                  return !(link.source === a || link.target === a);
                });
                links.forEach(function(link, i) {
                  link.edgeNum = i;
                });
                update();
              }
              else {
                if (!selected.empty()) {
                  deselectNode(selected);
                  var adj = isAdj(a, b);
                  if (a !== b) {
                    if (!adj) {
                      console.log('link');
                      links.push({ source: selected.datum(), target: d3.select(this).datum(), edgeNum: links.length });
                    }
                    else {
                      console.log('unlink');
                      links = _.filter(links, function(link) {
                        return !((link.source === a && link.target === b) || (link.source === b && link.target === a));
                      });
                    }
                    update();
                  }
                }
                if (selected.empty() || d3.event.shiftKey) {
                  console.log('select');
                  selectNode(d3.select(this));
                }
              }
            })
            .transition()
            .duration(200)
            .attr('r', 8);
    nodeObjs.exit()
            .transition()
            .duration(200)
            .attr('r', 0)
            .remove();

    force.start();

    $('#graph-prop-table tbody tr').each(function() {
      var name = $(this).attr('prop-name');
      console.log('updating prop ' + name);
      var val = graphProps[name].calc(nodes, links);
      graphProps[name].value = val;
      $(this).prop('__raw-value__', val);
      console.log(val);
      // TODO: change so that only sets value if shown, that way MathJax isn't calculating things that aren't shown
      if (graphProps[name].repr) {
        val = graphProps[name].repr(val);
      }
      else if (val === undefined) {
        val = '';
      }
      else if (val === null) {
        val = '?';
      }
      else if (isNaN(val)) {
        if (val instanceof Matrix) {
          val = '\\(' + val.toLatex() + '\\)';
        }
        else {
          val = '' + val;
        }
      }
      else if (!isFinite(val)) {
        val = '\\(' + (val === Number.NEGATIVE_INFINITY ? '-' : '') + '\\infty\\)';
      }
      else {
        val = val.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
      }
      $(this).find('td .value').text(val);
    });
    $('#graph-prop-table tbody td .value').filter(function() {
      return $(this).text() === 'true';
    }).css({ 'color': '#0a0' });
    $('#graph-prop-table tbody td .value').filter(function() {
      return $(this).text() === 'false';
    }).css({ 'color': '#a00' });
    updateTex();
  }

  d3.xml('img/erase.svg', 'image/svg+xml', function(xml) {
    svg.each(function() {
      d3.select(this.appendChild(xml.documentElement))
              .attr('class', 'erase-btn')
              .attr('width', '40')
              .attr('height', '40')
              .attr('x', '5')
              .attr('y', '5')
              .on('mouseover', function() {
                d3.select(this).select('svg circle')
                        .attr('fill', 'red');
              })
              .on('mouseout', function() {
                d3.select(this).select('svg circle')
                        .attr('fill', 'none');
              })
              .on('click', function() {
                clear();
                update();
              });
    });
  });

  var graphProps = {
    order: {
      desc: 'order',
      longDesc: 'The number of vertices of the graph.',
      math: '\\(\\left|V\\right|\\)',
      calc: function(vs, es) {
        return vs.length;
      }
    },
    size: {
      desc: 'size',
      longDesc: 'The number of edges of the graph.',
      math: '\\(\\left|E\\right|\\)',
      calc: function(vs, es) {
        return es.length;
      }
    },
    adjMatr: {
      desc: 'adjacency matrix',
      math: '\\(A\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        return Matrix.byFunc(vs.length, vs.length, function(i, j) {
          return +(_.any(es, function(e) {
            return (e.source === vs[i] && e.target === vs[j]) || (e.source === vs[j] && e.target === vs[i]);
          }));
        });
      }
    },
    degSeq: {
      desc: 'degree sequence',
      longDesc: 'The sequence of degrees of the vertices of the graph. The degree of a vertex is the number of edges it is connected to.',
      defaultHidden: true,
      calc: function(vs, es) {
        return _.sortBy(_.map(graphProps.adjMatr.value.values, function(row, i) {
          return {
            v: vs[i],
            deg: _.reduce(row, function(sum, adj) {
              return sum + adj;
            })
          };
        }), function(v) {
          return v.deg;
        });
      },
      repr: function(value) {
        return '(' + _.map(value, function(v) {
          return v.deg;
        }).join(', ') + ')';
      }
    },
    degMatr: {
      desc: 'degree matrix',
      math: '\\(D\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        return Matrix.byFunc(vs.length, vs.length, function(i, j) {
          if (i !== j)
            return 0;
          else
            return _.find(graphProps.degSeq.value, function(v) {
              return v.v === vs[i];
            }).deg;
        });
      }
    },
    laplMatr: {
      desc: 'Laplacian matrix',
      math: '\\(L\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        return graphProps.degMatr.value.sub(graphProps.adjMatr.value);
      }
    },
    inciMatr: {
      desc: 'incidence matrix',
      math: '\\(M\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        if (es.length === 0)
          return Matrix.zero(1, vs.length);
        return Matrix.byFunc(es.length, vs.length, function(e, v) {
          var i = es[e].source.vertNum;
          var j = es[e].target.vertNum;
          if (i < j) {
            var temp = i;
            i = j;
            j = temp;
          }
          if (v === i)
            return 1;
          if (v === j)
            return -1;
          return 0;
        });
      }
    },
    minDeg: {
      desc: 'min degree',
      math: '\\(\\delta\\)',
      calc: function(vs, es) {
        return vs.length === 0 ? Number.POSITIVE_INFINITY : graphProps.degSeq.value[0].deg;
      }
    },
    maxDeg: {
      desc: 'max degree',
      math: '\\(\\Delta\\)',
      calc: function(vs, es) {
        return vs.length === 0 ? Number.NEGATIVE_INFINITY : graphProps.degSeq.value[vs.length - 1].deg;
      }
    },
    regular: {
      desc: 'is regular',
      longDesc: 'Whether all vertices have the same degree or not.',
      calc: function(vs, es) {
        return graphProps.minDeg.value === graphProps.maxDeg.value;
      }
    },
    components: {
      desc: 'number of components',
      longDesc: 'The number of connected subgraphs.',
      calc: function(vs, es) {
        var count = 0;
        var seen = []; // kep track of if we've looked at this vertex already
        for (var i = 0; i < vs.length; i++) {
          seen.push(false);
        }
        var queue = []; // queue of adjacent vertices to look at
        var vs2 = _.range(vs.length);
        vs2.forEach(function(v) {
          if (!seen[v]) { // if haven't already seen this vertex
            count++; // increment component count
            seen[v] = true; // indicate we've seen this vertex
            queue.push(v);
            while (queue.length !== 0) { // go through all adjacent vertices
              var w = queue.shift();
              var row = graphProps.adjMatr.value.values[w]; // get adjacency row for w
              vs2.forEach(function(u) {
                if (row[u] === 1 && !seen[u]) { // if u is adj to w and we haven't seen u yet
                  seen[u] = true;
                  queue.push(u);
                }
              });
            }
          }
        });
        return count;
      }
    },
    connected: {
      desc: 'is connected',
      longDesc: 'Whether the number of components is one or not.',
      calc: function(vs, es) {
        return graphProps.components.value === 1;
      }
    },
    numTris: {
      desc: 'number of triangles',
      calc: function(vs, es) {
        return graphProps.adjMatr.value.pow(3).trace() / 6;
      }
    },
    bipartite: {
      desc: 'is bipartite',
      longDesc: 'Whether the vertex set can be partitioned into two disjoint subsets such that all edges are only between different subsets or not.',
      calc: function(vs, es) {
        if (!graphProps.connected.value || graphProps.numTris.value !== 0) // not connected or has triangles, not bipartite
          return { bipartite: false };
        var vs2 = _.range(vs.length);
        var color = []; // color array, 0 means uncolored, 1 and -1 are colors
        for (var i = 0; i < vs.length; i++) {
          color.push(0);
        }
        color[0] = 1; // color first vertex 1
        var queue = [0];
        while (queue.length !== 0) {
          var v = queue.shift();
          var row = graphProps.adjMatr.value.values[v]; // get adjacency row for v
          var uncolorable = _.any(vs2, function(u) {
            if (row[u] === 1) { // if u adj to v
              if (color[u] === 0) { // if u uncolored
                color[u] = -color[v]; // color opposite to v
                queue.push(u);
              }
              else if (color[u] === color[v]) { // if u is colored the same as v, not 2-colorable, not bipartite
                return true;
              }
            }
            return false;
          });
          if (uncolorable)
            return { bipartite: false };
        }
        var ret = _.countBy(color, function(c) {
          return c === 1 ? 'm' : 'n';
        });
        ret.bipartite = true;
        return ret;
      },
      repr: function(value) {
        return value.bipartite;
      }
    },
    spanningTrees: {
      desc: 'number of spanning trees',
      longDesc: 'The number of trees that contain all vertices of the graph and whose only edges are edges of the graph. A tree is a graph without any cycles.',
      calc: function(vs, es) {
        if (!graphProps.connected.value)
          return 0;
        if (vs.length === 1)
          return 1;
        if (es.length === vs.length * (vs.length - 1) / 2)
          return Math.pow(vs.length, vs.length - 2);
        if (graphProps.bipartite.value.bipartite) {
          var m = graphProps.bipartite.value.m;
          var n = graphProps.bipartite.value.n;
          if (es.length === m * n)
            return Math.pow(m, n - 1) * Math.pow(n, m - 1);
        }
        // TODO: broken for trees of order > 7?
        // TODO: not always exact, have to round, even then sometimes off (by a lot), maybe could do with fractions?
        return Math.floor(graphProps.laplMatr.value.cofactor(0, 0));
      }
    },
    tree: {
      desc: 'is a tree',
      longDesc: 'Whether the graph is connected and has no cycles or not.',
      calc: function(vs, es) {
        return graphProps.connected.value && graphProps.size.value === graphProps.order.value - 1;
      }
    }
    // TODO: num cycles, is tree, diameter
  };

  // TODO: show properties of selected vertex
  var vertProps = {
    deg: {
      desc: 'degree',
      math: '\\(\\mathrm{deg}\\)',
      value: function(v, vs, es) {
        return _.find(graphProps.degSeq.value, function(vert) {
          return vert.v === v;
        }).deg;
      }
    }
  };

  var graphLayouts = {
    radial: function(n) {
      return function(v) {
        var r = 100;
        var t = v * 2 * Math.PI / n - Math.PI / 2;
        return { x: r * Math.cos(t), y: r * Math.sin(t) };
      };
    },
    horizLines: function(n, topCount) {
      return function(v) {
        var sep = 100;
        if (v < topCount)
          return { x: v / (topCount - 1) * 2 * sep - sep, y: -sep };
        else
          return { x: (v - topCount) / (n - topCount - 1) * 2 * sep - sep, y: sep };
      };
    }
  };

  var commonGraphs = {
    complete: {
      desc: 'complete graph',
      math: '\\(K_n\\)',
      args: {
        'n': {
          min: 1,
          max: Number.POSITIVE_INFINITY,
          def: 5
        }
      },
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        for (var i = 0; i < n; i++) {
          for (var j = i + 1; j < n; j++) {
            es.push([i, j]);
          }
        }
        return { vs: vs, es: es, layout: graphLayouts.radial(n) };
      }
    },
    cycle: {
      desc: 'cycle graph',
      math: '\\(C_n\\)',
      args: {
        'n': {
          min: 1,
          def: 5
        }
      },
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        for (var i = 0; i < n; i++) {
          es.push([i, i + 1 >= n ? 0 : i + 1]);
        }
        return { vs: vs, es: es, layout: graphLayouts.radial(n) };
      }
    },
    path: {
      desc: 'path graph',
      math: '\\(P_n\\)',
      args: {
        'n': {
          min: 1,
          def: 4
        }
      },
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        for (var i = 0; i < n - 1; i++) {
          es.push([i, i + 1]);
        }
        return { vs: vs, es: es, layout: graphLayouts.radial(n) };
      }
    },
    completeBipartite: {
      desc: 'complete bipartite graph',
      math: '\\(K_{m,n}\\)',
      args: {
        'm': {
          min: 1,
          def: 4
        },
        'n': {
          min: 1,
          def: 3
        }
      },
      make: function(m, n) {
        var vs = _.range(n + m);
        var es = [];
        for (var i = 0; i < m; i++) {
          for (var j = 0; j < n; j++) {
            es.push([i, m + j]);
          }
        }
        return { vs: vs, es: es, layout: graphLayouts.horizLines(n + m, m) };
      }
    },
    hypercube: {
      desc: 'hypercube graph',
      math: '\\(Q_n\\)',
      args: {
        'n': {
          min: 0,
          def: 3
        }
      },
      make: function(n) {
        n = 1 << n;
        var vs = _.range(n);
        var es = [];
        for (var i = 0; i < n; i++) {
          for (var j = i + 1; j < n; j++) {
            var dist = 0;
            var val = i ^ j;
            while (val !== 0 && dist < 2) {
              dist++;
              val &= val - 1;
            }
            if (dist === 1)
              es.push([i, j]);
          }
        }
        return { vs: vs, es: es, layout: graphLayouts.radial(n) };
      }
    }
  };

  function isAdj(v1, v2) {
    return graphProps.adjMatr.value.get(v1.vertNum, v2.vertNum) === 1;
  }

  _.each(graphProps, function(prop, name) {
    // TODO: make longDesc some kind of fancy popup, also so can show math and links
    $('#graph-prop-table tbody').append('<tr prop-name="' + name + '"><td class="desc">' + prop.desc + (prop.longDesc ? ' <abbr title="' + prop.longDesc + '">(?)</abbr>' : '') + '</td><td class="math">' + (prop.math || '') + '</td><td class="value">' + (prop.defaultHidden ? '<span class="display-toggle">[show]</span><br>' : '') + '<span class="value"></span></td></tr>');
  });
  $('#graph-prop-table tbody td .display-toggle ~ .value').hide();
  $('#graph-prop-table tbody td .display-toggle').click(function(event) {
    var valueSpan = $(this).siblings('.value');
    if ($(this).text() === '[show]') {
      $(this).text('[hide]');
      valueSpan.show(500);
    } else {
      $(this).text('[show]');
      valueSpan.hide(500);
    }
  });
  _.each(commonGraphs, function(graph, name) {
    $('#common-graphs-table tbody').append('<tr graph-name="' + name + '"><td class="desc">' + graph.desc + '</td><td class="math">' + (graph.math || '') + '</td><td class="inputs"></td></tr>');
  });
  $('#common-graphs-table tbody tr').each(function() {
    var graphName = $(this).attr('graph-name');
    var graph = commonGraphs[graphName];
    var inputs = $(this).children('td.inputs');
    inputs.append('<form class="form-inline"></form>');
    inputs.submit(function() {
      console.log('make graph: ' + graph.desc);
      var args = [];
      $(this).find('input').each(function() {
        var argName = $(this).attr('id');
        argName = argName.substring(argName.lastIndexOf('-') + 1);
        var argVal = +$(this).val();
        args.push(argVal);
        console.log(argName + ': ' + argVal);
      });
      var g = graph.make.apply(null, args);
      console.log(g.vs.length + ' verts: ' + JSON.stringify(g.vs));
      console.log(g.es.length + ' edges: ' + JSON.stringify(g.es));
      console.log(g.layout);
      clear();
      g.vs.forEach(function(v) {
        var pos = g.layout(v);
        console.log(v + ' ' + JSON.stringify(pos));
        nodes.push({ x: pos.x + width / 2, y: pos.y + height / 2, vertNum: v });
      });
      g.es.forEach(function(e) {
        links.push({ source: nodes[e[0]], target: nodes[e[1]], edgeNum: e });
      });
      update();
      return false;
    });
    inputs = inputs.children('form');
    _.each(graph.args, function(arg, argName) {
      inputs.append('<div class="form-group" var-name="' + argName + '"><label for="graph-' + graphName + '-' + argName + '">' + argName + ':</label><input class="form-control" id="graph-' + graphName + '-' + argName + '" type="number"' + (arg.min === undefined ? '' : ' min="' + arg.min + '"') + (arg.max === undefined ? '' : ' max="' + arg.max + '"') + (arg.def === undefined ? '' : ' value="' + arg.def + '"') + '></input></div>');
    });
    inputs.append('<button type="submit" class="btn btn-default">Create</button>');
  });

  update();
});

$(document).on('click', 'a', function(event) {
  if ($(this).attr('href') === '#') {
    event.preventDefault();
    console.log('Default "#" link action prevented.');
  }
});
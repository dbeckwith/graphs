
/* global d3, _, MathJax, Matrix */

function updateTex() {
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

$(function() {

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
    // TODO: num cycles, is tree, diameter, vertex and edge connectivity, non-separable, girth, Eulerian
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
      verts: '\\(n\\)',
      edges: '\\(\\frac{n\\left(n-1\\right)}2\\)',
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
      verts: '\\(n\\)',
      edges: '\\(n\\)',
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
      verts: '\\(n\\)',
      edges: '\\(n-1\\)',
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
      verts: '\\(m+n\\)',
      edges: '\\(mn\\)',
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
      verts: '\\(2^n\\)',
      edges: '\\(2^{n-1}n\\)',
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

  var width, height;
  // TODO: move things when canvas resizes

  var svg = d3.select('#canvas');

  width = $('#canvas').width();
  height = $('#canvas').height();

  // data for nodes and links
  var nodes = [];
  var links = [];

  // node drag behavior
  var drag = d3.behavior.drag()
          .on('drag', function(d) {
            d.x = d3.event.x;
            d.y = d3.event.y;
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

  // sets the given node as selected
  function selectNode(node) {
    node.classed('node-selected', true)
            .transition()
            .duration(300)
            .attr('r', 13);
  }

  // sets the given node as deselected
  function deselectNode(node) {
    node.classed('node-selected', false)
            .transition()
            .duration(100)
            .attr('r', 8);
  }

  // clears all data
  function clear() {
    nodes.splice(0, nodes.length);
    links.splice(0, links.length);
    deselectNode(svg.selectAll('.node-selected'));
  }

  // setup draw layers
  var passiveLayer = svg.append('g').attr('class', 'passive-layer');
  var interactLayer = svg.append('g').attr('class', 'interact-layer');
  interactLayer.append('rect')
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

  // setup nodes and links
  var nodeObjs = interactLayer.append('g').selectAll('.node');
  var linkObjs = passiveLayer.append('g').selectAll('.link');

  // updates properties table
  function calcProps() {
    graphPropRow.each(function(prop) {
      console.log('updating prop ' + prop.name);
      var val = prop.value = prop.calc(nodes, links);
      console.log(val);
      if (prop.repr) {
        val = prop.repr(val);
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
        // TODO: val.toString() can sometimes produce scientific notation, don't use commas then
        val = val.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
      }
      var valueSpan = d3.select(this).select('td .value');
      valueSpan.text(val);
      if (valueSpan.text() === 'true')
        valueSpan.style('color', '#0a0');
      else if (valueSpan.text() === 'false')
        valueSpan.style('color', '#a00');
    });
    updateTex();
  }

  // updates view
  function update() {
    linkObjs = linkObjs.data(links);
    linkObjs.enter()
            .insert('line', '.node')
            .attr('class', 'link')
            .attr('opacity', 0)
            .transition()
            .duration(200)
            .attr('opacity', 1);
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
    nodeObjs
            .attr('cx', function(d) {
              return d.x;
            })
            .attr('cy', function(d) {
              return d.y;
            });
    nodeObjs.exit()
            .transition()
            .duration(200)
            .attr('r', 0)
            .remove();

    if (d3.select('#auto-calc-props-checkbox').property('checked'))
      calcProps();
  }

  // tests if two vertices are adjacent
  function isAdj(v1, v2) {
    return _.any(links, function(link) {
      return (link.source === v1 && link.target === v2) || (link.source === v2 && link.target === v1);
    });
  }

  // setup toolbar
  var toolbarButtonSize = 40;
  var toolbar = interactLayer.append('g').attr('class', 'toolbar')
          .selectAll('.toolbar-btn')
          .data([
            {
              img: 'img/help.svg',
              onclick: function() {
                $('#help').css({ 'background-color': 'hsla(240,100%,90%,0.9)' });
                $('body').scrollTo('#help', 500, function() {
                  $('#help').animate({ 'background-color': 'hsla(240,100%,90%,0.0)' }, 500);
                });
              },
              title: 'Help'
            },
            {
              img: 'img/erase.svg',
              onclick: function() {
                clear();
                update();
              },
              title: 'Clear'
            },
            {
              img: 'img/calc.svg',
              onclick: function() {
                calcProps();
              },
              title: 'Recalculate properties'
            }
          ])
          .enter()
          .append('g')
          .attr('class', 'toolbar-btn')
          .attr('transform', function(d, i) {
            return  'translate(' + (5 + (5 + toolbarButtonSize) * i) + ',5)';
          });
  toolbar.append('title')
          .text(function(d) {
            return d.title;
          });
  toolbar.append('circle')
          .attr('cx', toolbarButtonSize / 2)
          .attr('cy', toolbarButtonSize / 2)
          .attr('r', toolbarButtonSize / 2)
          .attr('fill', 'red')
          .attr('fill-opacity', '0.0')
          .attr('stroke', 'black')
          .attr('stroke-width', '2')
          .on('mouseover', function(d) {
            d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('fill-opacity', '1.0');
          })
          .on('mouseout', function() {
            d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('fill-opacity', '0.0');
          })
          .on('click', function(d, i) {
            d.onclick.apply(this, [d, i]);
          });
  toolbar.append('image')
          .attr('xlink:href', function(d) {
            return d.img;
          })
          .attr('width', toolbarButtonSize)
          .attr('height', toolbarButtonSize);


  // setup tables

  // TODO: make longDesc some kind of fancy popup, also so can show math and links
  var graphPropRow = d3.select('#graph-prop-table tbody').selectAll('tr')
          .data(_.map(graphProps, function(prop, name) {
            prop.name = name;
            return prop;
          }))
          .enter()
          .append('tr');
  graphPropRow.append('td')
          .html(function(prop) {
            var body = prop.desc;
            if (prop.longDesc)
              body += ' <abbr title="' + prop.longDesc + '">(?)</abbr>';
            return body;
          });
  graphPropRow.append('td')
          .text(function(prop) {
            return prop.math || '';
          });
  graphPropRow.append('td')
          .attr('class', 'value')
          .each(function(prop) {
            var cell = d3.select(this);
            var toggle = null;
            if (prop.defaultHidden)
              toggle = cell.append('span')
                      .attr('class', 'display-toggle')
                      .text('[show]');
            var valueSpan = cell.append('span')
                    .attr('class', 'value');
            if (prop.defaultHidden) {
              valueSpan.classed('tex2jax_ignore', true)
                      .style('display', 'none');
              toggle.on('click', function() {
                if (toggle.text() === '[show]') {
                  toggle.text('[hide]');
                  valueSpan.classed('tex2jax_ignore', false);
                  updateTex();
                  $(valueSpan[0]).show(500);
                }
                else {
                  toggle.text('[show]');
                  valueSpan.classed('tex2jax_ignore', true);
                  updateTex();
                  $(valueSpan[0]).hide(500);
                }
              });
            }
          });

  var commonGraphsRow = d3.select('#common-graphs-table tbody').selectAll('tr')
          .data(_.map(commonGraphs, function(prop, name) {
            prop.name = name;
            return prop;
          }))
          .enter()
          .append('tr');
  commonGraphsRow.append('td')
          .text(function(graph) {
            return graph.desc;
          });
  commonGraphsRow.append('td')
          .text(function(graph) {
            return graph.math || '';
          });
  commonGraphsRow.append('td')
          .text(function(graph) {
            return graph.verts;
          });
  commonGraphsRow.append('td')
          .text(function(graph) {
            return graph.edges;
          });
  commonGraphsRow.append('td')
          .attr('class', 'inputs')
          .each(function(graph) {
            var cell = d3.select(this);
            var form = cell.append('form')
                    .attr('class', 'form-inline');
            _.each(graph.args, function(arg, argName) {
              var div = form.append('div')
                      .attr('class', 'form-group');
              div.append('label')
                      .attr('for', 'graph-' + graph.name + '-' + argName)
                      .text(argName + ':');
              var input = div.append('input')
                      .attr('class', 'form-control')
                      .attr('id', 'graph-' + graph.name + '-' + argName)
                      .attr('type', 'number');
              if (arg.min !== undefined)
                input.attr('min', arg.min);
              if (arg.max !== undefined)
                input.attr('max', arg.max);
              if (arg.def !== undefined)
                input.attr('value', arg.def);
            });
            form.append('button')
                    .attr('type', 'submit')
                    .attr('class', 'btn btn-default')
                    .text('Create');
            form.on('submit', function() {
              console.log('make graph: ' + graph.name);
              var args = _.map(form.selectAll('input')[0], function(arg) {
                var argName = $(arg).attr('id');
                argName = argName.substring(argName.lastIndexOf('-') + 1);
                var argVal = +$(arg).val();
                console.log(argName + ': ' + argVal);
                return argVal;
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
              d3.event.preventDefault();
            });
          });

  update();
});

$(document).on('click', 'a', function(event) {
  if ($(this).attr('href') === '#') {
    event.preventDefault();
    console.log('Default "#" link action prevented.');
  }
});

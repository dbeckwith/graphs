
/* global d3, _, MathJax, Matrix */

function updateTex() {
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

// TODO: option to show labels on verts and edges
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
      longDesc: 'The matrix where the i,j<sup>th</sup> entry is a 1 if the i<sup>th</sup> and j<sup>th</sup> vertices are adjacent, or connected by an edge, and otherwise 0.',
      link: 'http://en.wikipedia.org/wiki/Adjacency_matrix',
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
      link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)#Degree_sequence',
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
      longDesc: 'The matrix where the i<sup>th</sup> diagonal entry is the degree of the i<sup>th</sup> vertex and all other entries are 0.',
      link: 'http://en.wikipedia.org/wiki/Degree_matrix',
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
      longDesc: 'The degree matrix of the graph minus the adjaceny matrix of the graph.',
      link: 'http://en.wikipedia.org/wiki/Laplacian_matrix',
      math: '\\(L\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        return graphProps.degMatr.value.sub(graphProps.adjMatr.value);
      }
    },
    inciMatr: {
      desc: 'incidence matrix',
      longDesc: 'The matrix where the i,j<sup>th</sup> entry is a 1 if the i<sup>th</sup> vertex and j<sup>th</sup> edge are indicent, or the vertex is part of the edge, and otherwise 0.',
      link: 'http://en.wikipedia.org/wiki/Incidence_matrix',
      math: '\\(M\\)',
      defaultHidden: true,
      calc: function(vs, es) {
        if (es.length === 0)
          return Matrix.zero(1, vs.length);
        return Matrix.byFunc(vs.length, es.length, function(v, e) {
          var i = es[e].source.vertNum;
          var j = es[e].target.vertNum;
          if (v === i || v === j)
            return 1;
          return 0;
        });
      }
    },
    minDeg: {
      desc: 'min degree',
      longDesc: 'The smallest degree of any vertex.',
      link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)',
      math: '\\(\\delta\\)',
      calc: function(vs, es) {
        return vs.length === 0 ? Number.POSITIVE_INFINITY : graphProps.degSeq.value[0].deg;
      }
    },
    maxDeg: {
      desc: 'max degree',
      longDesc: 'The largest degree of any vertex.',
      link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)',
      math: '\\(\\Delta\\)',
      calc: function(vs, es) {
        return vs.length === 0 ? Number.NEGATIVE_INFINITY : graphProps.degSeq.value[vs.length - 1].deg;
      }
    },
    regular: {
      desc: 'is regular',
      longDesc: 'Whether all vertices have the same degree or not.',
      link: 'http://en.wikipedia.org/wiki/Regular_graph',
      calc: function(vs, es) {
        return graphProps.minDeg.value === graphProps.maxDeg.value;
      }
    },
    complete: {
      desc: 'is complete',
      longDesc: 'Whether the graph has all possible edges or not.',
      link: 'http://en.wikipedia.org/wiki/Complete_graph',
      calc: function(vs, es) {
        return es.length === vs.length * (vs.length - 1) / 2;
      }
    },
    components: {
      desc: 'number of components',
      longDesc: 'The number of connected subgraphs.',
      link: 'http://en.wikipedia.org/wiki/Connected_component_(graph_theory)',
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
      link: 'http://en.wikipedia.org/wiki/Connectivity_(graph_theory)',
      calc: function(vs, es) {
        return graphProps.components.value === 1;
      }
    },
    diameter: {
      desc: 'diameter',
      longDesc: 'The longest length of the shortest path between any two vertices.',
      link: 'http://en.wikipedia.org/wiki/Distance_(graph_theory)',
      math: '\\(d\\)',
      calc: function(vs, es) {
        if (!graphProps.connected.value)
          return Number.NEGATIVE_INFINITY;
        var n = vs.length;
        var dists = Matrix.byFunc(n, n, function(i, j) {
          if (i === j)
            return 0;
          if (graphProps.adjMatr.value.get(i, j) === 1)
            return 1;
          return Number.POSITIVE_INFINITY;
        });
        for (var k = 0; k < n; k++)
          for (var i = 0; i < n; i++)
            for (var j = 0; j < n; j++) {
              var newDist = dists.get(i, k) + dists.get(k, j);
              if (dists.get(i, j) > newDist)
                dists.set(i, j, newDist);
            }
        var diam = Number.NEGATIVE_INFINITY;
        for (var i = 0; i < n; i++)
          for (var j = 0; j < n; j++)
            if (dists.get(i, j) > diam)
              diam = dists.get(i, j);
        return diam;
      }
    },
    numTris: {
      desc: 'number of triangles',
      longDesc: 'The number of disjoint cycles of length 3.',
      calc: function(vs, es) {
        return graphProps.adjMatr.value.pow(3).trace() / 6;
      }
    },
    bipartite: {
      desc: 'is bipartite',
      longDesc: 'Whether the vertex set can be partitioned into two disjoint subsets such that all edges are only between different subsets or not.',
      link: 'http://en.wikipedia.org/wiki/Bipartite_graph',
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
//    spanningTrees: {
//      desc: 'number of spanning trees',
//      longDesc: 'The number of trees that contain all vertices of the graph and whose only edges are edges of the graph. A tree is a graph without any cycles.',
//      link: 'http://en.wikipedia.org/wiki/Spanning_tree',
//      calc: function(vs, es) {
//        if (!graphProps.connected.value)
//          return 0;
//        if (vs.length === 1)
//          return 1;
//        if (graphProps.complete.value)
//          return Math.pow(vs.length, vs.length - 2);
//        if (graphProps.bipartite.value.bipartite) {
//          var m = graphProps.bipartite.value.m;
//          var n = graphProps.bipartite.value.n;
//          if (es.length === m * n)
//            return Math.pow(m, n - 1) * Math.pow(n, m - 1);
//        }
//        // TODO: broken for trees of order > 7?
//        // TODO: not always exact, have to round, even then sometimes off (by a lot), maybe could do with fractions?
//        return Math.floor(graphProps.laplMatr.value.cofactor(0, 0));
//      }
//    },
    tree: {
      desc: 'is a tree',
      longDesc: 'Whether the graph is connected and has no cycles or not.',
      link: 'http://en.wikipedia.org/wiki/Tree_(graph_theory)',
      calc: function(vs, es) {
        return graphProps.connected.value && graphProps.size.value === graphProps.order.value - 1;
      }
    },
    semiEulerian: {
      desc: 'is semi-Eulerian',
      longDesc: 'Whether the graph contains a walk which visits every edge exactly once or not.',
      link: 'http://en.wikipedia.org/wiki/Eulerian_path',
      calc: function(vs, es) {
        if (!graphProps.connected.value)
          return false;
        var oddCount = 0;
        return _.every(graphProps.degSeq.value, function(deg) {
          if (deg.deg % 2 === 1)
            oddCount++;
          if (oddCount > 2)
            return false;
          return true;
        }) && (oddCount === 2 || oddCount === 0);
      }
    },
    eulerian: {
      desc: 'is Eulerian',
      longDesc: 'Whether the graph contains a walk which visits every edge exactly once and returns to its origin or not.',
      link: 'http://en.wikipedia.org/wiki/Eulerian_path',
      calc: function(vs, es) {
        return graphProps.connected.value && ((graphProps.regular.value && graphProps.minDeg.value % 2 === 0) || _.every(graphProps.degSeq.value, function(deg) {
          return deg.deg % 2 === 0;
        }));
      }
    }
    // TODO: num cycles, vertex and edge connectivity, non-separable, girth, Hamiltonian, chromatic number and index, circuit rank,clique number, characteristic polynomial, chromatic polynomial, is perfect, is planar
  };

  // TODO: show properties of selected vertex (also is bridge vertex), maybe not worth showing
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
      var r = 100;
      if (n === 1)
        return function(v) {
          return { x: 0, y: 0 };
        };
      return function(v) {
        var t = v * 2 * Math.PI / n - Math.PI / 2 + (n % 2 === 0 ? Math.PI / n : 0);
        return { x: r * Math.cos(t), y: r * Math.sin(t) };
      };
    },
    multiRadial: function(n, counts, shifts, dists) {
      var scale = d3.scale.linear().domain([0, 2]).range([0, 100]);
      var coords = _.reduce(_.map(counts, function(count, i) {
        var coords = [];
        for (var j = 0; j < count; j++)
          coords.push({ r: dists ? dists[i] * 100 : scale(i), t: j * 2 * Math.PI / count - Math.PI / 2 + (count % 2 === 0 ? Math.PI / count : 0) + (shifts ? shifts[i] : 0) });
        return coords;
      }), function(all, coords) {
        return all.concat(coords);
      });
      if (coords.length !== n)
        throw new Error('counts do not add up to number of vertices');
      return function(v) {
        return { x: coords[v].r * Math.cos(coords[v].t), y: coords[v].r * Math.sin(coords[v].t) };
      };
    },
    horizLines: function(n, topCount) {
      var sep = 100;
      return function(v) {
        if (v < topCount)
          return { x: topCount === 1 ? 0 : v / (topCount - 1) * 2 * sep - sep, y: -sep };
        else
          return { x: n - topCount === 1 ? 0 : (v - topCount) / (n - topCount - 1) * 2 * sep - sep, y: sep };
      };
    }
  };

  var commonGraphs = {
    complete: {
      desc: 'complete graph',
      longDesc: 'The graph with n vertices that contains all possible edges.',
      link: 'http://en.wikipedia.org/wiki/Complete_graph',
      math: '\\(K_n\\)',
      args: {
        n: {
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
    empty: {
      desc: 'empty graph',
      longDesc: 'The graph with n vertices that contains no edges.',
      link: 'http://en.wikipedia.org/wiki/Null_graph#Edgeless_graph',
      math: '\\(\\overline K_n\\)',
      args: {
        n: {
          min: 1,
          max: Number.POSITIVE_INFINITY,
          def: 5
        }
      },
      verts: '\\(n\\)',
      edges: '\\(0\\)',
      make: function(n) {
        return { vs: _.range(n), es: [], layout: graphLayouts.radial(n) };
      }
    },
    cycle: {
      desc: 'cycle graph',
      longDesc: 'The graph with n vertices that consists of a single cycle.',
      link: 'http://en.wikipedia.org/wiki/Cycle_graph',
      math: '\\(C_n\\)',
      args: {
        n: {
          min: 1,
          def: 5
        }
      },
      verts: '\\(n\\)',
      edges: '\\(n\\)',
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        if (n === 2)
          es.push([0, 1]);
        if (n > 2)
          for (var i = 0; i < n; i++) {
            es.push([i, i + 1 >= n ? 0 : i + 1]);
          }
        return { vs: vs, es: es, layout: graphLayouts.radial(n) };
      }
    },
    path: {
      desc: 'path graph',
      longDesc: 'The graph with n vertices that consists of a single path.',
      link: 'http://en.wikipedia.org/wiki/Path_graph',
      math: '\\(P_n\\)',
      args: {
        n: {
          min: 1,
          def: 5
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
      longDesc: 'The bipartite graph with m vertices in one vertex partition and n in the other that has all possible edges between the two partitions.',
      link: 'http://en.wikipedia.org/wiki/Complete_bipartite_graph',
      math: '\\(K_{m,n}\\)',
      args: {
        m: {
          min: 1,
          def: 4
        },
        n: {
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
      longDesc: 'The graph that is a skeleton of an n-dimensional cube.',
      link: 'http://en.wikipedia.org/wiki/Hypercube_graph',
      math: '\\(Q_n\\)',
      args: {
        n: {
          min: 0,
          def: 3
        }
      },
      verts: '\\(2^n\\)',
      edges: '\\(2^{n-1}n\\)',
      make: function(n) {
        var nn = n;
        n = 1 << n;
        var vs = _.range(n);
        var es = [];
        for (var i = 0; i < n; i++) {
          for (var j = i + 1; j < n; j++) {
            var dist = 0;
            var val = i ^ j; // get bit difference between the numbers
            while (val !== 0 && dist < 2) { // count number of 1s
              dist++;
              val &= val - 1;
            }
            if (dist === 1)
              es.push([i, j]);
          }
        }
        var layout = null;
        if (nn <= 1) {
          layout = graphLayouts.radial(n);
        }
        else if (nn === 2) {
          vs = [0, 1, 3, 2];
          layout = graphLayouts.radial(n);
        }
        else if (nn === 3) {
          vs = [0, 1, 3, 2, 4, 5, 7, 6];
          layout = graphLayouts.multiRadial(n, [0, 4, 4], null, [0, Math.sqrt(2) / 4, 1]);
        }
        else if (nn === 4) {
          vs = [0, 13, 6, 1, 15, 2, 9, 14, 12, 4, 5, 7, 3, 11, 10, 8];
          layout = graphLayouts.multiRadial(n, [0, 8, 8], null, [0, 0.5, 0.5 * Math.sin(3 * Math.PI / 8) / Math.sin(Math.PI / 8)]);
        }
        else {
//          var counts = [0];
//          var shifts = [0];
//          for (var i = 0; i < n / 4; i++) {
//            counts.push(4);
//            shifts.push(Math.PI / 4 * i);
//          }
//          layout = graphLayouts.multiRadial(n, counts, shifts);
          layout = graphLayouts.radial(n);
        }
        return { vs: vs, es: es, layout: layout };
      }
    },
    star: {
      desc: 'star graph',
      longDesc: 'The graph with n vertices that has n-1 vertices connected to a single center vertex.',
      link: 'http://en.wikipedia.org/wiki/Star_(graph_theory)',
      math: '\\(S_n\\)',
      args: {
        n: {
          min: 1,
          def: 6
        }
      },
      verts: '\\(n\\)',
      edges: '\\(n-1\\)',
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        for (var i = 1; i < n; i++) {
          es.push([0, i]);
        }
        return { vs: vs, es: es, layout: graphLayouts.multiRadial(n, [1, 0, n - 1]) };
      }
    },
    wheel: {
      desc: 'wheel graph',
      longDesc: 'The graph with n vertices that has n-1 vertices connected in a cycle and to a single center vertex.',
      link: 'http://en.wikipedia.org/wiki/Wheel_graph',
      math: '\\(W_n\\)',
      args: {
        n: {
          min: 1,
          def: 6
        }
      },
      verts: '\\(n\\)',
      edges: '\\(2\\left(n-1\\right)\\)',
      make: function(n) {
        var vs = _.range(n);
        var es = [];
        for (var i = 1; i < n; i++) {
          es.push([0, i]);
        }
        if (i > 2)
          for (var i = 1; i < n; i++) {
            es.push([i, i + 1 >= n ? 1 : i + 1]);
          }
        return { vs: vs, es: es, layout: n === 3 ? graphLayouts.radial(n) : graphLayouts.multiRadial(n, [1, 0, n - 1]) };
      }
    },
    friendship: {
      desc: 'friendship graph',
      longDesc: 'The graph with n 3-cycles sharing only a single common vertex.',
      link: 'http://en.wikipedia.org/wiki/Friendship_graph',
      math: '\\(F_n\\)',
      args: {
        n: {
          min: 0,
          def: 3
        }
      },
      verts: '\\(2n+1\\)',
      edges: '\\(2n\\)',
      make: function(n) {
        var vs = [0];
        var es = [];
        for (var i = 0; i < n; i++) {
          var a, b;
          vs.push(a = vs.length);
          vs.push(b = vs.length);
          es.push([0, a]);
          es.push([0, b]);
          es.push([a, b]);
        }
        return { vs: vs, es: es, layout: n === 1 ? graphLayouts.radial(2 * n + 1) : graphLayouts.multiRadial(2 * n + 1, [1, 0, 2 * n + 1 - 1]) };
      }
    }
  };

  var staticGraphs = {
    peterson: {
      desc: 'Peterson graph',
      longDesc: 'The Petersen graph.',
      link: 'http://en.wikipedia.org/wiki/Petersen_graph',
      verts: 10,
      edges: [
        [0, 2],
        [0, 3],
        [0, 5],
        [1, 3],
        [1, 4],
        [1, 6],
        [2, 4],
        [2, 7],
        [3, 8],
        [4, 9],
        [5, 6],
        [6, 7],
        [7, 8],
        [8, 9],
        [9, 5]
      ],
      layout: {
        name: 'multiRadial',
        args: [10, [0, 5, 5]]
      }
    },
    durer: {
      desc: 'Dürer graph',
      longDesc: 'The Dürer graph',
      link: 'http://en.wikipedia.org/wiki/D%C3%BCrer_graph',
      verts: 12,
      edges: [
        [0, 2],
        [0, 4],
        [0, 6],
        [1, 3],
        [1, 5],
        [1, 7],
        [2, 4],
        [2, 8],
        [3, 5],
        [3, 9],
        [4, 10],
        [5, 11],
        [6, 7],
        [7, 8],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 6]
      ],
      layout: {
        name: 'multiRadial',
        args: [12, [0, 6, 6]]
      }
    },
    kantor: {
      desc: 'Möbius-Kantor graph',
      longDesc: 'The Möbius-Kantor graph',
      link: 'http://en.wikipedia.org/wiki/M%C3%B6bius%E2%80%93Kantor_graph',
      verts: 16,
      edges: [
        [0, 3],
        [0, 5],
        [0, 15],
        [1, 4],
        [1, 6],
        [1, 8],
        [2, 5],
        [2, 7],
        [2, 9],
        [3, 6],
        [3, 10],
        [4, 7],
        [4, 11],
        [5, 12],
        [6, 13],
        [7, 14],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 12],
        [12, 13],
        [13, 14],
        [14, 15],
        [15, 8]
      ],
      layout: {
        name: 'multiRadial',
        args: [16, [0, 8, 8], null, [0, 0.5, 0.5 * Math.sin(3 * Math.PI / 8) / Math.sin(Math.PI / 8)]]
      }
    },
    dodecahedron: {
      desc: 'dodecahedron Schlegel diagram',
      longDesc: 'The graph representing the skeleton of the dodecahedron.',
      link: 'http://en.wikipedia.org/wiki/Dodecahedron#Two-dimensional_symmetry_projections',
      verts: 20,
      edges: [
        [0, 1],
        [0, 4],
        [0, 5],
        [1, 2],
        [1, 6],
        [2, 3],
        [2, 7],
        [3, 4],
        [3, 8],
        [4, 9],
        [5, 10],
        [5, 14],
        [6, 10],
        [6, 11],
        [7, 11],
        [7, 12],
        [8, 12],
        [8, 13],
        [9, 13],
        [9, 14],
        [10, 15],
        [11, 16],
        [12, 17],
        [13, 18],
        [14, 19],
        [15, 16],
        [16, 17],
        [17, 18],
        [18, 19],
        [19, 15]
      ],
      layout: {
        name: 'multiRadial',
        args: [20, [0, 5, 5, 5, 5], [0, -Math.PI / 5, -Math.PI / 5, 0, 0]]
      }
    }
  };
  // add static graphs to common graphs
  _.each(staticGraphs, function(g, name) {
    commonGraphs[name] = {
      desc: g.desc,
      longDesc: g.longDesc,
      link: g.link,
      math: g.math,
      args: { },
      verts: '\\(' + g.verts + '\\)',
      edges: '\\(' + g.edges.length + '\\)',
      make: function() {
        return { vs: _.range(g.verts), es: g.edges, layout: graphLayouts[g.layout.name].apply(null, g.layout.args) };
      }
    };
  });

  var width, height;
  // TODO: move things when canvas resizes
  // TODO: consider canvas height more carefully, maybe should fill view vertically (use $(window).resize(handler) and window.innerHeight)

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
    console.log('select');
    node.classed('node-selected', true)
            .transition()
            .duration(300)
            .attr('r', 13);
  }

  // sets the given node as deselected
  function deselectNode() {
    console.log('deselect');
    svg.selectAll('.node-selected')
            .classed('node-selected', false)
            .transition()
            .duration(100)
            .attr('r', 8);
  }

  function hasSelection() {
    return !svg.selectAll('.node-selected').empty();
  }

  function addNode(node) {
    console.log('add node');
    nodes.push(node);
    linkNodes(node);
  }

  function linkNodes(node) {
    var selected = d3.selectAll('.node-selected');
    if (!selected.empty()) {
      var a = node, b = selected.datum();
      deselectNode();
      if (a !== b) {
        if (!isAdj(a, b)) {
          console.log('link');
          links.push({ source: b, target: a, edgeNum: links.length });
        }
        else {
          console.log('unlink');
          links = _.filter(links, function(link) {
            return !((link.source === a && link.target === b) || (link.source === b && link.target === a));
          });
          links.forEach(function(link, i) {
            link.edgeNum = i;
          });
        }
        return true;
      }
    }
    return false;
  }

  function removeNode(node) {
    console.log('remove');
    var selected = d3.selectAll('.node-selected');
    var a = node, b = selected.empty() ? null : selected.datum();
    if (a === b)
      deselectNode();
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
  }

  // clears all data
  function clear() {
    console.log('clear');
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
            if (svg.selectAll('.node-selected').empty()) {
              var m = d3.mouse(this);
              addNode({ x: m[0], y: m[1], vertNum: nodes.length });
              update();
            } else {
              deselectNode();
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
              var a = d3.select(this).datum();
              if (d3.event.ctrlKey) {
                removeNode(a);
                update();
              }
              else {
                var hadSelection = hasSelection();
                if (linkNodes(a))
                  update();
                if (!hadSelection || d3.event.shiftKey)
                  selectNode(d3.select(this));
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
            },
            {
              img: 'img/complement.svg',
              onclick: function() {
                var old_links = links.slice(0);
                links.splice(0, links.length);
                var e = 0;
                for (var i = 0; i < nodes.length; i++)
                  for (var j = i + 1; j < nodes.length; j++) {
                    var n1 = nodes[i];
                    var n2 = nodes[j];
                    if (!_.any(old_links, function(e) {
                      return (e.source === n1 && e.target === n2) || (e.source === n2 && e.target === n1);
                    })) {
                      links.push({ target: n1, source: n2, edgeNum: e++ });
                    }
                  }
                update();
              },
              title: 'Graph complement'
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

  var graphPropRow = d3.select('#graph-prop-table tbody').selectAll('tr')
          .data(_.map(graphProps, function(prop, name) {
            prop.name = name;
            return prop;
          }))
          .enter()
          .append('tr');
  graphPropRow.append('td')
          .attr('class', 'desc')
          .text(function(prop) {
            return prop.desc;
          })
          .filter(function(prop) {
            return prop.longDesc;
          })
          .append('a')
          .attr('tabindex', 0)
          .attr('role', 'button')
          .attr('class', 'help-popover btn btn-xs btn-default')
          .attr('data-toggle', 'popover')
          .attr('data-trigger', 'focus')
          .attr('data-placement', 'bottom')
          .attr('data-html', true)
          .attr('title', function(prop) {
            return prop.desc;
          })
          .attr('data-content', function(prop) {
            var s = prop.longDesc;
            if (prop.link)
              s += ' <a href=' + prop.link + ' target="_blank">[more]</a>';
            return s;
          })
          .text('?');
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
          })
          .filter(function(prop) {
            return prop.longDesc;
          })
          .append('a')
          .attr('tabindex', 0)
          .attr('role', 'button')
          .attr('class', 'help-popover btn btn-xs btn-default')
          .attr('data-toggle', 'popover')
          .attr('data-trigger', 'focus')
          .attr('data-placement', 'bottom')
          .attr('data-html', true)
          .attr('title', function(prop) {
            return prop.desc;
          })
          .attr('data-content', function(prop) {
            var s = prop.longDesc;
            if (prop.link)
              s += ' <a href=' + prop.link + ' target="_blank">[more]</a>';
            return s;
          })
          .text('?');
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
            var formGroup = form.selectAll('.form-group')
                    .data(_.map(graph.args, function(arg, name) {
                      arg.name = name;
                      return arg;
                    }))
                    .enter()
                    .append('div')
                    .attr('class', 'form-group');
            formGroup.append('label')
                    .text(function(arg) {
                      return '\\(' + arg.name + '\\):';
                    });
            formGroup.append('input')
                    .attr('class', 'form-control')
                    .attr('type', 'number')
                    .attr('min', function(arg) {
                      return arg.min;
                    })
                    .attr('max', function(arg) {
                      return arg.max;
                    })
                    .attr('defaultValue', function(arg) {
                      return arg.def;
                    })
                    .attr('value', function(arg) {
                      return arg.def;
                    })
                    .attr('required', true);
            form.append('button')
                    .attr('type', 'submit')
                    .attr('class', 'btn btn-default')
                    .text('Create');
            form.on('submit', function() {
              console.log('make graph: ' + graph.name);
              var args = [];
              formGroup.each(function(arg) {
                var val = +d3.select(this).select('input').property('value');
                console.log(arg.name + ': ' + val);
                args.push(val);
              });
              var g = graph.make.apply(null, args);
              console.log(g.vs.length + ' verts: ' + JSON.stringify(g.vs));
              console.log(g.es.length + ' edges: ' + JSON.stringify(g.es));
              clear();
              g.vs.forEach(function(v, i) {
                var pos = g.layout(i);
                console.log(v + ' ' + JSON.stringify(pos));
                nodes[v] = { x: pos.x + width / 2, y: pos.y + height / 2, vertNum: v };
              });
              g.es.forEach(function(e) {
                links.push({ source: nodes[e[0]], target: nodes[e[1]], edgeNum: e });
              });
              update();
              {
                var rect = $('#canvas')[0].getBoundingClientRect();
                console.log(JSON.stringify(rect));
                if (rect.top <= -width / 3)
                  $('body').scrollTo('#canvas', 500);
              }
              d3.event.preventDefault();
            });
          });

  update();

  // activate all popovers
  $('[data-toggle="popover"]').popover();

//  $('#file-load-input').change(function() {
//    var file = $(this)[0].files[0];
//    console.log(JSON.stringify(file));
//    if (file.type !== 'application/json')
//      return false;
//  });
});

$(document).on('click', 'a', function(event) {
  if ($(this).attr('href') === '#') {
    event.preventDefault();
    console.log('Default "#" link action prevented.');
  }
});

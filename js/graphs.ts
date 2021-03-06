/**
 * Created by djbsn_000 on 5/2/2015.
 */

///<reference path="libs/d3/d3.d.ts"/>
///<reference path="libs/jquery/jquery.d.ts"/>
///<reference path="libs/jquery-color/jquery.color.d.ts"/>
///<reference path="libs/jquery-scrollTo/jquery.scrollTo.d.ts"/>
///<reference path="libs/twitter-bootstrap/js/bootstrap.d.ts"/>
///<reference path="libs/underscore.js/underscore.d.ts"/>
///<reference path="libs/mathjax.d.ts"/>
///<reference path="Matrix.ts"/>

interface Point {
    x:number;
    y:number;
    dragstart?:number[];
}
interface PolarPoint {
    r:number;
    t:number;
}

interface Vertex extends Point {
    vertNum:number;
    __parent?:Vertex;
    __dist?:number;
}

interface Edge {
    source:Vertex;
    target:Vertex;
    edgeNum:number;
}

interface GraphProp<T> {
    desc:string;
    longDesc?:string;
    link?:string;
    math?:string;
    defaultHidden?:boolean;
    hidden?:boolean;
    calc:(vs:Vertex[], es:Edge[]) => T;
    repr?:(value:T) => string;
    value?:T;
    name?:string;
}

interface DegreeVertex {
    v:Vertex;
    deg:number;
}

interface Bipartite {
    bipartite:boolean;
    m?:number;
    n?:number;
}

interface GraphLayout {
    (v:number):Point;
}

interface GraphTemplate {
    vs:number[];
    es:number[][];
    layout:GraphLayout;
}

interface CommonGraphArg {
    min?:number;
    max?:number;
    def?:number;
    name?:string;
}

interface CommonGraph {
    desc:string;
    longDesc:string;
    link:string;
    math?:string;
    args:{[index:string]:CommonGraphArg};
    verts:string;
    edges:string;
    make:(...args) => GraphTemplate;
    name?:string;
}

interface StaticGraph {
    desc:string;
    longDesc:string;
    link:string;
    math?:string;
    verts:number;
    edges:number[][];
    layout: {
        name:string;
        args:any[];
    }
}

interface ToolbarButton {
    img:string;
    onclick:(button?:ToolbarButton, index?:number) => void;
    title:string;
    toggleable?:boolean;
    toggleState?:boolean;
    name?:string;
}

function updateTex() {
    MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
}

// TODO: option to show labels on edges
// TODO: should be able to rotate/scale view
$(function() {
    var graphProps:{[index:string]:GraphProp<any>} = {
        order: <GraphProp<number>>{
            desc: 'order',
            longDesc: 'The number of vertices of the graph.',
            math: '\\(\\left|V\\right|\\)',
            calc: (vs:Vertex[], es:Edge[]) => vs.length
        },
        size: <GraphProp<number>>{
            desc: 'size',
            longDesc: 'The number of edges of the graph.',
            math: '\\(\\left|E\\right|\\)',
            calc: (vs:Vertex[], es:Edge[]) => es.length
        },
        adjMatr: <GraphProp<Matrix>>{
            desc: 'adjacency matrix',
            longDesc: 'The matrix where the i,j<sup>th</sup> entry is a 1 if the i<sup>th</sup> and j<sup>th</sup> vertices are adjacent, or connected by an edge, and otherwise 0.',
            link: 'http://en.wikipedia.org/wiki/Adjacency_matrix',
            math: '\\(A\\)',
            defaultHidden: true,
            calc: (vs:Vertex[], es:Edge[]) => Matrix.byFunc(vs.length, vs.length,
                (i, j) => +(_.any(es,
                    (e:Edge) => (e.source === vs[i] && e.target === vs[j]) || (e.source === vs[j] && e.target === vs[i]))))
        },
        degSeq: <GraphProp<DegreeVertex[]>>{
            desc: 'degree sequence',
            longDesc: 'The sequence of degrees of the vertices of the graph. The degree of a vertex is the number of edges it is connected to.',
            link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)#Degree_sequence',
            defaultHidden: true,
            calc: (vs:Vertex[], es:Edge[]) => _.sortBy(_.map((<Matrix>graphProps['adjMatr'].value).values, (row, i) =><DegreeVertex>{
                v: vs[i],
                deg: _.reduce<number, number>(row, (sum, adj) => sum + adj, 0)
            }), (v:DegreeVertex) => v.deg),
            repr: (value:DegreeVertex[]) => '(' + _.map(value, (v:DegreeVertex) => v.deg).join(', ') + ')'
        },
        degMatr: <GraphProp<Matrix>>{
            desc: 'degree matrix',
            longDesc: 'The matrix where the i<sup>th</sup> diagonal entry is the degree of the i<sup>th</sup> vertex and all other entries are 0.',
            link: 'http://en.wikipedia.org/wiki/Degree_matrix',
            math: '\\(D\\)',
            defaultHidden: true,
            calc: (vs:Vertex[], es:Edge[]) => Matrix.byFunc(vs.length, vs.length,
                (i, j) =>  i !== j ? 0 : _.find(<DegreeVertex[]>graphProps['degSeq'].value,
                    (v:DegreeVertex) => v.v === vs[i]).deg)
        },
        laplMatr: <GraphProp<Matrix>>{
            desc: 'Laplacian matrix',
            longDesc: 'The degree matrix of the graph minus the adjaceny matrix of the graph.',
            link: 'http://en.wikipedia.org/wiki/Laplacian_matrix',
            math: '\\(L\\)',
            defaultHidden: true,
            calc: (vs:Vertex[], es:Edge[]) =>
                (<Matrix>graphProps['degMatr'].value).sub(<Matrix>graphProps['adjMatr'].value)
        },
        inciMatr: <GraphProp<Matrix>>{
            desc: 'incidence matrix',
            longDesc: 'The matrix where the i,j<sup>th</sup> entry is a 1 if the i<sup>th</sup> vertex and j<sup>th</sup> edge are indicent, or the vertex is part of the edge, and otherwise 0.',
            link: 'http://en.wikipedia.org/wiki/Incidence_matrix',
            math: '\\(M\\)',
            defaultHidden: true,
            calc: (vs:Vertex[], es:Edge[]) => es.length === 0 ? Matrix.zero(vs.length, 1) :
                Matrix.byFunc(vs.length, es.length,
                    (v, e) => +(v === es[e].source.vertNum || v === es[e].target.vertNum))
        },
        minDeg: <GraphProp<number>>{
            desc: 'min degree',
            longDesc: 'The smallest degree of any vertex.',
            link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)',
            math: '\\(\\delta\\)',
            calc: (vs:Vertex[], es:Edge[]) => vs.length === 0 ? Number.POSITIVE_INFINITY :
                (<DegreeVertex[]>graphProps['degSeq'].value)[0].deg
        },
        maxDeg: <GraphProp<number>>{
            desc: 'max degree',
            longDesc: 'The largest degree of any vertex.',
            link: 'http://en.wikipedia.org/wiki/Degree_(graph_theory)',
            math: '\\(\\Delta\\)',
            calc: (vs:Vertex[], es:Edge[]) => vs.length === 0 ? Number.NEGATIVE_INFINITY :
                (<DegreeVertex[]>graphProps['degSeq'].value)[vs.length - 1].deg
        },
        regular: <GraphProp<boolean>>{
            desc: 'is regular',
            longDesc: 'Whether all vertices have the same degree or not.',
            link: 'http://en.wikipedia.org/wiki/Regular_graph',
            calc: (vs:Vertex[], es:Edge[]) => <number>graphProps['minDeg'].value === <number>graphProps['maxDeg'].value
        },
        complete: <GraphProp<boolean>>{
            desc: 'is complete',
            longDesc: 'Whether the graph has all possible edges or not.',
            link: 'http://en.wikipedia.org/wiki/Complete_graph',
            calc: (vs:Vertex[], es:Edge[]) => es.length === vs.length * (vs.length - 1) / 2
        },
        components: <GraphProp<number>>{
            desc: 'number of components',
            longDesc: 'The number of connected subgraphs.',
            link: 'http://en.wikipedia.org/wiki/Connected_component_(graph_theory)',
            calc: (vs:Vertex[], es:Edge[]) => {
                var count = 0;
                var seen:boolean[] = []; // keep track of if we've looked at this vertex already
                for (var i = 0; i < vs.length; i++) {
                    seen.push(false);
                }
                var queue:number[] = []; // queue of adjacent vertices to look at
                var vs2 = _.range(vs.length);
                vs2.forEach((v) => {
                    if (!seen[v]) { // if haven't already seen this vertex
                        count++; // increment component count
                        seen[v] = true; // indicate we've seen this vertex
                        queue.push(v);
                        while (queue.length !== 0) { // go through all adjacent vertices
                            var w = queue.shift();
                            var row = (<Matrix>graphProps['adjMatr'].value).values[w]; // get adjacency row for w
                            vs2.forEach((u) => {
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
        connected: <GraphProp<boolean>>{
            desc: 'is connected',
            longDesc: 'Whether the number of components is one or not.',
            link: 'http://en.wikipedia.org/wiki/Connectivity_(graph_theory)',
            calc: (vs:Vertex[], es:Edge[]) => <number>graphProps['components'].value === 1
        },
        distMatr: <GraphProp<Matrix>>{
            desc: 'distance matrix',
            hidden: true,
            calc: (vs:Vertex[], es:Edge[]) => {
                var n = vs.length;
                var dists = Matrix.byFunc(n, n, (i, j) => {
                    if (i === j)
                        return 0;
                    if ((<Matrix>graphProps['adjMatr'].value).get(i, j) === 1)
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
                return dists;
            }
        },
        diameter: <GraphProp<number>>{
            desc: 'diameter',
            longDesc: 'The largest eccentricity of all vertices. The eccentricity of a vertex is the largest distance between it and any other vertex. The distance between two vertices is the length of the shortest path between them.',
            link: 'http://en.wikipedia.org/wiki/Distance_(graph_theory)',
            math: '\\(d\\)',
            calc: (vs:Vertex[], es:Edge[]) => {
                if (!<boolean>graphProps['connected'].value)
                    return 0;
                var diam = 0;
                for (var i = 0; i < vs.length; i++)
                    for (var j = i; j < vs.length; j++)
                        diam = Math.max(diam, (<Matrix>graphProps['distMatr'].value).get(i, j));
                return diam;
            }
        },
        radius: <GraphProp<number>>{
            desc: 'radius',
            longDesc: 'The smallest eccentricity of all vertices. The eccentricity of a vertex is the largest distance between it and any other vertex. The distance between two vertices is the length of the shortest path between them.',
            link: 'http://en.wikipedia.org/wiki/Distance_(graph_theory)',
            math: '\\(r\\)',
            calc: (vs:Vertex[], es:Edge[]) => {
                if (!(<boolean>graphProps['connected'].value))
                    return 0;
                var rad = Number.POSITIVE_INFINITY;
                for (var i = 0; i < vs.length; i++) {
                    var ecc = 0;
                    for (var j = 0; j < vs.length; j++)
                        if (i !== j)
                            ecc = Math.max(ecc, (<Matrix>graphProps['distMatr'].value).get(i, j));
                    rad = Math.min(rad, ecc);
                }
                return rad;
            }
        },
        numTris: <GraphProp<number>>{
            desc: 'number of triangles',
            longDesc: 'The number of disjoint cycles of length 3.',
            calc: (vs:Vertex[], es:Edge[]) => (<Matrix>graphProps['adjMatr'].value).pow(3).trace / 6
        },
        bipartite: <GraphProp<Bipartite>>{
            desc: 'is bipartite',
            longDesc: 'Whether the vertex set can be partitioned into two disjoint subsets such that all edges are only between different subsets or not.',
            link: 'http://en.wikipedia.org/wiki/Bipartite_graph',
            calc: (vs:Vertex[], es:Edge[]) => {
                if (!(<boolean>graphProps['connected'].value) || (<number>graphProps['numTris'].value) !== 0)
                // not connected or has triangles, not bipartite
                    return { bipartite: false };
                var color:number[] = []; // color array, 0 means uncolored, 1 and -1 are colors
                for (var i = 0; i < vs.length; i++) {
                    color.push(0);
                }
                color[0] = 1; // color first vertex 1
                var queue:number[] = [0];
                while (queue.length !== 0) {
                    var v = queue.shift();
                    var row = (<Matrix>graphProps['adjMatr'].value).getRow(v); // get adjacency row for v
                    for (var u = 0; u < vs.length; u++) {
                        if (row[u] === 1) { // if u adj to v
                            if (color[u] === 0) { // if u uncolored
                                color[u] = -color[v]; // color opposite to v
                                queue.push(u);
                            }
                            else if (color[u] === color[v]) {
                                // if u is colored the same as v, not 2-colorable, not bipartite
                                return { bipartite: false };
                            }
                        }
                    }
                }
                var counts = _.countBy(color, function(c) {
                    return c === 1 ? 'm' : 'n';
                });
                return {
                    bipartite: true,
                    m: counts['m'],
                    n: counts['n']
                };
            },
            repr: (value:Bipartite) => value.bipartite.toString()
        },
//    spanningTrees: {
//      desc: 'number of spanning trees',
//      longDesc: 'The number of trees that contain all vertices of the graph and whose only edges are edges of the graph. A tree is a graph without any cycles.',
//      link: 'http://en.wikipedia.org/wiki/Spanning_tree',
//      calc: function(vs, es) {
//        if (!(<>graphProps['connected'].value))
//          return 0;
//        if (vs.length === 1)
//          return 1;
//        if ((<>graphProps['complete'].value))
//          return Math.pow(vs.length, vs.length - 2);
//        if ((<>graphProps['bipartite'].value).bipartite) {
//          var m = (<>graphProps['bipartite'].value).m;
//          var n = (<>graphProps['bipartite'].value).n;
//          if (es.length === m * n)
//            return Math.pow(m, n - 1) * Math.pow(n, m - 1);
//        }
//        // TODO: broken for trees of order > 7?
//        // TODO: not always exact, have to round, even then sometimes off (by a lot), maybe could do with fractions?
//        return Math.floor((<>graphProps['laplMatr'].value).cofactor(0, 0));
//      }
//    },
        tree: <GraphProp<boolean>>{
            desc: 'is a tree',
            longDesc: 'Whether the graph is connected and has no cycles or not.',
            link: 'http://en.wikipedia.org/wiki/Tree_(graph_theory)',
            calc: (vs:Vertex[], es:Edge[]) => (<boolean>graphProps['connected'].value) && es.length === vs.length - 1
        },
        girth: <GraphProp<number>>{
            desc: 'girth',
            longDesc: 'The length of the shortest cycle.',
            link: 'http://en.wikipedia.org/wiki/Girth_(graph_theory)',
            math: '\\(g\\)',
            calc: (vs:Vertex[], es:Edge[]) => {
                if (!(<boolean>graphProps['connected'].value) || (<boolean>graphProps['tree'].value))
                    return Number.POSITIVE_INFINITY;
                if ((<number>graphProps['numTris'].value) > 0)
                    return 3;
                var g = Number.POSITIVE_INFINITY;
                // algorithm from http://webcourse.cs.technion.ac.il/234247/Winter2003-2004/ho/WCFiles/Girth.pdf

                for (var i = 0; i < vs.length; i++) {
                    var v:Vertex = vs[i];
                    var S:Vertex[] = [];
                    var R:Vertex[] = [v];
                    v.__parent = null;
                    v.__dist = 0;
                    while (R.length) {
                        var x = R.shift();
                        S.push(x);
                        _.filter(vs, (u:Vertex) => isAdj(u, x) && u !== x.__parent).forEach((y:Vertex) => {
                            if (S.indexOf(y) === -1) {
                                y.__parent = x;
                                y.__dist = x.__dist + 1;
                                R.push(y);
                            } else {
                                g = Math.min(g, x.__dist + y.__dist + 1);
                            }
                        });
                    }
                }
                vs.forEach((v) => {
                    delete v.__parent;
                    delete v.__dist;
                });

                return g;
            }
        },
        semiEulerian: <GraphProp<boolean>>{
            desc: 'is semi-Eulerian',
            longDesc: 'Whether the graph contains a walk which visits every edge exactly once or not.',
            link: 'http://en.wikipedia.org/wiki/Eulerian_path',
            calc: (vs:Vertex[], es:Edge[]) => {
                if (!(<boolean>graphProps['connected'].value))
                    return false;
                var oddCount = 0;
                return _.every((<DegreeVertex[]>graphProps['degSeq'].value), (deg:DegreeVertex) => {
                        if (deg.deg % 2 === 1)
                            oddCount++;
                        return oddCount <= 2;

                    }) && (oddCount === 2 || oddCount === 0);
            }
        },
        eulerian: <GraphProp<boolean>>{
            desc: 'is Eulerian',
            longDesc: 'Whether the graph contains a walk which visits every edge exactly once and returns to its origin or not.',
            link: 'http://en.wikipedia.org/wiki/Eulerian_path',
            calc: (vs:Vertex[], es:Edge[]) => (<boolean>graphProps['connected'].value) &&
            (((<boolean>graphProps['regular'].value) && (<number>graphProps['minDeg'].value) % 2 === 0) ||
            _.every((<DegreeVertex[]>graphProps['degSeq'].value), (deg:DegreeVertex) => deg.deg % 2 === 0))
        }
        // TODO: num cycles, vertex and edge connectivity, non-separable, Hamiltonian, chromatic number and index, circuit rank,clique number, characteristic polynomial, chromatic polynomial, is perfect, is planar
    };

    var graphLayouts = {
        radial: (n:number):GraphLayout => {
            var r = 100;
            if (n === 1)
                return (v:number) => <Point>{ x: 0, y: 0 };
            return (v:number) => {
                var t = v * 2 * Math.PI / n - Math.PI / 2 + (n % 2 === 0 ? Math.PI / n : 0);
                return <Point>{ x: r * Math.cos(t), y: r * Math.sin(t) };
            }
        },
        multiRadial: (n:number, counts:number[], shifts?:number[], dists?:number[]):GraphLayout => {
            var scale = d3.scale.linear().domain([0, 2]).range([0, 100]);
            var coords:PolarPoint[] = _.reduce<PolarPoint[], PolarPoint[]>(
                _.map<number, PolarPoint[]>(counts, (count, i) => {
                    var coords:PolarPoint[] = [];
                    for (var j = 0; j < count; j++)
                        coords.push({
                            r: dists ? dists[i] * 100 : scale(i),
                            t: j * 2 * Math.PI / count - Math.PI / 2 +
                            (count % 2 === 0 ? Math.PI / count : 0) +
                            (shifts ? shifts[i] : 0)
                        });
                    return coords;
                }), (all, coords) => all.concat(coords));
            if (coords.length !== n)
                throw new Error('counts do not add up to number of vertices');
            return (v:number) => <Point>{
                x: coords[v].r * Math.cos(coords[v].t),
                y: coords[v].r * Math.sin(coords[v].t)
            };
        },
        horizLines: (n:number, topCount:number) => {
            var sep = 100;
            return (v:number) => {
                if (v < topCount)
                    return <Point>{
                        x: topCount === 1 ? 0 : v / (topCount - 1) * 2 * sep - sep,
                        y: -sep
                    };
                else
                    return <Point>{
                        x: n - topCount === 1 ? 0 : (v - topCount) / (n - topCount - 1) * 2 * sep - sep,
                        y: sep
                    };
            }
        }
    };

    var commonGraphs:{[index:string]:CommonGraph} = {
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
            make: (n:number) => {
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: graphLayouts.radial(n) };
                for (var i = 0; i < n; i++) {
                    for (var j = i + 1; j < n; j++) {
                        gt.es.push([i, j]);
                    }
                }
                return gt;
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
            make: (n:number) => <GraphTemplate>{ vs: _.range(n), es: [], layout: graphLayouts.radial(n) }
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
            make: (n:number) => {
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: graphLayouts.radial(n) };
                if (n === 2)
                    gt.es.push([0, 1]);
                if (n > 2)
                    for (var i = 0; i < n; i++) {
                        gt.es.push([i, i + 1 >= n ? 0 : i + 1]);
                    }
                return gt;
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
            make: (n:number) => {
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: graphLayouts.radial(n) };
                for (var i = 0; i < n - 1; i++) {
                    gt.es.push([i, i + 1]);
                }
                return gt;
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
            make: (m:number, n:number) => {
                var gt:GraphTemplate = { vs: _.range(n + m), es: [], layout: graphLayouts.horizLines(n + m, m) };
                for (var i = 0; i < m; i++) {
                    for (var j = 0; j < n; j++) {
                        gt.es.push([i, m + j]);
                    }
                }
                return gt;
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
            make: (n:number) => {
                var nn = n;
                n = 1 << n;
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: null };
                for (var i = 0; i < n; i++) {
                    for (var j = i + 1; j < n; j++) {
                        var dist = 0;
                        var val = i ^ j; // get bit difference between the numbers
                        while (val !== 0 && dist < 2) { // count number of 1s
                            dist++;
                            val &= val - 1;
                        }
                        if (dist === 1)
                            gt.es.push([i, j]);
                    }
                }
                if (nn <= 1) {
                    gt.layout = graphLayouts.radial(n);
                }
                else if (nn === 2) {
                    gt.vs = [0, 1, 3, 2];
                    gt.layout = graphLayouts.radial(n);
                }
                else if (nn === 3) {
                    gt.vs = [0, 1, 3, 2, 4, 5, 7, 6];
                    gt.layout = graphLayouts.multiRadial(n, [0, 4, 4], null, [0, Math.sqrt(2) / 4, 1]);
                }
                else if (nn === 4) {
                    gt.vs = [0, 13, 6, 1, 15, 2, 9, 14, 12, 4, 5, 7, 3, 11, 10, 8];
                    gt.layout = graphLayouts.multiRadial(n, [0, 8, 8], null, [0, 0.5, 0.5 * Math.sin(3 * Math.PI / 8) / Math.sin(Math.PI / 8)]);
                }
                else {
//          var counts = [0];
//          var shifts = [0];
//          for (var i = 0; i < n / 4; i++) {
//            counts.push(4);
//            shifts.push(Math.PI / 4 * i);
//          }
//          gt.layout = graphLayouts.multiRadial(n, counts, shifts);
                    gt.layout = graphLayouts.radial(n);
                }
                return gt;
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
            make: (n:number) => {
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: graphLayouts.multiRadial(n, [1, 0, n - 1]) };
                for (var i = 1; i < n; i++) {
                    gt.es.push([0, i]);
                }
                return gt;
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
            make: (n:number) => {
                var gt:GraphTemplate = {
                    vs: _.range(n), es: [],
                    layout: n === 3 ? graphLayouts.radial(n) : graphLayouts.multiRadial(n, [1, 0, n - 1])
                };
                for (var i = 1; i < n; i++) {
                    gt.es.push([0, i]);
                }
                if (i > 2)
                    for (var i = 1; i < n; i++) {
                        gt.es.push([i, i + 1 >= n ? 1 : i + 1]);
                    }
                return gt;
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
            edges: '\\(3n\\)',
            make: (n:number) => {
                var gt:GraphTemplate = {
                    vs: [0],
                    es: [],
                    layout: n === 1 ? graphLayouts.radial(2 * n + 1) :
                        graphLayouts.multiRadial(2 * n + 1, [1, 0, 2 * n + 1 - 1])
                };
                for (var i = 0; i < n; i++) {
                    var a, b;
                    gt.vs.push(a = gt.vs.length);
                    gt.vs.push(b = gt.vs.length);
                    gt.es.push([0, a]);
                    gt.es.push([0, b]);
                    gt.es.push([a, b]);
                }
                return gt;
            }
        },
        book: {
            desc: 'book graph',
            longDesc: 'The graph with n 3-cycles sharing a common edge.',
            link: 'http://en.wikipedia.org/wiki/Book_%28graph_theory%29',
            math: '\\(B_n\\)',
            args: {
                n: {
                    min: 0,
                    def: 3
                }
            },
            verts: '\\(n+2\\)',
            edges: '\\(2n+1\\)',
            make: (n:number) => {
                var gt:GraphTemplate = { vs: [0, 1], es: [[0, 1]], layout: graphLayouts.horizLines(n + 2, 2) };
                for (var i = 0; i < n; i++) {
                    var a;
                    gt.vs.push(a = gt.vs.length);
                    gt.es.push([0, a]);
                    gt.es.push([1, a]);
                }
                return gt;
            }
        },
        johnson: {
            desc: 'Johnson graph',
            longDesc: 'The graph representing k-element subsets of n elements whch are related if they share exactly k-1 elements.',
            link: 'http://en.wikipedia.org/wiki/Johnson_graph',
            math: '\\(J\\!\\left(n,k\\right)\\)',
            args: {
                n: {
                    min: 1,
                    def: 5
                },
                k: {
                    min: 0,
                    def: 2
                }
            },
            verts: '\\(\\binom n k\\)',
            edges: '\\(\\frac{k\\left(n-k\\right)}2\\binom n k\\)',
            make: (n:number, k:number) => {
                if (k > n)
                    return <GraphTemplate>{ vs: [], es: [], layout: graphLayouts.radial(0) };
                var elements = _.range(k);
                var sets:number[][] = [];
                var done = false;
                while (!done) {
                    sets.push(elements.slice(0));
                    for (var i = k - 1; i >= 0; i--) {
                        if (i < k - 1)
                            elements[i + 1] = elements[i] + 2;
                        for (var j = i + 1; j + 1 < k; j++)
                            elements[j + 1] = elements[j] + 1;
                        if (++elements[i] > n - (k - i - 1) - 1) {
                            if (i === 0) {
                                done = true;
                                break;
                            }
                        }
                        else
                            break;
                    }
                }
                var gt:GraphTemplate = { vs: [], es: [], layout: null };
                sets.forEach(function(s, i) {
                    gt.vs.push(i);
                    for (var j = sets.length - 1; j > i; j--)
                        if (_.intersection(s, sets[j]).length === k - 1)
                            gt.es.push([i, j]);
                });
                gt.layout = graphLayouts.radial(gt.vs.length);
                return gt;
            }
        },
        kneser: {
            desc: 'Kneser graph',
            longDesc: 'The graph representing k-element subsets of n elements whch are related if they are disjoint.',
            link: 'http://en.wikipedia.org/wiki/Kneser_graph',
            math: '\\(K\\!\\left(n,k\\right)\\)',
            args: {
                n: {
                    min: 1,
                    def: 5
                },
                k: {
                    min: 0,
                    def: 2
                }
            },
            verts: '\\(\\binom n k\\)',
            edges: '\\(\\binom n k\\binom{n-k}k/2\\)',
            make: (n:number, k:number) => {
                if (k > n)
                    return <GraphTemplate>{ vs: [], es: [], layout: graphLayouts.radial(0) };
                var elements = _.range(k);
                var sets:number[][] = [];
                var done = false;
                while (!done) {
                    sets.push(elements.slice(0));
                    for (var i = k - 1; i >= 0; i--) {
                        if (i < k - 1)
                            elements[i + 1] = elements[i] + 2;
                        for (var j = i + 1; j + 1 < k; j++)
                            elements[j + 1] = elements[j] + 1;
                        if (++elements[i] > n - (k - i - 1) - 1) {
                            if (i === 0) {
                                done = true;
                                break;
                            }
                        }
                        else
                            break;
                    }
                }
                var gt:GraphTemplate = { vs: [], es: [], layout: null };
                sets.forEach(function(s, i) {
                    gt.vs.push(i);
                    for (var j = sets.length - 1; j > i; j--)
                        if (_.intersection(s, sets[j]).length === 0)
                            gt.es.push([i, j]);
                });
                gt.layout = graphLayouts.radial(gt.vs.length);
                return gt;
            }
        },
        turan: {
            desc: 'Turán graph',
            longDesc: 'The graph of order n that does not contain an (r+1)-clique and is of maximal size.',
            link: 'http://en.wikipedia.org/wiki/Tur%C3%A1n_graph',
            math: '\\(T\\!\\left(n,r\\right)\\)',
            args: {
                n: {
                    min: 1,
                    def: 6
                },
                r: {
                    min: 1,
                    def: 3
                }
            },
            verts: '\\(n\\)',
            edges: '\\(\\left\\lfloor\\frac{\\left(r-1\\right)n^2}{2r}\\right\\rfloor\\)',
            make: (n:number, r:number) => {
                if (r > n)
                    return <GraphTemplate>{ vs: [], es: [], layout: graphLayouts.radial(0) };
                var gt:GraphTemplate = { vs: _.range(n), es: [], layout: graphLayouts.radial(n) };
                var sets:number[][] = [];
                for (var i = 0, s = 0; i < n; s++) {
                    var setSize = s < n % r ? Math.floor(n / r) + 1 : Math.floor(n / r);
                    sets.push([]);
                    for (var j = 0; j < setSize; j++)
                        sets[s].push(j + i);
                    i += setSize;
                }
                for (var i = 0; i < sets.length; i++)
                    for (var j = i + 1; j < sets.length; j++)
                        sets[i].forEach((v1:number) => {
                            sets[j].forEach((v2:number) => {
                                gt.es.push([v1, v2]);
                            });
                        });
                return gt;
            }
        }
    };

    var staticGraphs:{[index:string]:StaticGraph} = {
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
        octahedron: {
            desc: 'octahedron Schlegel diagram',
            longDesc: 'The graph representing the skeleton of the octahedron.',
            link: 'http://en.wikipedia.org/wiki/Octahedron',
            verts: 6,
            edges: [
                [0, 1],
                [1, 2],
                [2, 0],
                [3, 0],
                [3, 2],
                [3, 4],
                [3, 5],
                [4, 0],
                [4, 1],
                [4, 5],
                [5, 1],
                [5, 2]
            ],
            layout: {
                name: 'multiRadial',
                args: [6, [0, 3, 0, 0, 3], [0, Math.PI / 3, 0, 0, 0]]
            }
        },
        icosahedron: {
            desc: 'icosahedron Schlegel diagram',
            longDesc: 'The graph representing the skeleton of the icosahedron.',
            link: 'http://en.wikipedia.org/wiki/Regular_icosahedron#Icosahedral_graph',
            verts: 12,
            edges: [
                [0, 1],
                [1, 2],
                [2, 0],
                [3, 0],
                [3, 2],
                [3, 4],
                [3, 8],
                [4, 0],
                [4, 5],
                [5, 0],
                [5, 1],
                [5, 6],
                [6, 1],
                [6, 7],
                [7, 1],
                [7, 2],
                [7, 8],
                [8, 2],
                [9, 3],
                [9, 4],
                [9, 8],
                [9, 10],
                [9, 11],
                [10, 4],
                [10, 5],
                [10, 6],
                [10, 11],
                [11, 6],
                [11, 7],
                [11, 8]
            ],
            layout: {
                name: 'multiRadial',
                args: [12, [0, 3, 6, 0, 0, 3], [0, Math.PI / 3, -Math.PI / 6, 0, 0, 0]]
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
    _.each(staticGraphs, (g:StaticGraph, name:string) => {
        commonGraphs[name] = {
            desc: g.desc,
            longDesc: g.longDesc,
            link: g.link,
            math: g.math,
            args: {},
            verts: '\\(' + g.verts + '\\)',
            edges: '\\(' + g.edges.length + '\\)',
            make: () => <GraphTemplate>{
                vs: _.range(g.verts),
                es: g.edges,
                layout: graphLayouts[g.layout.name].apply(null, g.layout.args)
            }
        };
    });

    var width:number;
    var height:number;
    // TODO: move things when canvas resizes
    // TODO: consider canvas height more carefully, maybe should fill view vertically (use $(window).resize(handler) and window.innerHeight)

    var svg = d3.select('#canvas');

    var $canvas = $('#canvas');
    width = $canvas.width();
    height = $canvas.height();

    // data for nodes and links
    var nodes:Vertex[] = [];
    var links:Edge[] = [];

    // sets the given node as selected
    function selectNode(node:D3.Selection) {
        console.log('select');
        node.classed('node-selected', true);
        node.select('circle')
            .transition()
            .duration(300)
            .attr('r', 13);
    }

    // sets the given node as deselected
    function deselectNode() {
        console.log('deselect');
        var node = svg.selectAll('.node-selected')
            .classed('node-selected', false);
        node.selectAll('circle')
            .transition()
            .duration(100)
            .attr('r', 8);
    }

    function hasSelection():boolean {
        return !svg.selectAll('.node-selected').empty();
    }

    function addNode(node:Vertex) {
        console.log('add node');
        nodes.push(node);
        linkNodes(node);
    }

    function linkNodes(node:Vertex):boolean {
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
                    links = _.filter(links, (link:Edge) =>
                        !((link.source === a && link.target === b) ||
                        (link.source === b && link.target === a)));
                    links.forEach((link, i) => {
                        link.edgeNum = i;
                    });
                }
                return true;
            }
        }
        return false;
    }

    function removeNode(node:Vertex) {
        console.log('remove');
        var selected = d3.selectAll('.node-selected');
        var a = node;
        if (!selected.empty() && a === selected.datum())
            deselectNode();
        nodes = _.filter(nodes, (n:Vertex) => n !== a);
        links = _.filter(links, (l:Edge) => !(l.source === a || l.target === a));
        nodes.forEach((node, i) => {
            node.vertNum = i;
        });
        links.forEach((link, i) => {
            link.edgeNum = i;
        });
    }

    // clears all data
    function clear() {
        console.log('clear');
        nodes.splice(0, nodes.length);
        links.splice(0, links.length);
        deselectNode();
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
    interactLayer.append('g').attr('class', 'node-group');
    passiveLayer.append('g').attr('class', 'link-group');

    // updates properties table
    function calcProps() {
        _.each(graphProps, (prop:GraphProp<any>) => {
            console.log('updating prop ' + prop.name);
            prop.value = prop.calc(nodes, links);
            console.log(prop.value);
        });
        graphPropRow.each(function(prop:GraphProp<any>) {
            var valStr = '';
            var val = prop.value;
            if (prop.repr) {
                valStr = prop.repr(val);
            }
            else if (val === undefined) {
                valStr = '';
            }
            else if (val === null) {
                valStr = '?';
            }
            else if (isNaN(val)) {
                if (val instanceof Matrix) {
                    valStr = '\\(' + (<Matrix>val).toLatex() + '\\)';
                }
                else {
                    valStr = '' + val;
                }
            }
            else if (!isFinite(val)) {
                valStr = '\\(' + (val === Number.NEGATIVE_INFINITY ? '-' : '') + '\\infty\\)';
            }
            else {
                // TODO: val.toString() can sometimes produce scientific notation, don't use commas then
                valStr = val.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
            }
            var valueSpan = d3.select(this).select('td .value');
            valueSpan.text(valStr);
            if (valueSpan.text() === 'true')
                valueSpan.style('color', '#0a0');
            else if (valueSpan.text() === 'false')
                valueSpan.style('color', '#a00');
        });
        updateTex();
    }

    // updates view
    function update() {
        updateObjData();
        updateObjs();
        if (d3.select('#auto-calc-props-checkbox').property('checked'))
            calcProps();
    }

    var nodeDrag = d3.behavior.drag()
        .on('dragstart', function(d:Point) {
            d.dragstart = d3.mouse(this);
        })
        .on('drag', function(d:Point) {
            var m = d3.mouse(this);
            d.x += m[0] - d.dragstart[0];
            d.y += m[1] - d.dragstart[1];
            updateObjs();
        })
        .on('dragend', function(d:Point) {
            delete d.dragstart;
        });

    // updates old and new nodes and links based on current data
    function updateObjData() {
        var linkObjs = svg.select('.link-group').selectAll('.link').data(links);
        var nodeObjs = svg.select('.node-group').selectAll('.node').data(nodes);

        nodeObjs.select('circle');
        nodeObjs.select('text');

        linkObjs.enter()
            .insert('line', '.node')
            .attr('class', 'link')
            .attr('opacity', 0)
            .transition()
            .duration(200)
            .attr('opacity', 1);

        var newNodes = nodeObjs.enter()
            .append('g')
            .attr('class', 'node')
            .call(nodeDrag);
        newNodes.append('circle')
            .attr('r', 0)
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
                        selectNode(d3.select(this.parentNode));
                }
            })
            .transition()
            .duration(200)
            .attr('r', 8);
        newNodes.append('text')
            .attr('text-anchor', 'middle')
            .attr('baseline-shift', '-3px')
            .attr('font-size', '10px')
            .attr('opacity', toolbarButtons['labels'].toggleState ? 1.0 : 0.0)
            .text(function(d:Vertex) {
                return d.vertNum + 1;
            });

        linkObjs.exit()
            .transition()
            .duration(200)
            .attr('opacity', 0)
            .remove();

        var oldNodes = nodeObjs.exit()
            .transition()
            .duration(200)
            .remove();
        oldNodes.selectAll('circle')
            .attr('r', 0);
    }

    // updates current nodes and links
    function updateObjs() {
        svg.selectAll('.link')
            .attr('x1', function(d:Edge) {
                return d.source.x;
            })
            .attr('y1', function(d:Edge) {
                return d.source.y;
            })
            .attr('x2', function(d:Edge) {
                return d.target.x;
            })
            .attr('y2', function(d:Edge) {
                return d.target.y;
            });
        svg.selectAll('.node')
            .attr('transform', function(d:Vertex) {
                return 'translate(' + d.x + ',' + d.y + ')';
            });
    }

    // tests if two vertices are adjacent
    function isAdj(v1:Vertex, v2:Vertex) {
        return _.any(links,
            (link:Edge) => (link.source === v1 && link.target === v2) || (link.source === v2 && link.target === v1));
    }

    // setup toolbar
    var toolbarButtonSize = 40;
    var toolbarButtons:{[index:string]:ToolbarButton} = {
        help: {
            img: 'img/help.svg',
            onclick: function() {
                $('#help').css({ 'background-color': 'hsla(240,100%,90%,0.9)' });
                $('body').scrollTo('#help', 500, function() {
                    $('#help').animate({ 'background-color': 'hsla(240,100%,90%,0.0)' }, 500);
                });
            },
            title: 'Help'
        },
        clear: {
            img: 'img/erase.svg',
            onclick: function() {
                clear();
                update();
            },
            title: 'Clear'
        },
        complement: {
            img: 'img/complement.svg',
            onclick: function() {
                var old_links = links.slice(0);
                links.splice(0, links.length);
                var e = 0;
                for (var i = 0; i < nodes.length; i++)
                    for (var j = i + 1; j < nodes.length; j++) {
                        var n1 = nodes[i];
                        var n2 = nodes[j];
                        if (!_.any(old_links, (e:Edge) =>
                            (e.source === n1 && e.target === n2) ||
                            (e.source === n2 && e.target === n1))) {
                            links.push({ target: n1, source: n2, edgeNum: e++ });
                        }
                    }
                update();
            },
            title: 'Graph complement'
        },
        linegraph: {
            img: 'img/linegraph.svg',
            onclick: function() {
                var old_links = links.slice(0);
                links.splice(0, links.length);
                nodes.splice(0, nodes.length);
                old_links.forEach((l, i) => {
                    nodes.push({ x: (l.target.x + l.source.x) / 2, y: (l.target.y + l.source.y) / 2, vertNum: i });
                });
                var e = 0;
                for (var i = 0; i < old_links.length; i++)
                    for (var j = i + 1; j < old_links.length; j++) {
                        var l1 = old_links[i];
                        var l2 = old_links[j];
                        if (l1.source === l2.source || l1.source === l2.target ||
                            l1.target === l2.source || l1.target === l2.target) {
                            links.push({ target: nodes[l1.edgeNum], source: nodes[l2.edgeNum], edgeNum: e++ });
                        }
                    }
                update();
            },
            title: 'Line graph'
        },
        labels: {
            img: 'img/numbers.svg',
            onclick: function(d:ToolbarButton) {
                svg.selectAll('.node').selectAll('text')
                    .transition()
                    .duration(500)
                    .attr('opacity', d.toggleState ? 1.0 : 0.0);
            },
            title: 'Show vertex labels',
            toggleable: true,
            toggleState: false
        },
        recalc: {
            img: 'img/calc.svg',
            onclick: function() {
                calcProps();
            },
            title: 'Recalculate properties'
        }
    };
    var toolbar = interactLayer.append('g').attr('class', 'toolbar')
        .selectAll('.toolbar-btn')
        .data(_.map(toolbarButtons, (button:ToolbarButton, name:string) => {
            button.name = name;
            return button;
        }))
        .enter()
        .append('g')
        .attr('class', 'toolbar-btn')
        .attr('transform', function(d:ToolbarButton, i) {
            return 'translate(' + (5 + (5 + toolbarButtonSize) * i) + ',5)';
        });
    toolbar.append('title')
        .text(function(d:ToolbarButton) {
            return d.title;
        });
    toolbar.append('circle')
        .attr('cx', toolbarButtonSize / 2)
        .attr('cy', toolbarButtonSize / 2)
        .attr('r', toolbarButtonSize / 2)
        .attr('fill', 'red')
        .attr('fill-opacity', function(d:ToolbarButton) {
            return d.toggleable ? (d.toggleState ? 1.0 : 0.0) : 0.0;
        })
        .attr('stroke', 'black')
        .attr('stroke-width', '2')
        .on('mouseover', function(d:ToolbarButton) {
            d3.select(this)
                .transition()
                .duration(100)
                .attr('fill-opacity', d.toggleable ? 0.5 : 1.0);
        })
        .on('mouseout', function(d:ToolbarButton) {
            d3.select(this)
                .transition()
                .duration(100)
                .attr('fill-opacity', d.toggleable ? (d.toggleState ? 1.0 : 0.0) : 0.0);
        })
        .on('click', function(d:ToolbarButton, i) {
            if (d.toggleable) {
                d.toggleState = !d.toggleState;
                d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('fill-opacity', d.toggleState ? 1.0 : 0.0);
            }
            d.onclick.apply(this, [d, i]);
        });
    toolbar.append('image')
        .attr('xlink:href', function(d:ToolbarButton) {
            return d.img;
        })
        .attr('width', toolbarButtonSize)
        .attr('height', toolbarButtonSize);


    // setup tables

    var graphPropRow = d3.select('#graph-prop-table tbody').selectAll('tr')
        .data(_.filter(_.map(graphProps, (prop:GraphProp<any>, name:string) => {
            prop.name = name;
            return prop;
        }), (prop:GraphProp<any>) => !prop.hidden))
        .enter()
        .append('tr');
    graphPropRow.append('td')
        .attr('class', 'desc')
        .text(function(prop:GraphProp<any>) {
            return prop.desc;
        })
        .filter(function(prop:GraphProp<any>) {
            return !!prop.longDesc;
        })
        .append('a')
        .attr('tabindex', 0)
        .attr('role', 'button')
        .attr('class', 'help-popover btn btn-xs btn-default')
        .attr('data-toggle', 'popover')
        .attr('data-trigger', 'focus')
        .attr('data-placement', 'bottom')
        .attr('data-html', true)
        .attr('title', function(prop:GraphProp<any>) {
            return prop.desc;
        })
        .attr('data-content', function(prop:GraphProp<any>) {
            var s = prop.longDesc;
            if (prop.link)
                s += ' <a href="' + prop.link + '" target="_blank">[more]</a>';
            return s;
        })
        .text('?');
    graphPropRow.append('td')
        .text(function(prop:GraphProp<any>) {
            return prop.math || '';
        });
    graphPropRow.append('td')
        .attr('class', 'value')
        .each(function(prop:GraphProp<any>) {
            var cell = d3.select(this);
            var toggle:D3.Selection = null;
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
        .data(_.map(commonGraphs, (graph:CommonGraph, name:string) => {
            graph.name = name;
            return graph;
        }))
        .enter()
        .append('tr');
    commonGraphsRow.append('td')
        .text(function(graph:CommonGraph) {
            return graph.desc;
        })
        .filter(function(graph:CommonGraph) {
            return !!graph.longDesc;
        })
        .append('a')
        .attr('tabindex', 0)
        .attr('role', 'button')
        .attr('class', 'help-popover btn btn-xs btn-default')
        .attr('data-toggle', 'popover')
        .attr('data-trigger', 'focus')
        .attr('data-placement', 'bottom')
        .attr('data-html', true)
        .attr('title', function(graph:CommonGraph) {
            return graph.desc;
        })
        .attr('data-content', function(graph:CommonGraph) {
            var s = graph.longDesc;
            if (graph.link)
                s += ' <a href="' + graph.link + '" target="_blank">[more]</a>';
            return s;
        })
        .text('?');
    commonGraphsRow.append('td')
        .text(function(graph:CommonGraph) {
            return graph.math || '';
        });
    commonGraphsRow.append('td')
        .text(function(graph:CommonGraph) {
            return graph.verts;
        });
    commonGraphsRow.append('td')
        .text(function(graph:CommonGraph) {
            return graph.edges;
        });
    commonGraphsRow.append('td')
        .attr('class', 'inputs')
        .each(function(graph:CommonGraph) {
            var cell = d3.select(this);
            var form = cell.append('form')
                .attr('class', 'form-inline');
            var formGroup = form.selectAll('.form-group')
                .data(_.map(graph.args, (arg:CommonGraphArg, name:string) => {
                    arg.name = name;
                    return arg;
                }))
                .enter()
                .append('div')
                .attr('class', 'form-group');
            formGroup.append('label')
                .text(function(arg:CommonGraphArg) {
                    return '\\(' + arg.name + '\\):';
                });
            formGroup.append('input')
                .attr('class', 'form-control')
                .attr('type', 'number')
                .attr('min', function(arg:CommonGraphArg) {
                    return arg.min;
                })
                .attr('max', function(arg:CommonGraphArg) {
                    return arg.max;
                })
                .attr('defaultValue', function(arg:CommonGraphArg) {
                    return arg.def;
                })
                .attr('value', function(arg:CommonGraphArg) {
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
                var g:GraphTemplate = graph.make.apply(null, args);
                console.log(g.vs.length + ' verts: ' + JSON.stringify(g.vs));
                console.log(g.es.length + ' edges: ' + JSON.stringify(g.es));
                clear();
                g.vs.forEach((v, i) => {
                    var pos = g.layout(i);
                    console.log(v + ' ' + JSON.stringify(pos));
                    nodes[v] = { x: pos.x + width / 2, y: pos.y + height / 2, vertNum: i };
                });
                g.es.forEach((e, i) => {
                    links.push({ source: nodes[e[0]], target: nodes[e[1]], edgeNum: i });
                });
                update();
                {
                    var rect = $('#canvas')[0].getBoundingClientRect();
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

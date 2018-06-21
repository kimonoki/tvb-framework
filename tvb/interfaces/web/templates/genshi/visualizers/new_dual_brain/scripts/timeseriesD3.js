/**
 * TheVirtualBrain-Framework Package. This package holds all Data Management, and
 * Web-UI helpful to run brain-simulations. To use it, you also need do download
 * TheVirtualBrain-Scientific Package (for simulators). See content of the
 * documentation-folder for more details. See also http://www.thevirtualbrain.org
 *
 * (c) 2012-2017, Baycrest Centre for Geriatric Care ("Baycrest") and others
 *
 * This program is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with this
 * program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/*

 tv.js should dump just a single public var named tv (T.VB V.isualizations)

 tv = {}

 with

 tv.ndar     array fun
 tv.plot     reusable plotting components
 tv.util     utility stuff

 */

/* global tv, d3 */

//added globals for time selection
var timeselection_interval = 0;
var timeselection = [];

// identify the initiator of the change of the time selection: brushing or movie timeline
var triggered_by_timeselection = true;
//store the unmapped selection value used to animate the time selection window
var selection_x = [];

tv = {};

tv.util = {

    // d3 style configurator. if this is slow, interp and eval source
    gen_access: function (obj, field) {
        return function (maybe) {
            if (maybe === undefined) {
                return obj["_" + field];
            } else {
                obj["_" + field] = maybe;
                return obj;
            }
        };
    },

    // helper to add usage notes to plots
    usage: function (root, heading, notes) {
        const p = root.append("p");
        p.classed("slice-info", true);
        p.append("h3").classed("instructions", true).text(heading);
        p.append("ul").selectAll("li").data(notes)
            .enter().append("li").classed("instructions", true).text(function (d) {
            return d;
        });
    },

    ord_nums: ["zeroeth", "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth",
        "eleventh", "twelfth", "thirteenth", "fourteenth", "fifteenth", "sixteenth", "seventeenth", "eighteenth", "nineteenth"],

    /* f is a templater/formatter cf. https://gist.github.com/984375 */
    fmt: function (f) { // fhe format specifier followed by any number of arguments

        var a = arguments; // store outer arguments
        return ("" + f) // force format specifier to String
            .replace( // replace tokens in format specifier
                /\{(?:(\d+)|(\w+))\}/g, // match {token} references
                function (s, // the matched string (ignored)
                          i, // an argument index
                          p // a property name
                ) {
                    return p && a[1] // if property name and first argument exist
                        ? a[1][p] // return property from first argument
                        : a[i]; // assume argument index and return i-th argument
                });
    },

    get_array_shape: function (baseURL, callback) {
        $.getJSON(baseURL + "/read_data_shape/False?kwd=0", callback);
    },

    get_array_slice: function (baseURL, slices, callback, channels, currentMode, currentStateVar) {
        var readDataURL = readDataChannelURL(baseURL, slices[0].lo, slices[0].hi,
            currentStateVar, currentMode, slices[0].di, JSON.stringify(channels));
        //NOTE: If we need to add slices for the other dimensions pass them as the 'specific_slices' parameter.
        //      Method called is from time_series.py.
        $.getJSON(readDataURL, callback);
    }
};

tv.ndar = function (data) {

    this.data = data;

    this.imap = function (f) {
        for (var i = 0; i < this.data.length; i++) {
            this.data[i] = f(this.data[i]);
        }
        return this;
    };

    this.map = function (f) {
        return (new tv.ndar(this.data.slice())).imap(f);
    };

    this.reduce = function (f, init) {
        for (var i = 0; i < this.data.length; i++) {
            init = f(init, this.data[i]);
        }
        return init;
    };

    this.max = function () {
        return this.reduce(function (l, r) {
            return l > r ? l : r;
        }, -1e300);
    };

    this.min = function () {
        return this.reduce(function (l, r) {
            return l < r ? l : r;
        }, 1e300);
    };

    this.sum = function () {
        return this.reduce(function (l, r) {
            return l + r;
        }, 0);
    };

    this.mean = function () {
        return this.sum() / this.length();
    };

    this.std = function () {
        var mean_sqr = this.map(function (x) {
                return x * x;
            }).mean(),
            mean = this.mean();
        return Math.sqrt(mean_sqr - mean * mean);
    };

    this.add = function (b) {
        return this.map(function (x) {
            return x + b;
        });
    };

    this.sub = function (b) {
        return this.add(-b);
    };

    this.mul = function (b) {
        return this.map(function (x) {
            return x * b;
        });
    };

    this.imul = function (b) {
        return this.imap(function (x) {
            return x * b;
        });
    };

    this.idiv = function (b) {
        return this.imul(1 / b);
    };

    this.div = function (b) {
        return this.mul(1 / b);
    };

    this.get = function (i) {
        return this.data[i];
    };

    this.set = function (i, val) {
        this.data[i] = val;
    };

    this.nd2lin = function (idx) {
        var l = 0;
        for (var i = 0; i < idx.length; i++) {
            l += this.strides[i] * idx[i];
        }
        return l;
    };

    this.length = function () {
        return this.data.length;
    };

    // return indices where condition is true
    this.where = function (f) {
        var indices = [];
        for (var i = 0; i < this.data.length; i++) {
            if (f(this.data[i], i)) {
                indices.push(i);
            }
        }
        return indices;
    };

    this.pretty_step = function (base) {
        return Math.pow(base, Math.floor(-1 + Math.log(this.max() - this.min()) / Math.log(base)));
    };

    this.pretty_ticks = function (base) {
        var d = this.pretty_step(base || 10), f = Math.floor;
        return tv.ndar.range(f(this.min() / d) * d, (f(this.max() / d) + 1) * d, d);
    };

    this.pretty_ticklabels = function (base) {
        return this.pretty_ticks(base).map(function (d) {
            return d.toPrecision(2);
        });
    };

    this.normalized = function () {
        var mn = this.min(), mx = this.max();
        return this.map(function (d) {
            return (d - mn) / (mx - mn);
        });
    };

    this.slice = function (lo, hi) {
        return tv.ndar.from(this.data.slice(lo, hi));
    };

};

tv.ndar.from = function (src) {
    return new tv.ndar(src);
};

tv.ndar.ndfrom = function (src) {
    var a = tv.ndar.from(src.data);
    a.shape = src.shape;
    a.strides = src.strides;
    return a;
};

tv.ndar.range = function (a, b, c) {
    var lo, hi, dx;

    if ((a || a === 0) && b) {
        if (c) {
            dx = c;
        }
        else {
            dx = 1;
        }
        lo = a;
        hi = b;
    } else {
        hi = a;
        lo = 0;
        dx = 1;
    }

    var end = Math.floor((hi - lo) / dx);
    var ar = new tv.ndar([]);
    for (var i = 0; i < end; i++) {
        ar.data[i] = dx * i + lo;
    }
    return ar;

};

tv.ndar.zeros = function (n) {
    return tv.ndar.range(n).imap(function () {
        return 0.0;
    });
};

tv.ndar.ones = function (n) {
    return tv.ndar.zeros(n).add(1.0);
};


tv.plot = {

    time_series: function () {

        var f = function (root) {

            f.p(f.p() || 0.1); // pad
            f.w(f.w() || 700);
            f.h(f.h() || 500);
            f.point_limit(f.point_limit() || 500);

            f.magic_fcs_amp_scl = 1;

            // make sure we got numbers not strings
            f.dt(+f.dt());
            f.t0(+f.t0());

            // Create the required UI elements.
            var svg = root.append("svg").attr("width", f.w()).attr("height", f.h());
            var rgp = svg.append("g").attr("transform", "scale(1, 1)");

            rgp.append("g").append("rect").attr("width", f.w()).attr("height", f.h()).classed("tv-fig-bg", true);

            f.status_line = svg.append("g").attr("transform", "translate(10, " + (f.h() - 10) + ")").append("text");

            // parts independent of data
            f.compute_layout();
            f.add_resizer(svg, rgp);
            f.do_scaffolding(rgp);

            // inversion of flow control in progress
            f.we_are_setup = false;
            f.render();
        }; // end function f()

        f.render = function () {
            f.status_line.text("waiting for data from server...");
            //console.log(f.baseURL(), f.current_slice())
            tv.util.get_array_slice(f.baseURL(), f.current_slice(), f.render_callback, f.channels(), f.mode(), f.state_var());
        };

        f.render_callback = function (data) {

            var kwd = kwd || {};

            f.status_line.text("handling data...");

            /* reformat data into normal ndar style */
            var flat = []
                , sl = f.current_slice()[0]
                , shape = [(sl.hi - sl.lo) / sl.di, f.shape()[2]]
                , strides = [f.shape()[2], 1];

            for (var i = 0; i < shape[0]; i++) {
                for (var j = 0; j < shape[1]; j++) {
                    flat.push(data[i][j]);
                }
            }

            var ts = [], t0 = f.t0(), dt = f.dt();

            for (var ii = 0; ii < shape[0]; ii++) {
                ts.push(t0 + dt * sl.lo + ii * dt * sl.di);
            }

            f.ts(tv.ndar.ndfrom({data: ts, shape: [shape[0]], strides: [1]}));
            f.ys(tv.ndar.ndfrom({data: flat, shape: shape, strides: strides}));

            f.status_line.text("examining data...");
            f.prepare_data();
            f.status_line.text("rendering data...");
            f.render_focus();

            if (!f.we_are_setup) {
                f.render_contexts();
                f.add_brushes();
                f.br_fcs_endfn(true); // no_render=true
                f.we_are_setup = true;
            }

            f.status_line.text("");
        };

        f.current_slice = function () {
            var dom = f.sc_fcs_x.domain()
                , lo = Math.floor((dom[0] - f.t0()) / f.dt())
                , hi = Math.floor((dom[1] - f.t0()) / f.dt())
                , di = Math.floor((hi - lo) / (2 * f.point_limit()));

            di = di === 0 ? 1 : di;

            if (lo > f.shape()[0]) {
                console.log("time_series.current_slice(): found lo>shape[0]: " + lo + ">" + f.shape()[0]);
                lo = f.shape()[0];
            }

            return [{lo: lo, hi: hi, di: di}];
        };

        // dimensions and placement of focus and context areas
        f.compute_layout = function () {
            // pad is only provisionally basis for dimensioning the context areas; later
            // we will need to have inner and outer pad
            f.pad = {x: (0 ? f.w() : f.h()) * f.p(), y: f.h() * f.p()};
            f.ul_ctx_y = {x: f.pad.x, y: f.pad.y};
            f.sz_ctx_y = {x: f.pad.x * 0.8, y: f.h() - 3 * f.pad.y - f.pad.y};
            f.ul_ctx_x = {x: f.pad.x, y: 2 * f.pad.y + f.sz_ctx_y.y};
            f.sz_ctx_x = {x: f.w() - 2 * f.pad.x, y: f.pad.y / 2};
            f.ul_fcs = {x: f.ul_ctx_x.x, y: f.ul_ctx_y.y};
            f.sz_fcs = {x: f.sz_ctx_x.x, y: f.sz_ctx_y.y};

        };

        // allows user to scale plot size dynamically
        // TODO refactor place in tv.util
        f.add_resizer = function (svg, rgp) {

            var resize_start;

            rgp.append("g").append("rect").classed("tv-resizer", true)
                .on("mouseover", function () {
                    rgp.attr("style", "cursor: se-resize");
                })
                .on("mouseout", function () {
                    rgp.attr("style", "");
                })
                .attr("x", f.w() - f.pad.x / 2).attr("y", f.h() - f.pad.y / 2)
                .attr("width", f.pad.x / 2).attr("height", f.pad.y / 2)
                .call(d3.drag().on("drag", function () {
                    var p1 = d3.mouse(svg.node())
                        , p2 = resize_start
                        , scl = {x: p1[0] / p2[0], y: p1[1] / p2[1]};
                    rgp.attr("transform", "scale(" + scl.x + ", " + scl.y + ")");
                    svg.attr("width", scl.x * f.w()).attr("height", scl.y * f.h());
                }).on("start", function () {
                    resize_start = d3.mouse(rgp.node());
                }));
        };

        // TODO migrate to tv.util
        var new_clip_path = function (el, id) {
            return el.append("defs").append("clipPath").attr("id", id);
        };

        f.mouse_scroll = function () {
            var ev = window.event
                , da = ev.detail ? ev.detail : ev.wheelDelta
                , sh = ev.shiftKey
                , dr = !!(da > 0);

            if (sh) {
                f.magic_fcs_amp_scl *= dr ? 1.2 : 1 / 1.2;
                // TODO scale transform instead via direct access...
                f.prepare_data();
                f.render_focus();
            } else {
                if (!(f.gp_br_fcs.node().__brush === null)) {
                    var dx = dr ? 1 : -1;
                    // stop scrolling if it is the end of the signals' list
                    if (f.dom_y[0] >= -1 && f.dom_y[1] <= f.channels().length) {
                        f.dom_y[0] += dx;
                        f.dom_y[1] += dx;
                    }
                    //lower bound
                    else if (f.dom_y[0] < -1) {
                        var delta = Math.abs(f.dom_y[0] - (-1));
                        f.dom_y[0] += delta;
                        f.dom_y[1] += delta;
                    }
                    //upper bound
                    else if (f.dom_y[1] > f.channels().length) {
                        var delta = Math.abs(f.channels().length - f.dom_y[1]);
                        f.dom_y[0] -= delta;
                        f.dom_y[1] -= delta;
                    }

                    //redraw the lines
                    var dom = f.dom_y;
                    var yscl = f.sz_fcs.y / (dom[1] - dom[0]) / 5;
                    f.sc_fcs_y.domain(dom).range([f.sz_ctx_y.y, 0]);
                    f.gp_ax_fcs_y.call(f.ax_fcs_y);
                    f.gp_lines.selectAll("g").attr("transform", function (d, i) {
                        return "translate(0, " + f.sc_fcs_y(i) + ")" + "scale (1, " + yscl + ")"
                    }).selectAll("path").attr("stroke-width", "" + (3 / yscl));
                    f.scale_focus_stroke();


                }


            }


        };

        f.signal_tick_labeler = function (tick_value) {
            return (tick_value % 1 === 0) ? f.labels()[tick_value] : "";
        };

        // setup groups, scales and axes for context and focus areas
        f.do_scaffolding = function (rgp) {

            // main groups for vertical and horizontal context areas and focus area
            f.gp_ctx_x = rgp.append("g").attr("transform", "translate(" + f.ul_ctx_x.x + ", " + f.ul_ctx_x.y + ")");
            f.gp_ctx_x.append("rect").attr("width", f.sz_ctx_x.x).attr("height", f.sz_ctx_x.y).classed("tv-data-bg", true);

            f.gp_fcs = rgp.append("g").attr("transform", "translate(" + f.ul_fcs.x + ", " + f.ul_fcs.y + ")");
            f.gp_fcs.on("mousewheel", f.mouse_scroll);
            f.gp_fcs.append("rect").attr("width", f.sz_fcs.x).attr("height", f.sz_fcs.y).classed("tv-data-bg", true);


            // the plotted time series in the focus and x ctx area are subject to a clipping region
            new_clip_path(rgp, "fig-lines-clip").append("rect").attr("width", f.sz_fcs.x).attr("height", f.sz_fcs.y);
            // new_clip_path(rgp, "fig-ctx-x-clip").append("rect").attr("width", f.sz_ctx_x.x).attr("height", f.sz_ctx_x.y);

            // group with clip path applied for the focus lines
            f.gp_lines = f.gp_fcs.append("g").attr("style", "clip-path: url(#fig-lines-clip)")
                .append("g").classed("line-plot", true);

            // scales for vertical and horizontal context, and the x and y axis of the focus area
            f.sc_ctx_y = d3.scaleLinear().domain([-1, f.shape()[2]]).range([f.sz_ctx_y.y, 0]);
            f.sc_ctx_x = d3.scaleLinear().domain([f.t0(), f.t0() + f.dt() * f.shape()[0]]).range([0, f.sz_ctx_x.x]);
            f.sc_fcs_x = d3.scaleLinear().domain([f.t0(), f.t0() + f.dt() * f.shape()[0]]).range([0, f.sz_fcs.x]);
            f.sc_fcs_y = d3.scaleLinear().domain([-1, f.shape()[2] + 1]).range([f.sz_fcs.y, 0]);


            f.dom_x = f.sc_ctx_x.domain();
            f.dom_y = f.sc_ctx_y.domain();

            // axes for each of the above scales
            f.ax_ctx_x = d3.axisBottom(f.sc_ctx_x);
            f.ax_fcs_x = d3.axisTop(f.sc_fcs_x);
            f.ax_fcs_y = d3.axisLeft(f.sc_fcs_y);

            f.ax_fcs_y.tickFormat(f.signal_tick_labeler);

            // groups for each of the above axes
            f.gp_ax_ctx_x = f.gp_ctx_x.append("g").classed("axis", true).call(f.ax_ctx_x)
                .attr("transform", "translate(0, " + f.sz_ctx_x.y + ")");
            f.gp_ax_fcs_x = f.gp_fcs.append("g").classed("axis", true).call(f.ax_fcs_x);
            f.gp_ax_fcs_y = f.gp_fcs.append("g").classed("axis", true).call(f.ax_fcs_y);

        };

        f.prepare_data = function () {

            var ts = f.ts();
            var ys = f.ys();
            var da_lines = [];
            var line_avg;
            var ys_std = ys.min();
            //To set this properly, we need to know:
            // nsig - how many signals on the screen?
            // std  - std of signals
            // pxav - vertical pixels available

            for (var sig_idx = 0; sig_idx < ys.shape[1]; sig_idx++) {

                da_lines[sig_idx] = [];
                for (var t_idx = 0; t_idx < ys.shape[0]; t_idx++) {
                    da_lines[sig_idx][t_idx] = ys.data[ys.strides[0] * t_idx + sig_idx];
                }

                line_avg = d3.mean(da_lines[sig_idx]);
                for (var tt_idx = 0; tt_idx < ys.shape[0]; tt_idx++) {
                    da_lines[sig_idx][tt_idx] = f.magic_fcs_amp_scl * (da_lines[sig_idx][tt_idx] - line_avg) / Math.abs(ys_std);
                    // multiply by -1 because the y axis points down
                    da_lines[sig_idx][tt_idx] *= -1;
                }


                da_lines[sig_idx] = {sig: da_lines[sig_idx], id: sig_idx};
            }

            // compute context data
            var da_x = []
                , da_xs = []
                , da_y = []
                , ys_mean = ys.mean()
                , ys_std = ys.std()
                , n_chan = ys.shape[1]
                , datum;

            // center an average signal
            for (var j = 0; j < ts.shape[0]; j++) {
                da_x[j] = 0;
                da_xs[j] = 0;
                for (var i = 0; i < n_chan; i++) {
                    datum = ys.data[j * n_chan + i];
                    da_x [j] += datum;
                    da_xs[j] += datum * datum;
                }
                da_xs[j] = Math.sqrt(da_xs[j] / n_chan - ((da_x[j] / n_chan) * (da_x[j] / n_chan)));
                da_x [j] = (da_x[j] / n_chan - ys_mean);
                // multiply by -1 because y axis points down
                da_x[j] *= -1;

                if ((isNaN(da_x[j])) || (isNaN(da_xs[j]))) {
                    console.log("encountered NaN in data: da_x[" + j + "] = " + da_x[j] + ", da_xs[" + j + "] = " + da_xs[j] + ".");
                }
            }

            // scale average signal by ptp
            var _dar = new tv.ndar(da_x);
            var da_max = _dar.max()
                , da_min = _dar.min()
                , da_ptp = da_max - da_min;

            for (var si = 0; si < da_x.length; si++) {
                da_x[si] = da_x[si] / da_ptp;
            }

            // center and scale the std line
            da_xs.min = tv.ndar.from(da_xs).min();
            for (var jj = 0; jj < da_xs.length; jj++) {
                da_xs[jj] -= da_xs.min;
                da_xs[jj] /= ys_std;
                // multiply by -1 because y axis points down
                da_xs[jj] *= -1;
            }

            // center and scale to std each signal
            for (var jjj = 0; jjj < n_chan; jjj++) {
                da_y[jjj] = [];
                // This computes a slice at the beginning of the signal to be displayed on the y axis
                // The signal might be shorter than the width hence the min
                for (var ii = 0; ii < Math.min(f.sz_ctx_y.x, ys.shape[0]); ii++) {
                    da_y[jjj][ii] = (ys.data[ii * n_chan + jjj] - ys_mean) / ys_std;
                    // multiply by -1 because y axis points down
                    da_y[jjj][ii] *= -1;
                }
            }

            f.da_lines = da_lines;
            f.da_x_dt = f.dt() * f.current_slice()[0].di;
            f.da_x = da_x;
            f.da_xs = [0, da_xs[da_xs.length - 1]].concat(da_xs, [0]); // filled area needs start == end
            f.da_y = da_y;
        };

        f.render_focus = function () {

            var ts = f.ts()
                , g = f.gp_lines.selectAll("g").data(f.da_lines, function (d) {
                return d.id;
            });


            if (!f.we_are_setup) {


                f.line_paths = g.enter()
                    .append("g")
                    .attr("transform", function (d, i) {
                        return "translate(0, " + f.sc_fcs_y(i) + ")";
                    })
                    .append("path")
                    .attr("vector-effect", "non-scaling-stroke");
            }


            f.line_paths.attr("d", function (d) {
                return d3.line()
                    .x(function (d, i) {
                        return f.sc_ctx_x(ts.data[i]);
                    })
                    .y(function (d) {
                        return d;
                    })
                    (d.sig);
            });


        };

        f.render_contexts = function () {

            // originally used to draw context lines and average


        };

        f.scale_focus_stroke = function () {
            var total = f.sz_fcs
                , xdom = f.sc_fcs_x.domain()
                , ydom = f.sc_fcs_y.domain()
                , dx = xdom[1] - xdom[0]
                , dy = ydom[1] - ydom[0]
                , area = dx * dy
                , area2 = total.x * total.y;

            //console.log(area / area2);
            if (window.navigator.userAgent.indexOf("Edge") > -1) {
                f.gp_lines.selectAll("g").selectAll("path").attr("stroke-width", "0.3px");//4*Math.sqrt(Math.abs(area / area2)))
            } else {
                f.gp_lines.selectAll("g").selectAll("path").attr("stroke-width", "1px");//4*Math.sqrt(Math.abs(area / area2)))
            }
        };


        f.add_brushes = function () {

            // horizontal context brush
            var br_ctx_x_fn = function () {

                    var event_selection_x = [];
                    // Different extent when it is:
                    //1.from the brush of 2D Focus Brush
                    if (d3.event.selection != null && d3.event.selection[0][0] != null) {
                        event_selection_x[0] = d3.event.selection[0][0];
                        event_selection_x[1] = d3.event.selection[1][0];
                    }
                    //2.from the end of focus brush
                    else if (d3.event.selection == null) {
                        event_selection_x = [f.sc_ctx_x.range()[0], f.sc_ctx_x.range()[1]];
                        f.dom_x = [f.t0(), f.t0() + f.dt() * f.shape()[0]];
                    }
                    //3.from itself
                    else {
                        event_selection_x = d3.event.selection;
                    }


                    var scale_brushed = d3.scaleLinear().domain(f.dom_x).range(f.sc_ctx_x.range());


                    //selection is now in coordinates and we have to map it using scales
                    event_selection_x = event_selection_x.map(scale_brushed.invert, scale_brushed);


                    dom = f.br_ctx_x === null ? f.sc_ctx_x.domain() : event_selection_x;

                    f.dom_x = dom;

                    sc = f.sc_fcs_x;
                    x_scaling = scale_brushed.domain()[1] / (dom[1] - dom[0]);
                    sc.domain(dom);
                    f.sc_ctx_x.domain(dom);
                    f.gp_ax_fcs_x.call(f.ax_fcs_x);
                    f.gp_ax_ctx_x.call(f.ax_ctx_x);


                    // TODO: This seems to cause problems with negative values and commenting it out does not seem to
                    // cause any additional problems. This could do with some double checking.
                    f.gp_lines.attr("transform", "translate(" + sc(0) + ", 0) scale(" + x_scaling + ", 1)");
                }

                // vertical changes
                , br_ctx_y_fn = function () {

                    var event_selection_y = [];

                    if (d3.event == null || d3.event.selection == null) {
                        event_selection_y = f.sc_ctx_y.range();
                        f.dom_y = [-1, f.shape()[2]];
                    }
                    else if (d3.event.selection != null && d3.event.selection[0][0] != null) {
                        event_selection_y[1] = d3.event.selection[0][1];
                        event_selection_y[0] = d3.event.selection[1][1];
                    }
                    else {
                        event_selection_y[0] = d3.event.selection[1];
                        event_selection_y[1] = d3.event.selection[0];
                    }

                    var scale_brushed = d3.scaleLinear().domain(f.dom_y).range(f.sc_ctx_y.range());


                    event_selection_y = event_selection_y.map(scale_brushed.invert, scale_brushed);
                    var dom = f.br_ctx_y === null ? f.sc_ctx_y.domain() : event_selection_y;
                    f.dom_y = dom;
                    var yscl = f.sz_fcs.y / (dom[1] - dom[0]) / 5;
                    f.sc_fcs_y.domain(dom).range([f.sz_ctx_y.y, 0]);
                    f.gp_ax_fcs_y.call(f.ax_fcs_y);
                    f.gp_lines.selectAll("g").attr("transform", function (d, i) {
                        return "translate(0, " + f.sc_fcs_y(i) + ")" + "scale (1, " + yscl + ")"
                    }).selectAll("path").attr("stroke-width", "" + (3 / yscl));

                };

            f.br_ctx_y_fn = br_ctx_y_fn;

            br_ctx_end = function () {

                //get the selected time range
                var event_selection_x = [];
                if (d3.event.selection != null) {
                    event_selection_x[0] = d3.event.selection[0];
                    event_selection_x[1] = d3.event.selection[1];
                    selection_x = event_selection_x;
                }
                var scale_brushed = d3.scaleLinear().domain(f.dom_x).range(f.sc_ctx_x.range());
                event_selection_x = event_selection_x.map(scale_brushed.invert, scale_brushed);
                dom = f.br_ctx_x === null ? f.sc_ctx_x.domain() : event_selection_x;
                timeselection = event_selection_x;

                // remove the last time's selection
                f.gp_ctx_x.selectAll(".selected-time").remove();

                //change the actual time point in the slider
                if (d3.event.selection != null) {
                    f.timeselection_update_fn(triggered_by_timeselection)
                }


            };

            // on end of focus brush
            // this is on f so that f can call it when everything else is done..
            f.br_fcs_endfn = function (no_render) {
                if (!d3.event || !d3.event.sourceEvent) {
                    br_ctx_y_fn();
                    f.scale_focus_stroke();
                    return;
                }
                br_ctx_x_fn();
                br_ctx_y_fn();
                f.gp_br_fcs.node().__brush.selection = null;
                f.gp_br_fcs.call(f.br_fcs);
                f.scale_focus_stroke();


            };


            f.br_fcs_startfn = function () {
                // we will use the left upper of the brush to do a tooltip

                //select a channel
                var event_selection_y = [];
                event_selection_y[1] = d3.event.selection[0][1];
                event_selection_y = event_selection_y.map(f.sc_ctx_y.invert);

                //choose the time point
                var event_selection_x = [];
                event_selection_x[1] = d3.event.selection[0][0];
                event_selection_x = event_selection_x.map(f.sc_ctx_x.invert);
                if (event_selection_x[1] < 0) {
                    event_selection_x[1] = 0 + f.da_x;
                }


                timerange = f.sc_fcs_x.domain()[1];
                channelID = parseInt(event_selection_y[1]);
                timepoint_length = f.da_lines[channelID].sig.length;

                timepoint = event_selection_x[1] / f.sc_fcs_x.domain()[1];
                timepoint = timepoint * timepoint_length;
                timepoint = parseInt(timepoint);

                valuearray = f.ys().data;
                channel_number = f.channels().length;
                channel_index = f.channels().indexOf(channelID);

                //print out the channel name(label) and value
                $("#info-channel").html(' ' + f.labels()[parseInt(event_selection_y[1])]);
                $("#info-time").html(" " + timepoint);
                $("#info-value").html(" " + valuearray[channel_number * timepoint + channel_index]);

            }


            // create brushes
            f.br_ctx_x = d3.brushX().extent([[f.sc_ctx_x.range()[0], 0], [f.sc_ctx_x.range()[1], f.sz_ctx_x.y]]).on("end", br_ctx_end);
            f.br_fcs = d3.brush().extent([[f.sc_fcs_x.range()[0], 0], [f.sc_fcs_x.range()[1], f.sz_fcs.y]])
                .on("end", f.br_fcs_endfn).on("start", f.br_fcs_startfn)
                .on("brush", f.br_fcs_brush);

            // add time selection brush group
            f.gp_br_ctx_x = f.gp_ctx_x.append("g");
            //add title for the time selection area
            f.timeselection_title = f.gp_br_ctx_x.append("text").text("Time Selection").attr("y", -10);
            f.gp_br_ctx_x.classed("brush", true).attr("class", "time-selection-brush").call(f.br_ctx_x).selectAll("rect").attr("height", f.sz_ctx_x.y);


            //add main focus brush group
            f.gp_br_fcs = f.gp_fcs.append("g").classed("brush", true).call(f.br_fcs);


        };


        //functions for the time selection window
        f.timeselection_update_fn = function (triggered) {

            //display the selected time range
            f.text = f.gp_ctx_x.append("text").attr("class", "selected-time").attr("id", "time-selection")
                .text("Selected Time Range: " + timeselection[0].toFixed(2) + "ms" + " to  " + timeselection[1].toFixed(2) + "ms");
            f.text_interval = f.gp_ctx_x.append("text").attr("class", "selected-time").attr("id", "time-selection-interval").text(" Interval:" + (timeselection[1] - timeselection[0]).toFixed(2)).attr("x", 100).attr("y", -10);

            if (triggered) {
                timeselection_interval = timeselection[1] - timeselection[0];

                //update the time in the input tag
                d3.select("#TimeNow").property('value', timeselection[0].toFixed(2));

                //update the time in the 3d viewer's time
                $('#slider').slider('value', timeselection[0].toFixed(2));
                loadFromTimeStep(parseInt(timeselection[0]));

            }

        };

        //move the time selection window with the slider
        f.timeselection_move_fn = function () {
            redrawSelection()
        };


        //increase and decease the interval by dt, need minus dt brought by the triggered change
        f.timeselection_interval_increase = function () {
            d3.select(f.gp_br_ctx_x.node()).call(f.br_ctx_x.move, [timeselection[0] - f.dt(), timeselection[1]].map(f.sc_ctx_x));
        }

        f.timeselection_interval_decrease = function () {
            d3.select(f.gp_br_ctx_x.node()).call(f.br_ctx_x.move, [timeselection[0]- f.dt(), timeselection[1] - 2*f.dt()].map(f.sc_ctx_x));

        }


        //need fix one additional step by any change
        function redrawSelection() {
            //>1 *timeStepsPerTick
            if (timeStepsPerTick > 1) {
                d3.select(f.gp_br_ctx_x.node()).call(f.br_ctx_x.move, [timeselection[0] + f.dt() * timeStepsPerTick, timeselection[1] + f.dt() * timeStepsPerTick].map(f.sc_ctx_x));
            }
            //<1  1/2 *0.33
            else if (timeStepsPerTick < 1) {
                d3.select(f.gp_br_ctx_x.node()).call(f.br_ctx_x.move, [timeselection[0] + f.dt() * 1 / (1 / timeStepsPerTick + 1), timeselection[1] + f.dt() * 1 / (1 / timeStepsPerTick + 1)].map(f.sc_ctx_x));
            }
            else if (timeStepsPerTick === 1) {
                d3.select(f.gp_br_ctx_x.node()).call(f.br_ctx_x.move, [timeselection[0] + f.dt(), timeselection[1] + f.dt()].map(f.sc_ctx_x));
            }
        }


        f.parameters = ["w", "h", "p", "baseURL", "preview", "labels", "shape",
            "t0", "dt", "ts", "ys", "point_limit", "channels", "mode", "state_var"];
        f.parameters.map(function (name) {
            f[name] = tv.util.gen_access(f, name);
        });

        return f;
    },


};

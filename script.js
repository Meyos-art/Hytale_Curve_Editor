// ── curve-sampling.js ──

// Numerical sampling of arbitrary functions. Error estimates are not analytic bounds.
(function (root) {
  'use strict';

  function sampleFunction(f, lo, hi, setting, value) {
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) throw Error('Invalid sampling range.');
    const point = x => {
      const y = f(x);
      if (!Number.isFinite(y)) throw Error('Non-finite value in the curve.');
      return [x, y];
    };
    if (setting === 'count') {
      if (!Number.isInteger(value) || value < 2 || value > 4097) throw Error('Choose an point count between 2 and 4097.');
      return Array.from({
        length: value
      }, (_, i) => point(i === value - 1 ? hi : lo + (hi - lo) * i / (value - 1)));
    }
    if (!Number.isFinite(value) || value < 1e-6 || value > 1) throw Error('Choose an precision between 0.000001 and 1.');
    const result = [point(lo)];
    function split(a, b, depth) {
      const probes = [.25, .5, .75].map(t => point(a[0] + (b[0] - a[0]) * t));
      const error = Math.max(...probes.map((p, i) => Math.abs(p[1] - (a[1] + (b[1] - a[1]) * (i + 1) / 4))));
      if (error > value / 2) {
        if (depth >= 40 || probes[1][0] === a[0] || probes[1][0] === b[0]) throw Error('Cannot reach this precision. Try a larger tolerance.');
        split(a, probes[1], depth + 1);
        split(probes[1], b, depth + 1);
      } else {
        if (result.length >= 4097) throw Error('Sampling needs more than 4097 points. Try a larger tolerance.');
        result.push(b);
      }
    }
    // Seed the full interval, then refine using multiple probes to catch inflections.
    for (let i = 0; i < 32; i++) split(result.at(-1), point(i === 31 ? hi : lo + (hi - lo) * (i + 1) / 32), 0);
    return result;
  }
  const api = {
    sampleFunction
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CurveSampling = api;
})(typeof globalThis === 'undefined' ? this : globalThis);


// ── curve-core.js ──

/* Local calculations; no network requests. */
(function (root) {
  'use strict';

  const {
    sampleFunction
  } = typeof module !== 'undefined' && module.exports ? require('./curve-sampling.js') : root.CurveSampling;
  const finite = (v, label) => {
    if (typeof v !== 'number' || !Number.isFinite(v)) throw Error(label + ' must be a finite number.');
    return v;
  };
  function compile(n) {
    if (!n || typeof n !== 'object') throw Error('Expected a JSON node.');
    const type = n.Type;
    const param = (k, d) => finite(n[k] ?? d, k);
    if (type === 'Constant') {
      const v = param('Value', 0);
      return () => v;
    }
    if (type === 'Manual') {
      if (n.TailMode && n.TailMode !== 'Constant') throw Error('Manual: only TailMode Constant is supported.');
      if (!Array.isArray(n.Points) || n.Points.length < 2) throw Error('Manual: at least two points are required.');
      const pts = n.Points.map(p => [finite(p.In ?? p.X ?? p.x ?? p[0], 'X coordinate'), finite(p.Out ?? p.Y ?? p.y ?? p[1], 'Y coordinate')]).sort((a, b) => a[0] - b[0]);
      if (pts.some((p, i) => i && p[0] === pts[i - 1][0])) throw Error('Two points have the same X coordinate.');
      return x => {
        if (x <= pts[0][0]) return pts[0][1];
        if (x >= pts.at(-1)[0]) return pts.at(-1)[1];
        let lo = 0,
          hi = pts.length - 1;
        while (hi - lo > 1) {
          const m = lo + hi >> 1;
          if (pts[m][0] <= x) lo = m;else hi = m;
        }
        const a = pts[lo],
          b = pts[hi];
        return a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0]);
      };
    }
    if (type === 'Not' || type === 'Inverter') {
      if (!n.Curve) throw Error(type + ': no curve connected in the Curve field.');
      const f = compile(n.Curve);
      return type === 'Not' ? x => 1 - f(x) : x => -f(x);
    }
    if (type === 'Sum' || type === 'Multiplier') {
      const children = n.Curves ?? n.Inputs;
      if (!Array.isArray(children) || !children.length) throw Error(type + ': a Curves or Inputs list is required.');
      const fs = children.map(compile);
      return type === 'Sum' ? x => fs.reduce((a, f) => a + f(x), 0) : x => fs.reduce((a, f) => a * f(x), 1);
    }
    if (type === 'DistanceExponential') {
      const R = param('Range', 1),
        p = param('Exponent', 1);
      if (R <= 0 || p <= 0) throw Error('Range and Exponent must be > 0 for this inverse.');
      return x => x < 0 ? 1 : x > R ? 0 : Math.pow(1 - x / R, p);
    }
    if (type === 'DistanceS') {
      const R = param('Range', 1),
        A = param('ExponentA', 1),
        B = param('ExponentB', 1),
        T = param('Transition', 1),
        s = param('TransitionSmooth', 1);
      if (R <= 0 || A <= 0 || B <= 0 || T < 0 || T > 1 || s < 0 || s > 1) throw Error('DistanceS: parameters out of range (Range and exponents must be > 0).');
      if (T === 0 && A !== B) throw Error('Zero transition with different exponents: discontinuity, no continuous Manual inverse.');
      const a = R * (1 - T) / 2,
        b = a + T * R;
      return x => {
        if (x >= R) return 0;
        const u = 1 - x / R,
          fa = Math.pow(u, A),
          fb = Math.pow(u, B);
        if (x < a) return fa;
        if (x > b) return fb;
        if (a === b) return (fa + fb) / 2;
        const t = (x - a) / (b - a),
          w = (1 - s) * t + s * (1 - Math.cos(Math.PI * t)) / 2;
        return (1 - w) * fa + w * fb;
      };
    }
    throw Error('Unsupported type: ' + String(type) + '. Paste the curve, not the entire graph.');
  }
  function prepare(node, c) {
    for (const k of ['h0', 'h1', 'scale', 'offset', 'gain', 'bias', 'n0', 'n1']) finite(c[k], k);
    if (c.h1 <= c.h0 || c.n1 <= c.n0 || c.scale === 0 || c.gain === 0) throw Error('Increasing intervals and nonzero scale and gain are required.');
    const f = compile(node),
      density = h => c.gain * f(c.scale * h + c.offset) + c.bias;
    const knots = manualKnots(node),
      steps = 16384,
      hs = [],
      ds = [];
    let sign = 0;
    const heights = Array.from({
      length: steps + 1
    }, (_, i) => c.h0 + (c.h1 - c.h0) * i / steps);
    if (knots) heights.push(...knots.map(x => (x - c.offset) / c.scale).filter(h => h > c.h0 && h < c.h1));
    heights.sort((a, b) => a - b);
    const uniqueHeights = [...new Set(heights)];
    for (let i = 0; i < uniqueHeights.length; i++) {
      const h = uniqueHeights[i],
        d = density(h);
      if (!Number.isFinite(d)) throw Error('Non-finite density at height ' + h);
      hs.push(h);
      ds.push(d);
      if (i) {
        const diff = d - ds[i - 1];
        if (Math.abs(diff) > 1e-12) {
          const s = Math.sign(diff);
          if (sign && s !== sign) throw Error('Non-monotonic curve: multiple possible heights. Choose a height range limited to one branch.');
          sign = s;
        }
      }
    }
    if (!sign) throw Error('Constant curve: no unique inverse.');
    const minD = Math.min(ds[0], ds.at(-1)),
      maxD = Math.max(ds[0], ds.at(-1));
    const terrain = c.mode !== 'inverse';
    const valid0 = terrain ? -maxD : minD,
      valid1 = terrain ? -minD : maxD;
    const lo = Math.max(c.n0, valid0),
      hi = Math.min(c.n1, valid1);
    if (hi <= lo) throw Error('No invertible interval in the requested input. Available domain: [' + valid0 + ', ' + valid1 + '].');
    const polarity = c.mode === 'sum' ? -1 : 1;
    const normalize = h => polarity * (c.out === 'signed' ? 2 * (h - c.h0) / (c.h1 - c.h0) - 1 : (h - c.h0) / (c.h1 - c.h0));
    const denormalize = y => c.h0 + (c.h1 - c.h0) * (c.out === 'signed' ? (y * polarity + 1) / 2 : y * polarity);
    const target = n => terrain ? -n : n;
    function inverse(n) {
      const d = target(n);
      if (d < minD - 1e-12 || d > maxD + 1e-12) return null;
      if (d === ds[0]) return c.h0;
      if (d === ds.at(-1)) return c.h1;
      let l = c.h0,
        r = c.h1;
      for (let i = 0; i < 60; i++) {
        const m = (l + r) / 2;
        if (density(m) < d === sign > 0) l = m;else r = m;
      }
      return (l + r) / 2;
    }
    // Plateaus produce intervals, not a single inverse. Reject rather than silently select one.
    for (let i = 1; i < ds.length; i++) if (ds[i] === ds[i - 1]) {
      const n = terrain ? -ds[i] : ds[i];
      if (n >= lo && n <= hi) throw Error('Plateau in the selected range: the inverse is not unique. Narrow the height bounds.');
    }
    const automatic = !knots && c.sampling === 'tolerance';
    const count = knots || automatic ? 2 : Number(c.count);
    if (!Number.isInteger(count) || count < 2 || count > 4097) throw Error('Point count: an integer between 2 and 4097.');
    const entries = knots ? [lo, hi, ...knots.map(x => (x - c.offset) / c.scale).filter(h => h > c.h0 && h < c.h1).map(h => terrain ? -density(h) : density(h)).filter(x => x > lo && x < hi)] : Array.from({
      length: count
    }, (_, i) => lo + (hi - lo) * i / (count - 1));
    const points = [...new Set(entries)].sort((a, b) => a - b).map(x => [x, Number(normalize(inverse(x)).toPrecision(15))]);
    if (automatic) points.splice(0, points.length, ...sampleFunction(x => normalize(inverse(x)), lo, hi, 'tolerance', Number(c.tolerance)));
    let error = 0;
    for (let i = 0; i < points.length - 1; i++) for (let j = 1; j <= 8; j++) {
      const t = j / 9,
        x = points[i][0] * (1 - t) + points[i + 1][0] * t,
        y = points[i][1] * (1 - t) + points[i + 1][1] * t;
      error = Math.max(error, Math.abs(density(denormalize(y)) - target(x)));
    }
    const inputLo = c.normalizeInput ? c.out === 'signed' ? -1 : 0 : lo;
    const inputHi = c.normalizeInput ? 1 : hi;
    const toRawInput = q => c.normalizeInput ? lo + (q - inputLo) / (inputHi - inputLo) * (hi - lo) : q;
    const fromRawInput = n => c.normalizeInput ? inputLo + (n - lo) / (hi - lo) * (inputHi - inputLo) : n;
    const exportedPoints = c.normalizeInput ? points.map(([x, y]) => [Number(fromRawInput(x).toPrecision(15)), y]) : points;
    return {
      points: exportedPoints,
      error,
      exact: !!knots,
      lo: inputLo,
      hi: inputHi,
      rawLo: lo,
      rawHi: hi,
      toRawInput,
      fromRawInput,
      clipped: lo > c.n0 || hi < c.n1,
      inverse: q => inverse(toRawInput(q)),
      rawInverse: inverse,
      normalize,
      density,
      hs,
      ds,
      valid0,
      valid1
    };
  }
  function manualKnots(n) {
    if (n.Type === 'Manual') return n.Points.map(p => p.In ?? p.X ?? p.x ?? p[0]);
    if (['Not', 'Inverter'].includes(n.Type) && n.Curve) return manualKnots(n.Curve);
    return null;
  }
  function suggestedOutputRange(n) {
    while (n && ['Not', 'Inverter'].includes(n.Type)) n = n.Curve;
    if (n?.Type === 'Manual' && Array.isArray(n.Points)) return n.Points.some(p => (p.Out ?? p.Y ?? p.y ?? p[1]) < 0) ? 'signed' : 'unit';
    return 'unit';
  }
  function noiseInterval(node, c) {
    if (c.out === 'signed') return {
      n0: -1,
      n1: 1
    };
    const f = compile(node);
    const a = -(c.gain * f(c.scale * c.h0 + c.offset) + c.bias);
    const b = -(c.gain * f(c.scale * c.h1 + c.offset) + c.bias);
    return Math.max(a, b) <= 0 && Math.min(a, b) < 0 ? {
      n0: -1,
      n1: 0
    } : {
      n0: 0,
      n1: 1
    };
  }
  function infer(n, mode = 'terrain') {
    const f = compile(n),
      chain = [];
    let leaf = n;
    while (['Not', 'Inverter'].includes(leaf.Type)) {
      chain.push(leaf.Type);
      leaf = leaf.Curve;
    }
    let h0, h1;
    if (leaf.Type === 'Manual') {
      const xs = manualKnots(leaf);
      h0 = Math.min(...xs);
      h1 = Math.max(...xs);
    } else if (['DistanceS', 'DistanceExponential'].includes(leaf.Type)) {
      h0 = 0;
      h1 = leaf.Range ?? 1;
    } else throw Error('Cannot infer the range for ' + leaf.Type + '. Disable automatic mode to set the bounds.');
    chain.push(leaf.Type);
    const a = f(h0),
      b = f(h1),
      factor = mode === 'inverse' ? 1 : -1;
    return {
      h0,
      h1,
      n0: Math.min(factor * a, factor * b),
      n1: Math.max(factor * a, factor * b),
      scale: 1,
      offset: 0,
      gain: 1,
      bias: 0,
      chain: chain.reverse().join(' → ')
    };
  }
  function parseClipboard(text) {
    try {
      return JSON.parse(text);
    } catch (original) {
      const values = [];
      let start = -1,
        depth = 0,
        quoted = false,
        escaped = false;
      for (let i = 0; i < text.length; i++) {
        const ch = text[i];
        if (quoted) {
          if (escaped) escaped = false;else if (ch === '\\') escaped = true;else if (ch === '"') quoted = false;
          continue;
        }
        if (ch === '"') {
          quoted = true;
          continue;
        }
        if (ch === '{' || ch === '[') {
          if (depth === 0) start = i;
          depth++;
        } else if (ch === '}' || ch === ']') {
          depth--;
          if (depth < 0) throw original;
          if (depth === 0) {
            values.push(JSON.parse(text.slice(start, i + 1)));
            start = -1;
          }
        } else if (depth === 0 && !/\s/.test(ch)) throw original;
      }
      if (depth || quoted || !values.length) throw original;
      return values;
    }
  }
  function extract(value) {
    const results = [];
    function walk(n, path) {
      if (!n || typeof n !== 'object') return;
      if (n.Type === 'CurveMapper' && n.Curve) {
        results.push({
          curve: n.Curve,
          path: path + '.Curve',
          context: 'CurveMapper input: ' + (n.Inputs?.map(v => v.Type).join(', ') || 'not provided') + '. Slider/BaseHeight offsets are not applied automatically.'
        });
        return;
      }
      if (n.Type && ((n.$NodeId || '').endsWith('.Curve') || ['DistanceS', 'DistanceExponential', 'Manual'].includes(n.Type) || ['Not', 'Inverter', 'Sum', 'Multiplier'].includes(n.Type) && (n.Curve || n.Curves))) {
        results.push({
          curve: n,
          path,
          context: 'Curve extracted; other graph branches are not included.'
        });
        return;
      }
      for (const [k, v] of Object.entries(n)) if (v && typeof v === 'object') walk(v, path + (Array.isArray(n) ? '[' + k + ']' : '.' + k));
    }
    walk(value, '$');
    return results;
  }
  function exportManual(points, schema = 'hytale') {
    // Fill columns in curve order. Within each row, emit left nodes last so they overlap right nodes.
    const ordered = schema === 'hytale' ? [...points].sort((a, b) => a[0] - b[0]) : points;
    const Points = ordered.map(([x, y], i) => schema === 'hytale' ? {
      $NodeId: 'CurvePoint.Curve',
      $Position: {
        $x: 968 + Math.floor(i / 10) * 300,
        $y: i % 10 * 160
      },
      In: x,
      Out: y
    } : schema === 'inout' ? {
      In: x,
      Out: y
    } : schema === 'xy' ? {
      X: x,
      Y: y
    } : {
      X: x,
      Y: y,
      Z: 0
    });
    if (schema === 'hytale') Points.sort((a, b) => a.$Position.$y - b.$Position.$y || b.$Position.$x - a.$Position.$x);
    return schema === 'hytale' ? {
      $NodeId: 'Manual.Curve',
      $Position: {
        $x: 0,
        $y: 0
      },
      Type: 'Manual',
      TailMode: 'Constant',
      Points
    } : {
      Type: 'Manual',
      Points
    };
  }
  const api = {
    compile,
    prepare,
    parseClipboard,
    extract,
    exportManual,
    infer,
    suggestedOutputRange,
    noiseInterval
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CurveInverse = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);


// ── curve-editor-math.js ──

// Pure Bezier geometry: no DOM, events, or editor session state.
(function (root) {
  'use strict';

  const {
    sampleFunction
  } = typeof module !== 'undefined' && module.exports ? require('./curve-sampling.js') : root.CurveSampling;
  const clamp = (n, a, b) => Math.max(a, Math.min(b, n));
  function knots(points) {
    return points.map(([x, y], i) => ({
      x,
      y,
      linked: true,
      l: {
        x: x - (x - (points[i - 1]?.[0] ?? x)) / 3,
        y: y - (y - (points[i - 1]?.[1] ?? y)) / 3
      },
      r: {
        x: x + ((points[i + 1]?.[0] ?? x) - x) / 3,
        y: y + ((points[i + 1]?.[1] ?? y) - y) / 3
      }
    }));
  }
  function constrain(ps) {
    ps.forEach((p, i) => {
      // A handle may use the whole adjacent segment. The two controls of a
      // segment still keep their X order so a Manual Curve export is valid.
      p.l.x = clamp(p.l.x, i ? ps[i - 1].r.x : p.x, p.x);
      p.r.x = clamp(p.r.x, p.x, i < ps.length - 1 ? ps[i + 1].l.x : p.x);
    });
  }
  function smooth(ps) {
    const slopes = ps.slice(1).map((p, i) => (p.y - ps[i].y) / (p.x - ps[i].x));
    ps.forEach((p, i) => {
      let m = slopes[i === ps.length - 1 ? i - 1 : i];
      if (i > 0 && i < ps.length - 1) {
        const a = slopes[i - 1],
          b = slopes[i],
          h0 = p.x - ps[i - 1].x,
          h1 = ps[i + 1].x - p.x;
        m = a * b <= 0 ? 0 : 3 * (h0 + h1) / ((2 * h1 + h0) / a + (h1 + 2 * h0) / b);
      }
      const dl = i ? (p.x - ps[i - 1].x) / 3 : 0,
        dr = i < ps.length - 1 ? (ps[i + 1].x - p.x) / 3 : 0;
      p.l = {
        x: p.x - dl,
        y: p.y - m * dl
      };
      p.r = {
        x: p.x + dr,
        y: p.y + m * dr
      };
    });
  }
  function sample(ps, mode, tolerance = .001) {
    if (ps.length < 2 || ps.some((p, i) => !Number.isFinite(p.x) || !Number.isFinite(p.y) || i && p.x <= ps[i - 1].x)) throw Error('At least two points with strictly increasing In values are required.');
    if (mode === 'linear') return ps.map(p => [p.x, p.y]);
    if (!(tolerance >= 1e-6 && tolerance <= 1)) throw Error('Precision must be between 0.000001 and 1.');
    const result = [[ps[0].x, ps[0].y]],
      mid = (a, b) => ({
        x: (a.x + b.x) / 2,
        y: (a.y + b.y) / 2
      });
    function split(a, b, c, d, depth) {
      if (result.length > 32768) throw Error('Too many points: increase the Bezier tolerance.');
      const line = p => a.y + (d.y - a.y) * (p.x - a.x) / (d.x - a.x);
      if (Math.max(Math.abs(b.y - line(b)), Math.abs(c.y - line(c))) <= tolerance) {
        result.push([d.x, d.y]);
        return;
      }
      if (depth >= 24) throw Error('Cannot reach the requested precision for this segment.');
      const ab = mid(a, b),
        bc = mid(b, c),
        cd = mid(c, d),
        abc = mid(ab, bc),
        bcd = mid(bc, cd),
        center = mid(abc, bcd);
      split(a, ab, abc, center, depth + 1);
      split(center, bcd, cd, d, depth + 1);
    }
    for (let i = 0; i < ps.length - 1; i++) {
      const a = ps[i],
        d = ps[i + 1];
      if ([a.r.x, a.r.y, d.l.x, d.l.y].some(n => !Number.isFinite(n)) || a.r.x < a.x || a.r.x > d.l.x || d.l.x > d.x) throw Error('Handles must follow the order of In values.');
      split(a, a.r, d.l, d, 0);
    }
    return result;
  }
  // Split the segment with the largest vertical error bound first.
  // A cubic lies in the convex hull of its controls, so this also detects
  // S-shaped segments whose midpoint happens to lie on their chord.
  function sampleCount(ps, mode, count) {
    sample(ps, 'linear');
    if (!Number.isInteger(count) || count < ps.length || count > 4097) throw Error('Choose an integer point count between ' + ps.length + ' and 4097 (control points are preserved).');
    const heap = [],
      mid = (a, b) => ({
        x: (a.x + b.x) / 2,
        y: (a.y + b.y) / 2
      });
    const higher = (a, b) => a.error > b.error || a.error === b.error && a.width > b.width;
    function segment(a, b, c, d) {
      if ([b.x, b.y, c.x, c.y].some(n => !Number.isFinite(n)) || b.x < a.x || b.x > c.x || c.x > d.x) throw Error('Handles must follow the order of In values.');
      const line = p => a.y + (d.y - a.y) * (p.x - a.x) / (d.x - a.x);
      return {
        a,
        b,
        c,
        d,
        width: d.x - a.x,
        error: Math.max(Math.abs(b.y - line(b)), Math.abs(c.y - line(c)))
      };
    }
    function push(s) {
      heap.push(s);
      let i = heap.length - 1;
      while (i) {
        const p = i - 1 >> 1;
        if (!higher(heap[i], heap[p])) break;
        [heap[i], heap[p]] = [heap[p], heap[i]];
        i = p;
      }
    }
    function pop() {
      const best = heap[0],
        last = heap.pop();
      if (heap.length) {
        heap[0] = last;
        let i = 0;
        for (;;) {
          let j = 2 * i + 1;
          if (j >= heap.length) break;
          if (j + 1 < heap.length && higher(heap[j + 1], heap[j])) j++;
          if (!higher(heap[j], heap[i])) break;
          [heap[i], heap[j]] = [heap[j], heap[i]];
          i = j;
        }
      }
      return best;
    }
    for (let i = 0; i < ps.length - 1; i++) {
      const a = ps[i],
        d = ps[i + 1];
      push(mode === 'bezier' ? segment(a, a.r, d.l, d) : segment(a, {
        x: a.x + (d.x - a.x) / 3,
        y: a.y + (d.y - a.y) / 3
      }, {
        x: a.x + 2 * (d.x - a.x) / 3,
        y: a.y + 2 * (d.y - a.y) / 3
      }, d));
    }
    while (heap.length < count - 1) {
      const {
          a,
          b,
          c,
          d
        } = pop(),
        ab = mid(a, b),
        bc = mid(b, c),
        cd = mid(c, d),
        abc = mid(ab, bc),
        bcd = mid(bc, cd),
        p = mid(abc, bcd);
      if (!(p.x > a.x && p.x < d.x)) throw Error('Numerical precision limit: reduce the point count.');
      push(segment(a, ab, abc, p));
      push(segment(p, bcd, cd, d));
    }
    const errorBound = heap[0].error;
    heap.sort((a, b) => a.a.x - b.a.x);
    return {
      points: [[ps[0].x, ps[0].y], ...heap.map(s => [s.d.x, s.d.y])],
      errorBound
    };
  }
  function validateBounds(b) {
    if (['x0', 'x1', 'y0', 'y1'].some(k => !Number.isFinite(b[k])) || b.x0 >= b.x1 || b.y0 >= b.y1) throw Error('Bounds must be finite, with each minimum below its maximum.');
  }
  function fitCurve(ps, b) {
    validateBounds(b);
    const x0 = ps[0].x,
      x1 = ps.at(-1).x,
      ys = ps.flatMap(p => [p.y, p.l.y, p.r.y]),
      y0 = Math.min(...ys),
      y1 = Math.max(...ys);
    for (const p of ps) for (const v of [p, p.l, p.r]) {
      v.x = b.x0 + (v.x - x0) / (x1 - x0) * (b.x1 - b.x0);
      v.y = y1 === y0 ? (b.y0 + b.y1) / 2 : b.y0 + (v.y - y0) / (y1 - y0) * (b.y1 - b.y0);
    }
    ps[0].x = b.x0;
    ps.at(-1).x = b.x1;
    constrain(ps);
  }
  function applyBounds(ps, b, lockX, clampY) {
    validateBounds(b);
    if (lockX && (ps[0].x !== b.x0 || ps.at(-1).x !== b.x1)) {
      const x0 = ps[0].x,
        x1 = ps.at(-1).x;
      for (const p of ps) for (const v of [p, p.l, p.r]) v.x = b.x0 + (v.x - x0) / (x1 - x0) * (b.x1 - b.x0);
      ps[0].x = b.x0;
      ps.at(-1).x = b.x1;
    }
    if (clampY) for (const p of ps) for (const v of [p, p.l, p.r]) v.y = clamp(v.y, b.y0, b.y1);
    constrain(ps);
  }
  function alignHandles(ps, i, side = 'r', b = null) {
    const p = ps[i],
      other = side === 'r' ? 'l' : 'r',
      sign = side === 'r' ? 1 : -1;
    let dx = (p[side].x - p.x) * sign,
      dy = (p[side].y - p.y) * sign;
    if (Math.hypot(dx, dy) < 1e-14) {
      dx = (p.x - p[other].x) * sign;
      dy = (p.y - p[other].y) * sign;
    }
    if (dx < 0) {
      dx = 0;
    }
    if (Math.hypot(dx, dy) < 1e-14) {
      dx = 1;
      dy = 0;
    }
    if (b && i > 0 && i < ps.length - 1 && (p.y <= b.y0 || p.y >= b.y1)) dy = 0;
    if (Math.hypot(dx, dy) < 1e-14) dx = 1;
    const length = Math.hypot(dx, dy),
      ux = dx / length,
      uy = dy / length;
    for (const key of ['l', 'r']) {
      const dir = key === 'r' ? 1 : -1,
        neighbor = ps[i + dir];
      let distance = Math.hypot(p[key].x - p.x, p[key].y - p.y);
      if (!neighbor) {
        p[key] = {
          x: p.x,
          y: p.y
        };
        continue;
      }
      // Do not impose a half-segment limit. The shared segment constraint
      // above resolves the relation with the opposite control point.
      if (ux > 0) distance = Math.min(distance, Math.abs(neighbor.x - p.x) / ux);
      const vy = uy * dir;
      if (b && vy !== 0) distance = Math.min(distance, Math.max(0, (vy > 0 ? b.y1 - p.y : b.y0 - p.y) / vy));
      p[key] = {
        x: p.x + dir * ux * distance,
        y: p.y + dir * uy * distance
      };
    }
    // If the moved handle reaches far into a segment, the other control of
    // that segment is its only horizontal limit.
    constrain(ps);
  }
  const api = {
    knots,
    constrain,
    smooth,
    sample,
    sampleCount,
    fitCurve,
    applyBounds,
    alignHandles,
    sampleFunction,
    validateBounds
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CurveEditorMath = api;
})(typeof globalThis === 'undefined' ? this : globalThis);


// ── ui-shared.js ──

// Small shared UI services. The app remains dependency-free and CodePen-compatible.
(function (root) {
  'use strict';

  function createHistory(limit = 10) {
    const past = [];
    const future = [];
    const snapshot = value => JSON.parse(JSON.stringify(value));
    function push(value) {
      past.push(snapshot(value));
      if (past.length > limit) past.shift();
    }
    return {
      record(value) {
        push(value);
        future.length = 0;
      },
      undo(current) {
        if (!past.length) return null;
        future.push(snapshot(current));
        return past.pop();
      },
      redo(current) {
        if (!future.length) return null;
        push(current);
        return future.pop();
      }
    };
  }
  async function copyText(text) {
    try {
      await navigator.clipboard.writeText(text);
      return;
    } catch {
      // Supports the standalone file when Clipboard API access is unavailable.
      const focused = document.activeElement;
      const field = document.createElement('textarea');
      field.value = text;
      field.className = 'clipboard-fallback';
      document.body.append(field);
      field.select();
      let copied;
      try {
        copied = document.execCommand('copy');
      } finally {
        field.remove();
        focused?.focus({
          preventScroll: true
        });
      }
      if (!copied) throw Error('The browser blocked copying. Allow clipboard access.');
    }
  }
  function syncSamplingFields(select, count, tolerance) {
    const fixed = select.value === 'count';
    count.closest('label').hidden = !fixed;
    tolerance.closest('label').hidden = fixed;
    count.disabled = !fixed;
    tolerance.disabled = fixed;
    return fixed;
  }
  function createConfirmation(doc) {
    const get = id => doc.getElementById(id);
    const dialog = get('curve-confirm');
    let pending = null;
    function dismiss() {
      pending = null;
      dialog.close();
    }
    get('curve-confirm-cancel').onclick = dismiss;
    dialog.addEventListener('cancel', () => {
      pending = null;
    });
    get('curve-confirm-accept').onclick = () => {
      const action = pending;
      dismiss();
      action?.();
    };
    return function confirm(title, description, label, action) {
      for (const [id, text] of [['curve-confirm-title', title], ['curve-confirm-description', description]]) {
        const nodes = text.split(/(all points|every handle)/g).map(part => {
          if (part !== 'all points' && part !== 'every handle') return doc.createTextNode(part);
          const strong = doc.createElement('strong');
          strong.textContent = part;
          return strong;
        });
        get(id).replaceChildren(...nodes);
      }
      get('curve-confirm-accept').textContent = label;
      pending = action;
      dialog.showModal();
      get('curve-confirm-cancel').focus();
    };
  }
  const api = {
    createHistory,
    copyText,
    syncSamplingFields,
    createConfirmation
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CurveUI = api;
})(typeof globalThis === 'undefined' ? this : globalThis);


// ── app.js ──

(function () {
  'use strict';

  const {
    copyText,
    syncSamplingFields
  } = CurveUI;
  const $ = id => document.getElementById(id);
  $('inverter-help-toggle').onclick = () => {
    const panel = $('inverter-help');
    panel.hidden = !panel.hidden;
    $('inverter-help-toggle').setAttribute('aria-expanded', String(!panel.hidden));
  };
  const numberFields = ['h0', 'h1', 'scale', 'offset', 'gain', 'bias'];
  const range = () => $('range-signed').checked ? 'signed' : 'unit';
  const fmt = n => Number(n.toPrecision(7)).toString();
  let model = null,
    sourcePlot = null,
    exportText = '',
    timer = null,
    noiseBounds = {
      n0: -1,
      n1: 1
    };
  function chosen() {
    const found = CurveInverse.extract(CurveInverse.parseClipboard($('source').value));
    if (found.length === 1) return found[0].curve;
    if (found.length > 1) throw Error('Multiple independent curves: paste only the chain to invert, including its Curve children.');
    throw Error('No recognized curve in this JSON.');
  }
  function config(node) {
    const c = {
      mode: 'terrain',
      out: range(),
      n0: range() === 'signed' ? -1 : 0,
      n1: 1
    };
    for (const k of numberFields) {
      if ($(k).value.trim() === '') throw Error('Empty field: ' + k);
      c[k] = Number($(k).value);
      if (!Number.isFinite(c[k])) throw Error('Invalid number: ' + k);
    }
    c.count = Number($('count').value);
    c.sampling = $('inverse-output-mode').value;
    c.tolerance = Number($('inverse-tolerance').value);
    noiseBounds = CurveInverse.noiseInterval(node, c);
    return Object.assign(c, noiseBounds);
  }
  function makeSourcePlot(node, c) {
    if (c.h1 <= c.h0) throw Error('Maximum height must exceed minimum height.');
    const f = CurveInverse.compile(node),
      hs = Array.from({
        length: 1025
      }, (_, i) => c.h0 + (c.h1 - c.h0) * i / 1024);
    function knots(n) {
      if (!n || typeof n !== 'object') return;
      if (n.Type === 'Manual') for (const p of n.Points) {
        const h = ((p.In ?? p.X ?? p.x ?? p[0]) - c.offset) / c.scale;
        if (h > c.h0 && h < c.h1) hs.push(h);
      }
      if (n.Curve) knots(n.Curve);
      for (const child of n.Curves ?? n.Inputs ?? []) knots(child);
    }
    knots(node);
    hs.sort((a, b) => a - b);
    const points = [...new Set(hs)].map(h => [h, c.gain * f(c.scale * h + c.offset) + c.bias]);
    if (points.some(p => !Number.isFinite(p[1]))) throw Error('The curve produces a non-finite density within this range.');
    let y0 = Math.min(...points.map(p => p[1])),
      y1 = Math.max(...points.map(p => p[1]));
    if (y0 === y1) {
      y0 -= .5;
      y1 += .5;
    }
    return {
      points,
      x0: c.h0,
      x1: c.h1,
      y0,
      y1
    };
  }
  function chart(id, data, color) {
    const canvas = $(id),
      w = canvas.clientWidth,
      h = 290,
      dpr = window.devicePixelRatio || 1;
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, w, h);
    if (!data) {
      ctx.fillStyle = '#a0abba';
      ctx.font = '13px system-ui';
      ctx.textAlign = 'center';
      ctx.fillText('Waiting for nodes', w / 2, h / 2);
      return;
    }
    const {
        points,
        x0,
        x1,
        y0,
        y1
      } = data,
      L = 56,
      R = w - 20,
      T = 24,
      B = h - 38;
    const X = x => L + (x - x0) / (x1 - x0) * (R - L),
      Y = y => B - (y - y0) / (y1 - y0) * (B - T);
    ctx.font = '12px system-ui';
    const ticks = w < 360 ? 2 : 4;
    for (let i = 0; i <= ticks; i++) {
      const t = i / ticks,
        x = x0 + (x1 - x0) * t,
        y = y0 + (y1 - y0) * t;
      ctx.strokeStyle = '#37414f';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(X(x), T);
      ctx.lineTo(X(x), B);
      ctx.moveTo(L, Y(y));
      ctx.lineTo(R, Y(y));
      ctx.stroke();
      ctx.fillStyle = '#a0abba';
      ctx.textAlign = i === 0 ? 'left' : i === ticks ? 'right' : 'center';
      ctx.fillText(Number(x.toPrecision(4)).toString(), X(x), B + 23);
      ctx.textAlign = 'right';
      ctx.fillText(Number(y.toPrecision(4)).toString(), L - 9, Y(y) + 4);
    }
    ctx.save();
    ctx.beginPath();
    ctx.rect(L, T, R - L, B - T);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    points.forEach(([x, y], i) => i ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y)));
    ctx.stroke();
    if (points.length <= 20) for (const [x, y] of points) {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(X(x), Y(y), 3, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.restore();
  }
  function draw() {
    chart('source-plot', sourcePlot, '#a0abba');
    const lo = range() === 'signed' ? -1 : 0;
    chart('plot', model ? {
      points: model.points,
      x0: noiseBounds.n0,
      x1: noiseBounds.n1,
      y0: lo,
      y1: 1
    } : null, '#7ca5cc');
  }
  function probe() {
    if (!model) {
      $('probeResult').textContent = '';
      return;
    }
    const n = Number($('probe').value),
      lo = noiseBounds.n0,
      hi = noiseBounds.n1;
    if ($('probe').value === '' || n < lo || n > hi) {
      $('probeResult').textContent = 'Choose a value in [' + lo + ', ' + hi + '].';
      return;
    }
    const h = model.inverse(n);
    $('probeResult').textContent = h === null ? 'No surface for this noise value within the height range.' : 'Height ' + fmt(h) + ' → Out ' + fmt(model.normalize(h));
  }
  function reset() {
    model = null;
    sourcePlot = null;
    exportText = '';
    $('copy').disabled = true;
    $('copy-status').textContent = '';
    $('probeResult').textContent = '';
    draw();
  }
  function run() {
    clearTimeout(timer);
    reset();
    if (!$('source').value.trim()) {
      $('status').className = '';
      $('status').textContent = 'Paste your nodes to calculate their inverse.';
      return;
    }
    try {
      const node = chosen();
      let leaf = node;
      while (['Not', 'Inverter'].includes(leaf.Type) && leaf.Curve) leaf = leaf.Curve;
      $('inverse-output-settings').hidden = leaf.Type === 'Manual';
      syncSamplingFields($('inverse-output-mode'), $('count'), $('inverse-tolerance'));
      if ($('auto').checked) {
        const inferred = CurveInverse.infer(node, 'terrain');
        for (const k of numberFields) $(k).value = inferred[k];
        $('context').textContent = inferred.chain;
      }
      const c = config(node);
      sourcePlot = makeSourcePlot(node, c);
      draw();
      const span = c.h1 - c.h0;
      $('probe').min = c.n0;
      $('probe').max = c.n1;
      model = CurveInverse.prepare(node, c);
      exportText = JSON.stringify(CurveInverse.exportManual(model.points, 'hytale'), null, 2);
      $('copy').disabled = false;
      $('status').className = 'inverse-result';
      const height = c.out === 'signed' ? 'Surface height = ' + fmt(span / 2) + ' × (G + 1) + ' + fmt(c.h0) : 'Surface height = ' + fmt(span) + ' × G + ' + fmt(c.h0);
      $('status').textContent = 'Noise [' + c.n0 + ', ' + c.n1 + '] · ' + model.points.length + ' points · ' + height + ' · measured residual ' + model.error.toExponential(2);
      if (model.clipped) {
        const warning = document.createElement('span');
        warning.className = 'warn';
        warning.textContent = 'Surface available only for In ∈ [' + fmt(model.lo) + ', ' + fmt(model.hi) + ']. The export covers only this interval; its constant extension outside this range is not a solution.';
        $('status').append(warning);
      }
      draw();
      probe();
    } catch (e) {
      $('status').className = 'error';
      $('status').textContent = e.message;
    }
  }
  function schedule() {
    reset();
    $('context').textContent = '';
    $('status').className = '';
    $('status').textContent = 'Reading nodes…';
    clearTimeout(timer);
    timer = setTimeout(run, 300);
  }
  function autoFields() {
    for (const k of numberFields) $(k).disabled = $('auto').checked;
  }
  $('auto').onchange = () => {
    autoFields();
    schedule();
  };
  autoFields();
  for (const id of ['source', ...numberFields, 'count', 'inverse-output-mode', 'inverse-tolerance', 'range-unit', 'range-signed']) $(id).addEventListener('input', schedule);
  $('probe').oninput = probe;
  $('copy').onclick = async () => {
    if (!exportText) return;
    try {
      await copyText(exportText);
      $('copy-status').textContent = 'Nodes copied.';
    } catch (e) {
      $('copy-status').textContent = e.message;
    }
  };
  document.addEventListener('curve:viewchange', draw);
  const observer = new ResizeObserver(draw);
  observer.observe($('source-plot'));
  observer.observe($('plot'));
  draw();
})();


// ── curve-editor.js ──

// Curve editor controller: session state, rendering, and user interactions.
(function (root) {
  'use strict';

  if (typeof document === 'undefined') {
    if (typeof module !== 'undefined' && module.exports) module.exports = require('./curve-editor-math.js');
    return;
  }
  const {
    knots,
    constrain,
    smooth,
    sample,
    sampleCount,
    fitCurve,
    applyBounds,
    alignHandles,
    sampleFunction,
    validateBounds
  } = root.CurveEditorMath;
  const {
    copyText,
    createHistory,
    createConfirmation,
    syncSamplingFields
  } = root.CurveUI;
  const clone = value => JSON.parse(JSON.stringify(value));
  const clamp = (n, a, b) => Math.max(a, Math.min(b, n));
  const $ = id => document.getElementById(id),
    canvas = $('editor-canvas'),
    ctx = canvas.getContext('2d');
  let ps = knots([[-1, -1], [1, 1]]),
    mode = 'linear',
    selected = 0,
    bounds = null,
    drag = null,
    geometry = null,
    outputModel = null,
    squareView = false;
  let limits = {
      x0: -1,
      x1: 1,
      y0: -1,
      y1: 1
    },
    lockX = true,
    clampY = true;
  const num = n => Number(n.toPrecision(9));
  const state = () => ({
    ps: clone(ps),
    mode,
    selected,
    limits: {
      ...limits
    },
    lockX,
    clampY
  });
  const history = createHistory(10);
  function pushHistory(snapshot) {
    history.record(snapshot);
  }
  function remember() {
    pushHistory(state());
  }
  function restore(s) {
    ps = s.ps;
    mode = s.mode;
    selected = s.selected;
    limits = s.limits;
    lockX = s.lockX;
    clampY = s.clampY;
    syncBounds();
    fit();
    update();
  }
  function syncBounds() {
    for (const k of ['x0', 'x1', 'y0', 'y1']) $('bound-' + k).value = limits[k];
    $('bounds-preset').value = limits.x0 === -1 && limits.y0 === -1 && limits.x1 === 1 && limits.y1 === 1 ? 'signed' : limits.x0 === 0 && limits.y0 === 0 && limits.x1 === 1 && limits.y1 === 1 ? 'unit' : 'custom';
    $('bound-custom').hidden = $('bounds-preset').value !== 'custom';
    $('lock-endpoints').checked = lockX;
    $('clamp-output').checked = clampY;
  }
  function readBounds() {
    const b = {};
    for (const k of ['x0', 'x1', 'y0', 'y1']) {
      const field = $('bound-' + k);
      b[k] = field.value.trim() === '' ? NaN : Number(field.value);
    }
    validateBounds(b);
    return b;
  }
  function message(text, error = false) {
    $('editor-status').textContent = text;
    $('editor-status').className = error ? 'error' : 'hint';
  }
  function fit() {
    const ys = ps.flatMap(p => mode === 'bezier' ? [p.y, p.l.y, p.r.y] : [p.y]);
    const x0 = Math.min(ps[0].x, limits.x0),
      x1 = Math.max(ps.at(-1).x, limits.x1),
      y0 = Math.min(...ys, limits.y0),
      y1 = Math.max(...ys, limits.y1),
      dy = y1 - y0 || 1,
      dx = x1 - x0;
    bounds = {
      x0: x0 - dx * .08,
      x1: x1 + dx * .08,
      y0: y0 - dy * .14,
      y1: y1 + dy * .14
    };
  }
  function syncAspectButton() {
    const button = $('editor-aspect');
    button.textContent = squareView ? '▣' : '▭';
    button.setAttribute('aria-pressed', String(squareView));
    button.setAttribute('aria-label', squareView ? 'Square scale: on' : 'Rectangular view: on');
    button.title = squareView ? 'Square scale: on — switch to rectangular view' : 'Rectangular view: on — switch to square scale';
  }
  function makeOutput() {
    if (mode === 'linear') return {
      points: ps.map(p => [p.x, p.y]),
      errorBound: 0
    };
    if ($('editor-output-mode').value === 'count') return sampleCount(ps, mode, Number($('editor-count').value));
    const tolerance = Number($('editor-tolerance').value);
    return {
      points: sample(ps, mode, tolerance),
      errorBound: mode === 'linear' ? 0 : tolerance
    };
  }
  function update(editingKey = null) {
    document.querySelector('.export-settings').hidden = mode === 'linear';
    applyBounds(ps, limits, lockX, clampY);
    ps.forEach((p, i) => {
      if (p.linked) alignHandles(ps, i, 'r', clampY ? limits : null);
    });
    // Alignment can move a handle after bounds are applied. Restore the
    // left-to-right control order required by a Manual Curve export.
    constrain(ps);
    const endpoint = selected === 0 || selected === ps.length - 1;
    $('editor-x').disabled = lockX && endpoint;
    $('editor-x').title = lockX && endpoint ? 'Endpoint X is locked to the range boundary.' : '';
    $('bounds-status').textContent = '';
    $('editor-mode').value = mode;
    $('editor-handles').hidden = mode !== 'bezier';
    $('editor-bezier-tools').hidden = mode !== 'bezier';
    for (const id of ['editor-smooth', 'editor-vector', 'editor-auto', 'editor-tolerance']) $(id).disabled = mode !== 'bezier';
    const precisionOption = $('editor-output-mode').querySelector('option[value="tolerance"]');
    precisionOption.hidden = mode !== 'bezier';
    precisionOption.disabled = mode !== 'bezier';
    if (mode !== 'bezier') $('editor-output-mode').value = 'count';
    const fixed = syncSamplingFields($('editor-output-mode'), $('editor-count'), $('editor-tolerance'));
    $('editor-count').min = ps.length;
    $('editor-tolerance').disabled = fixed || mode === 'linear';
    $('editor-output-help').textContent = fixed ? 'Choose the total number of exported points. Endpoints and control points are included.' + (mode === 'bezier' ? ' More points are placed where the curve needs more detail.' : '') : 'Choose the maximum approximation error. Smaller values create more points and follow the Bezier curve more closely.';
    $('editor-selection').replaceChildren(...ps.map((p, i) => new Option('Point ' + (i + 1), String(i))));
    $('editor-selection').value = selected;
    $('editor-handle-mode').value = ps[selected].linked ? 'aligned' : 'free';
    $('editor-handle-mode').disabled = mode !== 'bezier';
    const p = ps[selected];
    for (const [id, v] of Object.entries({
      x: p.x,
      y: p.y,
      lx: p.l.x,
      ly: p.l.y,
      rx: p.r.x,
      ry: p.r.y
    })) if (id !== editingKey) $('editor-' + id).value = num(v);
    for (const k of ['lx', 'ly']) $('editor-' + k).disabled = selected === 0;
    for (const k of ['rx', 'ry']) $('editor-' + k).disabled = selected === ps.length - 1;
    $('editor-delete').disabled = ps.length <= 2 || lockX && endpoint;
    outputModel = null;
    try {
      outputModel = makeOutput();
      $('editor-copy').disabled = false;
      message(outputModel.points.length + ' exported points · Out error bounded by ' + outputModel.errorBound.toExponential(2) + ' (excluding floating-point rounding).');
    } catch (e) {
      $('editor-copy').disabled = true;
      message(e.message, true);
    }
    draw();
    drawPreview();
  }
  function drawPreview() {
    if ($('editor-view').hidden) return;
    const cv = $('editor-preview'),
      g = cv.getContext('2d'),
      w = cv.clientWidth,
      h = 220,
      dpr = window.devicePixelRatio || 1;
    cv.width = Math.round(w * dpr);
    cv.height = h * dpr;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    if (!outputModel) {
      g.fillStyle = '#b8b8b8';
      g.font = '12px system-ui';
      g.fillText('Preview unavailable: check the settings.', 15, 35);
      return;
    }
    const vals = ps.flatMap(p => mode === 'bezier' ? [p.y, p.l.y, p.r.y] : [p.y]);
    let ymin = Math.min(...vals),
      ymax = Math.max(...vals);
    const pad = (ymax - ymin || 1) * .08;
    ymin -= pad;
    ymax += pad;
    const x0 = ps[0].x,
      x1 = ps.at(-1).x,
      L = 45,
      R = w - 14,
      T = 14,
      B = h - 28;
    const X = x => L + (x - x0) / (x1 - x0) * (R - L),
      Y = y => B - (y - ymin) / (ymax - ymin) * (B - T);
    g.font = '10px system-ui';
    for (let i = 0; i <= 4; i++) {
      const t = i / 4,
        x = x0 + t * (x1 - x0),
        y = ymin + t * (ymax - ymin);
      g.strokeStyle = '#414141';
      g.lineWidth = 1;
      g.beginPath();
      g.moveTo(X(x), T);
      g.lineTo(X(x), B);
      g.moveTo(L, Y(y));
      g.lineTo(R, Y(y));
      g.stroke();
      g.fillStyle = '#aaa';
      g.textAlign = 'right';
      g.fillText(Number(y.toPrecision(3)), L - 5, Y(y) + 3);
    }
    g.textAlign = 'left';
    g.fillText(num(x0), L, h - 8);
    g.textAlign = 'right';
    g.fillText(num(x1), R, h - 8);
    g.save();
    g.beginPath();
    g.rect(L - 3, T - 3, R - L + 6, B - T + 6);
    g.clip();
    if ($('editor-reference').checked) {
      g.strokeStyle = '#d3d3d3';
      g.lineWidth = 2;
      g.setLineDash([4, 4]);
      g.beginPath();
      g.moveTo(X(ps[0].x), Y(ps[0].y));
      for (let i = 1; i < ps.length; i++) {
        const a = ps[i - 1],
          b = ps[i];
        if (mode === 'bezier') g.bezierCurveTo(X(a.r.x), Y(a.r.y), X(b.l.x), Y(b.l.y), X(b.x), Y(b.y));else g.lineTo(X(b.x), Y(b.y));
      }
      g.stroke();
      g.setLineDash([]);
    }
    g.strokeStyle = '#e4b477';
    g.lineWidth = 1.5;
    g.beginPath();
    outputModel.points.forEach(([x, y], i) => i ? g.lineTo(X(x), Y(y)) : g.moveTo(X(x), Y(y)));
    g.stroke();
    g.fillStyle = '#e4b477';
    for (const [x, y] of outputModel.points) {
      g.beginPath();
      g.arc(X(x), Y(y), 2.5, 0, Math.PI * 2);
      g.fill();
    }
    g.restore();
  }
  function draw() {
    if ($('editor-view').hidden) return;
    if (!bounds) fit();
    canvas.classList.toggle('is-rectangular', !squareView);
    const w = canvas.clientWidth,
      h = squareView ? w : 440,
      dpr = window.devicePixelRatio || 1;
    canvas.style.height = h + 'px';
    canvas.width = Math.round(w * dpr);
    canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const leftGutter = 58,
      rightGutter = 22,
      topGutter = 24,
      bottomGutter = 42,
      {
        x0,
        x1,
        y0,
        y1
      } = bounds;
    let L, R, T, B, viewX0, viewX1, viewY0, viewY1;
    if (squareView) {
      const size = Math.min(w - leftGutter - rightGutter, h - topGutter - bottomGutter),
        visibleSpan = Math.max(x1 - x0, y1 - y0);
      L = leftGutter + (w - leftGutter - rightGutter - size) / 2;
      R = L + size;
      T = topGutter + (h - topGutter - bottomGutter - size) / 2;
      B = T + size;
      viewX0 = (x0 + x1 - visibleSpan) / 2;
      viewX1 = viewX0 + visibleSpan;
      viewY0 = (y0 + y1 - visibleSpan) / 2;
      viewY1 = viewY0 + visibleSpan;
    } else {
      L = leftGutter;
      R = w - rightGutter;
      T = topGutter;
      B = h - bottomGutter;
      viewX0 = x0;
      viewX1 = x1;
      viewY0 = y0;
      viewY1 = y1;
    }
    const X = x => L + (x - viewX0) / (viewX1 - viewX0) * (R - L),
      Y = y => B - (y - viewY0) / (viewY1 - viewY0) * (B - T);
    geometry = {
      X,
      Y,
      from: (x, y) => ({
        x: viewX0 + (x - L) / (R - L) * (viewX1 - viewX0),
        y: viewY0 + (B - y) / (B - T) * (viewY1 - viewY0)
      }),
      L,
      R,
      T,
      B,
      x0: viewX0,
      x1: viewX1,
      y0: viewY0,
      y1: viewY1
    };
    ctx.font = '12px system-ui';
    for (let i = 0; i <= 20; i++) {
      const x = viewX0 + (viewX1 - viewX0) * i / 20,
        y = viewY0 + (viewY1 - viewY0) * i / 20;
      ctx.strokeStyle = i % 5 === 0 ? '#494949' : '#303030';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(X(x), T);
      ctx.lineTo(X(x), B);
      ctx.moveTo(L, Y(y));
      ctx.lineTo(R, Y(y));
      ctx.stroke();
    }
    for (let i = 0; i <= 4; i++) {
      const x = viewX0 + (viewX1 - viewX0) * i / 4,
        y = viewY0 + (viewY1 - viewY0) * i / 4;
      ctx.fillStyle = '#b3b3b3';
      ctx.textAlign = 'center';
      ctx.fillText(Number(x.toPrecision(4)), X(x), B + 24);
      ctx.textAlign = 'right';
      ctx.fillText(Number(y.toPrecision(4)), L - 8, Y(y) + 4);
    }
    ctx.save();
    ctx.beginPath();
    ctx.rect(L, T, R - L, B - T);
    ctx.clip();
    ctx.fillStyle = '#91bded10';
    ctx.fillRect(X(limits.x0), Y(limits.y1), X(limits.x1) - X(limits.x0), Y(limits.y0) - Y(limits.y1));
    ctx.strokeStyle = '#8aa8c3';
    ctx.lineWidth = 1;
    ctx.setLineDash([5, 4]);
    ctx.strokeRect(X(limits.x0), Y(limits.y1), X(limits.x1) - X(limits.x0), Y(limits.y0) - Y(limits.y1));
    ctx.setLineDash([]);
    ctx.strokeStyle = '#ededed';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(X(ps[0].x), Y(ps[0].y));
    for (let i = 1; i < ps.length; i++) {
      const a = ps[i - 1],
        b = ps[i];
      if (mode === 'bezier') ctx.bezierCurveTo(X(a.r.x), Y(a.r.y), X(b.l.x), Y(b.l.y), X(b.x), Y(b.y));else ctx.lineTo(X(b.x), Y(b.y));
    }
    ctx.stroke();
    const dot = (p, r, color) => {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(X(p.x), Y(p.y), r, 0, 2 * Math.PI);
      ctx.fill();
      ctx.strokeStyle = '#191d24';
      ctx.lineWidth = 2;
      ctx.stroke();
    };
    if (mode === 'bezier') {
      const p = ps[selected];
      for (const key of ['l', 'r']) {
        if (key === 'l' && selected === 0 || key === 'r' && selected === ps.length - 1) continue;
        ctx.strokeStyle = '#cdb790';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(X(p.x), Y(p.y));
        ctx.lineTo(X(p[key].x), Y(p[key].y));
        ctx.stroke();
        dot(p[key], 5, '#e4b477');
      }
    }
    ps.forEach((p, i) => dot(p, i === selected ? 6 : 4, i === selected ? '#e4b477' : '#a9a9a9'));
    if (lockX) for (const p of [ps[0], ps.at(-1)]) {
      const x = X(p.x) + 10,
        y = Y(p.y) - 14;
      ctx.strokeStyle = '#b8cde0';
      ctx.lineWidth = 1.4;
      ctx.strokeRect(x, y + 5, 9, 7);
      ctx.beginPath();
      ctx.arc(x + 4.5, y + 5, 3, Math.PI, 0);
      ctx.stroke();
    }
    ctx.restore();
  }
  function movePoint(i, x, y) {
    const p = ps[i],
      span = ps.at(-1).x - ps[0].x,
      eps = Math.max(span * 1e-7, 1e-10);
    if (lockX && (i === 0 || i === ps.length - 1)) x = p.x;
    if (clampY) y = clamp(y, limits.y0, limits.y1);
    x = clamp(x, i ? ps[i - 1].x + eps : -Infinity, i < ps.length - 1 ? ps[i + 1].x - eps : Infinity);
    const dx = x - p.x,
      dy = y - p.y;
    p.x = x;
    p.y = y;
    p.l.x += dx;
    p.r.x += dx;
    p.l.y += dy;
    p.r.y += dy;
    constrain(ps);
  }
  function add(x, y) {
    if (lockX && (x <= limits.x0 || x >= limits.x1)) {
      message('Add points inside the locked X range.', true);
      return;
    }
    if (clampY) y = clamp(y, limits.y0, limits.y1);
    if (ps.some(p => Math.abs(p.x - x) < (ps.at(-1).x - ps[0].x) * 1e-7)) {
      message('Choose an In value different from the existing points.', true);
      return;
    }
    remember();
    const p = {
      x,
      y,
      linked: true,
      l: {
        x,
        y
      },
      r: {
        x,
        y
      }
    };
    ps.push(p);
    ps.sort((a, b) => a.x - b.x);
    selected = ps.indexOf(p);
    const prev = ps[selected - 1],
      next = ps[selected + 1];
    if (prev) p.l = {
      x: x - (x - prev.x) / 3,
      y: y - (y - prev.y) / 3
    };
    if (next) p.r = {
      x: x + (next.x - x) / 3,
      y: y + (next.y - y) / 3
    };
    constrain(ps);
    update();
  }
  function remove() {
    if (ps.length <= 2 || lockX && (selected === 0 || selected === ps.length - 1)) return;
    remember();
    ps.splice(selected, 1);
    selected = Math.min(selected, ps.length - 1);
    constrain(ps);
    update();
  }
  function pointer(e) {
    const rect = canvas.getBoundingClientRect();
    return {
      x: e.clientX - rect.left,
      y: e.clientY - rect.top
    };
  }
  function snapPoint(i, v) {
    const step = Number($('editor-snap').value),
      round = n => Number((Math.round(n / step) * step).toPrecision(12));
    let x = round(v.x),
      y = round(v.y);
    // Preserve point ordering rather than nudging a snapped point off the grid.
    if (i && x <= ps[i - 1].x || i < ps.length - 1 && x >= ps[i + 1].x) x = ps[i].x;
    return {
      x,
      y
    };
  }
  canvas.addEventListener('pointerdown', e => {
    if (!geometry || drag) return;
    const at = pointer(e);
    if (e.button === 1) {
      e.preventDefault();
      if (at.x >= geometry.L && at.x <= geometry.R && at.y >= geometry.T && at.y <= geometry.B) {
        drag = {
          k: 'pan',
          start: at,
          bounds: {
            ...bounds
          }
        };
        canvas.setPointerCapture(e.pointerId);
        canvas.style.cursor = 'grabbing';
      }
      return;
    }
    if (e.button !== 0) return;
    const {
        X,
        Y
      } = geometry,
      p = ps[selected];
    let hit = null;
    if (mode === 'bezier') for (const k of ['l', 'r']) {
      if (k === 'l' && selected === 0 || k === 'r' && selected === ps.length - 1) continue;
      if (Math.hypot(X(p[k].x) - at.x, Y(p[k].y) - at.y) < 11) hit = {
        i: selected,
        k
      };
    }
    if (!hit) ps.forEach((p, i) => {
      if (Math.hypot(X(p.x) - at.x, Y(p.y) - at.y) < 11) hit = {
        i,
        k: 'point'
      };
    });
    if (!hit) return;
    e.preventDefault();
    canvas.focus({
      preventScroll: true
    });
    selected = hit.i;
    drag = {
      ...hit,
      before: state(),
      moved: false
    };
    canvas.setPointerCapture(e.pointerId);
    if (e.shiftKey && hit.k === 'point') {
      const v = snapPoint(hit.i, ps[hit.i]);
      if (v.x !== ps[hit.i].x || v.y !== ps[hit.i].y) {
        pushHistory(drag.before);
        drag.moved = true;
        movePoint(hit.i, v.x, v.y);
      }
    }
    update();
  });
  canvas.addEventListener('pointermove', e => {
    if (!drag) return;
    const at = pointer(e);
    if (drag.k === 'pan') {
      const b = drag.bounds,
        dx = (at.x - drag.start.x) / (geometry.R - geometry.L) * (geometry.x1 - geometry.x0),
        dy = (at.y - drag.start.y) / (geometry.B - geometry.T) * (geometry.y1 - geometry.y0);
      bounds = {
        x0: b.x0 - dx,
        x1: b.x1 - dx,
        y0: b.y0 + dy,
        y1: b.y1 + dy
      };
      draw();
      return;
    }
    let v = geometry.from(clamp(at.x, geometry.L, geometry.R), clamp(at.y, geometry.T, geometry.B));
    if (e.shiftKey && drag.k === 'point') v = snapPoint(drag.i, v);
    if (!drag.moved) {
      pushHistory(drag.before);
      drag.moved = true;
    }
    if (drag.k === 'point') movePoint(drag.i, v.x, v.y);else {
      ps[drag.i][drag.k] = v;
      if (ps[drag.i].linked) alignHandles(ps, drag.i, drag.k, clampY ? limits : null);else constrain(ps);
    }
    update();
  });
  for (const name of ['pointerup', 'pointercancel', 'lostpointercapture']) canvas.addEventListener(name, () => {
    drag = null;
    canvas.style.cursor = '';
  });
  canvas.addEventListener('auxclick', e => {
    if (e.button === 1) e.preventDefault();
  });
  canvas.addEventListener('dblclick', e => {
    const at = pointer(e);
    if (!geometry || at.x < geometry.L || at.x > geometry.R || at.y < geometry.T || at.y > geometry.B) return;
    const p = geometry.from(at.x, at.y);
    add(p.x, p.y);
  });
  canvas.addEventListener('contextmenu', e => {
    e.preventDefault();
    const at = pointer(e);
    if (!geometry || at.x < geometry.L || at.x > geometry.R || at.y < geometry.T || at.y > geometry.B) return;
    const p = geometry.from(at.x, at.y);
    add(p.x, p.y);
  });
  canvas.addEventListener('keydown', e => {
    if (!e.ctrlKey && !e.metaKey && !e.altKey && (e.key === 'Delete' || e.key === 'Backspace' || e.key.toLowerCase() === 'x')) {
      e.preventDefault();
      remove();
    }
  });
  $('editor-delete').onclick = remove;
  $('editor-mode').onchange = () => {
    remember();
    mode = $('editor-mode').value;
    update();
  };
  const confirmCurveAction = createConfirmation(document);
  $('editor-smooth').onclick = () => confirmCurveAction('Smooth all points?', 'This recalculates every handle and changes the curve shape. Your current handle adjustments will be replaced.', 'Smooth all', () => {
    remember();
    smooth(ps);
    fit();
    update();
  });
  $('editor-vector').onclick = () => {
    remember();
    const next = knots(ps.map(p => [p.x, p.y]));
    ps[selected].linked = false;
    ps[selected].l = next[selected].l;
    ps[selected].r = next[selected].r;
    update();
  };
  function setHandleMode(linked) {
    remember();
    ps[selected].linked = linked;
    if (linked) {
      const next = clone(ps);
      smooth(next);
      ps[selected].l = next[selected].l;
      ps[selected].r = next[selected].r;
    }
    update();
  }
  $('editor-handle-mode').onchange = () => setHandleMode($('editor-handle-mode').value === 'aligned');
  document.addEventListener('keydown', e => {
    if ($('curve-confirm').open || e.defaultPrevented || e.isComposing || e.repeat || e.ctrlKey || e.metaKey || e.altKey || e.key.toLowerCase() !== 'v' || $('editor-view').hidden || mode !== 'bezier' || drag) return;
    if (e.target.closest('input,textarea,select,[contenteditable]:not([contenteditable="false"])')) return;
    e.preventDefault();
    setHandleMode(!ps[selected].linked);
  });
  $('editor-auto').onclick = () => {
    remember();
    const next = clone(ps);
    smooth(next);
    ps[selected].l = next[selected].l;
    ps[selected].r = next[selected].r;
    fit();
    update();
  };
  function zoom(factor) {
    const mx = (bounds.x0 + bounds.x1) / 2,
      my = (bounds.y0 + bounds.y1) / 2,
      dx = (bounds.x1 - bounds.x0) * factor / 2,
      dy = (bounds.y1 - bounds.y0) * factor / 2;
    if (dx < 1e-10 || dy < 1e-10 || dx > 1e12 || dy > 1e12) return;
    bounds = {
      x0: mx - dx,
      x1: mx + dx,
      y0: my - dy,
      y1: my + dy
    };
    draw();
  }
  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    zoom(e.deltaY > 0 ? 1.15 : 1 / 1.15);
  }, {
    passive: false
  });
  $('editor-fit').onclick = () => {
    fit();
    draw();
  };
  $('editor-aspect').onclick = () => {
    squareView = !squareView;
    syncAspectButton();
    draw();
  };
  $('editor-selection').onchange = () => {
    selected = Number($('editor-selection').value);
    update();
  };
  function commitBounds(rescale = false) {
    try {
      const next = readBounds();
      if (!rescale && Object.keys(next).every(k => next[k] === limits[k])) {
        update();
        return;
      }
      remember();
      limits = next;
      if (rescale) fitCurve(ps, limits);
      applyBounds(ps, limits, lockX, clampY);
      fit();
      update();
    } catch (e) {
      $('bounds-status').textContent = e.message;
    }
  }
  $('bounds-preset').onchange = () => {
    const choice = $('bounds-preset').value;
    $('bound-custom').hidden = choice !== 'custom';
    if (choice !== 'custom') {
      for (const k of ['x0', 'y0']) $('bound-' + k).value = choice === 'unit' ? 0 : -1;
      for (const k of ['x1', 'y1']) $('bound-' + k).value = 1;
      commitBounds();
    }
  };
  for (const k of ['x0', 'x1', 'y0', 'y1']) {
    const field = $('bound-' + k);
    field.onchange = () => commitBounds();
    field.addEventListener('keydown', e => {
      if (e.key === 'Enter') {
        e.preventDefault();
        commitBounds();
      }
    });
  }
  $('bounds-fit').onclick = () => confirmCurveAction('Fit curve to bounds?', 'This rescales all points and handles to fit the selected bounds. Their coordinates will change.', 'Fit curve', () => commitBounds(true));
  $('lock-endpoints').onchange = () => {
    remember();
    lockX = $('lock-endpoints').checked;
    applyBounds(ps, limits, lockX, clampY);
    fit();
    update();
  };
  $('clamp-output').onchange = () => {
    remember();
    clampY = $('clamp-output').checked;
    applyBounds(ps, limits, lockX, clampY);
    fit();
    update();
  };
  for (const key of ['x', 'y', 'lx', 'ly', 'rx', 'ry']) {
    const field = $('editor-' + key);
    field.oninput = () => {
      const v = Number(field.value);
      if (field.value.trim() === '' || !Number.isFinite(v)) {
        message('Enter a finite number.', true);
        return;
      }
      if (key === 'x' && (selected && v <= ps[selected - 1].x || selected < ps.length - 1 && v >= ps[selected + 1].x)) {
        message('The In value must stay between neighboring points.', true);
        return;
      }
      remember();
      const p = ps[selected];
      if (key === 'x' || key === 'y') movePoint(selected, key === 'x' ? v : p.x, key === 'y' ? v : p.y);else {
        p[key[0]][key[1]] = v;
        if (p.linked) alignHandles(ps, selected, key[0], clampY ? limits : null);else constrain(ps);
      }
      fit();
      update(key);
    };
    field.onblur = () => update();
  }
  $('editor-tolerance').oninput = update;
  $('editor-count').oninput = () => update();
  $('editor-output-mode').onchange = () => update();
  $('editor-reference').onchange = drawPreview;
  document.addEventListener('keydown', e => {
    if ($('curve-confirm').open || e.defaultPrevented || e.isComposing || $('editor-view').hidden || !(e.ctrlKey || e.metaKey) || e.altKey || e.key.toLowerCase() !== 'z') return;
    // Keep native text undo in the JSON import field and other text editors.
    if (e.target.closest('textarea,[contenteditable="true"],input[type="text"]')) return;
    e.preventDefault();
    if (drag) return;
    const previous = e.shiftKey ? history.redo(state()) : history.undo(state());
    if (previous) restore(previous);
  });
  function syncImportOptions() {
    let nonlinear = false;
    try {
      const found = CurveInverse.extract(CurveInverse.parseClipboard($('editor-source').value));
      if (found.length === 1) {
        let leaf = found[0].curve;
        while ((leaf.Type === 'Not' || leaf.Type === 'Inverter') && leaf.Curve) leaf = leaf.Curve;
        nonlinear = leaf.Type !== 'Manual';
      }
    } catch {}
    $('import-sampling').hidden = !nonlinear;
    syncSamplingFields($('import-mode'), $('import-count'), $('import-tolerance'));
  }
  $('editor-source').addEventListener('input', syncImportOptions);
  $('import-mode').onchange = syncImportOptions;
  syncImportOptions();
  $('editor-import').onclick = () => {
    try {
      const found = CurveInverse.extract(CurveInverse.parseClipboard($('editor-source').value));
      if (found.length !== 1) throw Error('Paste a single curve chain.');
      const node = found[0].curve,
        c = CurveInverse.infer(node, 'terrain'),
        f = CurveInverse.compile(node);
      let leaf = node;
      while ((leaf.Type === 'Not' || leaf.Type === 'Inverter') && leaf.Curve) leaf = leaf.Curve;
      const points = leaf.Type === 'Manual' ? leaf.Points.map(p => p.In ?? p.X ?? p.x ?? p[0]).sort((a, b) => a - b).map(x => [x, f(x)]) : sampleFunction(f, c.h0, c.h1, $('import-mode').value, Number($('import-mode').value === 'count' ? $('import-count').value : $('import-tolerance').value));
      if (points.some(p => p.some(v => !Number.isFinite(v)))) throw Error('Non-finite value in the curve.');
      const next = knots(points);
      sample(next, 'linear');
      remember();
      ps = next;
      selected = 0;
      mode = 'linear';
      $('editor-count').value = Math.max(Number($('editor-count').value) || 32, ps.length);
      const incompatibleX = ps[0].x !== limits.x0 || ps.at(-1).x !== limits.x1,
        incompatibleY = ps.some(p => p.y < limits.y0 || p.y > limits.y1);
      if (incompatibleX) lockX = false;
      if (incompatibleY) clampY = false;
      syncBounds();
      fit();
      update();
      $('editor-import-status').className = 'hint';
      $('editor-import-status').textContent = (leaf.Type === 'Manual' ? 'Curve loaded with its vertices.' : 'Curve loaded: approximated with ' + points.length + ' points.' + ($('import-mode').value === 'tolerance' ? ' Precision is estimated numerically, not a guaranteed error bound.' : '')) + (incompatibleX || incompatibleY ? ' Incompatible locks were disabled to preserve the import. Use Fit curve to normalize it, then enable the locks.' : '');
    } catch (e) {
      $('editor-import-status').className = 'error';
      $('editor-import-status').textContent = e.message;
    }
  };
  $('editor-copy').onclick = async () => {
    try {
      if (!outputModel) throw Error('Correct the output settings before copying.');
      const text = JSON.stringify(CurveInverse.exportManual(outputModel.points, 'hytale'), null, 2);
      await copyText(text);
      message('Nodes copied · Manual Curve ' + (mode === 'bezier' ? 'with sampled Bezier segments.' : 'linear.'));
    } catch (e) {
      message(e.message, true);
    }
  };
  document.addEventListener('curve:viewchange', () => {
    draw();
    drawPreview();
  });
  new ResizeObserver(draw).observe(canvas);
  new ResizeObserver(drawPreview).observe($('editor-preview'));
  syncBounds();
  syncAspectButton();
  fit();
  update();
})(typeof globalThis === 'undefined' ? this : globalThis);


// ── app-shell.js ──

// Navigation is independent of both tools; switching never resets their sessions.
(function () {
  'use strict';

  const get = id => document.getElementById(id);
  function showView(editor) {
    get('editor-view').hidden = !editor;
    get('inverter-view').hidden = editor;
    get('tab-editor').setAttribute('aria-pressed', String(editor));
    get('tab-inverter').setAttribute('aria-pressed', String(!editor));
    document.querySelector('.view-switch').classList.toggle('is-inverter', !editor);
    document.querySelector('header h1').textContent = editor ? 'Curve editor' : 'Curve Inverter';
    document.querySelector('header p').textContent = editor ? 'Create and adjust your curves, then copy the Hytale nodes.' : 'Invert DistanceS, DistanceExponential, or monotonic Manual curves for your 2D noise.';
    document.dispatchEvent(new CustomEvent('curve:viewchange', {
      detail: {
        editor
      }
    }));
  }
  get('tab-editor').onclick = () => showView(true);
  get('tab-inverter').onclick = () => showView(false);
})();


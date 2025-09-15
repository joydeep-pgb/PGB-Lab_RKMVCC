#!/usr/bin/env python3
"""
Hybrid Web App:
Find 19–20 nt subsequences that start/end with G or C
AND pass a GC% threshold.
"""

from flask import Flask, request, render_template_string

app = Flask(__name__)

# ================= Core search logic =================
def gc_fraction(seq):
    """Return GC fraction (0–1)."""
    s = seq.upper()
    g = s.count("G")
    c = s.count("C")
    a = s.count("A")
    t = s.count("T")
    denom = a+t+g+c
    return (g+c)/denom if denom > 0 else 0

def find_gc_start_end(seq, min_len=19, max_len=20, gc_threshold=0.0):
    """Primer Finder"""
    seq = ''.join(ch for ch in seq if ch.isalpha()).upper()
    n = len(seq)
    hits = []

    for L in range(min_len, max_len+1):
        for i in range(0, n-L+1):
            sub = seq[i:i+L]
            if sub[0] in ("G","C") and sub[-1] in ("G","C"):
                gc = gc_fraction(sub)
                if gc >= gc_threshold:
                    hits.append({
                        "start": i+1,          # 1-based
                        "end": i+L,            # 1-based
                        "length": L,
                        "gc_percent": round(gc*100,2),
                        "sequence": sub
                    })
    return hits


# ================= Web UI =================
HTML_TEMPLATE = """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Primer Finder</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <style>
    body { background: #f9fafb; }
    .container { max-width: 900px; margin-top: 2em; }
    textarea { font-family: monospace; }
    table code { font-size: 0.9em; }
  </style>
</head>
<body>
<div class="container">
  <h2 class="mb-4 text-primary">Primer Finder</h2>
  <div class="card shadow-sm p-4 mb-4">
    <form method="POST" enctype="multipart/form-data">
      <div class="mb-3">
        <label class="form-label">Paste DNA Sequence</label>
        <textarea class="form-control" name="sequence" rows="6">{{sequence}}</textarea>
      </div>
      <div class="mb-3">
        <label class="form-label">Or upload FASTA file</label>
        <input type="file" class="form-control" name="fasta">
      </div>
      <div class="mb-3">
        <label class="form-label">GC% threshold (0–1)</label>
        <input type="number" step="0.01" class="form-control" name="gc_threshold" value="{{gc_threshold}}">
      </div>
      <button type="submit" class="btn btn-primary">Search</button>
    </form>
  </div>

  {% if hits is not none %}
    {% if hits|length > 0 %}
      <div class="card shadow-sm p-4">
        <h5 class="mb-3">Found {{hits|length}} sequences:</h5>
        <div class="table-responsive">
          <table class="table table-striped table-hover align-middle">
            <thead class="table-dark">
              <tr>
                <th>Start</th><th>End</th><th>Length</th><th>GC%</th><th>Sequence</th><th></th>
              </tr>
            </thead>
            <tbody>
              {% for h in hits %}
                <tr>
                  <td>{{h.start}}</td>
                  <td>{{h.end}}</td>
                  <td>{{h.length}}</td>
                  <td>{{h.gc_percent}}</td>
                  <td><code>{{h.sequence}}</code></td>
                  <td>
                    <button class="btn btn-sm btn-outline-secondary" onclick="copySeq('{{h.sequence}}')">Copy</button>
                  </td>
                </tr>
              {% endfor %}
            </tbody>
          </table>
        </div>
      </div>
    {% else %}
      <div class="alert alert-warning">No valid sequences found with current parameters.</div>
    {% endif %}
  {% endif %}
</div>

<script>
function copySeq(seq) {
  navigator.clipboard.writeText(seq).then(() => {
    alert("Copied: " + seq);
  });
}
</script>
</body>
</html>
"""

@app.route("/", methods=["GET","POST"])
def index():
    seq = ""
    hits = None
    gc_threshold = 0.0

    if request.method=="POST":
        seq = request.form.get("sequence","").strip()
        f = request.files.get("fasta")
        if f and f.filename:
            content=f.read().decode("utf-8")
            seq=''.join(l.strip() for l in content.splitlines() if not l.startswith(">"))

        gc_threshold = float(request.form.get("gc_threshold", 0.0))

        if seq:
            hits = find_gc_start_end(seq, gc_threshold=gc_threshold)

    return render_template_string(HTML_TEMPLATE,
                                  sequence=seq,
                                  gc_threshold=gc_threshold,
                                  hits=hits)

if __name__=="__main__":
    app.run(debug=True)

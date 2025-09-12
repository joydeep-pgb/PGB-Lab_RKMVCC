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
    return (g+c)/denom if denom>0 else 0

def find_gc_start_end(seq, min_len=19, max_len=20, gc_threshold=0.0):
    """
    Find subsequences of length 19–20 nt that start/end with G/C
    AND have GC% >= gc_threshold.
    """
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
<html>
<head>
  <title>GC-start/end 19–20 nt Finder with GC% Filter</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 2em; }
    textarea { width: 100%; height: 200px; }
    table { border-collapse: collapse; margin-top: 1em; width: 100%; }
    th, td { border: 1px solid #ccc; padding: 6px 10px; text-align: left; }
    code { background: #f5f5f5; padding: 2px 4px; }
  </style>
</head>
<body>
  <h2>Find 19–20 nt Sequences (Start & End with G/C + GC% Threshold)</h2>
  <form method="POST" enctype="multipart/form-data">
    <label>Paste DNA Sequence:</label><br>
    <textarea name="sequence">{{sequence}}</textarea><br><br>
    OR upload FASTA file: <input type="file" name="fasta"><br><br>
    GC% threshold (0–1): <input type="number" step="0.01" name="gc_threshold" value="{{gc_threshold}}">
    <br><br>
    <button type="submit">Search</button>
  </form>

  {% if hits is not none %}
    {% if hits|length > 0 %}
      <h3>Found {{hits|length}} sequences:</h3>
      <table>
        <tr><th>Start</th><th>End</th><th>Length</th><th>GC%</th><th>Sequence</th></tr>
        {% for h in hits %}
          <tr>
            <td>{{h.start}}</td>
            <td>{{h.end}}</td>
            <td>{{h.length}}</td>
            <td>{{h.gc_percent}}</td>
            <td><code>{{h.sequence}}</code></td>
          </tr>
        {% endfor %}
      </table>
    {% else %}
      <p><b>No valid sequences found with current parameters.</b></p>
    {% endif %}
  {% endif %}
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
        # file upload
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

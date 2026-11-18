#!/usr/bin/env python3
"""make_witnesses.py -- writes data/witnesses_738.json: the certificate that closure.js produces for C_738 from the five
starting values of the note (the text the page shows under "Show the certificate").  Runs closure.js in headless Firefox,
like test_closure.py.  usage: python3 make_witnesses.py"""
import os, json, test_closure as T
HERE = os.path.dirname(os.path.abspath(__file__))
def cert_text():
    case = {"coef": " ".join(map(str, T.C738)), "seeds": "0, -1/2, 1/2, -3/2, 1/3", "K": 16}
    page = """<!doctype html><meta charset="utf-8"><title>certificate</title><script>%s</script><script>
var c = %s, out; try { var st = Closure.start(c.coef, c.seeds, { K: c.K }); while (!st.step(1000)) {} out = { text: Closure.certText(st) }; } catch (e) { out = { error: String(e) }; }
var x = new XMLHttpRequest(); x.open('POST', '/result', false); x.send(JSON.stringify(out));
</script>"""
    r = T.run_page(page % (open(os.path.join(HERE, "closure.js"), encoding="utf-8").read(), json.dumps(case)))
    if "error" in r: raise SystemExit("closure.js failed: " + r["error"])
    return r["text"]
if __name__ == "__main__":
    t = cert_text(); open(os.path.join(HERE, "data", "witnesses_738.json"), "w", encoding="utf-8").write(t)
    print("data/witnesses_738.json: %d bytes, %d witnesses" % (len(t.encode("utf-8")), len(json.loads(t)["witnesses"])))

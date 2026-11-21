"""CI: the in-house solver must match Clarabel on the same problems.
Fuel within 0.01%, and it must never call a problem feasible that the reference proves infeasible."""
import json, subprocess, sys
subprocess.run(['node', 'bench/inner_js.js', '40', '/tmp/inner.json', '2500'], check=True)
subprocess.run([sys.executable, 'bench/inner_ref.py', '/tmp/inner.json', '/tmp/inner_cmp.json'], check=True, stderr=subprocess.DEVNULL)
d = json.load(open('/tmp/inner_cmp.json'))
unsafe = [r['i'] for r in d if r['js']['feasible'] and r['CLARABEL']['status'] != 'optimal']
gaps = [abs(r['js']['fuel'] - r['CLARABEL']['fuel']) / r['CLARABEL']['fuel'] for r in d if r['js']['feasible'] and r['CLARABEL']['status'] == 'optimal']
print(f'{len(gaps)} matched, worst fuel gap {max(gaps) * 100:.5f}%, false feasible {unsafe}')
sys.exit(0 if not unsafe and max(gaps) < 1e-4 else 1)

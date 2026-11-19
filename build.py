# Inlines the solver and flight model into the page: python3 build.py
import pathlib
s = pathlib.Path('src')
page = (s / 'index.src.html').read_text().replace('/*SOLVER*/', (s / 'solver.js').read_text()).replace('/*SIM*/', (s / 'sim.js').read_text())
pathlib.Path('index.html').write_text(page); print('index.html', len(page), 'bytes')

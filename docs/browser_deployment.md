# Browser calculation and GitHub Pages

Live application: https://bioinfo-ccnu.github.io/PocketSeekR/

The static website uses an ES module Web Worker to initialize Pyodide, load NumPy/SciPy/PyYAML and import a ZIP of the same `src/pocketseekr` package used by the CLI. The user's PDB text is passed to the worker and parsed in its in-memory filesystem. Prediction calls `pocketseekr.browser.predict_text`, which calls the native public API. No ligand coordinates are consumed. Cancellation terminates the worker; the next calculation starts a fresh runtime.

`web/runtime-lock.json` pins the runtime, numerical wheels and 3Dmol viewer versions and SHA256 hashes. `scripts/build_site.py` verifies downloaded assets and the upstream package dependency lock, then copies the required files into `outputs/site`. All application and runtime requests are same-origin static GETs. No API server or external runtime CDN is needed after publication. Browser cache behavior is controlled by the hosting service; offline operation is not promised.

The generated site is approximately 30.2 MiB, including a 3D viewer, a real RNA-only example and the algorithm PDF. Runtime downloads occur when calculation starts, not on initial page load. Browser computation uses one worker; large RNAs can be expensive in memory and time. The browser input limit is 20 MB, and ordinary PDB export field limits also apply. NMR inputs use the first model, following the native parser.

## Interface

The page uses a compact scientific-tool layout with a blue navigation bar, RNA input fieldset, interactive structure viewer, instructions and download resources. Its visual style is inspired by [SigSegmenter](http://bioinfo.isyslab.info/sigsegmenter/); the RNA banner and all PocketSeekR text are specific to this tool. Navigation example buttons load the same RNA-only input without starting calculation automatically. Result-table buttons highlight the corresponding RNA crop; the table reports both regional quality and the diversity-adjusted selection gain.

## Local build and preview

Run from the repository root with Python 3.10+:

```bash
python -m pip install -e '.[test]'
python scripts/build_site.py
python -m http.server 8000 --directory outputs/site
```

Open `http://localhost:8000/`. Serving over HTTP is required; do not open `index.html` as a `file://` URL. Generated assets and download cache are ignored by Git.

## Verification

```bash
python -m pytest -q
npm ci
npx playwright install chromium
npm run test:browser
```

On Linux, use `npx playwright install --with-deps chromium` to install system browser dependencies. An existing Chrome can be selected with `CHROME_PATH=/absolute/path/to/chrome`; a non-default Python can be selected with `PYTHON=/absolute/path/to/python`.

The test serves the artifact under `/PocketSeekR/`, exercises upload, cancellation/restart, viewer and exports, compares two real-structure reports to native Python, checks mobile overflow, rejects missing RNA, and verifies that all observed requests are same-origin GETs. Screenshots and comparison reports are saved to `outputs/browser-verification/`. To verify a deployment, run `SITE_URL=https://bioinfo-ccnu.github.io/PocketSeekR/ npm run test:browser` against the same source revision.

## Publication

The Pages workflow builds and tests on pull requests but publishes only from `main` or a manual workflow run. The repository's Pages source is **GitHub Actions**. The build job produces the static artifact after browser verification; the deploy job uses the `github-pages` environment with Pages write and OIDC permissions. See `.github/workflows/pages.yml`.

Changes to the Python package and web source are deployed together. Numerical behavior changes must update verification records and document their scientific implications. This deployment does not introduce a new benchmark or change the method's developmental evaluation scope.

Official references: [Pyodide deployment](https://pyodide.org/en/stable/usage/downloading-and-deploying.html), [module workers](https://pyodide.org/en/stable/usage/webworker.html), [3Dmol viewer API](https://3dmol.csb.pitt.edu/doc/GLViewer.html), [GitHub Pages workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).
